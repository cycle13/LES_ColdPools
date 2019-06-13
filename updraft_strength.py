import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import netCDF4 as nc
import argparse
import json as simplejson
import os

label_size = 12
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'


def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    # parser.add_argument("dTh", type=int)
    # parser.add_argument("--zparams", nargs='+', type=int)
    # parser.add_argument('--rparams', nargs='+', type=int)
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    # parser.add_argument("--k0")
    args = parser.parse_args()
    set_input_parameters(args)

    path_out_figs = os.path.join(path, 'figs_vorticity')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)

    ''' updraft generation from mass conservation: see importance of vorticity '''

    # a) read in CP radius r_torus(t) and CP spreading velocity u_torus(t) from tracer statistics
    # b) read in updrafts from angular average w(t,r,z)
    # c) compute vertical gradient of updrafts (as a function of t,r,z)
    # d) compute residual if only linear velocity
    # e) read in vorticity and compare to residual (??? how to best compute vorticity for mature CP?

    # # used files:
    # CP radius and spreading vel:      path/tracer_k0/output/coldpool_tracer_out.txt
    # >> use: plot_tracer_analysis_all.py>get_radius_vel

    # radial velocity:        path/fields_v_rad/v_rad.nc  >> v_rad
    # vertical velocity:      path/fields/t0.nc >> w
    #
    # radial + vertical velocity:      path/data_analsysis/stats_radial_average.nc >> v_rad, w



    krange = [0,1,2]
    nk = len(krange)
    k0 = 0
    nt = len(times)
    cp_id = 1  # circle ID that is used for statistics

    # a) read in CP radius r_torus(t)=r_av and CP spreading velocity u_torus(t)=U_rad_av from tracer statistics
    id = os.path.basename(path[:-1])
    fullpath_in = os.path.join(path, 'tracer_k' + str(k0), 'output')
    n_tracers = get_number_tracers(fullpath_in)
    n_cps = get_number_cps(fullpath_in)
    print('number of CPs: ', n_cps)
    print('number of tracers per CP: ', n_tracers)
    print ''

    dist_av = np.zeros((nt, nk))
    # r_av = np.zeros((nt, nk))
    # drdt_av = np.zeros((nt, nk))
    U_rad_av = np.zeros((nt, nk))
    for it, t0 in enumerate(times):
        print('---t0: ' + str(t0) + '---', it)
        dist_av[it, k0], U_rad_av[it, k0] = get_radius_vel(fullpath_in, it, cp_id, n_tracers, n_cps)
    r_av = dist_av * dx[0]
    del dist_av


    # b) read in updrafts w(t,r,z) and radial velocity v_rad(t,r,z) from angular average
    fullpath_in = os.path.join(path, 'data_analysis', 'stats_radial_averaged.nc')
    print fullpath_in
    rootgrp = nc.Dataset(fullpath_in, 'r')
    w_av = rootgrp.groups['stats'].variables['w'][:,:,:]                # nt, nr, nz
    v_rad_av = rootgrp.groups['stats'].variables['v_rad'][:,:,:]        # nt, nr, nz
    radius_rad_av = rootgrp.groups['stats'].variables['r'][:]           # nt, nr, nz
    time_rad_av = rootgrp.groups['timeseries'].variables['time'][:]     # nt
    rootgrp.close()

    k0_tracer = 0
    vel_at_rim = np.zeros((nt,nk), dtype=np.double)
    for it,t0 in enumerate(times):
        for k0 in krange:
            ir_tracer = np.where(radius_rad_av == np.int(np.round(r_av[it, k0_tracer], -2)))[0][0]
            vel_at_rim[it,k0] = v_rad_av[it,ir_tracer,k0]
    # test plot v_rad_av vs. r_av and U_rad_av
    plot_vel_at_rim(r_av, U_rad_av, radius_rad_av, v_rad_av, w_av, vel_at_rim, time_rad_av, k0_tracer, krange,
                        path_out_figs)




    # c) compute vertical gradient of updrafts (as a function of t,r,z)
    w_grad_av = np.zeros(shape=w_av.shape, dtype=np.double)
    dzi = 1. / dx[2]
    for k in range(1, w_av.shape[2]):
        w_grad_av[:, :, k] = dzi * (w_av[:, :, k] - w_av[:, :, k - 1])
    print''
    print('time, rad av: ', time_rad_av)
    print('times:        ', times)
    print''


    # d) read in width of CP torus ring
    fullpath_in = os.path.join(path, 'data_analysis', 'stats_radial_averaged_rimwidth.nc')
    rootgrp = nc.Dataset(fullpath_in, 'r')
    r_wmin = rootgrp.groups['rim_width'].variables['r_wmin'][:, :]      # r_wmin[nt, nz]
    r_wcenter = rootgrp.groups['rim_width'].variables['r_wcenter'][:, :]
    r_wmax = rootgrp.groups['rim_width'].variables['r_wmax'][:, :]
    R = 0.5*(r_wmax - r_wmin)
    rootgrp.close()

    fig_name = 'rim_width_test_k' + str(k0) + '.png'
    fig, axis = plt.subplots(2, 1, figsize=(9, 10), sharex='all')
    plt.tight_layout()
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    # # e) compute residual if only linear velocity
    # # !!! need to define width of CP torus
    # dUdr = np.zeros((nt, nk))
    # for it, t0 in enumerate(times):
    #     dUdr[it, k0] = U_rad_av[it,k0]*( 1./r_av[it,k0] + R[it,k0]/r_av[it,k0]**2 )
    #
    #
    # # f) compute residual
    # omega_res = np.zeros()
    # for it, t0 in enumerate(times_stats):
    #     omega_res[it,:] = w_grad_av[it, r_wcenter[it,k0],:] - dUdr[it,:]        # [t,r,z]
    #
    # # g) compute velocity form angular velocity
    # fullpath_in = os.path.join(path, 'data_analysis', 'stats_radial_averaged.nc')
    # rootgrp = nc.Dataset(fullpath_in, 'r')
    # times_stats = rootgrp.groups['timeseries'].variables['time'][:]
    # omega_plus =rootgrp.groups['rim_width'].variables['omega_plus'][:, :]  # r_wmin[nt, nz]
    # omega_min =rootgrp.groups['rim_width'].variables['omega_minus'][:, :]  # r_wmin[nt, nz]
    # rootgrp.close()
    #
    #
    # print(len(times_stats))
    # print(omega_res.shape)
    # print(r_wmin.shape)
    # fig_name = 'updraft_strength_t'+str(t0)+'.png'
    # fig, axis = plt.subplots(1,2)
    # ax = axis[0]
    # ax.plot(times_stats, omega_plus[:,k0], label='omega_plus')
    # ax.plot(times_stats, omega_min[:,k0], label='omega_minus')
    # ax.plot(times_stats, omega_res[:,k0], label='omega_res')
    # ax.legend()
    # ax = axis[1]
    # ax.plot(times, dUdr[:,k0], label='dUdr')
    # ax.plot(times, w_grad_av[:, r_wcenter, k0])
    # ax.plot(times, dUdr[:,k0]+omega_plus[:,k0], label='dUdr')
    # ax.legend()
    # fig.savefig(os.path.join(path_out_figs, fig_name))
    # plt.close(fig)

    return

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def plot_vel_at_rim(r_av, U_rad_av, radius_rad_av, v_rad_av, w_av, vel_at_rim, time_rad_av,
                    k0_tracer, krange, path_out_figs):
    rmax_plot = 9e3
    nt = len(times)
    irmax = np.where(radius_rad_av == rmax_plot)[0][0]
    for k0 in krange:
        fig_name = 'v_rad_test_k'+str(k0)+'.png'
        print path_out_figs
        fig, axis = plt.subplots(2,1, figsize=(9,10), sharex='all')
        ax0 = axis[0]
        ax1 = axis[1]
        # for it,t0 in enumerate(time_rad_av[1::2]):
        for it,t0 in enumerate(times[1::2]):
            count_color = 2 * np.double(it) / len(time_rad_av)
            ir_tracer = np.where(radius_rad_av == np.int(np.round(r_av[2*it+1,k0_tracer],-2)))[0][0]
            ax0.plot(radius_rad_av[:irmax], v_rad_av[2*it+1,:irmax,k0], color=cm.jet(count_color), label='t='+str(t0)+'s')
            ax1.plot(radius_rad_av[:irmax], w_av[2*it+1,:irmax,k0], color=cm.jet(count_color), label='t='+str(t0)+'s')
            ax0.plot([r_av[2*it+1,k0_tracer], r_av[2*it+1,k0_tracer]], [-10,10],
                     ':', color='0.25', linewidth=1)
            ax1.plot([r_av[2*it+1,k0_tracer], r_av[2*it+1,k0_tracer]], [-10,10],
                     ':', color='0.25', linewidth=1)
            ax0.plot(radius_rad_av[ir_tracer],v_rad_av[2*it+1,ir_tracer,k0], 'ko', markersize=5)
            ax1.plot(radius_rad_av[ir_tracer],w_av[2*it+1,ir_tracer,k0], 'ko', markersize=5)
        ax0.set_ylim(np.amin(v_rad_av[:,:irmax,k0]), np.amax(v_rad_av[:,:irmax,k0]))
        ax1.set_ylim(np.minimum(np.amin(w_av[:,:irmax,k0]), -np.amax(w_av[:,:irmax,k0])), np.amax(w_av[:,:irmax,k0]))
        ax0.legend(loc=1,ncol=2)
        ax1.legend(loc=4,ncol=2)
        ax0.set_title('radial velocity: v_rad')
        ax1.set_title('vertical velocity: w')
        # ax0.set_xlabel('radius r  [m]')
        ax1.set_xlabel('radius r  [m]')
        ax0.set_ylabel('v_rad  [m/s]')
        ax1.set_ylabel('w    [m/s]')
        plt.tight_layout()
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)

        fig_name = 'v_rad_vel_k'+str(k0)+'.png'
        fig, axis = plt.subplots(1, 3, figsize=(18, 5))
        ax0 = axis[0]
        ax1 = axis[1]
        ax2 = axis[2]
        for it, t0 in enumerate(times[1::2]):
            count_color = 2 * np.double(it) / len(time_rad_av)
            ir_tracer = np.where(radius_rad_av == np.int(np.round(r_av[2 * it + 1, k0_tracer], -2)))[0][0]
            ax0.plot(radius_rad_av[:irmax], v_rad_av[2 * it + 1, :irmax, k0], color=cm.jet(count_color), label='t='+str(t0)+'s')
            ax0.plot([r_av[2 * it + 1, k0_tracer], r_av[2 * it + 1, k0_tracer]], [-10, 10],
                     ':', color='0.25', linewidth=1)
            # ax0.plot([r_av[2 * it + 1, k0_tracer], r_av[2 * it + 1, k0_tracer]],
            #          [np.amin(v_rad_av[:, :irmax, k0]), np.amax(v_rad_av[:, :irmax, k0])], '-k', linewidth=1)
            ax0.plot(radius_rad_av[ir_tracer], v_rad_av[2 * it + 1, ir_tracer, k0], 'ko', markersize=5)
        ax1.plot(times, vel_at_rim[:,k0], label='radial velocity at tracer position')
        ax1.plot(times, np.amax(v_rad_av[:,:irmax,k0], axis=1), label='max(radial velocity)')
        ax1.plot(times, U_rad_av[:,k0_tracer], label='rim vel from tracer')
        ax2.plot(times, vel_at_rim[:,k0]/np.amax(v_rad_av[:,:irmax,k0], axis=1), label='vel at rim')
        ax2.plot(times, np.amax(v_rad_av[:,:irmax,k0], axis=1)/np.amax(v_rad_av[:,:irmax,k0], axis=1), label='vel at rim')
        ax2.plot(times, U_rad_av[:,k0_tracer]/np.amax(v_rad_av[:,:irmax,k0], axis=1), label='vel at rim')
        # ax1.plot(times, np.amax(v_rad_av[:,:irmax,k0], axis=1), label='max vel')
        ax2.fill_between(times, 0.6*np.ones(nt), 0.4*np.ones(nt), alpha=0.2, color='0.3')
        ax0.legend(loc=1, ncol=2)
        ax1.legend()
        ax0.set_ylim(np.amin(v_rad_av[:,:irmax,k0]), np.amax(v_rad_av[:,:irmax,k0]))
        ax1.set_xlim(times[0], times[-1])
        ax2.set_xlim(times[0], times[-1])
        ax1.grid()
        ax2.grid()
        ax0.set_title('radial velocity: v_rad')
        ax1.set_title('radial velocity at rim')
        ax2.set_title('(radial velocity at rim)/max(v_rad)')
        ax0.set_xlabel('radius r  [m]')
        ax1.set_xlabel('time  [s]')
        ax2.set_xlabel('time  [s]')
        ax0.set_ylabel('v_r  [m/s]')
        ax1.set_ylabel('w    [m/s]')
        plt.tight_layout()
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)

    return

# ----------------------------------------------------------------------
# ---------------------------- TRACER STATISTICS -----------------------

def get_number_tracers(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in + '/coldpool_tracer_out.txt', 'r')
    lines = f.readlines()
    count = 0
    # while CP age is 0 and CP ID is 1
    cp_age = int(lines[count].split()[0])
    cp_ID = int(lines[count].split()[3])
    print('cp_age', cp_age)
    while (cp_age == 1 and cp_ID == 1):
        count += 1
        cp_age = int(lines[count].split()[0])
        cp_ID = int(lines[count].split()[3])
    n_tracers = count
    f.close()

    return n_tracers



def get_number_cps(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in + '/coldpool_tracer_out.txt', 'r')
    lines = f.readlines()
    # count = 0
    # # while CP age is 0 and CP ID is 1
    # while (int(lines[count].split()[0]) == 1):
    #     count += 1
    # cp_number = int(lines[count-1].split()[3])
    cp_number = int(lines[-1].split()[3])

    f.close()

    return cp_number

def get_radius_vel(fullpath_in, t0, cp_id, n_tracers, n_cps):
    # print('in', fullpath_in)
    f = open(fullpath_in + '/coldpool_tracer_out.txt', 'r')
    # f = open(DIR+EXPID+'/'+child+'/output/irt_tracks_output_pure_sort.txt', 'r')
    lines = f.readlines()
    count = 0
    dist = []
    vel = []

    count = t0 * n_cps * n_tracers + (cp_id - 1)*n_tracers
    # while CP age is 0 and CP ID is cp_id
    timestep = int(lines[count].split()[0])
    cp_ID = int(lines[count].split()[3])
    # print(timestep, cp_ID)
    while (timestep-1 == t0 and int(lines[count].split()[3])==cp_id):
        columns = lines[count].split()
        dist.append(float(columns[8]))
        # vel.append(np.sqrt(float(columns[10])**2 + float(columns[11])**2))
        vel.append(float(columns[12]))
        count += 1
        timestep = int(lines[count].split()[0])
    f.close()
    r_av = np.average(dist)
    vel_av = np.average(vel)

    return r_av, vel_av
# ----------------------------------------------------------------------

def set_input_parameters(args):
    print('--- set input parameters ---')
    global path, path_fields, case_name
    path = args.path
    path_fields = os.path.join(path, 'fields')
    case_name = args.casename

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dV, gw
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = np.zeros(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx[0] * dx[1] * dx[2]



    global tmin, tmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100

    ''' time range '''
    global times
    times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
             and tmin <= np.int(name[:-3]) <= tmax]
    times.sort()
    print('tmin, tmax: ', tmin, tmax)
    print('times: ', times)
    print('')
    return
# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()