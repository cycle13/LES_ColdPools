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
    parser.add_argument("--rmax")
    parser.add_argument("--path_tracers")
    args = parser.parse_args()
    set_input_parameters(args)

    path_out_figs = os.path.join(path, 'figs_rim')
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

    # FIELDS
    # r_tracers_av:             tracer radius
    # U_rad_av:                 tracer velocity
    # time_av:                  times at which v_rad_av computed
    # r_av[nr]                  radius at which v_rad_av computed
    # v_rad_av[nt,nr,nz]        average radial velocity from 3D fields
    # w_av[nt,nr,nz]            average vertical velocit from 3D fields
    # v_rad_av_at_rim:          v_rad_av[t,r==r_tracers_av(t),z]

    krange = [0,1,2]
    nk = len(krange)
    k0 = 0
    nt = len(times)
    if args.rmax:
        irmax = np.int(np.int(args.rmax)/dx[0])
    else:
        irmax = np.int(10e3/dx[0])
    cp_id = 1  # circle ID that is used for statistics
    print('krange: '+str(krange))
    print ''


    # a) read in CP radius r_torus(t)=r_av and CP spreading velocity u_torus(t)=U_rad_av from tracer statistics
    if args.path_tracers:
        fullpath_in = os.path.join(path, args.path_tracers, 'output')
    else:
        fullpath_in = os.path.join(path, 'tracer_k' + str(k0), 'output')
    ID = os.path.basename(path[:-1])
    n_tracers = get_number_tracers(fullpath_in)
    n_cps = get_number_cps(fullpath_in)
    print('number of CPs: ', n_cps)
    print('number of tracers per CP: ', n_tracers)
    print ''

    dist_av = np.zeros((nt, nk))
    U_rad_av = np.zeros((nt, nk))
    for it, t0 in enumerate(times):
        print('---t0: ' + str(t0) + '---', it)
        dist_av[it, k0], U_rad_av[it, k0] = get_radius_vel(fullpath_in, it, cp_id, n_tracers, n_cps)
    r_tracers_av = dist_av * dx[0]
    del dist_av


    # b) read in updrafts w(t,r,z) and radial velocity v_rad(t,r,z) from angular average
    fullpath_in = os.path.join(path, 'data_analysis', 'stats_radial_averaged.nc')
    print fullpath_in
    rootgrp = nc.Dataset(fullpath_in, 'r')
    w_av = rootgrp.groups['stats'].variables['w'][:,:,:]                # nt, nr, nz
    v_rad_av = rootgrp.groups['stats'].variables['v_rad'][:,:,:]        # nt, nr, nz
    v_tan_av = rootgrp.groups['stats'].variables['v_tan'][:,:,:]        # nt, nr, nz
    r_av = rootgrp.groups['stats'].variables['r'][:]                    # nr
    time_av = rootgrp.groups['timeseries'].variables['time'][:]         # nt
    rootgrp.close()

    k0_tracer = 0
    v_rad_av_at_rim = np.zeros((nt,nk), dtype=np.double)
    print('v_rad_av, v_rad_av_at_rim', v_rad_av.shape, v_rad_av_at_rim.shape, nt, times)
    for it,t0 in enumerate(times):
        for k0 in krange:
            ir_tracer = np.where(r_av == np.int(np.round(r_tracers_av[it, k0_tracer], -2)))[0][0]
            v_rad_av_at_rim[it,k0] = v_rad_av[it,ir_tracer,k0]
    # test plot v_rad_av vs. r_av and U_rad_av
    plot_vel_at_rim(r_tracers_av, U_rad_av, r_av, v_rad_av, w_av, v_rad_av_at_rim,
                    time_av, k0_tracer, krange, irmax,
                    path_out_figs)



    # c) compute vertical gradient of updrafts:  w_grad_av[t,r,z]
    #       >> w_grad_av is at levels z (like u_rad, u, v  and scalars)
    #       - read in zrange from Stats-file to get right position for gradient
    #           >> Positions in PyCLES: half-levels = box centres (scalars, u, v); full levels = box sides (w)
    #           in Stats-file output: z_half >> z, z >> z_full (same for rho0, etc)
    # >> scalars, u, v:     z[0] = dz/2, z[1] = dz + dz/2, ..., z[k] = (k+1/2)*dz
    # >> w:                 z_full[0] = dz, z[1] = 2*dz, ..., z[k] = (k+1)*dz
    statsfile = nc.Dataset(os.path.join(path, 'stats', 'Stats.' + case_name + '.nc'), 'r')
    z_full = statsfile.groups['reference'].variables['z_full'][:]  # height of w
    z = statsfile.groups['reference'].variables['z'][:]  # height of scalars, u, v
    statsfile.close()
    dzi = 1. / dx[2]
    w_grad_av = np.zeros(shape=w_av.shape, dtype=np.double)
    w_grad_av[:,:,0] = dzi * w_av[:,:,0]        # w=0 at surface (k=-1)
    for k in range(1, w_av.shape[2]):
        w_grad_av[:, :, k] = dzi * (w_av[:, :, k] - w_av[:, :, k - 1])
    print''

    fig_name = 'w_grad_test.png'
    fig, axis = plt.subplots(2, len(krange), figsize=(5*len(krange)+1, 10), sharex='all', sharey='all')
    for k in krange:
        ax = axis[0,k]
        cf = ax.contourf(r_av[:irmax], times[:nt], w_av[:nt,:irmax,k])
        plt.colorbar(cf, ax=ax)
        ax.set_title('w (k='+str(k)+', z='+str(z[k])+')')
        ax.set_xlabel('x')
        ax.set_ylabel('time')
        ax = axis[1,k]
        cf = ax.contourf(r_av[:irmax], times[:nt], w_grad_av[:nt,:irmax,k])
        plt.colorbar(cf, ax=ax)
        ax.set_title('dw/dz (k='+str(k)+', z='+str(z[k])+')')
        ax.set_xlabel('x')
        ax.set_ylabel('time')
    plt.tight_layout()
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)


    # d) read in width of CP torus ring
    fullpath_in = os.path.join(path, 'data_analysis', 'stats_radial_averaged_rimwidth.nc')
    rootgrp = nc.Dataset(fullpath_in, 'r')
    r_wmin = rootgrp.groups['rim_width'].variables['r_wmin'][:, :]      # r_wmin[nt, nz]
    r_wcenter = rootgrp.groups['rim_width'].variables['r_wcenter'][:, :]
    r_wmax = rootgrp.groups['rim_width'].variables['r_wmax'][:, :]
    R = 0.5*(r_wmax - r_wmin)
    rootgrp.close()

    # fig_name = 'model_urad_w.png'
    # k0_plot = 0
    # t0 = 1800
    # plot_model_figure(r_tracers_av, U_rad_av, r_av,
    #                   v_rad_av, w_av,
    #                   r_wmin, r_wcenter, r_wmax,
    #                   times, t0, k0_plot, fig_name, path_out_figs)


    # # plot different radii from tracers and rim_width.py
    # k0 = 0
    # plot_rim_width(r_tracers_av, U_rad_av, r_av, v_rad_av, w_av,
    #                r_wmax, r_wmin, r_wcenter,
    #                k0, k0_tracer, time_av, irmax, path_out_figs)



    # v_rot=omega*R defined as difference between max(v_rad)=v_rad(r_c) and v_rad(r_tracer)
    v_rot = np.zeros(nt)
    v_rad_max = np.amax(v_rad_av[:, :, 0], axis=1)
    for it,t0 in enumerate(times):
        ir_tracer = np.where(r_av == np.int(np.round(r_tracers_av[it, k0_tracer], -2)))[0][0]
        v_rot[it] = v_rad_max[it] - v_rad_av[it, ir_tracer, 0]
    # radius of vortex ring defined as R=rw_max - rw_center (or R*=0.5*(rw_max - rw_min))
    R_ = r_wmax[:,0] - r_wcenter[:,0]
    R = 0.5*(r_wmax - r_wmin)
    # iR_ = np.double(R)/dx[2]
    ir = np.zeros(nt)
    w_rot = np.zeros(nt)
    w_rot_ = np.zeros(nt)
    w_rot_max = np.zeros(nt)
    print('eeeeeeeeeeeeeeeeeeee times: ', times)
    for it,t0 in enumerate(times):
        print('---t0: ' + str(t0) + '---', it)
        ir_tracer = np.where(r_av == np.int(np.round(r_tracers_av[it, k0_tracer], -2)))[0][0]
        if R[it,0] >= 0:
            ir[it] = np.int(np.double(R[it,0])/dx[0])
            print it, R[it,0], ir[it]
            w_rot[it] = w_av[it, ir_tracer, ir[it]]
        else:
            print R[it,:]
            R[it,:] = 0.0
            w_rot[it] = w_av[it, ir_tracer, 0]
            # w_rot_[it] = w_av[it, ir_tracer, 0]
        w_rot_[it] = np.amax(w_av[it, ir_tracer, :])
        w_rot_max[it] = np.amax(w_av[it, :, :])

    fig_name = 'model_urot.png'
    fig, axis = plt.subplots(1, 2, figsize=(10, 4), sharex='all')
    ax = axis[0]
    ax.plot(times, v_rot, label='v_rot')
    ax.plot(times, w_rot, label='w_av[r=r_tracer, z=R]')
    ax.legend()
    ax.set_xlabel('time')
    ax.grid()
    ax = axis[1]
    ax.plot(times, v_rot, label='v_rot')
    ax.plot(times, w_rot, label='w_av[r=r_tracer, z=R]')
    ax.plot(times, w_rot_, label='max(w_av[r=r_tracer])')
    ax.plot(times, w_rot_max, label='max(w_av)')
    ax.legend()
    ax.grid()
    ax.set_xlabel('time')
    plt.suptitle('v_rot = max(v_rad)-v_rad[r_tracer,z=0] = omega*R')
    # ax.plot(times, w_rot_)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)




    # # e) compute residual if only linear velocity (assume that U maximum
    # # !!! need to define width of CP torus
    # dUdr = np.zeros((nt, nk))
    # # for it, t0 in enumerate(times):
    # #     dUdr[it, k0] = U_rad_av[it,k0]*( 1./r_av[it,k0] + R[it,k0]/r_av[it,k0]**2 )
    # for it, t0 in enumerate(times):
    #     dUdr[it,k0] = 0.5 * U_rad_av[it,k0] / (r_tracers_av[it,k0_tracer] - r_wcenter[it,k0])
    #
    # fig_name = 'dUdr_test_' + str(k0) + '.png'
    # nt = len(times)
    # fig, axis = plt.subplots(2, 1, figsize=(9, 10), sharex='all')
    # ax0 = axis[0]
    # ax1 = axis[1]
    # ax0.plot(time_av, dUdr[:,k0], label='dUdr')
    # for it,t0 in enumerate(times):
    #     ax1.plot(r_av[:irmax], w_grad_av[it,:irmax,k0])
    # ax0.set_title('r(w=max) & tracer radius')
    # ax1.set_title('r(w=min), r(w=0) & tracer radius')
    # ax1.set_xlabel('radius r  [m]')
    # ax0.set_ylabel('v_rad  [m/s]')
    # ax1.set_ylabel('w  [m/s]')
    # plt.tight_layout()
    # fig.savefig(os.path.join(path_out_figs, fig_name))
    # plt.close(fig)
    #
    #
    #
    # # # f) compute residual
    # # omega_res = np.zeros()
    # # for it, t0 in enumerate(times_stats):
    # #     omega_res[it,:] = w_grad_av[it, r_wcenter[it,k0],:] - dUdr[it,:]        # [t,r,z]
    # #
    # # # g) compute velocity form angular velocity
    # # fullpath_in = os.path.join(path, 'data_analysis', 'stats_radial_averaged.nc')
    # # rootgrp = nc.Dataset(fullpath_in, 'r')
    # # times_stats = rootgrp.groups['timeseries'].variables['time'][:]
    # # omega_plus =rootgrp.groups['rim_width'].variables['omega_plus'][:, :]  # r_wmin[nt, nz]
    # # omega_min =rootgrp.groups['rim_width'].variables['omega_minus'][:, :]  # r_wmin[nt, nz]
    # # rootgrp.close()
    # #
    # #
    # # print(len(times_stats))
    # # print(omega_res.shape)
    # # print(r_wmin.shape)
    # # fig_name = 'updraft_strength_t'+str(t0)+'.png'
    # # fig, axis = plt.subplots(1,2)
    # # ax = axis[0]
    # # ax.plot(times_stats, omega_plus[:,k0], label='omega_plus')
    # # ax.plot(times_stats, omega_min[:,k0], label='omega_minus')
    # # ax.plot(times_stats, omega_res[:,k0], label='omega_res')
    # # ax.legend()
    # # ax = axis[1]
    # # ax.plot(times, dUdr[:,k0], label='dUdr')
    # # ax.plot(times, w_grad_av[:, r_wcenter, k0])
    # # ax.plot(times, dUdr[:,k0]+omega_plus[:,k0], label='dUdr')
    # # ax.legend()
    # # fig.savefig(os.path.join(path_out_figs, fig_name))
    # # plt.close(fig)

    return

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def plot_model_figure(r_tracers_av, U_rad_av, r_av,
                      v_rad_av, w_av,
                      r_wmin, r_wcenter, r_wmax,
                      times, t0, k0, fig_name, path_out_figs):
    irmax = np.int(7e3/dx[0])
    it = np.int(t0/dt_fields)
    k0_tracer = 0
    ir_tracer = np.where(r_av == np.int(np.round(r_tracers_av[it, k0_tracer], -2)))[0][0]
    ir_max = np.where(r_av == np.int(np.round(r_wmax[it, k0], -2)))[0][0]
    ir_min = np.where(r_av == np.int(np.round(r_wmin[it, k0], -2)))[0][0]
    ir_c = np.where(r_av == np.int(np.round(r_wcenter[it, k0], -2)))[0][0]

    romps = np.zeros(r_av.shape)
    romps[:ir_tracer+1] = v_rad_av[it, ir_tracer, k0] / r_tracers_av[it,k0] * r_av[:ir_tracer+1]

    fig, ax = plt.subplots(1, 1, figsize=(10, 6), sharex='all')
    ax.plot([0, r_av[-1]], [0,0], 'k', linewidth=1)
    ax.plot([r_av[ir_tracer], r_av[ir_tracer]], [-1,5], 'k', linewidth=1)
    ax.plot([r_av[ir_min], r_av[ir_min]], [-1,5], 'k', linewidth=1)
    ax.plot([r_av[ir_c], r_av[ir_c]], [-1,5], 'k', linewidth=1)
    ax.plot(r_av, v_rad_av[it, :, k0], color='k', label='v_rad', linewidth=3)
    ax.plot(r_av, w_av[it, :, k0], color='0.5', label='w', linewidth=3)

    ax.plot(r_av, romps, 'b', label='linear model (Romps)')

    ax.plot(r_av[ir_tracer], v_rad_av[it, ir_tracer, k0],'s', color='0.5', markersize=9, markeredgewidth=0., label='r tracer')
    ax.plot(r_av[ir_tracer], w_av[it, ir_tracer, k0],'s', color='0.5', markersize=9, markeredgewidth=0., )
    ax.plot(r_wmax[it, k0], v_rad_av[it, ir_max, k0], 'ko', markersize=6, label='r(w=max)')
    ax.plot(r_wmax[it, k0], w_av[it, ir_max, k0], 'ko', markersize=6)
    ax.plot(r_wmin[it, k0], v_rad_av[it, ir_min, k0], 'b^', markersize=9, label='r(w=min)')
    ax.plot(r_wmin[it, k0], w_av[it, ir_min, k0], 'b^', markersize=9)
    ax.plot(r_wcenter[it, k0], v_rad_av[it, ir_c, k0], 'kv', markersize=9, label='r(w=0)')
    ax.plot(r_wcenter[it, k0], w_av[it, ir_c, k0], 'kv', markersize=9)
    ax.set_xlabel('radius r [m]')
    ax.set_ylabel('radial / vertical velocity [m/s]')
    ax.set_xlim(0, r_av[irmax])
    if dx[0] == 100:
        ax.set_ylim(-0.6, 3.)
    elif dx[0] == 50:
        ax.set_ylim(-0.6, 4)
    else:
        ax.set_ylim(-0.6, 4)
    ax.legend(loc='best')
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return


def plot_vel_at_rim(r_av, U_rad_av, radius_rad_av,
                    v_rad_av, w_av, vel_at_rim, time_rad_av,
                    k0_tracer, krange, irmax, path_out_figs):
    nt = len(times)
    for k0 in krange:
        fig_name = 'v_rad_test_k'+str(k0)+'.png'
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
            if it == 0:
                ax0.plot(radius_rad_av[ir_tracer],U_rad_av[2*it+1,k0], 'o', color='0.5', markersize=6, label='tracer velocity')
                ax0.plot(radius_rad_av[ir_tracer],v_rad_av[2*it+1,ir_tracer,k0], 'ko', markersize=6, label='tracer radius')
                ax1.plot(radius_rad_av[ir_tracer],w_av[2*it+1,ir_tracer,k0], 'ko', markersize=6, label='tracer radius')
            else:
                ax0.plot(radius_rad_av[ir_tracer],U_rad_av[2*it+1,k0], 'o', color='0.5', markersize=6)
                ax0.plot(radius_rad_av[ir_tracer],v_rad_av[2*it+1,ir_tracer,k0], 'ko', markersize=6)
                ax1.plot(radius_rad_av[ir_tracer],w_av[2*it+1,ir_tracer,k0], 'ko', markersize=6)
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
        ax1.plot(times, np.amax(v_rad_av[:nt,:irmax,k0], axis=1), label='max(radial velocity)')
        ax1.plot(times, U_rad_av[:nt,k0_tracer], label='rim vel from tracer')
        ax2.plot(times, vel_at_rim[:nt,k0]/np.amax(v_rad_av[:nt,:irmax,k0], axis=1), label='vel at rim')
        ax2.plot(times, np.amax(v_rad_av[:nt,:irmax,k0], axis=1)/np.amax(v_rad_av[:nt,:irmax,k0], axis=1), label='vel at rim')
        ax2.plot(times, U_rad_av[:nt,k0_tracer]/np.amax(v_rad_av[:nt,:irmax,k0], axis=1), label='vel at rim')
        ax1.plot(times, np.amax(v_rad_av[:nt,:irmax,k0], axis=1), label='max vel')
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


def plot_rim_width(r_tracers_av, U_rad_av, r_av, v_rad_av, w_av,
                   r_wmax, r_wmin, r_wcenter,
                   k0, k0_tracer, time_rad_av, irmax, path_out_figs):
    nt = len(times)
    fig_name = 'rim_width_test_k' + str(k0) + '.png'
    fig, axis = plt.subplots(2, 2, figsize=(20, 10), sharex='all')
    ax0 = axis[0, 0]
    ax1 = axis[1, 0]
    ax2 = axis[0, 1]
    ax3 = axis[1, 1]
    for it, t0 in enumerate(times[1::2]):
        count_color = 2 * np.double(it) / len(time_rad_av)
        ax0.plot(r_av[:irmax], v_rad_av[2 * it + 1, :irmax, k0], color=cm.jet(count_color))
        ax1.plot(r_av[:irmax], w_av[2 * it + 1, :irmax, k0], color=cm.jet(count_color))
        ax1.plot(r_av[:irmax], w_av[2 * it + 1, :irmax, k0], color=cm.jet(count_color))
        ir = np.where(r_av == np.int(np.round(r_wmax[2 * it + 1, k0], -2)))[0][0]
        ir_tracer = np.where(r_av == np.int(np.round(r_tracers_av[2 * it + 1, k0_tracer], -2)))[0][0]
        if it == 0:
            ax0.plot(r_av[ir_tracer], v_rad_av[2 * it + 1, ir_tracer, k0],
                     'ok', markersize=6, label='r tracer')
            ax1.plot(r_av[ir_tracer], w_av[2 * it + 1, ir_tracer, k0], 'ok', markersize=6, label='r tracer')
            ax0.plot(r_wmax[2 * it + 1, k0], v_rad_av[2 * it + 1, ir, k0], 'o', color='0.5', markersize=6,
                     label='r(w=max)')
            ax1.plot(r_wmax[2 * it + 1, k0], w_av[2 * it + 1, ir, k0], 'o', color='0.5', markersize=6, label='r(w=max)')
        else:
            ax0.plot(r_av[ir_tracer], v_rad_av[2 * it + 1, ir_tracer, k0], 'ko', markersize=10)
            ax1.plot(r_av[ir_tracer], w_av[2 * it + 1, ir_tracer, k0], 'ko', markersize=10)
            ax0.plot(r_wmax[2 * it + 1, k0], v_rad_av[2 * it + 1, ir, k0], 'o', color='0.5', markersize=6)
            ax1.plot(r_wmax[2 * it + 1, k0], w_av[2 * it + 1, ir, k0], 'o', color='0.5', markersize=6)

        ax2.plot(r_av[:irmax], v_rad_av[2 * it + 1, :irmax, k0], color=cm.jet(count_color))
        ax3.plot(r_av[:irmax], w_av[2 * it + 1, :irmax, k0], color=cm.jet(count_color))
        ir = np.where(r_av == np.int(np.round(r_wmin[2 * it + 1, k0], -2)))[0][0]
        ir_c = np.where(r_av == np.int(np.round(r_wcenter[2 * it + 1, k0], -2)))[0][0]
        ir_tracer = np.where(r_av == np.int(np.round(r_tracers_av[2 * it + 1, k0_tracer], -2)))[0][0]
        if it == 0:
            ax2.plot(r_av[ir_tracer], v_rad_av[2 * it + 1, ir_tracer, k0], 'ko', markersize=10, label='r tracer')
            ax3.plot(r_av[ir_tracer], w_av[2 * it + 1, ir_tracer, k0], 'ko', markersize=10, label='r tracer')
            ax2.plot(r_wmin[2 * it + 1, k0], v_rad_av[2 * it + 1, ir, k0], 'o', color='0.5', markersize=6, label='r(w=min)')
            ax3.plot(r_wmin[2 * it + 1, k0], w_av[2 * it + 1, ir, k0], 'o', color='0.5', markersize=6, label='r(w=min)')
            ax2.plot(r_wcenter[2 * it + 1, k0], v_rad_av[2 * it + 1, ir_c, k0], 'kd', markersize=6, label='r(w=0)')
            ax3.plot(r_wcenter[2 * it + 1, k0], w_av[2 * it + 1, ir_c, k0], 'kd', markersize=6, label='r(w=0)')
        else:
            ax2.plot(r_av[ir_tracer], v_rad_av[2 * it + 1, ir_tracer, k0], 'ko', markersize=10)
            ax3.plot(r_av[ir_tracer], w_av[2 * it + 1, ir_tracer, k0], 'ko', markersize=10)
            ax2.plot(r_wmin[2 * it + 1, k0], v_rad_av[2 * it + 1, ir, k0], 'o', color='0.5', markersize=6)
            ax3.plot(r_wmin[2 * it + 1, k0], w_av[2 * it + 1, ir, k0], 'o', color='0.5', markersize=6)
            ax2.plot(r_wcenter[2 * it + 1, k0], v_rad_av[2 * it + 1, ir_c, k0], 'kd', markersize=6)
            ax3.plot(r_wcenter[2 * it + 1, k0], w_av[2 * it + 1, ir_c, k0], 'kd', markersize=6)

    ax0.legend()
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax0.set_xlim(0,r_av[irmax])
    ax1.set_xlim(0,r_av[irmax])
    ax2.set_xlim(0,r_av[irmax])
    ax3.set_xlim(0,r_av[irmax])
    ax0.set_title('radial velocity with r(w=max) & tracer radius')
    ax1.set_title('vertical velocity with r(w=min), r(w=0) & tracer radius')
    ax2.set_title('radial velocity with r(w=max) & tracer radius')
    ax3.set_title('vertical velocity with r(w=min), r(w=0) & tracer radius')
    ax1.set_xlabel('radius r  [m]')
    ax3.set_xlabel('radius r  [m]')
    ax0.set_ylabel('v_rad  [m/s]')
    ax1.set_ylabel('w  [m/s]')
    plt.tight_layout()
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)



    var_list = ['w', 'v_rad']
    # var_list = ['v_rad']
    for var_name in var_list:
        fig_name = 'rim_width_test_'+var_name+'.png'
        if var_name == 'w':
            var = w_av
        elif var_name == 'v_rad':
            var = v_rad_av
        fig, axis = plt.subplots(4, 1, figsize=(10, 15), sharex='all')
        plt.tight_layout()
        if dx[0] == 100:
            krange = [0, 1, 2, 3]
        elif dx[0] == 50:
            krange = [1, 3, 5, 7]
            krange = [0, 1, 2, 3]
        elif dx[0] == 25:
            krange = [2, 6, 10, 14]
        for k,k_ in enumerate(krange):
            ax = axis[k]
            for it, t0 in enumerate(times[1::2]):
                count_color = 2 * np.double(it) / len(time_rad_av)
                ir_max = np.where(r_av == np.int(np.round(r_wmax[2*it + 1, k0], -2)))[0][0]
                ir_min = np.where(r_av == np.int(np.round(r_wmin[2*it + 1, k0], -2)))[0][0]
                ir_c = np.where(r_av == np.int(np.round(r_wcenter[2*it + 1, k0], -2)))[0][0]
                ir_tracer = np.where(r_av == np.int(np.round(r_tracers_av[2 * it + 1, k0_tracer], -2)))[0][0]
                ax.plot(r_av[:irmax], var[2*it+1,:irmax,k], color=cm.jet(count_color), label='t='+str(t0))

                if it == 0:
                    ax.plot(r_av[ir_tracer], var[2*it+1, ir_tracer, k], 'o', color='0.5', markersize=10, label='r tracer (k=0)')
                    ax.plot(r_wmax[2*it+1, k0], var[2*it+1, ir_max, k], 'ko', markersize=6, label='r(w=max)')
                    ax.plot(r_wmin[2*it+1, k0], var[2*it+1, ir_min, k], 'kd', markersize=6, label='r(w=min)')
                    ax.plot(r_wcenter[2*it+1, k0], var[2*it+1, ir_c, k], 'kx', markersize=6, label='r(w=0)')
                else:
                    ax.plot(r_av[ir_tracer], var[2*it+1, ir_tracer, k], 'o', color='0.5', markersize=10)
                    ax.plot(r_wmax[2*it+1, k0], var[2*it+1, ir_max, k], 'ko', markersize=6)
                    ax.plot(r_wmin[2*it+1, k0], var[2*it+1, ir_min, k], 'kd', markersize=6)
                    ax.plot(r_wcenter[2*it+1, k0], var[2*it+1, ir_c, k], 'kx', markersize=6)
            # ax.legend()
            ax.set_title('z='+str((k_+1)*dx[2])+ ' (k='+str(k_)+')')
            ax.set_xlim(0, r_av[irmax])
        axis[-1].legend(loc='best', ncol=4)
        if var_name == 'w':
            if dx[0] == 100:
                axis[0].set_ylim(-1.5,1.6)
                axis[1].set_ylim(-1.5,1.6)
                axis[2].set_ylim(-1,1.2)
                axis[3].set_ylim(-0.7,.9)
            else:
                axis[0].set_ylim(-1.5,1.5)
                axis[1].set_ylim(-2,2.)
                axis[2].set_ylim(-2,2.2)
                axis[3].set_ylim(-2,2,)
        elif var_name == 'v_rad':
            if dx[0] == 100:
                axis[0].set_ylim(-.75, 7)
                axis[1].set_ylim(-.5, 4)
                axis[2].set_ylim(-.5, .5)
                axis[3].set_ylim(-1.5, .75)
            else:
                axis[0].set_ylim(-.75,8)
                axis[1].set_ylim(-.5,5)
                axis[2].set_ylim(-.5,2.5)
                axis[3].set_ylim(-1.5,.75)
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

    global dt_fields
    dt_fields = nml['fields_io']['frequency']


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
