import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os
import scipy
from scipy import optimize

execfile('settings.py')


label_size = 10
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.labelsize'] = 15
# plt.rcParams['text.usetex'] = 'true'

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    case_name = 'ColdPoolDry_single_3D'
    path_root = '/nbi/ac/conv4/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/single_3D_noise/run5_PE_scaling_dx100m/'
    path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single'
    print('path figures: '  + path_out_figs)
    dTh, z_params, r_params = set_input_parameters(args, path_root, case_name)

    # reference case: dTh3_z1000_r1000
    rstar_ref = 1000
    zstar_ref = 1000
    id_ref = 'dTh3_z' + str(zstar_ref) + '_r' + str(rstar_ref)
    path_ref = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/' + id_ref
    dt_fields = 100
    cp_id = 2  # circle ID that is used for statistics

    id0 = 'dTh' + str(dTh) + '_z' + str(z_params[0]) + '_r' + str(r_params[0])
    fullpath_in = os.path.join(path_root, id0, 'tracer_k' + str(k0), 'output')
    n_tracers = get_number_tracers(fullpath_in)
    n_cps = get_number_cps(fullpath_in)
    print('number of CPs: ', n_cps)
    print('number of tracers per CP: ', n_tracers)
    print ''


    colorlist5 = ['maroon', 'indianred', 'orange', 'darkcyan', 'navy']

    # --------------------------------------
    ''' (a) read in data from tracer output (text-file)'''
    dist_av = np.zeros((n_params, nt))
    r_av = np.zeros((n_params, nt))  # absolute radius
    r_av_abs = np.zeros((n_params, nt))  # radius minus initial radius (r(t) - r(t=0))
    drdt_av = np.zeros((n_params, nt))
    U_rad_av = np.zeros((n_params, nt))

    print('--- reading in tracer data')
    for istar in range(n_params):
        zstar = z_params[0]
        rstar = r_params[istar]
        if rstar == 1000:
            id = id_ref
            fullpath_in = os.path.join(path_ref, 'tracer_k' + str(k0), 'output')
        else:
            id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
            fullpath_in = os.path.join(path_root, id, 'tracer_k' + str(k0), 'output')
        print('id', id)
        # read_in_txtfile(fullpath_in)
        for it, t0 in enumerate(times):
            print('---t0: ' + str(t0) + '---', it)
            dist_av[istar, it], U_rad_av[istar, it] = get_radius_vel(fullpath_in, it, cp_id, n_tracers, n_cps)
        r_av = dist_av * dx[0]
        r_av_abs[istar, :] = r_av[istar, :] - rstar
        r_av_abs[istar, :] = r_av[istar, :] - r_av[istar, 0]
        print('rstar: ', dx[0], rstar, r_av[istar, 0], r_av_abs[istar, 0])
        for it, t0 in enumerate(times[1:]):
            drdt_av[:, it] = 1. / dt_fields * (r_av[:, it] - r_av[:, it - 1])
        print ''
    print ''

    # --------------------------------------
    fig_name = 'R_U_scaling.png'

    tmin_r = 8
    tmin_u = 6
    t0_r = np.double(tmin_r * dt_fields)
    t0_u = np.double(tmin_u * dt_fields)

    R0 = r_av[0, tmin_r]
    #dR_log = np.log(r_av[:, tmin_r]) - np.log(R0)
    dR_log = np.log(r_av[:, tmin_r])
    #dR = R0 / r_av[:, tmin_r]
    dR = 1. / r_av[:, tmin_r]
    r_av_ens_mean = np.average(r_av[:, :], axis=0)
    #dR_ens_mean = R0 / r_av_ens_mean[tmin_r]
    dR_ens_mean = 1. / r_av_ens_mean[tmin_r]
    def fit_r(R0, m, times, t0):
        #return np.log(R0) + m * (np.log(times / t0))
        return m * (np.log(times / t0))

    U0 = drdt_av[0,tmin_r]
    #dU_log = np.log(drdt_av[:, tmin_u]) - np.log(U0)
    dU_log = np.log(drdt_av[:, tmin_r])
    #dU = U0 / drdt_av[:, tmin_u]
    dU = 1. / drdt_av[:, tmin_r]
    U_ens_mean = np.average(drdt_av[:, :], axis=0)
    #dU_ens_mean = U0 / U_ens_mean[tmin_u]
    dU_ens_mean = 1. / U_ens_mean[tmin_r]


    fig, axes = plt.subplots(2, 2, sharex='none', figsize=(18, 12))

    ax0 = axes[0,0]
    ax2 = axes[0,1]
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar])
        ax0.plot(times[1:-1], r_av[istar, 1:-1], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[1:-1] / t0_r), np.log(r_av[istar, 1:-1]) - dR_log[istar], '-', color='0.5')
    #ax2.plot(np.log(times[3:-1]/t0_r), np.log(r_av_ens_mean[3:-1] * dR_ens_mean), '-', color='k', linewidth=3, label=lbl)
    #ax2.plot(np.log(times[3:-1]/t0_r), np.log(r_av_ens_mean[3:-1] * dR_ens_mean), '-', color='k', linewidth=3,
    #         label=lbl)
    ax2.plot(np.log(times[3:-1]/t0_r), np.log(r_av_ens_mean[3:-1]) - np.log(dR_ens_mean), '-', color='k', linewidth=3,
             label=lbl)
    m = 0.6
    ax2.plot(np.log(times[3:-1]/t0_r), m * (np.log(times[3:-1] / t0_r)), '-r', label='m='+str(m))
    ax0.set_xlabel('t  [s]')
    ax0.set_ylabel('R  [m]')
    ax2.set_xlabel('log(t/t0)  [-]')
    ax2.set_ylabel('log(R/R0)  [-]')
    ax0.legend(loc='best')
    ax2.set_xlim(-1., 1.5)

    ax0 = axes[1, 0]
    ax2 = axes[1, 1]
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar])
        ax0.plot(times[1:-1], drdt_av[istar, 1:-1], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[3:-1] / t0_r), np.log(drdt_av[istar, 3:-1]) - dU_log[istar], '-', color='0.5')
    ax2.plot(np.log(times[3:-1] / t0_r), np.log(U_ens_mean[3:-1] * dU_ens_mean), '-', color='k', linewidth=3, label=lbl)
    times_ = np.append(times, np.arange(times[-1]+dt_fields, times[-1]+20*dt_fields, dt_fields))
    for m in np.arange(0.5,1,0.1):
        #ax2.plot(np.log(times_[7:] / t0_u), (np.log(U0)) - m*(np.log((times_[7:])/t0_u)-0.5), '-r', linewidth=1, label='m=-'+str(m))
        ax2.plot(np.log(times[3:-1] / t0_r), -m*(np.log((times[3:-1])/t0_r)), '-r', linewidth=1, label='m=-'+str(m))
    ax0.set_xlabel('t  [s]')
    ax0.set_ylabel('U  [m]')
    ax2.set_xlabel('log(t/t0)  [-]')
    ax2.set_ylabel('log(U/U0)  [-]')
    ax2.set_xlim(-1., 1.5)
    #ax2.set_ylim(-.7, 1.3)
    # ax0.plot([t0_u, t0_u], [0,6], 'k')

    plt.subplots_adjust(bottom=0.05, right=.95, left=0.06, top=0.95, wspace=0.25, hspace=0.2)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)


    return



# ----------------------------------------------------------------------

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
    cp_number = int(lines[-1].split()[3])

    f.close()

    return cp_number

# ----------------------------------------------------------------------

def set_input_parameters(args, path_root, case_name):
    print('--- set input parameters ---')

    global n_params
    dTh = 5
    z_params = [1000]  # run5, 6
    r_params = [500, 1000, 1100, 1600, 2300]
    n_params = len(r_params)
    print('dTh: ', dTh)
    print('z*: ', z_params)
    print('r*: ', r_params)
    print('n_params: ', n_params)

    id0 = 'dTh' + str(dTh) + '_z' + str(z_params[0]) + '_r' + str(r_params[0])
    nml = simplejson.loads(open(os.path.join(path_root, id0, case_name + '.in')).read())
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

    global times, nt, k0
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = np.int(100)
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = np.int(3600)
    times = np.arange(tmin, tmax + 100, 100)
    times.sort()

    nt = len(times)
    k0 = 0
    # krange = [0]
    # nk = len(krange)

    print('times', times)
    print('nt: '+str(nt))
    print('k0: ' + str(k0))
    print ''

    return dTh, z_params, r_params

# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()