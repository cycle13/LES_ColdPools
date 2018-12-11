import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 12

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path_root")
    parser.add_argument("dTh", type=int)
    parser.add_argument("--zparams", nargs='+', type=int)
    parser.add_argument('--rparams', nargs='+', type=int)
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--perc")
    args = parser.parse_args()

    # percentile for threshold
    if args.perc:
        perc = args.perc
    else:
        perc = 98  # tested for triple 3D, dTh=10K, t=100-400s

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')
    cm_fall = plt.cm.get_cmap('winter')
    cm_summer = plt.cm.get_cmap('spring')



    dTh, z_params, r_params = set_input_parameters(args)

    # test file
    id0 = 'dTh' + str(dTh) + '_z' + str(z_params[0]) + '_r' + str(r_params[0])
    filename = 'rimstats_perc' + str(perc) + 'th.nc'
    fullpath_in = os.path.join(path_root, id0, 'fields_CP_rim', filename)
    print('fullpath:', fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    ts_grp = rootgrp.groups['timeseries']
    r_av = ts_grp.variables['r_av'][:, :]
    [nt, nz] = r_av.shape
    ts_grp = rootgrp.groups['time']
    times = ts_grp.variables['time']
    z_grp = rootgrp.groups['profiles']
    krange = z_grp.variables['krange'][:]
    k_dumped = z_grp.variables['k_dumped'][:]
    krange = krange*k_dumped
    rootgrp.close()


    n_params = len(z_params)
    r_av = np.zeros((n_params, nt, nz))
    U_av = np.zeros((n_params, nt, nz))
    dU_av = np.zeros((n_params, nt, nz))
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, sharex='all', figsize=(18, 5))
    for istar in range(n_params):
        zstar = z_params[istar]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        print('id', id)
        # filename = 'rimstats_perc' + str(perc) + 'th.nc'
        fullpath_in = os.path.join(path_root, id, 'fields_CP_rim', filename)
        rootgrp = nc.Dataset(fullpath_in, 'r')
        ts_grp = rootgrp.groups['timeseries']
        r_av[istar, :,:] = ts_grp.variables['r_av'][:,:]
        U_av[istar, :,:] = ts_grp.variables['U_av'][:,:]
        dU_av[istar, :,:] = ts_grp.variables['dU_av'][:,:]
        rootgrp.close()
        ax0.plot(times, r_av[istar, :, 0], '-o', label=id)
        ax1.plot(times, U_av[istar, :, 0], '-o', label=id)
        ax2.plot(times, dU_av[istar, :, 0], '-o', label=id)

    ax0.legend()
    ax1.legend()
    ax2.legend()
    ax0.set_title('r_av')
    ax1.set_title('U_av')
    ax2.set_title('dU_av')
    fig.suptitle('CP rim (dTh=' + str(dTh)+')')
    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, 'CP_rim_dTh' + str(dTh) + '.png'))
    plt.close(fig)


    return


# ----------------------------------------------------------------------

def set_input_parameters(args):
    print('--- set input parameters ---')
    global case_name
    global path_root, path_out_figs
    global times

    path_root = args.path_root
    path_out_figs = os.path.join(path_root, 'figs_CP_height')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)

    dTh = args.dTh
    z_params = args.zparams
    r_params = args.rparams
    print('dTh: ', dTh)
    print('z*: ', z_params)
    print('r*: ', r_params)

    case_name = args.casename
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

    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = np.int(100)
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = np.int(100)
    times = np.arange(tmin, tmax + 100, 100)
    # times = [np.int(name[:-3]) for name in files]
    times.sort()
    print('times', times)


    return dTh, z_params, r_params

# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()