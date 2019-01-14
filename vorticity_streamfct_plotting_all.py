import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import netCDF4 as nc
import argparse
import json as simplejson
import os

# compute potential temperature by integrating over anomaly
#   PE = \int dz g * (th_anomaly(z) - th_env(z)) * z

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35
# plt.rcParams['savefig.edgecolor'] = 'white'


def main():
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path_root")
    parser.add_argument("dTh", type=int)
    parser.add_argument("--zparams", nargs='+', type=int)
    parser.add_argument('--rparams', nargs='+', type=int)
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')

    nml, dTh, z_params, r_params, tmin, tmax = set_input_parameters(args)
    i0_center, j0_center = define_geometry(case_name, nml)

    ng = len(z_params)
    kmax = np.amax(z_params) + 2000. / dx[2]

    print ' '
    fig, axes = plt.subplots(1, 4, figsize=(20,4))
    for istar in range(ng):
        zstar = z_params[istar]
        rstar = r_params[istar]
        irstar = np.int(np.round(rstar / dx[0]))
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)

        path = os.path.join(path_root, id)
        path_in = os.path.join(path_root, id, 'figs_vorticity')
        path_fields = os.path.join(path_root, id, 'fields')
        path_out = os.path.join(path_root, 'figs_vorticity')
        if not os.path.exists(path_out):
            os.mkdir(path_out)

        filename = 'Stats_rotational.nc'
        rootgrp = nc.Dataset(os.path.join(path_in, filename), 'r+', format='NETCDF4')
        ts_grp = rootgrp.groups['timeseries']
        times = ts_grp.variables['time'][:]
        max = ts_grp.variables['vort_yz_max'][:]
        min = ts_grp.variables['vort_yz_min'][:]
        sum = ts_grp.variables['vort_yz_sum'][:]
        env = ts_grp.variables['vort_yz_env'][:]
        rootgrp.close()

        ax = axes[0]
        ax.plot(times, max, '-o', markeredgecolor='w', label=id)
        ax.set_title('max(vort_yz)')
        ax.set_xlabel('time  [s]')
        ax = axes[1]
        ax.plot(times, min, '-o', markeredgecolor='w', label=id)
        ax.set_title('min(vort_yz)')
        ax.set_xlabel('time  [s]')
        ax = axes[2]
        ax.plot(times, sum, '-o', markeredgecolor='w', label=id)
        ax.set_title('sum_ijk(vort_yz)')
        ax.set_xlabel('time  [s]')
        ax = axes[3]
        ax.plot(times, env, '-o', markeredgecolor='w', label=id)
        ax.set_title('environmental vort_yz')
        ax.set_xlabel('time  [s]')

    for ax in axes:
        ax.legend(loc='best')
    # plt.lines.set_markeredgewidth(0.0)
    plt.tight_layout
    plt.savefig(os.path.join(path_out, 'dTh_' + str(dTh) + '_vort_yz_domain.png'))
    plt.close(fig)


    print ' '

    fig, axes = plt.subplots(1, 4, figsize=(20, 4))
    for istar in range(ng):
        zstar = z_params[istar]
        rstar = r_params[istar]
        irstar = np.int(np.round(rstar / dx[0]))
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)

        path = os.path.join(path_root, id)
        path_in = os.path.join(path_root, id, 'figs_vorticity')
        path_fields = os.path.join(path_root, id, 'fields')
        path_out = os.path.join(path_root, 'figs_vorticity')
        if not os.path.exists(path_out):
            os.mkdir(path_out)

        filename = 'Stats_rotational.nc'
        rootgrp = nc.Dataset(os.path.join(path_in, filename), 'r+', format='NETCDF4')
        ts_grp = rootgrp.groups['timeseries']
        times = ts_grp.variables['time'][:]
        max = ts_grp.variables['vort_yz_max'][:]
        min = ts_grp.variables['vort_yz_min'][:]
        sum = ts_grp.variables['vort_yz_sum'][:]
        env = ts_grp.variables['vort_yz_env'][:]
        rootgrp.close()

        ax = axes[0]
        ax.plot(times, max, '-o', markeredgecolor='w', label=id)
        ax.set_title('max(vort_yz)')
        ax.set_xlabel('time  [s]')
        ax = axes[1]
        ax.plot(times, min, '-o', markeredgecolor='w', label=id)
        ax.set_title('min(vort_yz)')
        ax.set_xlabel('time  [s]')
        ax = axes[2]
        ax.plot(times, sum, '-o', markeredgecolor='w', label=id)
        ax.set_title('sum_ijk(vort_yz)')
        ax.set_xlabel('time  [s]')
        ax = axes[3]
        ax.plot(times, env, '-o', markeredgecolor='w', label=id)
        ax.set_title('environmental vort_yz')
        ax.set_xlabel('time  [s]')

    for ax in axes:
        ax.legend(loc='best')
    # plt.lines.set_markeredgewidth(0.0)
    plt.tight_layout
    plt.savefig(os.path.join(path_out, 'dTh_' + str(dTh) + '_vort_yz_domain.png'))
    plt.close(fig)

    # --------------------------------------
    print ' '
    print('--- plotting r=1km ---')
    dTh_params = [2, 3, 4]
    z_params = [1225, 1000, 870]
    r_params = z_params
    n_params = len(dTh_params)
    fig, axes = plt.subplots(1, 4, figsize=(20, 4))
    for istar in range(n_params):
        dTh = dTh_params[istar]
        zstar = z_params[istar]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        print('id', id)
        path = os.path.join(path_root, id)
        path_in = os.path.join(path_root, id, 'figs_vorticity')
        path_out = os.path.join(path_root, 'figs_vorticity')
        if not os.path.exists(path_out):
            os.mkdir(path_out)
        print''
        rootgrp = nc.Dataset(os.path.join(path_in, filename), 'r+', format='NETCDF4')
        ts_grp = rootgrp.groups['timeseries']
        times = ts_grp.variables['time'][:]
        max = ts_grp.variables['vort_yz_max'][:]
        min = ts_grp.variables['vort_yz_min'][:]
        sum = ts_grp.variables['vort_yz_sum'][:]
        env = ts_grp.variables['vort_yz_env'][:]
        rootgrp.close()

        ax = axes[0]
        ax.plot(times, max, '-o', markeredgecolor='w', label=id)
        ax.set_title('max(vort_yz)')
        ax.set_xlabel('time  [s]')
        ax = axes[1]
        ax.plot(times, min, '-o', markeredgecolor='w', label=id)
        ax.set_title('min(vort_yz)')
        ax.set_xlabel('time  [s]')
        ax = axes[2]
        ax.plot(times, sum, '-o', markeredgecolor='w', label=id)
        ax.set_title('sum_ijk(vort_yz)')
        ax.set_xlabel('time  [s]')
        ax = axes[3]
        ax.plot(times, env, '-o', markeredgecolor='w', label=id)
        ax.set_title('environmental vort_yz')
        ax.set_xlabel('time  [s]')


    for ax in axes:
        ax.legend(loc='best')
        # plt.lines.set_markeredgewidth(0.0)
    plt.tight_layout
    figname = 'vort_yz_domain_r1km.png'
    plt.savefig(os.path.join(path_out, figname))

    plt.close(fig)

    return


# ----------------------------------------------------------------------

def set_input_parameters(args):
    print('--- set input parameters ---')
    global case_name
    global path_root
    global timerange, nt

    path_root = args.path_root
    path_out_figs = os.path.join(path_root, 'figs_CP_height')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)

    dTh = args.dTh
    z_params = args.zparams
    r_params = args.rparams
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


    ''' determine file range '''
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = tmin
    timerange = np.arange(tmin, tmax + 100, 100)
    nt = len(timerange)
    # times = [np.int(name[:-3]) for name in files]
    # times.sort()
    print('timerange', timerange)

    return nml, dTh, z_params, r_params, tmin, tmax




def define_geometry(case_name, nml):
    print('--- define geometry ---')
    global x_half, y_half, z_half
    global ic_arr, jc_arr
    # global rstar, irstar, zstar, kstar

    x_half = np.empty((nx), dtype=np.double, order='c')
    y_half = np.empty((ny), dtype=np.double, order='c')
    z_half = np.empty((nz), dtype=np.double, order='c')
    count = 0
    for i in xrange(nx):
        x_half[count] = (i + 0.5) * dx[0]
        count += 1
    count = 0
    for j in xrange(ny):
        y_half[count] = (j + 0.5) * dx[1]
        count += 1
    count = 0
    for i in xrange(nz):
        z_half[count] = (i + 0.5) * dx[2]
        count += 1

    # set coordinates for plots
    if case_name == 'ColdPoolDry_single_3D':
        try:
            ic = nml['init']['ic']
            jc = nml['init']['jc']
            # print('(ic,jc) from nml')
        except:
            ic = np.int(nx/2)
            jc = np.int(ny/2)
            # print('(ic,jc) NOT from nml')
        ic_arr = np.zeros(1)
        jc_arr = np.zeros(1)
        ic_arr[0] = ic
        jc_arr[0] = jc
    elif case_name == 'ColdPoolDry_double_2D':
        # try:
        #     rstar = nml['init']['r']
        # except:
        #     rstar = 5000.0  # half of the width of initial cold-pools [m]
        # irstar = np.int(np.round(rstar / dx[0]))
        # zstar = nml['init']['h']
        # kstar = np.int(np.round(zstar / dx[2]))
        # isep = 4 * irstar
        ic1 = np.int(nx / 3)  # np.int(Gr.dims.ng[0] / 3)
        ic2 = ic1 + isep
        jc1 = np.int(ny / 2)
        jc2 = jc1
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_double_3D':
        # try:
        #     rstar = nml['init']['r']
        # except:
        #     rstar = 5000.0  # half of the width of initial cold-pools [m]
        # irstar = np.int(np.round(rstar / dx[0]))
        # zstar = nml['init']['h']
        # kstar = np.int(np.round(zstar / dx[2]))
        # isep = 4 * irstar
        # jsep = 0
        ic1 = np.int(np.round((nx + 2 * gw) / 3)) - gw
        jc1 = np.int(np.round((ny + 2 * gw) / 2)) - gw
        ic2 = ic1 + isep
        jc2 = jc1 + jsep
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_triple_3D':
        # try:
        #     rstar = nml['init']['r']
        # except:
        #     rstar = 5000.0  # half of the width of initial cold-pools [m]
        # irstar = np.int(np.round(rstar / dx[0]))
        # zstar = nml['init']['h']
        # kstar = np.int(np.round(zstar / dx[2]))
        d = np.int(np.round(ny / 2))
        dhalf = np.int(np.round(ny / 4))
        a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        ic1 = np.int(np.round(a / 2))  # + gw
        ic2 = ic1
        ic3 = ic1 + np.int(np.round(a))
        jc1 = np.int(np.round(d / 2))  # + gw
        jc2 = jc1 + d
        jc3 = jc1 + np.int(np.round(d / 2))
        ic_arr = [ic1, ic2, ic3]
        jc_arr = [jc1, jc2, jc3]

        isep = dhalf

    ''' plotting parameters '''
    if case_name == 'ColdPoolDry_single_3D':
        i0_coll = 10
        j0_coll = 10
        i0_center = ic_arr[0]
        j0_center = jc_arr[0]
    elif case_name == 'ColdPoolDry_double_3D':
        i0_coll = 0.5 * (ic_arr[0] + ic_arr[1])
        i0_center = ic_arr[0]
        j0_coll = jc_arr[0]
        j0_center = jc_arr[0]
        # domain boundaries for plotting
    elif case_name == 'ColdPoolDry_triple_3D':
        i0_coll = ic_arr[2]
        i0_center = ic_arr[0]
        j0_coll = jc_arr[2]
        j0_center = jc_arr[0]
        # domain boundaries for plotting

    return i0_center, j0_center

# ----------------------------------


# ----------------------------------

if __name__ == '__main__':
    main()
