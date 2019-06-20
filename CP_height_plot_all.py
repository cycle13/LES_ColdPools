import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as patch
import netCDF4 as nc
import argparse
import json as simplejson
import os
import time


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35

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
    parser.add_argument("--s_crit")
    args = parser.parse_args()

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')
    cm_fall = plt.cm.get_cmap('winter')
    cm_summer = plt.cm.get_cmap('spring')

    nml, dTh, z_params, r_params = set_input_parameters(args)
    i0_center, j0_center, xmin_plt, xmax_plt, ymin_plt, ymax_plt = define_geometry(case_name, nml)

    ''' threshold for entropy '''
    if args.s_crit:
        s_crit = args.scrit
    else:
        s_crit = 5e-1
    print('threshold for ds=s-s_bg: ' + str(s_crit) + 'J/K')
    print('')

    # plot_CP_height_dTh(z_params, r_params, dTh, s_crit)
    plot_CP_height_rad_av_dTh(z_params, r_params, dTh, s_crit)

    # plot_w_max_height_dTh(z_params, r_params, dTh, s_crit)



    return


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

def plot_CP_height_dTh(z_params, r_params, dTh, s_crit):
    print('start')
    ng = len(z_params)
    kmax = np.amax(z_params) + 2000./dx[2]
    fig, axis = plt.subplots(1, 4, sharex='all', figsize=(18, 5))
    ax0 = axis[0]
    ax1 = axis[1]
    ax2 = axis[2]
    ax3 = axis[3]
    for istar in range(ng):
        zstar = z_params[istar]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        print('id', id)
        filename = 'CP_height_' + id + '_sth' + str(s_crit) + '.nc'
        fullpath_in = os.path.join(path_root, id, 'data_analysis', filename)
        print(fullpath_in)
        rootgrp = nc.Dataset(fullpath_in, 'r')
        time_CP_height = rootgrp.groups['timeseries'].variables['time'][:]
        CP_height_max = rootgrp.groups['timeseries'].variables['CP_height_max'][:]
        CP_height_grad = rootgrp.groups['timeseries'].variables['CP_height_gradient_max'][:]
        w_max = rootgrp.groups['fields_2D'].variables['w_max_2d'][:,:,:]
        w_max_height = rootgrp.groups['fields_2D'].variables['w_max_height_2d'][:,:,:]
        rootgrp.close()
        print('max', CP_height_max)
        print('max', time_CP_height)
        w_max_height_max = np.amax(np.amax(w_max_height, axis=2), axis=1)
        w_max_max = np.amax(np.amax(w_max, axis=2), axis=1)
        ax0.plot(time_CP_height, CP_height_grad, '-o', label=id)
        ax1.plot(time_CP_height, CP_height_max, '-o', label=id)
        ax2.plot(time_CP_height, w_max_height_max, '-o', label=id)
        ax3.plot(time_CP_height, w_max_max, '-o', label=id)
    for i in range(4):
        axis[i].legend(loc='best')
        axis[i].set_xlabel('time  [s]')
    ax0.set_title('max(CP_height_grad)')
    ax1.set_title('max(CP_height)')
    ax2.set_title('max(height of w_max[i,j])')
    ax3.set_title('max(w_max[ij])')
    fig.suptitle('dTh=' + str(dTh))
    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, 'CP_height_dTh' + str(dTh) + '.png'))
    plt.close(fig)

    return


def plot_CP_height_rad_av_dTh(z_params, r_params, dTh, s_crit):
    ng = len(z_params)
    fig, axis = plt.subplots(1, 3, sharex='none', figsize=(18, 5))
    ax0 = axis[0]
    ax1 = axis[1]
    ax2 = axis[2]
    # ax3 = axis[3]
    colorlist = ['black', 'firebrick', 'navy', 'darkgreen', 'darkcyan']

    for istar in range(ng):
        zstar = z_params[istar]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        print('id', id)
        filename = 'CP_height_' + id + '_sth' + str(s_crit) + '.nc'
        fullpath_in = os.path.join(path_root, id, 'data_analysis', filename)
        rootgrp = nc.Dataset(fullpath_in, 'r')
        time_CP_height = rootgrp.groups['timeseries'].variables['time'][:]
        CP_height_max = rootgrp.groups['timeseries'].variables['CP_height_max'][:]
        CP_height_radav = rootgrp.groups['stats'].variables['CP_height_rad'][:,:]
        r_range = rootgrp.groups['stats'].variables['r'][:]
        rootgrp.close()

        ax0.set_title('max_ijk(CP_height)')
        ax0.plot(time_CP_height, CP_height_max, '-o', color=colorlist[istar], label=id)
        ax0.set_xlabel('time')
        ax0.set_ylabel('max(CP_height)')
        # ax0.set_xlim(tmin, tmax)
        ax0.legend()

        ax1.set_title('max_r(CP_height(r))')
        ax1.set_xlabel('time')
        ax1.set_ylabel('max(CP_height(r))')
        ax1.plot(time_CP_height, np.amax(CP_height_radav, axis=1), '-o', color=colorlist[istar], label=id)
        ax1.legend()


        ax2.set_title('CP_height(r,t)')
        ax2.set_xlabel('radius')
        ax2.set_ylabel('CP_height')


        ta = np.where(time_CP_height == tmin)[0][0]
        tb = np.where(time_CP_height == tmax)[0][0]
        rmax_plot = 9e3
        irmax = np.where(r_range == rmax_plot)[0][0]
        count_color = np.double(istar)/ng
        if len(times) >= 2:
            for it, t0 in enumerate(time_CP_height[1::2]):
                if it >= ta and it <= tb:
                    # count_color = 2 * np.double(it) / len(times)
                    ax2.plot(r_range[:irmax], CP_height_radav[2*it+1, :irmax],
                             color=colorlist[istar])#color=cm.jet(count_color))
        else:
            for it, t0 in enumerate(time_CP_height):
                if it >= ta and it <= tb:
                    # count_color = np.double(it) / len(times)
                    # ax2.plot([0, r_range[irmax]], [0.,0.], 'k')
                    ax2.plot(r_range[:irmax], CP_height_radav[it, :irmax],
                             color=colorlist[istar])  #color=cm.jet(count_color))
    #  ax2.legend()
    ax2.legend(loc='upper center', bbox_to_anchor=(1.2, 1.),
                    fancybox=True, shadow=True, ncol=1, fontsize=10)
    fig.suptitle('dTh=' + str(dTh), fontsize=18)
    fig.tight_layout()
    plt.subplots_adjust(bottom=0.12, right=.9, left=0.04, top=0.8, wspace=0.15)
    fig.savefig(os.path.join(path_out_figs, 'CP_height_radav_dTh' + str(dTh) + '.png'))
    plt.close(fig)
    return


def plot_w_max_height_dTh(z_params, r_params, dTh, s_crit):
    ng = len(z_params)
    fig, axis = plt.subplots(2, 2, sharex='all', figsize=(10, 10))
    ax0 = axis[0, 0]
    ax1 = axis[0, 1]
    ax2 = axis[1, 0]
    ax3 = axis[1, 1]
    for istar in range(ng):
        zstar = z_params[istar]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        print('id', id)
        filename = 'CP_height_' + id + '_sth' + str(s_crit) + '.nc'
        fullpath_in = os.path.join(path_root, id, 'data_analysis', filename)
        print(fullpath_in)
        rootgrp = nc.Dataset(fullpath_in, 'r')
        time = rootgrp.groups['timeseries'].variables['time'][:]
        w_max = rootgrp.groups['fields_2D'].variables['w_max_2d'][:, :, :]
        w_max_height = rootgrp.groups['fields_2D'].variables['w_max_height_2d'][:, :, :]
        rootgrp.close()
        print('max', time)
        w_max_height_max = np.amax(np.amax(w_max_height, axis=2), axis=1)
        w_max_height_av = np.mean(np.mean(w_max_height, axis=2), axis=1)
        w_max_max = np.amax(np.amax(w_max, axis=2), axis=1)
        w_max_av = np.mean(np.mean(w_max, axis=1), axis=1)
        ax0.plot(time, w_max_max, '-o', label=id)
        ax1.plot(time, w_max_av, '-o', label=id)
        ax2.plot(time, w_max_height_max, '-o', label=id)
        ax3.plot(time, w_max_height_av, '-o', label=id)

    for i in range(4):
        axis[0, np.mod(i, 2)].legend(loc='best')
        axis[1, np.mod(i, 2)].legend(loc='best')
        axis[1, np.mod(i, 2)].set_xlabel('time  [s]')
    axis[0, 0].set_ylabel('w  [m/s]')
    axis[1, 0].set_ylabel('z  [m]')
    ax0.set_title('max(w_max[ij])')
    ax1.set_title('mean(w_max[ij])')
    ax2.set_title('max(height of w_max[i,j])')
    ax3.set_title('mean(height of w_max[i,j])')
    fig.suptitle('dTh=' + str(dTh))
    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, 'w_max_height_dTh' + str(dTh) + '.png'))
    plt.close(fig)
    return


# ----------------------------------------------------------------------
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
    print('')
    print('path figs out: ')
    print('   ' + path_out_figs)
    print('')

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
    global tmin, tmax
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

    return nml, dTh, z_params, r_params




def define_geometry(case_name, nml):
    print('--- define geometry ---')
    global x_half, y_half, z_half
    global ic_arr, jc_arr
    global rstar, irstar, zstar, kstar

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
        rstar = nml['init']['r']
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
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
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        isep = 4 * irstar
        ic1 = np.int(nx / 3)  # np.int(Gr.dims.ng[0] / 3)
        ic2 = ic1 + isep
        jc1 = np.int(ny / 2)
        jc2 = jc1
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_double_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        kstar = np.int(np.round(zstar / dx[2]))
        isep = 4 * irstar
        jsep = 0
        ic1 = np.int(np.round((nx + 2 * gw) / 3)) - gw
        jc1 = np.int(np.round((ny + 2 * gw) / 2)) - gw
        ic2 = ic1 + isep
        jc2 = jc1 + jsep
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_triple_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        kstar = np.int(np.round(zstar / dx[2]))
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
        xmin_plt = 100
        xmax_plt = nx - xmin_plt
        ymin_plt = xmin_plt
        ymax_plt = xmax_plt
    elif case_name == 'ColdPoolDry_double_3D':
        i0_coll = 0.5 * (ic_arr[0] + ic_arr[1])
        i0_center = ic_arr[0]
        j0_coll = jc_arr[0]
        j0_center = jc_arr[0]
        # domain boundaries for plotting
        xmin_plt = 30
        xmax_plt = 230
        ymin_plt = xmin_plt
        ymax_plt = xmax_plt
    elif case_name == 'ColdPoolDry_triple_3D':
        i0_coll = ic_arr[2]
        i0_center = ic_arr[0]
        j0_coll = jc_arr[2]
        j0_center = jc_arr[0]
        # domain boundaries for plotting
        xmin_plt = 0
        xmax_plt = nx
        ymin_plt = xmin_plt
        ymax_plt = xmax_plt

    return i0_center, j0_center, xmin_plt, xmax_plt, ymin_plt, ymax_plt

# _____________________________________________________________________
# ----------------------------------
# ----------------------------------
def read_in_netcdf_fields(variable_name, fullpath_in):
    print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    # shape = var.shape
    # data = np.ndarray(shape = var.shape)
    data = var[:,:,:]
    rootgrp.close()
    return data

# _____________________________________________________________________
if __name__ == '__main__':
    main()