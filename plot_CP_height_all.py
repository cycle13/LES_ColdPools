import numpy as np
import matplotlib.pyplot as plt
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
    print('')

    print('start')
    ng = len(z_params)
    kmax = np.amax(z_params) + 2000./dx[2]
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(1, 4, sharex='all', figsize=(18, 5))
    for istar in range(ng):
        zstar = z_params[istar]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        print('id', id)
        filename = 'CP_height_' + id + '_sth' + str(s_crit) + '.nc'
        fullpath_in = os.path.join(path_root, id, filename)
        print(fullpath_in)
        rootgrp = nc.Dataset(fullpath_in, 'r')
        time = rootgrp.groups['timeseries'].variables['time'][:]
        CP_height_max = rootgrp.groups['timeseries'].variables['CP_height'][:]
        CP_height_grad = rootgrp.groups['timeseries'].variables['CP_height_gradient'][:]
        w_max = rootgrp.groups['fields_2D'].variables['w_max'][:,:,:]
        w_max_height = rootgrp.groups['fields_2D'].variables['w_max_height'][:,:,:]
        rootgrp.close()
        print('max', CP_height_max)
        print('max', time)
        w_max_height_max = np.amax(np.amax(w_max_height, axis=2), axis=1)
        w_max_max = np.amax(np.amax(w_max, axis=2), axis=1)
        ax0.plot(time, CP_height_grad, '-o', label=id)
        ax1.plot(time, CP_height_max, '-o', label=id)
        ax2.plot(time, w_max_height_max, '-o', label=id)
        ax3.plot(time, w_max_max, '-o', label=id)

    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax0.set_title('max(CP_height_grad)')
    ax1.set_title('max(CP_height)')
    ax2.set_title('max(height of w_max[i,j])')
    ax3.set_title('max(w_max[ij])')
    fig.suptitle('dTh=' + str(dTh))
    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, 'CP_height_dTh' + str(dTh) + '.png'))
    plt.close(fig)







    return


# ----------------------------------------------------------------------
def plot_contourf_test_yz(s, smin, smax, CP_top, kmax, t0):
    ic1 = ic_arr[0]
    levels = np.arange(smin, smax, 0.1)

    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    ax1 = axes[0]
    a = ax1.contourf(s[:, :, 1].T, levels=levels)
    ax1.plot([ic1, ic1], [0, ny], 'k-')
    plt.colorbar(a, ax=ax1)
    ax2 = axes[1]
    ax2.contourf(s[ic1, :, :kmax].T, levels=levels)
    ax2.plot(CP_top[ic1, :], 'w-')
    fig.suptitle('t=' + str(t0) + 's')
    fig.tight_layout()
    fig.savefig(os.path.join(path_out, 'CP_height_yz_t' + str(t0) + '.png'))
    plt.close(fig)
    return


def plot_contourf_xy(CP_top, w_max, t0):
    wmax = np.ceil(np.maximum(np.amax(w_max[0,:,:]), np.abs(np.amin(w_max[0,:,:]))))
    if wmax > 0.0:
        levels = np.arange(-wmax, wmax+1, 1)
    else:
        levels = np.arange(-1,2,1)
    # print('')
    # print levels

    fig, axes = plt.subplots(1,3, figsize=(16,5), sharey=True)
    ax1 = axes[0]
    a = ax1.contourf(dx[2]*CP_top.T)
    plt.colorbar(a, ax=ax1)
    ax1.set_title('CP height  (max='+str(dx[2]*np.amax(CP_top))+'m)')
    ax1.set_xlabel('x  (dx=' + str(dx[0]) + 'm)')
    ax1.set_ylabel('y  (dy=' + str(dx[1]) + 'm)')

    ax2 = axes[1]
    b = ax2.contourf(w_max[0,:,:].T, levels=levels, cmap=cm_bwr)
    plt.colorbar(b, ax=ax2)
    ax2.set_title('max(w), (max='+str(np.round(np.amax(w_max[0,:,:]),2))+'m/s)')
    ax2.set_xlabel('x  (dx='+str(dx[0])+')')

    ax3 = axes[2]
    b = ax3.contourf(dx[2]*w_max[1, :, :].T)
    plt.colorbar(b, ax=ax3)
    ax3.set_title('height of max(w), (max='+str(dx[2]*np.int(np.amax(w_max[1,:,:])))+'m)')
    ax3.set_xlabel('x  (dx='+str(dx[0])+')')
    fig.suptitle('t='+str(t0)+'s')
    fig.tight_layout()
    fig.savefig(os.path.join(path_out, 'CP_height_t'+str(t0)+'.png'))
    plt.close(fig)

    return



def plot_x_crosssections(s0, CP_top, w_max, jp1, jp2):
    cm = plt.cm.get_cmap('coolwarm')

    f1, axes = plt.subplots(3, 1, figsize=(6, 12))
    ax1 = axes[0]
    ax1.set_title('entropy')
    ax1.set_xlabel('x  (dx=' + str(dx[0]) + 'm)')
    ax1.set_ylabel('y  (dy=' + str(dx[1]) + 'm)')
    ax2 = axes[1]
    ax2.set_title('x-crossection through collision point')
    ax2.set_ylabel('CP height [m]')
    ax3 = axes[2]
    ax3.set_title('x-crossection through center of coldpool #1')
    ax3.set_xlabel('y  (dy='+str(dx[1])+'m)')
    ax3.set_ylabel('CP height [m]')
    ax1.imshow(s0[:,:,1].T)
    ax1.set_xlim(0,nx)
    ax1.set_ylim(0,ny)
    ax1.plot([0,nx],[jp1,jp1], 'k', linewidth=2, label='collision center')
    ax1.plot([0,nx],[jp2,jp2], 'k--', linewidth=2, label='center CP #1')
    ax1.legend()
    for it,t0 in enumerate(times):
        count_color = np.double(it) / len(times)
        ax2.plot(CP_top[it, :, jp1], color=cm(count_color), label='t=' + str(t0))
        ax3.plot(CP_top[it, :, jp2], color=cm(count_color), label='t=' + str(t0))
    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
               fancybox=True, shadow=True, ncol=5)
    plt.suptitle('CP height', fontsize=21)
    plt.savefig(os.path.join(path_out_figs, 'CP_height_xz_plane.png'))
    plt.close()
    return



def plot_y_crosssections(s0, CP_top, w_max, ip1, ip2):
    cm = plt.cm.get_cmap('coolwarm')

    f1, axes = plt.subplots(3, 1, figsize=(6, 12))
    ax1 = axes[0]
    ax1.set_title('entropy')
    ax1.set_xlabel('x  (dx=' + str(dx[0]) + 'm)')
    ax1.set_ylabel('y  (dy=' + str(dx[1]) + 'm)')
    ax2 = axes[1]
    ax2.set_title('y-crossection through 2 CP-collision point')
    ax2.set_ylabel('CP height [m]')
    ax3 = axes[2]
    ax3.set_title('y-crossection through center of coldpool #3')
    ax3.set_xlabel('y  (dy='+str(dx[1])+'m)')
    ax3.set_ylabel('CP height [m]')
    ax1.imshow(s0[:,:,1].T)
    ax1.set_xlim(0,nx)
    ax1.set_ylim(0,ny)
    ax1.plot([ip1,ip1],[0,ny], 'k', linewidth=2, label='2-CP collision')
    ax1.plot([ip2,ip2],[0,ny], 'k--', linewidth=2, label='center of CP #3')
    ax1.legend()
    for it,t0 in enumerate(times):
        count_color = np.double(it) / len(times)
        ax2.plot(CP_top[it, ip1, :], color=cm(count_color), label='t=' + str(t0))
        ax3.plot(CP_top[it, ip2, :], color=cm(count_color), label='t=' + str(t0))
    # # ax3.legend( )
    # # Put a legend below current axis
    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
               fancybox=True, shadow=True, ncol=5)
    plt.suptitle('CP height', fontsize=21)
    plt.savefig(os.path.join(path_out, 'CP_height_yz_plane.png'))
    plt.close()
    return




def plot_x_crosssections_CPtop_w(s0, CP_top, w_max, jp1, jp2):
    cm = plt.cm.get_cmap('coolwarm')

    f1, axes = plt.subplots(3, 1, figsize=(6, 12))
    ax1 = axes[0]
    ax1.set_title('entropy')
    ax1.set_xlabel('x  (dx='+str(dx[0])+'m)')
    ax1.set_ylabel('y  (dy='+str(dx[1])+'m)')
    ax2 = axes[1]
    ax2.set_title('CP height')
    ax2.set_ylabel('CP height [m]')
    ax3 = axes[2]
    ax3.set_title('max w')
    ax3.set_xlabel('y  (dy='+str(dx[1])+'m)')
    ax3.set_ylabel('max(w) [m/s]')
    ax1.imshow(s0[:,:,1].T)
    ax1.set_xlim(0,nx)
    ax1.set_ylim(0,ny)
    ax1.plot([0,nx],[jp1,jp1], 'k', linewidth=2, label='jp1')
    for it,t0 in enumerate(times):
        count_color = np.double(it) / len(times)
        ax2.plot(CP_top[it, :, jp1], color=cm(count_color), label='t=' + str(t0))
        ax3.plot(w_max[0, it, :, jp1], color=cm(count_color), label='t=' + str(t0))
    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
               fancybox=True, shadow=True, ncol=5)
    plt.suptitle('x-crossection through collision point', fontsize=21)
    plt.savefig(os.path.join(path_out, 'CP_height_wmax_xz_plane.png'))
    plt.close()
    return



def plot_y_crosssections_CPtop_w(s0, CP_top, w_max, ip1, ip2):
    cm = plt.cm.get_cmap('coolwarm')

    f1, axes = plt.subplots(3, 1, figsize=(6, 12))
    ax1 = axes[0]
    ax1.set_title('entropy')
    ax1.set_xlabel('x  (dx='+str(dx[0])+'m)')
    ax1.set_ylabel('y  (dy='+str(dx[1])+'m)')
    ax2 = axes[1]
    ax2.set_title('CP height')
    ax2.set_ylabel('CP height [m]')
    ax3 = axes[2]
    ax3.set_title('max w')
    ax3.set_xlabel('y  (dy='+str(dx[1])+'m)')
    ax3.set_ylabel('max(w) [m/s]')
    ax1.imshow(s0[:,:,1].T)
    ax1.set_xlim(0,nx)
    ax1.set_ylim(0,ny)
    ax1.plot([ip1,ip1],[0,ny], 'k', linewidth=2, label='ip1')
    for it,t0 in enumerate(times):
        count_color = np.double(it) / len(times)
        ax2.plot(CP_top[it, ip1, :], color=cm(count_color), label='t=' + str(t0))
        ax3.plot(w_max[0, it, ip1, :], color=cm(count_color), label='t=' + str(t0))
    # # ax3.legend( )
    # # Put a legend below current axis
    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
               fancybox=True, shadow=True, ncol=5)
    plt.suptitle('y-crossection through 2 CP-collision point', fontsize=21)
    plt.savefig(os.path.join(path_out, 'CP_height_wmax_yz_plane.png'))
    plt.close()
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
            print('(ic,jc) from nml')
        except:
            ic = np.int(nx/2)
            jc = np.int(ny/2)
            print('(ic,jc) NOT from nml')
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
        # zstar = nml['init']['h']
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


# ----------------------------------
def create_output_file(id, nx, ny, times, path_out):
    # output for each CP:
    # - min, max (timeseries)
    # - CP height (field; max=timeseries)
    # - (ok) CP rim (field)
    nt = len(times)

    filename = 'CP_height_' + id + '.nc'
    rootgrp = nc.Dataset(os.path.join(path_out, filename), 'w', format='NETCDF4')

    ts_grp = rootgrp.createGroup('timeseries')
    ts_grp.createDimension('nt', nt)
    var = ts_grp.createVariable('time', 'f8', ('nt'))
    var.units = "s"
    var[:] = times

    var = ts_grp.createVariable('CP_height', 'f8', ('nt'))
    var.units = "m"
    var = ts_grp.createVariable('CP_height_gradient', 'f8', ('nt'))
    var.units = "m"

    field_grp = rootgrp.createGroup('fields_2D')
    field_grp.createDimension('nt', nt)
    field_grp.createDimension('nx', nx)
    field_grp.createDimension('ny', ny)
    var = field_grp.createVariable('w_max', 'f8', ('nt', 'nx', 'ny'))
    var.units = "m/s"
    var = field_grp.createVariable('w_max_height', 'f8', ('nt', 'nx', 'ny'))
    var.units = "m"

    rootgrp.close()
    return


def dump_output_file(CP_top_max, w_max, it, id, path_out):
    filename = 'CP_height_' + id + '.nc'
    print('dumping', os.path.join(path_out, filename))
    rootgrp = nc.Dataset(os.path.join(path_out, filename), 'r+', format='NETCDF4')

    ts_grp = rootgrp.groups['timeseries']
    var = ts_grp.variables['CP_height']
    var[it] = CP_top_max[0, it]*dx[2]
    var = ts_grp.variables['CP_height_gradient']
    var[it] = CP_top_max[1, it]*dx[2]

    field_grp = rootgrp.groups['fields_2D']
    var = field_grp.variables['w_max']
    var[it,:,:] = w_max[0,:,:]
    var = field_grp.variables['w_max_height']
    var[it,:,:] = w_max[1,:,:]

    rootgrp.close()
    return
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


if __name__ == '__main__':
    main()