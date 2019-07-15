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
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35

def main():

    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    case, times, files, krange = set_input_parameters(args)
    icoll, jcoll, ic_arr, jc_arr, x_half, y_half, z_half = define_geometry(case_name, files)

    ''' --- auxiliary arrays (since no Grid.pyx) ---'''
    # test file:
    var = read_in_netcdf_fields('s', os.path.join(path_fields, files[0]))
    plot_configuration(var, icoll, jcoll, ic_arr, jc_arr)

    #
    # ''' (a) compute horizontal convergence''''
    #

    ''' (b) plot profiles for w, s for all time-steps '''
    var_list = ['w', 's']
    # var_list = ['s']
    for k0 in krange:
        print('---- k0: '+str(k0) + '----')
        if case == 'double':
            # -- 3D --
            i0 = ic_arr[0] + isep
            plot_y_crosssections(var_list, k0, i0, ic_arr[0], jc_arr[0], ic_arr[1], jc_arr[1], files, path_out, path_fields)

            # -- 2D --
            # j0 = np.int(jc1 + (ic2 - ic1) / 4)
            # plot_xz_crosssections(var_list, j0, k0, ic_arr[0], jc_arr[0], ic_arr[1], jc_arr[1], files, path_out, path_fields)
        if case == 'triple':
            j0 = jcoll
            imin = (icoll - 200)
            imax = (icoll + 200)
            delta = r_int
            plot_xz_crosssections(var_list, j0, k0,
                                  ic_arr[0], jc_arr[0], ic_arr[2], jc_arr[2],
                                  icoll, jcoll,
                                  x_half, y_half, imin, imax, delta,
                                  files, path_out, path_fields)

    ''' (c) plot profiles at different levels for collision time '''
    krange_ = [1, 4, 10, 15]
    if case == 'single':
        jmin = jc_arr[0] - 80
        jmax = jc_arr[0] + 80
        plot_yz_crosssections_multilevel('s', ic_arr[0], krange_, jmin, jmax,  jc_arr, files, path_out, path_fields)
    if case == 'double':
        # -- 2D --
        plot_yz_crosssections_multilevel('w', icoll, krange, jmin, jmax,
                                         jc_arr, files, path_out, path_fields)
    elif case == 'triple':
        print 'triple'
        # -- 3D --
        j0 = jcoll
        imin = (icoll - 200)
        imax = (icoll + 200)
        # plot_yz_crosssections_multilevel('w', i0, krange_, jmin, jmax,   jc_arr, files, path_out, path_fields)
        plot_xz_crosssections_multilevel('w', j0, krange_, imin, imax, icoll, jcoll, ic_arr, jc_arr,
                                         files, path_out, path_fields)
    # print ''

    return


# ----------------------------------

def plot_configuration(var, icoll, jcoll, ic_arr, jc_arr):
    fig, (ax1, ax2)= plt.subplots(1,2)
    ax1.imshow(var[:, :, 0].T, origin='lower')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    plt.tight_layout()
    ax1.set_xlim(0,nx)
    ax1.set_ylim(0,ny)
    ax2.contourf(var[:, :, 0].T, origin='lower')
    ax2.set_xlabel('x')
    for ic in ic_arr:
        ax1.plot([ic, ic], [0, ny-1], 'w', linewidth=1)
        ax2.plot([ic, ic], [0, ny-1], 'w', linewidth=1)
    for jc in jc_arr:
        ax1.plot([0, nx], [jc,jc], 'w', linewidth=1)
        ax2.plot([0, nx], [jc,jc], 'w', linewidth=1)
    ax1.plot([icoll,icoll], [0,ny], 'w', linewidth=2)
    ax2.plot([icoll,icoll], [0,ny], 'w', linewidth=2)
    ax1.plot([0, nx], [jcoll,jcoll], 'w', linewidth=2)
    ax2.plot([0, nx], [jcoll,jcoll], 'w', linewidth=2)

    plt.tight_layout()
    ax2.set_xlim(icoll-200,icoll+200)
    ax2.set_ylim(jcoll-200,jcoll+200)
    ax2.set_aspect('equal')
    plt.savefig(os.path.join(path_out, 'initial.png'))
    plt.close()

    return

# ----------------------------------

def plot_xz_crosssections(var_list, j0, k0, ic1, jc1, ic2, jc2,
                          icoll, jcoll, x_arr, y_arr, imin, imax, delta, files, path_out, path_fields):
    # cm = plt.cm.get_cmap('bwr')
    cm = plt.cm.get_cmap('coolwarm')
    # cm = plt.cm.get_cmap('winter')
    # cm = plt.cm.get_cmap('viridis')
    # cm = plt.cm.get_cmap('viridis_r')

    print('xz-crosssection at j0='+str(j0)+', jc2='+str(jc2))

    for var_name in var_list:
        var_min1 = 9999.9       # min/max along axis through CP 1
        var_max1 = -9999.9
        var_min2 = 9999.9       # min/max along axis through CP 2 (collision)
        var_max2 = -9999.9
        var_min3 = 9999.9       # min/max along axis through CP 2 (2 CP collision)
        var_max3 = -9999.9

        fig, axis = plt.subplots(3, 1, figsize=(12, 12))
        ax1 = axis[0]
        ax2 = axis[1]
        ax3 = axis[2]
        s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
        ax1.imshow(s0[:, :, k0].T, origin="lower")
        ax1.plot([0, nx-1], [jc1, jc1], 'k', linewidth=2)
        ax1.plot([0, nx-1], [j0, j0], 'k', linewidth=2)
        ax1.plot([ic1, ic1], [0, ny-1], 'k--', linewidth=1)
        ax1.plot([ic2, ic2], [0, ny-1], 'k--', linewidth=1)
        ax1.plot([icoll, icoll], [0, ny-1], 'k:', linewidth=2)
        ax1.set_xlim([imin,imax])
        ax1.set_ylim([imin,imax])
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        count = 0
        for file in files:
            var = read_in_netcdf_fields(var_name, os.path.join(path_fields,file))
            if var_name == 's':
                var_min1 = np.minimum(var_min1,np.maximum(6000,np.amin(var[:,jc1,k0])))
                var_min2 = np.minimum(var_min2,np.maximum(6000,np.amin(var[ic1-delta:ic1+delta-3,jc2,k0])))
                var_min3 = np.minimum(var_min3,np.maximum(6000,np.amin(var[ic1+delta-3:,jc2,k0])))
            else:
                var_min1 = np.minimum(var_min1, np.amin(var[:, jc1, k0]))
                var_min2 = np.minimum(var_min2, np.amin(var[ic1-delta:ic1+delta-3, jc2, k0]))
                var_min3 = np.minimum(var_min3, np.amin(var[ic1+delta-3:, jc2, k0]))
            var_max1 = np.maximum(var_max1, np.amax(var[:,jc1,k0]))
            var_max2 = np.maximum(var_max2, np.amax(var[ic1-delta:ic1+delta-3,jc2,k0]))
            var_max3 = np.maximum(var_max3, np.amax(var[ic1+delta-3:,jc2,k0]))
            count_color = np.double(count) / len(files)
            ax2.plot(x_arr,var[:,jc1,k0], color=cm(count_color), label='t='+str(file[:-3]))
            ax3.plot(x_arr,var[:,j0,k0], color=cm(count_color), label='t='+str(file[:-3]))
            ax2.set_xlim(imin*dx[0], imax*dx[0])
            ax3.set_xlim(imin*dx[0], imax*dx[0])
            count += 1

        max = np.maximum(var_max2, var_max3)
        ax3.plot([dx[0]*icoll, dx[0]*icoll], [var_min1, max], ':', linewidth=2)
        ax3.fill_betweenx(np.linspace(var_min1, max, 2), dx[0]*(ic1-delta)*np.ones(2), dx[0]*(ic1+delta-2)*np.ones(2),
                          facecolor='green', alpha=0.5, linewidth=0)
        ax3.fill_betweenx(np.linspace(var_min1, max, 2), dx[0]*(ic1+delta-2)*np.ones(2), dx[0]*nx*np.ones(2),
                          facecolor='gray', linewidth=0)
        for i,ax in enumerate(axis[1:]):
            ax.plot([dx[0]*ic1, dx[0]*ic1], [var_min1, max], 'k--', linewidth=1.5)
            ax.plot([dx[0]*ic2, dx[0]*ic2], [var_min1, max], 'k--', linewidth=1.5)
            ax.set_xlabel('y')
            ax.set_ylabel(var_name)
            ax.grid()
        ax2.set_title('through center of coldpool #1 (max='+str(np.round(var_max1,1))+')')
        ax3.set_title('through collision point (max(2CP)='+str(np.round(var_max2,1))
                      + ', max(3CP)='+str(np.round(var_max3,1))+')')
        ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
                   fancybox=True, shadow=False, ncol=6, fontsize=10)
        plt.suptitle(var_name + ' (z='+str(k0*dx[2])+'m, k='+str(k0)+')', fontsize=21  )
        plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.95, hspace=0.26)
        plt.savefig(os.path.join(path_out, 'xz_plane_' + var_name + '_z' + str(np.int(k0*dx[2])) + 'm.png'))
        plt.close()
    return



def plot_y_crosssections(var_list, k0, i0, ic1, jc1, ic2, jc2, files, path_out, path_fields):
    cm = plt.cm.get_cmap('coolwarm')
    # cm = plt.cm.get_cmap('viridis')

    var_min1 = 9999.9
    var_max1 = -9999.9
    for var_name in var_list:
        f1 = plt.figure(figsize=(6, 12))
        ax1 = plt.subplot(3, 1, 1)
        ax1.set_title('entropy')
        ax2 = plt.subplot(3, 1, 2)
        ax2.set_title('y-crossection through center of coldpool #1')
        ax3 = plt.subplot(3, 1, 3)
        ax3.set_title('y-crossection through collision point')
        s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
        ax1.imshow(s0[:,:,k0].T)
        ax1.plot([ic1,ic1],[0,2*jc1], 'k', linewidth=2)
        ax1.plot([i0,i0],[0,2*jc1], 'k', linewidth=2)

        for file in files:
            var = read_in_netcdf_fields(var_name, os.path.join(path_fields, file))
            var_min1 = np.minimum(var_min1, np.amin(var[:, jc1, k0]))
            var_max1 = np.maximum(var_max1, np.amax(var[:, jc1, k0]))
            count_color = np.double(file[:-3]) / np.double(files[-1][:-3])
            ax2.plot(var[ic1, :, k0], color=cm(count_color), label='t=' + str(file[:-3]))
            ax3.plot(var[i0, :, k0], color=cm(count_color), label='t=' + str(file[:-3]))

        plt.subplot(3, 1, 1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.subplot(3, 1, 2)
        plt.xlabel('y')
        plt.ylabel(var_name)
        plt.subplot(3, 1, 3)
        plt.xlabel('y')
        plt.ylabel(var_name)
        #
        # ax3.legend( )
        # Put a legend below current axis
        ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                   fancybox=True, shadow=True, ncol=5)
        plt.suptitle(var_name + '  (k=' + str(k0) + ')', fontsize=21)
        plt.savefig(os.path.join(path_out, 'yz_plane_' + var_name + '_k'+str(k0)+'.png'))
        plt.close()
    return



def plot_xz_crosssections_multilevel(var_name, j0, krange, imin, imax,
                                     icoll, jcoll, ic_arr, jc_arr,
                                     files, path_out, path_fields):
    cm = plt.cm.get_cmap('coolwarm')
    cm_cw = plt.cm.get_cmap('coolwarm')

    nline = 5

    fig, axis = plt.subplots(nline, 1, figsize=(6, 18))
    ax1 = axis[0]
    s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
    ax1.imshow(s0[:, :, 1].T, origin='lower', cmap=cm_cw)
    ax1.plot([0, nx], [j0, j0], 'k', linewidth=2)
    ax1.plot([icoll, icoll], [0, ny], 'k:', linewidth=2)
    ax1.set_xlim(imin, imax)
    ax1.set_ylim(imin, imax)

    max = np.zeros(len(krange))
    var_min = 9999.9
    var_max = -9999.9
    for file in files:
        var = read_in_netcdf_fields(var_name, os.path.join(path_fields, file))
        var_min = np.minimum(var_min, np.amin(var))
        var_max = np.maximum(var_max, np.amax(var))
        count_color = np.double(file[:-3]) / np.double(files[-1][:-3])
        for i,ax in enumerate(axis[1:]):
            ax.plot(var[:, j0, krange[i]], color=cm(count_color), label='t=' + str(file[:-3]))
            ax.plot([icoll, icoll], [var_min, var_max], ':', linewidth=2)
            k = krange[i]
            max[i] = np.maximum(max[i], np.amax(var[:,j0,k]))
    for i,ax in enumerate(axis[1:]):
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_xlim(imin, imax)
        k = krange[i]
        ax.set_title('k=' + str(k) + ', z=' + str(dx[2] * k) + 'm (max=' + str(np.round(max[i], 1)) + ')')
        ax.set_xlabel('y')
        ax.set_ylabel(var_name)
    axis[-1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
                   fancybox=True, shadow=True, ncol=5)
    plt.suptitle(var_name + ' multilevel')
    plt.subplots_adjust(bottom=0.11, right=.95, left=0.1, top=0.95, hspace=0.3)
    plt.savefig(os.path.join(path_out, 'xz_multilevel_'+var_name+'.png'))
    plt.close(fig)
    return



def plot_yz_crosssections_multilevel(var_name, i0, krange, jmin, jmax,  jc_arr, files, path_out, path_fields):
    cm = plt.cm.get_cmap('coolwarm')
    jc1 = jc_arr[0]
    nline = 5

    fig, axis = plt.subplots(nline, 1, figsize=(6, 18))
    ax1 = axis[0]
    s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
    ax1.imshow(s0[:, :, 1].T, origin='lower')
    ax1.plot([i0, i0], [0, 4 * jc1], 'k', linewidth=2)
    ax1.set_xlim(jmin,jmax)
    ax1.set_ylim(jmin,jmax)

    max = np.zeros(len(krange))
    for file in files:
        var = read_in_netcdf_fields(var_name, os.path.join(path_fields, file))
        count_color = np.double(file[:-3]) / np.double(files[-1][:-3])
        for i,ax in enumerate(axis[1:]):
            ax.plot(var[i0, :, krange[i]], color=cm(count_color), label='t=' + str(file[:-3]))
            k = krange[i]
            max[i] = np.maximum(max[i], np.amax(var[i0,:,k]))

    for i,ax in enumerate(axis[1:]):
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_ylim(jmin, jmax)
        k = krange[i]
        ax.set_title('k=' + str(k) + ', z=' + str(dx[2] * k) + 'm (max=' + str(np.round(max[i], 1)) + ')')
        ax.set_xlabel('y')
        ax.set_ylabel(var_name)
    axis[-1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
                    fancybox=True, shadow=True, ncol=5)
    plt.suptitle(var_name + ' multilevel')
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.1, top=0.95, hspace=0.3)
    plt.savefig(os.path.join(path_out, 'yz_multilevel_'+var_name+'.png'))
    plt.close()
    return




# ----------------------------------
def set_input_parameters(args):
    print ''' setting parameters '''
    global path_in, path_fields, path_out
    path_in = args.path
    if os.path.exists(os.path.join(path_in, 'fields')):
        path_fields = os.path.join(path_in, 'fields')
    path_out = os.path.join(path_in, 'figs_crosssections')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    print('path: ', path_in)
    print('path out: ', path_out)
    print('')

    global case_name
    case_name = args.casename
    case = case_name[12:-3]
    print('')
    print('casename: ' + case_name)
    print('case: ' + case)
    print('')


    ''' determine file range '''
    if args.tmin:
        time_min = np.int(args.tmin)
    else:
        time_min = np.int(100)
    if args.tmax:
        time_max = np.int(args.tmax)
    else:
        time_max = np.int(10000)
    times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
             and np.int(name[:-3]) >= time_min and np.int(name[:-3]) <= time_max]
    times.sort()
    print('times', times)
    files = [str(t) + '.nc' for t in times]
    # print('files', files)
    print('')

    if args.kmin:
        kmin = np.int(args.kmin)
    else:
        kmin = 0
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = kmin
    krange = np.arange(kmin, kmax + 1)
    print('krange: ', krange)
    print ''

    return case, times, files, krange

# ----------------------------------
def define_geometry(case_name, files):
    print 'define geometry'
    global nx, ny, nz, dx, gw
    nml = simplejson.loads(open(os.path.join(path_in, case_name + '.in')).read())
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = np.ndarray(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']

    # set coordinates for plots
    # (a) double 3D
    global isep

    if case_name == 'ColdPoolDry_single_3D':
        rstar = nml['init']['r']
        try:
            ic = nml['init']['ic']
            jc = nml['init']['jc']
        except:
            ic = np.int(nx/2)
            jc = np.int(ny/2)
        ic_arr = np.zeros(1)
        jc_arr = np.zeros(1)
        ic_arr[0] = ic
        jc_arr[0] = jc
    # (b) double 2D
    elif case_name == 'ColdPoolDry_double_2D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
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
        irstar = np.int(np.round(rstar / dx))
        isep = 4 * irstar
        jsep = 0
        # ic1 = np.int(np.round((nx + 2 * gw) / 3)) - gw
        ic1 = np.int(np.round((nx - gw) / 3)) + 1
        ic2 = ic1 + isep
        jc1 = np.int(np.round((ny + 2 * gw) / 2)) - gw
        jc2 = jc1 + jsep
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
        icoll = 0.5*(ic1+ic2)
        jcoll = 0.5*(jc1+jc2)
    elif case_name == 'ColdPoolDry_triple_3D':
        global d, r_int
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        d = nml['init']['d']
        i_d = np.int(np.round(d / dx[0]))
        a = np.int(np.round(i_d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        idhalf = np.int(np.round(i_d / 2))
        r_int = np.int(np.round(np.sqrt(3.) / 6 * i_d))  # radius of inscribed circle
        # point of 3-CP collision (ic, jc)
        ic = np.int(np.round(nx / 2))
        jc = np.int(np.round(ny / 2))
        ic1 = ic - r_int
        ic2 = ic1
        print('------', ic, a, r_int)
        ic3 = ic + (a - r_int)
        jc1 = jc - idhalf
        jc2 = jc + idhalf
        jc3 = jc
        ic_arr = [ic1, ic2, ic3]
        jc_arr = [jc1, jc2, jc3]
        icoll = ic
        jcoll = jc

    print('ic, jc: ', ic_arr, jc_arr)

    ''' --- auxiliary arrays (since no Grid.pyx) ---'''
    global nx_, ny_, nz_
    # test file:
    var = read_in_netcdf_fields('u', os.path.join(path_fields, files[0]))
    [nx_, ny_, nz_] = var.shape
    x_half = np.empty((nx_), dtype=np.double, order='c')
    y_half = np.empty((ny_), dtype=np.double, order='c')
    z_half = np.empty((nz_), dtype=np.double, order='c')
    count = 0
    for i in xrange(nx_):
        x_half[count] = (i + 0.5) * dx[0]
        count += 1
    count = 0
    for j in xrange(ny_):
        y_half[count] = (j + 0.5) * dx[1]
        count += 1
    count = 0
    for i in xrange(nz_):
        z_half[count] = (i + 0.5) * dx[2]
        count += 1

    return icoll, jcoll, ic_arr, jc_arr, x_half, y_half, z_half

# ----------------------------------
def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    # shape = var.shape
    # data = np.ndarray(shape = var.shape)
    data = var[:,:,:]
    rootgrp.close()
    return data



if __name__ == '__main__':
    main()

