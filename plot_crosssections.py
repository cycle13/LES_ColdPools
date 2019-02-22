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

    global path_in, path_fields, path_out
    case_name = args.casename
    case = case_name[12:-3]
    print('')
    print('casename: ' + case_name)
    print('case: '+ case)
    print('')

    path_in = args.path
    if os.path.exists(os.path.join(path_in, 'fields')):
        path_fields = os.path.join(path_in, 'fields')
    elif os.path.exists(os.path.join(path_in, 'fields_k120')):
        path_fields = os.path.join(path_in, 'fields_k120')
    path_out = os.path.join(path_in, 'figs_crosssections')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    print('path: ', path_in)
    print('path out: ', path_out)
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
             and np.int(name[:-3])>=time_min and np.int(name[:-3])<=time_max]
    times.sort()
    print('times', times)
    files = [str(t)+'.nc' for t in times]
    print('files', files)

    x_half, y_half, z_half = define_geometry(case_name, files)

    ''' --- auxiliary arrays (since no Grid.pyx) ---'''
    # test file:
    var = read_in_netcdf_fields('s', os.path.join(path_fields, files[0]))
    plot_configuration(var)

    # [nx_, ny_, nz_] = var.shape
    # x_half = np.empty((nx_), dtype=np.double, order='c')
    # y_half = np.empty((ny_), dtype=np.double, order='c')
    # z_half = np.empty((nz_), dtype=np.double, order='c')
    # count = 0
    # for i in xrange(nx_):
    #     x_half[count] = (i + 0.5) * dx
    #     count += 1
    # count = 0
    # for j in xrange(ny_):
    #     y_half[count] = (j + 0.5) * dy
    #     count += 1
    # count = 0
    # for i in xrange(nz_):
    #     z_half[count] = (i + 0.5) * dz
    #     count += 1
    #
    #
    # ''' (a) compute horizontal convergence''''
    #
    #
    ''' (b) plot profiles for w, s for all time-steps '''
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

    zrange = np.arange(0, 1)
    var_list = ['w', 's']
    # var_list = ['s']
    for k0 in krange:
        print('---- k0: '+str(k0) + '----')
        # if case == 'double':
        #     # -- 2D --
        #     i0 = ic1 + isep
        #     i_collision = np.int(ic1 + np.round((ic2 - ic1) / 2))  # double 3D
        #     # plot_y_crosssections(var_list, k0, i0, ic1, jc1, ic2, jc2, files, path_out, path_fields)
        #
        #     # -- 2D --
        #     # j0 = np.int(jc1 + (ic2 - ic1) / 4)
        #     # plot_xz_crosssections(var_list, j0, k0, ic1, jc1, ic2, jc2, files, path_out, path_fields)
        if case == 'triple':
            # -- 3D --
            j0 = np.int(jc_arr[0] + isep)
            plot_xz_crosssections(var_list, j0, k0,
                                  ic_arr[0], jc_arr[0], ic_arr[2], jc_arr[2],
                                  x_half, y_half, files, path_out, path_fields)

    #
    # ''' (c) plot profiles at different levels for collision time '''
    # # krange_ = [0, 2, 5,10]
    # krange_ = [1, 3, 4, 6]
    # krange_ = [0, 1, 2, 3]
    # # krange_ = [7, 8, 9, 10]
    # # krange_ = [0, 1, 2, 3]
    # # krange_ = [4, 5, 6, 7]
    # if case == 'single':
    #     jmin = jc_arr[0] - 80
    #     jmax = jc_arr[0] + 80
    #     plot_yz_crosssections_multilevel('s', ic_arr[0], krange_, jmin, jmax,  jc_arr, files, path_out, path_fields)
    # if case == 'double':
    #     # -- 2D --
    #     i_collision = np.int(1/2*(ic_arr[0]+ic_arr[1]))
    #     plot_yz_crosssections_multilevel('w', i_collision, krange, jmin, jmax,  jc_arr, files, path_out, path_fields)
    # elif case == 'triple':
    #     # -- 3D --
    #     jc1 = jc_arr[0]
    #     j0 = np.int(jc1 + isep)
    #     # # plot_xz_crosssections_multilevel('w', j0, krange, ic_arr, jc_arr, files, path_out, path_fields)
    #     ic1 = ic_arr[0]
    #     ic2 = ic_arr[1]
    #     ic3 = ic_arr[2]
    #     i0 = ic3 + np.int(gw/2)
    #     i0 = 0.5*(ic1+ic2)
    #     plot_yz_crosssections_multilevel('w', i0, krange_, jmin, jmax,   jc_arr, files, path_out, path_fields)
    # print ''

    return


# ----------------------------------

def plot_configuration(var):
    fig, (ax1, ax2)= plt.subplots(1,2)
    ax1.imshow(var[:, :, 0].T, origin='lower')
    ax1.plot([ic_arr[0], ic_arr[0]], [0, ny-1], 'w', linewidth=1)
    ax1.plot([0, nx], [jc_arr[0], jc_arr[0]], 'w', linewidth=1)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    plt.tight_layout()
    ax1.set_xlim(0,nx)
    ax1.set_ylim(0,ny)
    ax2.contourf(var[:, :, 0].T, origin='lower')
    ax2.plot([ic_arr[0], ic_arr[0]], [0, ny-1], 'w', linewidth=1)
    ax2.plot([0, nz], [jc_arr[0], jc_arr[0]], 'w', linewidth=1)
    ax2.set_xlabel('x')
    plt.tight_layout()
    ax2.set_xlim(70,130)
    ax2.set_ylim(120,180)
    print('gw', gw, ic_arr[0])
    ax2.set_aspect('equal')
    plt.savefig(os.path.join(path_out, 'initial.png'))
    plt.close()

    return

# ----------------------------------

def plot_xz_crosssections(var_list, j0, k0, ic1, jc1, ic2, jc2,
                          x_arr, y_arr, files, path_out, path_fields):
    # cm = plt.cm.get_cmap('bwr')
    cm = plt.cm.get_cmap('coolwarm')
    # cm = plt.cm.get_cmap('winter')
    # cm = plt.cm.get_cmap('viridis')
    # cm = plt.cm.get_cmap('viridis_r')

    # x_arr = dx*np.arange(0,nx)


    for var_name in var_list:
        var_min1 = 9999.9
        var_max1 = -9999.9
        var_min2 = 9999.9
        var_max2 = -9999.9

        f1 = plt.figure(figsize=(12,12))
        # f1 = plt.figure(figsize=(6,12))
        ax1 = plt.subplot(3, 1, 1)
        ax2 = plt.subplot(3, 1, 2)
        ax3 = plt.subplot(3, 1, 3)
        s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
        ax1.imshow(s0[:, :, k0].T, origin="lower")
        ax1.plot([0, nx-1], [jc1, jc1], 'k', linewidth=2)
        ax1.plot([0, nx-1], [j0, j0], 'k', linewidth=2)
        ax1.plot([ic1, ic1], [0, ny-1], 'k--', linewidth=1)
        ax1.plot([ic2, ic2], [0, ny-1], 'k--', linewidth=1)
        ax1.set_xlim([0,nx])
        ax1.set_ylim([0,ny])
        count = 0
        for file in files:
            var = read_in_netcdf_fields(var_name, os.path.join(path_fields,file))
            if var_name == 's':
                # print('s',np.amin(var[:,jc1,k0]), var_min)
                var_min1 = np.minimum(var_min1,np.maximum(6000,np.amin(var[:,jc1,k0])))
                var_min2 = np.minimum(var_min2,np.maximum(6000,np.amin(var[:,jc2,k0])))
            else:
                var_min1 = np.minimum(var_min1, np.amin(var[:, jc1, k0]))
                var_min2 = np.minimum(var_min2, np.amin(var[:, jc2, k0]))
            var_max1 = np.maximum(var_max1, np.amax(var[:,jc1,k0]))
            var_max2 = np.maximum(var_max2, np.amax(var[:,jc2,k0]))
            # count_color = np.double(file[:-3]) / np.double(files[-1][:-3])
            count_color = np.double(count) / len(files)
            ax2.plot(x_arr,var[:,jc1,k0], color=cm(count_color), label='t='+str(file[:-3]))
            ax3.plot(x_arr,var[:,j0,k0], color=cm(count_color), label='t='+str(file[:-3]))
            count += 1

        # print('minmax', var_min1, var_max1)
        ax2.plot([dx*ic1, dx*ic1], [var_min1, var_max1], 'k--', linewidth=1.5)
        ax2.plot([dx*ic2, dx*ic2], [var_min1, var_max1], 'k--', linewidth=1.5)
        ax3.plot([dx*ic1, dx*ic1], [var_min2, var_max2], 'k--', linewidth=1.5)
        ax3.plot([dx*ic2, dx*ic2], [var_min2, var_max2], 'k--', linewidth=1.5)

        plt.subplot(3, 1, 1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.subplot(3,1,2)
        ax2.set_title('through center of coldpool #1 (max='+str(np.round(var_max1,1))+')')
        plt.grid()
        plt.xlabel('y')
        plt.ylabel(var_name)
        plt.subplot(3,1,3)
        ax3.set_title('through collision point (max='+str(np.round(var_max2,1))+')')
        plt.grid()
        plt.xlabel('y')
        plt.ylabel(var_name)

        # Put a legend below current axis
        ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                   fancybox=True, shadow=True, ncol=5)
        plt.suptitle(var_name + ' z='+str(k0*dz)+'m  (k='+str(k0)+')', fontsize=21  )
        plt.savefig(os.path.join(path_out, 'xz_plane_' + var_name + '_z' + str(np.int(k0*dz)) + 'm.png'))
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



def plot_xz_crosssections_multilevel(var_name, j0, krange, ic_arr, jc_arr, files, path_out, path_fields):
    cm = plt.cm.get_cmap('coolwarm')
    # cm = plt.cm.get_cmap('coolwarm')
    ic1 = ic_arr[0]
    ic2 = ic_arr[1]
    jc1 = jc_arr[0]
    jc2 = jc_arr[1]

    nline = 5

    plt.figure(figsize=(6, 18))

    ax1 = plt.subplot(nline, 1, 1)
    ax2 = plt.subplot(nline, 1, 2)
    ax3 = plt.subplot(nline, 1, 3)
    ax4 = plt.subplot(nline, 1, 4)
    ax5 = plt.subplot(nline, 1, 5)
    s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
    ax1.imshow(s0[:, :, 1].T, origin='lower')
    ax1.plot([0, nx], [j0, j0], 'k', linewidth=2)

    for file in files:
        var = read_in_netcdf_fields(var_name, os.path.join(path_fields, file))
        # var_min1 = np.minimum(var_min1, np.amin(var[:, jc1, k0]))
        # var_max1 = np.maximum(var_max1, np.amax(var[:, jc1, k0]))
        count_color = np.double(file[:-3]) / np.double(files[-1][:-3])

        ax2.plot(var[:, j0, krange[0]], color=cm(count_color), label='t=' + str(file[:-3]))
        ax3.plot(var[:, j0, krange[1]], color=cm(count_color), label='t=' + str(file[:-3]))
        ax4.plot(var[:, j0, krange[2]], color=cm(count_color), label='t=' + str(file[:-3]))
        ax5.plot(var[:, j0, krange[3]], color=cm(count_color), label='t=' + str(file[:-3]))

    plt.subplot(nline, 1, 1)
    plt.xlabel('x')
    plt.ylabel('y')
    for n,k in enumerate(krange):
        plt.subplot(nline, 1, n+2)
        plt.title('k='+str(k)+', z='+str(dz*k)+'m')
        plt.ylabel(var_name)
        plt.xlabel('y')
    plt.subplot(nline,1,nline)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                   fancybox=True, shadow=True, ncol=5)
    plt.suptitle('w multilevel')
    plt.savefig(os.path.join(path_out, 'yz_multilevel_'+var_name+'.png'))
    plt.close()
    return



def plot_yz_crosssections_multilevel(var_name, i0, krange, jmin, jmax,  jc_arr, files, path_out, path_fields):
    cm = plt.cm.get_cmap('coolwarm')
    # ic1 = ic_arr[0]
    # ic2 = ic_arr[1]
    jc1 = jc_arr[0]
    # jc2 = jc_arr[1]

    nline = 5

    plt.figure(figsize=(6, 18))

    ax1 = plt.subplot(nline, 1, 1)
    ax2 = plt.subplot(nline, 1, 2)
    ax3 = plt.subplot(nline, 1, 3)
    ax4 = plt.subplot(nline, 1, 4)
    ax5 = plt.subplot(nline, 1, 5)
    s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
    ax1.imshow(s0[:, :, 1].T, origin='lower')
    # ax1.plot([ic1, ic1], [0, 2 * jc1], 'k', linewidth=2)
    ax1.plot([i0, i0], [0, 4 * jc1], 'k', linewidth=2)
    ax1.set_xlim(jmin,jmax)
    ax1.set_ylim(jmin,jmax)

    max = np.zeros(len(krange))
    for file in files:
        var = read_in_netcdf_fields(var_name, os.path.join(path_fields, file))
        # var_min1 = np.minimum(var_min1, np.amin(var[:, jc1, k0]))
        # var_max1 = np.maximum(var_max1, np.amax(var[:, jc1, k0]))
        count_color = np.double(file[:-3]) / np.double(files[-1][:-3])

        ax2.plot(var[i0, jmin:jmax, krange[0]], color=cm(count_color), label='t=' + str(file[:-3]))
        ax3.plot(var[i0, jmin:jmax, krange[1]], color=cm(count_color), label='t=' + str(file[:-3]))
        ax4.plot(var[i0, jmin:jmax, krange[2]], color=cm(count_color), label='t=' + str(file[:-3]))
        ax5.plot(var[i0, jmin:jmax, krange[3]], color=cm(count_color), label='t=' + str(file[:-3]))
        for n,k in enumerate(krange):
            max[n] = np.maximum(max[n], np.amax(var[i0,:,k]))

    plt.subplot(nline, 1, 1)
    plt.xlabel('x')
    plt.ylabel('y')
    for n,k in enumerate(krange):
        plt.subplot(nline, 1, n+2)
        plt.title('k='+str(k)+', z='+str(dz*k)+'m, max=' + str(np.round(max[n], 1)))
        plt.ylabel(var_name)
        plt.xlabel('y')
    plt.subplot(nline,1,nline)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                   fancybox=True, shadow=True, ncol=5)
    plt.suptitle(var_name + ' multilevel')
    plt.savefig(os.path.join(path_out, 'yz_multilevel_'+var_name+'.png'))
    plt.close()
    return




# ----------------------------------
def define_geometry(case_name, files):
    print 'define geometry'
    global nx, ny, nz, dx, dy, dz, gw
    nml = simplejson.loads(open(os.path.join(path_in, case_name + '.in')).read())
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    gw = nml['grid']['gw']

    # set coordinates for plots
    # (a) double 3D
    global ic_arr, jc_arr
    global isep

    if case_name == 'ColdPoolDry_single_3D':
        rstar = nml['init']['r']
        irstar = np.int(np.round(rstar / dx))
        zstar = nml['init']['h']
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
        irstar = np.int(np.round(rstar / dx))
        # zstar = nml['init']['h']
        isep = 4 * irstar
        jsep = 0
        # ic1 = np.int(np.round((nx + 2 * gw) / 3)) - gw
        ic1 = np.int(np.round((nx - gw) / 3)) + 1
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
        irstar = np.int(np.round(rstar / dx))
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
        x_half[count] = (i + 0.5) * dx
        count += 1
    count = 0
    for j in xrange(ny_):
        y_half[count] = (j + 0.5) * dy
        count += 1
    count = 0
    for i in xrange(nz_):
        z_half[count] = (i + 0.5) * dz
        count += 1

    return x_half, y_half, z_half

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

