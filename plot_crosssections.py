import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os


def main():
    print('hello')

    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path")
    args = parser.parse_args()

    case_name = args.casename
    path_in = args.path
    path_fields = os.path.join(path_in, 'fields')
    path_out = os.path.join(path_in, 'profiles')
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    # out_path = os.path.join(path, 'Updrafts_figures_Labels')
    # if not os.path.exists(out_path):
    #     os.mkdir(out_path)

    # global nx, ny, nz, dx, dy, dz
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
    rstar = 5000.0  # half of the width of initial cold-pools [m]
    irstar = np.int(np.round( rstar / dx ))
    # zstar = nml['init']['h']
    isep = 4 * irstar
    jsep = 0
    ic1 = np.int(np.round( ( nx + 2*gw) / 3)) - gw
    jc1 = np.int(np.round( (ny + 2*gw) / 2)) - gw
    ic2 = ic1 + isep
    jc2 = jc1 + jsep
    i0 = np.int(ic1 + np.round( (ic2-ic1)/2 ))
    # xc1 = Gr.x_half[ic1]  # center of cold-pool
    # yc1 = Gr.y_half[jc1]  # center of cold-pool
    # xc2 = Gr.x_half[ic1]  # center of cold-pool
    # yc2 = Gr.y_half[jc2]  # center of cold-pool

    files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc']
    files.sort(key=len)
    print('files', files)

    # compute horizontal convergence


    # plot profiles for w, s
    zrange = np.arange(0, 1)
    var_list = ['w', 's']
    # var_list = ['s']
    k0 = 0
    plot_y_crosssections(var_list, k0, i0, ic1, jc1, ic2, jc2, files, path_out, path_fields)
    plot_x_crosssections(var_list, k0, i0, ic1, jc1, ic2, jc2, files, path_out, path_fields)


    return


# ----------------------------------

def plot_x_crosssections(var_list, k0, i0, ic1, jc1, ic2, jc2, files, path_out, path_fields):
    cm = plt.cm.get_cmap('coolwarm')
    # cm = plt.cm.get_cmap('viridis')
    # cm = plt.cm.get_cmap('viridis_r')
    j0 = np.int(jc1 + (ic2 - ic1) / 4)

    for var_name in var_list:
        var_min1 = 9999.9
        var_max1 = -9999.9
        f1 = plt.figure(figsize=(6,12))
        ax1 = plt.subplot(3, 1, 1)
        ax2 = plt.subplot(3, 1, 2)
        ax3 = plt.subplot(3, 1, 3)
        s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
        ax1.imshow(s0[:, :, k0].T)
        ax1.plot([0, 3*ic1], [jc1, jc1], 'k', linewidth=2)
        ax1.plot([0, 3*ic1], [j0, j0], 'k', linewidth=2)
        ax1.plot([ic1, ic1], [0, 2 * jc1], 'k--', linewidth=1)
        ax1.plot([ic2, ic2], [0, 2 * jc1], 'k--', linewidth=1)
        for file in files:
            var = read_in_netcdf_fields(var_name, os.path.join(path_fields,file))
            if var_name == 's':
                print('s',np.amin(var[:,jc1,k0]), var_min1)
                var_min1 = np.minimum(var_min1,np.maximum(6000,np.amin(var[:,jc1,k0])))
            else:
                var_min1 = np.minimum(var_min1, np.amin(var[:, jc1, k0]))
            var_max1 = np.maximum(var_max1, np.amax(var[:,jc1,k0]))
            count_color = np.double(file[:-3]) / np.double(files[-1][:-3])
            ax2.plot(var[:,jc1,k0], color=cm(count_color), label='t='+str(file[:-3]))
            ax3.plot(var[:,j0,k0], color=cm(count_color), label='t='+str(file[:-3]))

        print('minmax', var_max1, var_min1)
        ax2.plot([ic1, ic1], [var_min1, var_max1], 'k--', linewidth=1)
        ax2.plot([ic2, ic2], [var_min1, var_max1], 'k--', linewidth=1)
        ax3.plot([ic1, ic1], [var_min1, var_max1], 'k--', linewidth=1)
        ax3.plot([ic2, ic2], [var_min1, var_max1], 'k--', linewidth=1)
        plt.subplot(3, 1, 1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.subplot(3,1,2)
        plt.grid()
        plt.title('y-crossection through center of coldpool #1')
        plt.xlabel('y')
        plt.ylabel(var_name)
        plt.subplot(3,1,3)
        plt.grid()
        plt.title('y-crossection through collision point')
        plt.xlabel('y')
        plt.ylabel(var_name)

        # Put a legend below current axis
        ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                   fancybox=True, shadow=True, ncol=5)
        plt.suptitle(var_name + '  (k='+str(k0)+')', fontsize=21  )
        plt.savefig(os.path.join(path_out, 'xz_plane_' + var_name + '_k' + str(k0) + '.png'))
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
        ax2 = plt.subplot(3, 1, 2)
        ax3 = plt.subplot(3, 1, 3)
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
        plt.title('entropy')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.subplot(3, 1, 2)
        plt.title('y-crossection through center of coldpool #1')
        plt.xlabel('y')
        plt.ylabel(var_name)
        plt.subplot(3, 1, 3)
        plt.title('y-crossection through collision point')
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

