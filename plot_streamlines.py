import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import netCDF4 as nc
import argparse
import json as simplejson
import os

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35$

def main():


    print('hello')

    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path")
    args = parser.parse_args()

    case_name = args.casename
    path_in = args.path
    if os.path.exists(os.path.join(path_in, 'fields')):
        path_fields = os.path.join(path_in, 'fields')
    elif os.path.exists(os.path.join(path_in, 'fields_k120')):
        path_fields = os.path.join(path_in, 'fields_k120')
    path_out = os.path.join(path_in, 'profiles')
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    global nx, ny, nz, dx, dy, dz
    nml = simplejson.loads(open(os.path.join(path_in, case_name + '.in')).read())
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    gw = nml['grid']['gw']
    # nx = 100

    # auxiliary since no Grid.pyx
    x_half = np.empty((nx), dtype=np.double, order='c')
    y_half = np.empty((ny), dtype=np.double, order='c')
    count = 0
    for i in xrange(nx):
        x_half[count] = (i + 0.5) * dx
        count += 1
    count = 0
    for j in xrange(ny):
        y_half[count] = (j + 0.5) * dy
        count += 1

    files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc']
    files.sort(key=len)
    print('files', files)

    nt = 0
    nx_ = 200
    vel = np.ndarray((2,nx_,ny,nz))
    file = files[nt]
    t0 = file[:-3]
    print('t', t0)
    vel[0,:,:,:] = read_in_netcdf_fields('u', os.path.join(path_fields, file))
    vel[1,:,:,:] = read_in_netcdf_fields('v', os.path.join(path_fields, file))
    vel_ = vel[:,:nx,:,:]
    speed = np.sqrt(vel[0,:] * vel[0,:] + vel[1,:] * vel[1,:])
    w = read_in_netcdf_fields('w', os.path.join(path_fields, file))


    k0 = 1
    # plot_streamplot_xy(w, vel_, x_half, y_half, k0, t0)
    plot_streamplot_xy_varythickness(w, vel, x_half, y_half, speed, k0, t0)


    # example_streamplot()

    return




def plot_streamplot_xy(w, vel, x_arr, y_arr, k0, t0):
    # x_max = 100
    x_max = nx
    cm = plt.cm.get_cmap('bwr')
    plt.figure()
    # plt.subplot(1,2,1)
    # plt.imshow(w[:,:,k0].T)
    # plt.contourf(x_arr, y_arr, w[:,:,k0].T, cmap = cm)
    # !!! only transposed field gives x- and y-coordinate in right order...!!!
    plt.contourf(x_arr, y_arr, w[:x_max,:,k0].T, cmap = cm)
    plt.streamplot(x_arr, y_arr, vel[0,:x_max,:,k0].T, vel[1,:x_max,:,k0].T,
                   density=1.5, linewidth=1.5)
    #
    # plt.xlabel('x [m]   (dx='+str(dx)+')')
    # plt.ylabel('y [m]   (dy='+str(dy)+')')
    # plt.title('t='+str(t0) + ', z='+str(dz*k0))
    plt.savefig('./streamlines_xy_t'+str(t0)+'_k'+str(k0)+'.png')
    plt.close()
    return




def plot_streamplot_xy_varythickness(w, vel, x_arr, y_arr, speed, k0, t0):

    plt.figure()
    lw = 5 * speed[:,:,k0] / speed[:,:,k0].max()
    plt.contourf(x_arr, y_arr, w[:,:,k0].T, alpha=0.5)
    plt.streamplot(x_arr, y_arr, vel[0,:,:,k0].T ,vel[1,:,:,k0].T,
                   color='k', density=1.5, linewidth=lw[:,:])
    plt.xlabel('x [m]   (dx='+str(dx)+')')
    plt.ylabel('y [m]   (dy='+str(dy)+')')
    plt.title('t='+str(t0) + ', z='+str(dz*k0))
    plt.savefig('./streamlines_xy_t'+str(t0)+'_k'+str(k0)+'_lw.png')
    plt.close()
    return


def example_streamplot():
    w = 3
    Y, X = np.mgrid[-w:w:100j, -w:w:100j]
    U = -1 - X ** 2 + Y
    V = 1 + X - Y ** 2
    speed = np.sqrt(U * U + V * V)

    fig = plt.figure(figsize=(7, 9))
    gs = gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[1, 1, 2])

    #  Varying density along a streamline
    ax0 = fig.add_subplot(gs[0, 0])
    ax0.streamplot(X, Y, U, V, density=[0.5, 1])
    ax0.set_title('Varying Density')

    # # Varying color along a streamline
    # ax1 = fig.add_subplot(gs[0, 1])
    # strm = ax1.streamplot(X, Y, U, V, color=U, linewidth=2, cmap='autumn')
    # fig.colorbar(strm.lines)
    # ax1.set_title('Varying Color')
    #
    # #  Varying line width along a streamline
    # ax2 = fig.add_subplot(gs[1, 0])
    # lw = 5 * speed / speed.max()
    # ax2.streamplot(X, Y, U, V, density=0.6, color='k', linewidth=lw)
    # ax2.set_title('Varying Line Width')
    #
    # # Controlling the starting points of the streamlines
    # seed_points = np.array([[-2, -1, 0, 1, 2, -1], [-2, -1, 0, 1, 2, 2]])
    #
    # ax3 = fig.add_subplot(gs[1, 1])
    # strm = ax3.streamplot(X, Y, U, V, color=U, linewidth=2,
    #                       cmap='autumn', start_points=seed_points.T)
    # fig.colorbar(strm.lines)
    # ax3.set_title('Controlling Starting Points')
    #
    # # Displaying the starting points with blue symbols.
    # ax3.plot(seed_points[0], seed_points[1], 'bo')
    # ax3.axis((-w, w, -w, w))
    #
    # # Create a mask
    # mask = np.zeros(U.shape, dtype=bool)
    # mask[40:60, 40:60] = True
    # U[:20, :20] = np.nan
    # U = np.ma.array(U, mask=mask)
    #
    # ax4 = fig.add_subplot(gs[2:, :])
    # ax4.streamplot(X, Y, U, V, color='r')
    # ax4.set_title('Streamplot with Masking')
    #
    # ax4.imshow(~mask, extent=(-w, w, -w, w), alpha=0.5,
    #            interpolation='nearest', cmap='gray', aspect='auto')
    # ax4.set_aspect('equal')
    #
    # plt.tight_layout()
    plt.savefig('./streamlines_example.png')
    plt.close()
    # plt.show()
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
    # return

if __name__ == '__main__':
    main()