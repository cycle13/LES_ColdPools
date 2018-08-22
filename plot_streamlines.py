import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import netCDF4 as nc
import argparse
import json as simplejson
import os

label_size = 12
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 21
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
plt.rcParams['figure.titlesize'] = 35

def main():



    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--timemin")
    parser.add_argument("--timemax")
    parser.add_argument("--k0")
    args = parser.parse_args()

    case_name = args.casename
    path_in = args.path
    if os.path.exists(os.path.join(path_in, 'fields')):
        path_fields = os.path.join(path_in, 'fields')
    elif os.path.exists(os.path.join(path_in, 'fields_k120')):
        path_fields = os.path.join(path_in, 'fields_k120')
    path_out = os.path.join(path_in, 'streamlines')
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


    ''' determine file range '''
    files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc']
    files.sort(key=len)
    print('files', files)
    print('len', len(files[0]))
    print('')
    if args.timemin:
        time_min = np.int(args.timemin)
    else:
        try:
            time_min = np.int(files[0][:-3])
        except:
            time_min = 100
    if args.timemax:
        time_max = np.int(args.timemax)
    else:
        try:
            time_max = np.int(files[-1][:-3])
        except:
            time_max = 3000

    if len(files[0])<7:     # 100.nc
        files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc'
             and np.int(name[:-3])>=time_min and np.int(name[:-3])<=time_max]
        times = [np.int(name[:-3]) for name in files]
    else:
        # print(name, name[:-8])
        files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc'
                 and np.int(name[:-8]) >= time_min and np.int(name[:-8]) <= time_max]
        times = [np.int(name[:-8]) for name in files]
    files.sort(key=len)
    times.sort()
    print('')
    print('files', files)
    print('')
    print('times', times)
    print('')

    ''' --- set locations ---'''
    global ic_arr, jc_arr

    if case_name == 'ColdPoolDry_double_2D':
        rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        # zstar = nml['init']['h']
        isep = 4 * irstar
        ic1 = np.int(nx / 3)
        ic2 = ic1 + isep
        jc1 = np.int(ny / 2)
        jc2 = jc1
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_double_3D':
        rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        # zstar = nml['init']['h']
        isep = 4 * irstar
        jsep = 0
        ic1 = np.int(np.round((nx + 2 * gw) / 3)) - gw
        jc1 = np.int(np.round((ny + 2 * gw) / 2)) - gw
        ic2 = ic1 + isep
        jc2 = jc1 + jsep
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_triple_3D':
        d = np.int(np.round(ny / 2))
        # dhalf = np.int(np.round(ny / 4))
        a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        ic1 = np.int(np.round(a / 2)) + gw
        ic2 = ic1
        ic3 = ic1 + np.int(np.round(a))
        jc1 = np.int(np.round(d / 2) + gw)
        jc2 = jc1 + d
        jc3 = jc1 + np.int(np.round(d / 2))


    ''' --- auxiliary arrays (since no Grid.pyx) ---'''
    # test file:
    var = read_in_netcdf_fields('u', os.path.join(path_fields, files[0]))
    [nx_,ny_,nz_] = var.shape

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

    ''' --- call plotting functions --- '''
    if args.k0:
        k0 = np.int(args.k0)
    else:
        k0 = 1
    nt = 0
    vel = np.ndarray((2,nx_,ny_,nz_))
    for it, file in enumerate(files):
        t0 = times[it]
        print('t', t0)
        vel[0,:,:,:] = read_in_netcdf_fields('u', os.path.join(path_fields, file))
        vel[1,:,:,:] = read_in_netcdf_fields('v', os.path.join(path_fields, file))
        speed_h = np.sqrt(vel[0,:] * vel[0,:] + vel[1,:] * vel[1,:])
        w = read_in_netcdf_fields('w', os.path.join(path_fields, file))
        speed_yz = np.sqrt(vel[1, :] * vel[1, :] + w * w)
        speed_xz = np.sqrt(vel[0, :] * vel[0, :] + w * w)


        ''' (a) xy-level '''
        i0 = np.int(np.round( ic1 + np.double(isep) / 2 ))
        di = 15
        j0 = jc1
        dj = np.int(2*irstar)
        # plot_streamplot_xy_collision(w, vel, speed_h, x_half, y_half, i0, di, j0, dj, k0, t0, path_out)
        # plot_streamplot_xy_varythickness(w, vel, x_half, y_half, speed_h, k0, t0, path_out)

        ''' (b) yz-plane at collision point'''
        i0 = np.int(np.round( ic1 + np.double(isep) / 2 ))  # at collision point
        # plot_streamplot_yz(w, vel, speed_yz, y_half, z_half, i0, t0, path_out, True)

        ''' (c) yz-plane at center of cold pool #1'''
        i0 = ic1    # through center of cold pool #1
        # plot_streamplot_yz(w, vel, speed_yz, y_half, z_half, i0, t0, path_out, True)

        ''' (d) xz-plane at center of cold pool #1'''
        j0 = jc1
        # plot_streamplot_xz(w, vel, speed_xz, x_half, z_half, j0, t0, path_out, True)

    return


def plot_streamplot_xz(w, vel, speed, x_arr, z_arr, j0, t0, path_out, vary=False):
    print(path_out)
    if t0 <= 100:
        plt.figure()
        plt.contourf(w[:, :, 1].T)
        plt.plot([0, nx],[j0, j0], 'k')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig(os.path.join(path_out, 'streamlines_xz_check_i'+str(j0)+'.png'))
        plt.close()

    cm = plt.cm.get_cmap('bwr')
    w_ = w[:, j0, :]
    wmax = np.maximum(np.abs(np.amin(w_)), np.abs(np.amax(w_)))
    levels = np.linspace(-wmax, wmax, 1e3)

    plt.figure(figsize=(12, 10))
    ax = plt.contourf(x_arr, z_arr, w_.T,cmap=cm, levels=levels)
    plt.colorbar(ax)

    if vary:
        lw = 5 * speed[:, j0, :] / speed[:, j0, :].max()
        plt.streamplot(x_arr, z_arr, vel[0,:,j0,:].T, w_.T,
                       color='k', density = 1.5, linewidth=lw[:,:].T)
    else:
        plt.streamplot(x_arr, z_arr, vel[0, :, j0, :].T, w_.T,
                       color='k', density=1.5, linewidth=2)

    plt.xlabel('x [m]   (dx=' + str(dx) + ')')
    plt.ylabel('z [m]   (dz=' + str(dy) + ')')
    plt.title('t=' + str(t0) + ' x=' + str(j0 * dy) + 'm')
    plt.savefig(os.path.join(path_out, 'streamlines_xz_t' + str(t0) + '_i' + str(j0) + '.png'))
    plt.close()

    return



def plot_streamplot_yz(w, vel, speed, y_arr, z_arr, i0, t0, path_out, vary=True):
    if t0 <= 100:
        plt.figure()
        plt.contourf(w[:,:,1].T)
        plt.plot([i0,i0],[0,ny],'k')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig(os.path.join(path_out, 'streamlines_yz_check_i'+str(i0)+'.png'))
        plt.close()

    cm = plt.cm.get_cmap('bwr')
    w_ = w[i0, :,:]
    # levels = np.linspace(np.amin(w_), np.amax(w_), 1e3)
    wmax = np.maximum(np.abs(np.amin(w_)), np.abs(np.amax(w_)))
    levels = np.linspace(-wmax, wmax, 1e3)

    plt.figure(figsize=(12, 10))
    ax = plt.contourf(y_arr, z_arr, w_.T,cmap=cm, levels=levels)

    plt.colorbar(ax)
    if not vary:
        plt.streamplot(y_arr, z_arr, vel[1,i0,:,:].T, w_.T, color='k', density=1.5, linewidth=1.5)
    elif vary:
        #  Varying line width along a streamline
        lw = 5 * speed[i0,:,:] / speed[i0,:,:].max()
        plt.streamplot(y_arr, z_arr, vel[1, i0, :, :].T, w_.T, color='k',
                    density=1.5,
                    linewidth = lw[:, :].T)

    plt.xlabel('y [m]   (dy=' + str(dy) + ')')
    plt.ylabel('z [m]   (dz=' + str(dy) + ')')
    plt.title('t=' + str(t0)+ ' x='+str(i0*dx)+'m')
    plt.savefig(os.path.join(path_out, 'streamlines_yz_i'+str(i0)+'_t'+str(t0)+'.png'))
    plt.close()
    return



def plot_streamplot_xy_collision(w, vel, speed, x_arr, y_arr, ic, di, jc, dj, k0, t0, path_out):
    cm_cont = plt.cm.get_cmap('bwr')
    cm_lines = plt.cm.get_cmap('winter')
    # cm_lines = plt.cm.get_cmap('autumn')

    w_ = w[ic-di:ic+di+1, jc-dj:jc+dj+1, k0]
    wmax = np.maximum(np.abs(np.amin(w_)), np.abs(np.amax(w_)))
    levels = np.linspace(-wmax, wmax, 1e3)

    #  Varying line width along a streamline
    lw = 5 * speed[ic-di:ic+di+1, jc-dj:jc+dj+1, k0] / speed[ic-di:ic+di+1, jc-dj:jc+dj+1, k0].max()
    plt.figure(figsize=(10, 12))
    ax = plt.contourf(x_arr[ic - di:ic + di + 1], y_arr[jc - dj:jc + dj + 1], w_.T,
                      cmap=cm_cont, levels=levels)
    plt.colorbar(ax)
    plt.streamplot(x_arr[ic - di:ic + di + 1], y_arr[jc - dj:jc + dj + 1],
                   vel[0, ic - di:ic + di + 1, jc - dj:jc + dj + 1, k0].T,
                   vel[1, ic - di:ic + di + 1, jc - dj:jc + dj + 1, k0].T,
                   density=1.5, linewidth=lw[:,:].T, color='k')

    plt.xlabel('x [m]   (dx=' + str(dx) + ')')
    plt.ylabel('y [m]   (dy=' + str(dy) + ')')
    plt.title('t=' + str(t0) + ', z=' + str(dz * k0))
    plt.savefig(os.path.join(path_out, 'streamlines_xy_collision_lw_t' + str(t0) + '_k' + str(k0) + '.png'))
    plt.close()

    # Varying color along a streamline
    fig = plt.figure(figsize=(12, 12))
    lw = 5 * speed[ic - di:ic + di + 1, jc - dj:jc + dj + 1, k0] / \
         speed[ic - di:ic + di + 1, jc - dj:jc + dj + 1, k0].max()
    ax = plt.contourf(x_arr[ic - di:ic + di + 1], y_arr[jc - dj:jc + dj + 1], w_.T,
                      cmap=cm_cont, levels=levels, alpha=1.)
    plt.contourf(x_arr[ic - di:ic + di + 1], y_arr[jc - dj:jc + dj + 1], np.zeros((2*di+1,2*dj+1)).T,
                      colors='w', alpha=0.5)
    plt.colorbar(ax, shrink=0.5)
    u = vel[0,ic-di:ic+di+1,jc-dj:jc+dj+1,k0].T
    v = vel[1,ic-di:ic+di+1,jc-dj:jc+dj+1,k0].T
    strm = plt.streamplot(x_arr[ic-di:ic+di+1], y_arr[jc-dj:jc+dj+1], u, v,
                          color=u, cmap=cm_lines,
                          linewidth=lw[:,:].T)
    fig.colorbar(strm.lines, shrink=0.5)
    plt.title('Varying Color')
    plt.savefig(os.path.join(path_out, 'streamlines_xy_collision_col_t' + str(t0) + '_k' + str(k0) + '.png'))
    plt.close()
    return



def plot_streamplot_xy_varythickness(w, vel, x_arr, y_arr, speed, k0, t0, path_out):

    plt.figure(figsize=(12,10))
    lw = 5 * speed[:,:,k0] / speed[:,:,k0].max()
    ax = plt.contourf(x_arr, y_arr, w[:,:,k0].T, alpha=0.5)
    plt.colorbar(ax)
    plt.streamplot(x_arr, y_arr, vel[0,:,:,k0].T ,vel[1,:,:,k0].T,
                   color='k', density=1.5, linewidth=lw[:,:].T)
    plt.xlabel('x [m]   (dx='+str(dx)+')')
    plt.ylabel('y [m]   (dy='+str(dy)+')')
    plt.title('t='+str(t0) + ', z='+str(dz*k0))
    plt.savefig(os.path.join(path_out, 'streamlines_xy_lw_t'+str(t0)+'_k'+str(k0)+'.png'))
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