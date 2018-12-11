import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import netCDF4 as nc
import argparse
import json as simplejson
import os

global tick_size, label_size
tick_size = 12
label_size = 21
plt.rcParams['xtick.labelsize'] = tick_size
plt.rcParams['ytick.labelsize'] = tick_size
plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = label_size
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35

def main():

    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    global cm_bwr
    cm_bwr = plt.cm.get_cmap('bwr')

    # set_input_parameters(args)
    files, times, nml = set_input_parameters(args)
    x_half, y_half, z_half = define_geometry(case_name, nml, files)


    # print('')
    ''' --- call plotting functions --- '''
    print(path_out)
    var_list = ['w', 'temperature'] # 's'
    cont_var_name = 'w'
    vel_h = np.ndarray((2,nx_,ny_,nz_))

    for it, file in enumerate(files):
        t0 = times[it]
        print('--- time: ', t0, '('+str(it), file+') ---')
        vel_h[0,:,:,:] = read_in_netcdf_fields('u', os.path.join(path_fields, file))
        vel_h[1,:,:,:] = read_in_netcdf_fields('v', os.path.join(path_fields, file))
        w = read_in_netcdf_fields('w', os.path.join(path_fields, file))
        if cont_var_name != 'w':
            cont_var = read_in_netcdf_fields(cont_var_name, os.path.join(path_fields, file))
        else:
            cont_var = w
        speed_h = np.sqrt(vel_h[0, :] * vel_h[0, :] + vel_h[1, :] * vel_h[1, :])
        speed_yz = np.sqrt(vel_h[1, :] * vel_h[1, :] + w * w)
        speed_xz = np.sqrt(vel_h[0, :] * vel_h[0, :] + w * w)


        # --- 1D ---
        ic1 = ic_arr[0]
        jc1 = jc_arr[0]
        i0 = ic1 + irstar
        j0 = jc1 + irstar
        di = np.int(3 * irstar)
        dj = di


        ''' (a) xy-plane '''
        for k0 in krange:
            plot_streamplot_xy_varythickness(cont_var_name, cont_var, vel_h, x_half, y_half, speed_h, k0, t0, path_out)

        #     plot_streamplot_xy_collision(cont_var_name, cont_var, vel, speed_h, x_half, y_half,
        #                                  i0, di, j0, dj, k0, t0, path_out)
    #
    #     # ''' (b) yz-plane at collision point'''
    #     # # --- 2D ---
    #     # i0 = np.int(np.round( ic1 + np.double(isep) / 2 ))  # at collision point
    #     # plot_streamplot_yz(cont_var_name, cont_var, w, vel, speed_yz, y_half, z_half, i0, t0, path_out, True)
    #     #
    #     # ''' (c) yz-plane at center of cold pool #1'''
        # i0 = ic1    # through center of cold pool #1
        # kmax = 40
        # jmin = 20
        # jmax = ny-jmin
        # plot_streamplot_yz(cont_var_name, cont_var, w, vel, speed_yz, y_half, z_half,
        #                    i0, jmin, jmax, kmax, t0, path_out, True)

    #     # ''' (d) xz-plane at center of cold pool #1'''
    #     # # --- 2D ---
    #     # j0 = jc1
    #     # # --- 3D ---
    #     # # j0 = jc3
    #     # # plot_streamplot_xz(cont_var_name, cont_var, w, vel, speed_xz, x_half, z_half, j0, t0, path_out, True)

    return


def plot_streamplot_xz(cont_var_name, cont_var, w, vel, speed, x_arr, z_arr, j0, t0, path_out, vary=False):
    # print(path_out)
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



def plot_streamplot_yz(cont_var_name, cont_var, w, vel, speed,
                       y_arr, z_arr, i0, jmin, jmax, kmax, t0, path_out, vary=True):
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
    # # levels = np.linspace(np.amin(w_), np.amax(w_), 1e3)
    # wmax = np.maximum(np.abs(np.amin(w_)), np.abs(np.amax(w_)))
    # levels = np.linspace(-wmax, wmax, 1e3)
    var_ = cont_var[i0,:,:]
    if cont_var_name == 'w':
        varmax = np.maximum(np.abs(np.amin(var_)), np.abs(np.amax(var_)))
        levels = np.linspace(-varmax, varmax, 1e3)
    else:
        varmin = np.amin(var_)
        varmax = np.amax(var_)
        levels = np.linspace(varmin, varmax, 1e3)

    plt.figure(figsize=(12, 10))
    ax = plt.contourf(y_arr[jmin:jmax], z_arr[:kmax], var_[jmin:jmax,:kmax].T,cmap=cm, levels=levels)

    plt.colorbar(ax)
    if not vary:
        plt.streamplot(y_arr[jmin:jmax], z_arr[:kmax], vel[1,i0,jmin:jmax,:kmax].T, w_[:,:kmax].T,
                       color='k', density=1.5, linewidth=1.5)
    elif vary:
        #  Varying line width along a streamline
        lw = 5 * speed[i0,:,:] / speed[i0,:,:].max()
        plt.streamplot(y_arr[jmin:jmax], z_arr[:kmax], vel[1, i0, jmin:jmax, :kmax].T, w_[jmin:jmax,:kmax].T,
                       color='k',
                    density=1.5,
                    linewidth = lw[jmin:jmax, :kmax].T)

    plt.xlabel('y [m]   (dy=' + str(dy) + ')')
    plt.ylabel('z [m]   (dz=' + str(dy) + ')')
    plt.title('t=' + str(t0)+ ', x='+str(i0*dx)+'m')
    plt.savefig(os.path.join(path_out, 'streamlines_yz_i'+str(i0)+'_t'+str(t0)+'.png'))
    plt.close()
    return





def plot_streamplot_xy_varythickness(cont_var_name, cont_var, vel, x_arr, y_arr, speed, k0, t0, path_out):
    cm = plt.cm.get_cmap('bwr')
    cm_lines = plt.cm.get_cmap('winter')

    cont_var_ = cont_var[:, :, k0]
    if cont_var_name == 'w':
        wmax = np.maximum(np.abs(np.amin(cont_var_)), np.abs(np.amax(cont_var_)))
        levels = np.linspace(-wmax, wmax, 1e3)
        del wmax
    else:
        levels = np.linspace(np.amin(cont_var_), np.amax(cont_var_))
    del cont_var_

    fig, ax = plt.subplots(figsize=(16,10))
    ax.set_aspect('equal')    # ax.set_aspect(1.0)
    if np.abs(speed[:,:,k0].max()) > 0.0:
        lw = 5 * speed[:,:,k0] / speed[:,:,k0].max()
    else:
        lw = 2 * np.ones(shape=speed[:,:,k0].shape)
    ax1 = plt.contourf(x_arr, y_arr, cont_var[:,:,k0].T, levels=levels, alpha=1., cmap = cm)
    plt.colorbar(ax1, shrink=0.5)
    # # plt.streamplot(x_arr, y_arr, vel[0,:,:,k0].T ,vel[1,:,:,k0].T,
    # #                color='k', density=1.5, linewidth=lw[:,:].T)
    strm = plt.streamplot(x_arr, y_arr, vel[0,:,:,k0].T ,vel[1,:,:,k0].T,
                   color=vel[0,:,:,k0], cmap=cm_lines, density=1.5, linewidth=lw[:,:].T)
    plt.colorbar(strm.lines, shrink=0.5)
    plt.xlabel('x [m]   (dx='+str(dx)+')')
    plt.ylabel('y [m]   (dy='+str(dy)+')')
    plt.title('t='+str(t0) + ', z='+str(dz*k0), fontsize=label_size)

    plt.savefig(os.path.join(path_out, 'streamlines_xy_lw_t'+str(t0)+'_k'+str(k0)+'.png'))
    plt.close()
    return




# ----------------------------------
def set_input_parameters(args):
    global path_in, path_out, path_fields
    path_in = args.path
    if os.path.exists(os.path.join(path_in, 'fields')):
        path_fields = os.path.join(path_in, 'fields')
    elif os.path.exists(os.path.join(path_in, 'fields_k120')):
        path_fields = os.path.join(path_in, 'fields_k120')
    path_out = os.path.join(path_in, 'streamlines')
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    global case_name
    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path_in, case_name + '.in')).read())
    global nx, ny, nz, dx, dy, dz
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    gw = nml['grid']['gw']
    # nx = 100

    files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc']
    files.sort(key=len)
    if args.tmin:
        time_min = np.int(args.tmin)
    else:
        try:
            time_min = np.int(files[0][:-3])
        except:
            time_min = 100
    if args.tmax:
        time_max = np.int(args.tmax)
    else:
        try:
            time_max = np.int(files[-1][:-3])
        except:
            time_max = 3000

    global krange
    if args.kmin:
        kmin = np.int(args.kmin)
    else:
        kmin = 1
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = 1
    krange = np.arange(kmin, kmax+1)


    ''' determine file range '''
    if len(files[0]) < 7:  # 100.nc
        files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc'
                 and np.int(name[:-3]) >= time_min and np.int(name[:-3]) <= time_max]
        times = [np.int(name[:-3]) for name in files]
        times.sort()
        # print(type(times), times)
        for it,t0 in enumerate(times):
            files[it] = str(t0)+'.nc'
    else:  # 100_k120.nc
        files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc'
                 and np.int(name[:-8]) >= time_min and np.int(name[:-8]) <= time_max]
        times = [np.int(name[:-8]) for name in files]
        times = times.sort()
        for it,t0 in enumerate(times):
            files[it] = str(t0)+files[0][3:]

    print('')
    print('files', files)
    print('len', len(files[0]))
    print('times', times)
    print('krange', krange)
    print('')

    return files, times, nml



def define_geometry(case_name, nml, files):
    # print('--- define geometry ---')
    global ic_arr, jc_arr
    global ic, jc, isep
    global rstar, irstar
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
    elif case_name == 'ColdPoolDry_single_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        # zstar = nml['init']['h']
        dTh = nml['init']['dTh']
        marg_i = np.int(500. / np.round(dx))  # width of margin
        marg = marg_i * dx  # width of margin
        ic = np.int(nx / 2)
        jc = np.int(ny / 2)
        # xc = Gr.x_half[ic + Gr.dims.gw]  # center of cold-pool
        # yc = Gr.y_half[jc + Gr.dims.gw]  # center of cold-pool
        ic_arr = [ic]
        jc_arr = [jc]
    elif case_name == 'ColdPoolDry_double_3D':
        try:
            rstar = nml['init']['r']
        except:
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
    # return

if __name__ == '__main__':
    main()