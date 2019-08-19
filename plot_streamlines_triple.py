import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import netCDF4 as nc
import argparse
import json as simplejson
import os

global tick_size, label_size
tick_size = 15
label_size = 21
plt.rcParams['xtick.labelsize'] = tick_size
plt.rcParams['ytick.labelsize'] = tick_size
plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = label_size
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['xtick.minor.width'] = 1.
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['ytick.minor.width'] = 1.
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
    parser.add_argument("--hor")
    parser.add_argument("--vert")
    args = parser.parse_args()

    global cm_bwr
    cm_bwr = plt.cm.get_cmap('bwr')

    # set_input_parameters(args)
    files, times, nml = set_input_parameters(args)
    i0_coll, j0_coll, ic_arr, jc_arr, x_half, y_half, z_half = define_geometry(case_name, nml, files)
    ID = os.path.split(path_in)[1]
    print ('id: ', ID)


    print('')
    ''' --- call plotting functions --- '''
    print(path_out)
    var_list = ['w', 'temperature', 's', 'temperature_anomaly']
    var_list = ['w']

    # # --- 2D ---
    # ic1 = ic_arr[0]
    # jc1 = jc_arr[0]
    # i0 = np.int(np.round( ic1 + np.double(isep) / 2 ))
    # di = 15
    # j0 = jc1
    # dj = np.int(2*irstar)
    # --- 3D ---
    di = np.int(i0_coll*2/3)
    dj = (jc_arr[1]-jc_arr[0])

    for it, file in enumerate(files):
        t0 = times[it]
        print('--- time: ', t0, '('+str(it), file+') ---')
        vel_h = np.ndarray((2, nx_, ny_, nz_))
        vel_h[0,:,:,:] = read_in_netcdf_fields('u', os.path.join(path_fields, file))
        vel_h[1,:,:,:] = read_in_netcdf_fields('v', os.path.join(path_fields, file))
        w = read_in_netcdf_fields('w', os.path.join(path_fields, file))
        speed_h = np.sqrt(vel_h[0, :] * vel_h[0, :] + vel_h[1, :] * vel_h[1, :])
        speed_yz = np.sqrt(vel_h[1, :] * vel_h[1, :] + w * w)
        speed_xz = np.sqrt(vel_h[0, :] * vel_h[0, :] + w * w)
        for cont_var_name in var_list:
            print('contourfigure for ' + cont_var_name)

            if cont_var_name == 'w':
                cont_var = w
            elif cont_var_name == 'temperature_anomaly':
                cont_var = read_in_netcdf_fields('temperature', os.path.join(path_fields, file))
                t_mean = np.average(np.average(cont_var[:, :, :kmax+1], axis=0), axis=0)
                cont_var = cont_var[:, :, :kmax] - t_mean
            else:
                cont_var = read_in_netcdf_fields(cont_var_name, os.path.join(path_fields, file))




            if args.hor in ['True','true']:
                ''' (a) xy-plane '''
                for k0 in krange:
                    print('-- k0=' + str(k0) + ' --')
                    plot_streamplot_xy_collision(cont_var_name, cont_var[:,:,:kmax+1], vel_h, speed_h,
                                                 x_half, y_half, ic_arr, jc_arr,
                                                 i0_coll, di, j0_coll, dj, k0, t0, path_out)
                    # plot_streamplot_xy_varythickness(cont_var_name, cont_var, vel_h,
                    #                                  x_half, y_half, speed_h, k0, t0, path_out)

            if args.vert in ['True', 'true']:
                ''' (b) xz-plane at collision point'''
                # --- 2D ---
                # i0 = np.int(np.round( ic1 + np.double(isep) / 2 ))  # at collision point
                # --- 3D ---
                imin = 50
                imax = nx - 50
                plot_streamplot_xz(cont_var_name, cont_var[:,:,:kmax+1], w, vel_h, speed_xz,
                                   x_half, z_half, j0_coll, imin, imax, kmax, t0, path_out, ID, True)

                ''' (c) yz-plane at center of cold pool #1'''
                i0 = ic_arr[0]    # through center of cold pool #1
                jmin = 50
                jmax = 350
                plot_streamplot_yz(cont_var_name, cont_var[:,:,:kmax+1], w, vel_h, speed_yz,
                                   y_half, z_half, i0, jmin, jmax, kmax, t0, path_out, True)

    #         ''' (d) xz-plane at center of cold pool #1'''
    #         imin = 10
    #         imax = nx - 40
    #         # # --- 2D ---
    #         # j0 = jc1
    #         # # --- 3D ---
    #         j0 = jc_arr[2]
    #         # if cont_var_name == 'temperature_anomaly':
    #         #     t_mean = np.average(cont_var[:, j0, :kmax], axis=0)
    #         #     cont_var = cont_var[:, j0, :kmax] - t_mean
    #         # else:
    #         #     cont_var = cont_var[:, j0, :kmax]
    #         plot_streamplot_xz(cont_var_name, cont_var, w, vel_h, speed_xz,
    #                            x_half, z_half, j0, imin, imax, kmax, t0, path_out, ID, True)

    return


def plot_streamplot_xz(cont_var_name, cont_var, w, vel, speed, x_arr, z_arr, j0, imin, imax, kmax,
                       t0, path_out, ID, vary=False):
    print('plot streamplot xz')

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
    var_ = cont_var[:, j0, :]
    if cont_var_name == 'w':
        max = np.ceil(np.maximum(np.abs(np.amin(var_)), np.abs(np.amax(var_))))
        # max = 4
        min = -max
    elif cont_var_name == 'temperature_anomaly':
        if ID[4] == '_':
            dTh = np.int(ID[3:4])
        else:
            dTh = np.int(ID[3:5])
        max = dTh
        min = -max
    else:
        min = np.floor(np.amin(var_))
        max = np.ceil(np.amax(var_))
    levels = np.linspace(min, max, 1e3)

    fig, ax = plt.subplots(figsize=(12, 9))
    cf = plt.contourf(x_arr[imin:imax], z_arr[:kmax], var_[imin:imax,:kmax].T,cmap=cm,
                      levels=levels, extend='both')
    cbar = plt.colorbar(cf, shrink=0.6, ticks=np.arange(min, max+1, 1), aspect=12)
    # cbar.ax.set_yticklabels(['0', '1', '2', '>3'])
    if cont_var_name == 'temperature_anomaly':
        # cbar.ax.set_yticklabels(np.arange(min,max+1,1))
        # cbar.set_label('T - <T>', rotation=90)
        cbar.set_label('T - <T>  [K]')
    else:
        cbar.set_label(cont_var_name, rotation=90)
    if vary:
        lw = 5 * speed[:, j0, :kmax] / speed[:, j0, :kmax].max()
        plt.streamplot(x_arr[imin:imax], z_arr[:kmax], vel[0,imin:imax,j0,:kmax].T, w_[imin:imax,:kmax].T,
                       color='k', density = 1.5, linewidth=lw[imin:imax,:].T)
    else:
        plt.streamplot(x_arr[imin:imax], z_arr[:kmax], vel[0, imin:imax, j0, :kmax].T, w_.T,
                       color='k', density=1.5, linewidth=2)

    x_ticks = ax.get_xticks() * 1e-3
    y_ticks = ax.get_yticks() * 1e-3
    ax.set_xticklabels(x_ticks)
    ax.set_yticklabels(y_ticks)
    plt.xlabel('y [km]   (dy=' + str(dy) + ')')
    plt.ylabel('z [km]   (dz=' + str(dy) + ')')
    plt.title('crosssection through center of CP, (t=' + str(t0) + 's)', fontsize=18)#' y=' + str(j0 * dy) + 'm')

    fig_name = cont_var_name + '_streamlines_xz_t' + str(t0) + '_j' + str(j0) + '.png'
    plt.savefig(os.path.join(path_out, fig_name))
    plt.close()

    return



def plot_streamplot_yz(cont_var_name, cont_var, w, vel, speed,
                       y_arr, z_arr, i0, jmin, jmax, kmax, t0, path_out, vary=True):
    print('plot streamplot yz')
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
    var_ = cont_var[i0,:,:]
    if cont_var_name == 'w':
        varmax = np.maximum(np.abs(np.amin(var_)), np.abs(np.amax(var_)))
        # varmax = 4
        varmin = -varmax
    else:
        varmin = np.amin(var_)
        varmax = np.amax(var_)
    levels = np.linspace(varmin, varmax, 1e3)

    fig, ax = plt.subplots(figsize=(12, 9))
    cf = plt.contourf(y_arr[jmin:jmax], z_arr[:kmax], var_[jmin:jmax,:kmax].T,cmap=cm,
                      levels=levels, extend='both')
    cb = plt.colorbar(cf, shrink=0.6, ticks=np.arange(np.floor(varmin), np.floor(varmax) + 1, 1), aspect=12)
    cb.ax.tick_params(width=1, size=4)  # , labelsize=6)
    if not vary:
        plt.streamplot(y_arr[jmin:jmax], z_arr[:kmax], vel[1,i0,jmin:jmax,:kmax].T, w_[:,:kmax].T,
                       color='k', density=1.5, linewidth=1.5)
    elif vary and speed[i0,:,:].max()>0.:
        #  Varying line width along a streamline
        lw = 5 * speed[i0,:,:] / speed[i0,:,:].max()
        plt.streamplot(y_arr[jmin:jmax], z_arr[:kmax], vel[1, i0, jmin:jmax, :kmax].T, w_[jmin:jmax,:kmax].T,
                       color='k',
                    density=1.5,
                    linewidth = lw[jmin:jmax, :kmax].T)

    x_ticks = ax.get_xticks() * 1e-3
    y_ticks = ax.get_yticks() * 1e-3
    ax.set_xticklabels(x_ticks)
    ax.set_yticklabels(y_ticks)
    plt.xlabel('y [km]   (dy=' + str(dy) + ')')
    plt.ylabel('z [km]   (dz=' + str(dy) + ')')
    plt.title('t=' + str(t0)+ ', x='+str(i0*dx)+'m')
    fullpath_out = os.path.join(path_out, 'streamlines_yz_i'+str(i0)+'_t'+str(t0)+'.png')
    plt.savefig(fullpath_out)
    plt.close(fig)
    return





def plot_streamplot_xy_collision(cont_var_name, cont_var, vel, speed,
                                 x_arr, y_arr, ic_arr, jc_arr, ic, di, jc, dj, k0, t0, path_out):
    print('plot streamplot xy (t0='+str(t0)+', k0='+str(k0)+')')
    # cm_cont = plt.cm.get_cmap('bwr')
    cm_cont = plt.cm.get_cmap('gnuplot')
    cm_lines = plt.cm.get_cmap('winter')
    # cm_lines = plt.cm.get_cmap('autumn')

    var_ = cont_var[ic-di:ic+di+1, jc-dj:jc+dj+1, :]
    x_arr_ = x_arr[ic-di:ic+di+1]
    y_arr_ = y_arr[jc-dj:jc+dj+1]
    if cont_var_name == 'w':
        cm_cont = plt.cm.get_cmap('bwr')
        varmax = np.maximum(np.abs(np.amin(var_)), np.abs(np.amax(var_)))
        varmin = -varmax
    elif cont_var_name == 'temperature':
        # varmin = np.amin(var_[:,:,k0])
        # varmax = np.amax(var_[:,:,k0])
        varmin = 296
        varmax = 298.5
    else:
        varmin = np.amin(var_)
        varmax = np.amax(var_)
    levels = np.linspace(varmin, varmax, 1e3)


    #  Varying line width along a streamline
    if np.abs(speed[ic-di:ic+di+1, jc-dj:jc+dj+1, k0].max()) > 0.0:
        lw = 5 * speed[ic-di:ic+di+1, jc-dj:jc+dj+1, k0] / speed[ic-di:ic+di+1, jc-dj:jc+dj+1, k0].max()
    else:
        lw = 5 * np.ones(shape=speed[ic-di:ic+di+1, jc-dj:jc+dj+1, k0].shape)

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_aspect('equal')  # ax.set_aspect(1.0)

    cf = plt.contourf(x_arr_, y_arr_, var_[:,:,k0].T,
                      cmap=cm_cont, levels=levels, extend='min')

    if cont_var_name == 's':
        cb = plt.colorbar(cf, shrink=0.6, ticks=np.arange(np.floor(varmin), np.floor(varmax)+1, 2), aspect=12)
        cf = plt.contourf(x_arr_, y_arr_, var_[:, :, k0].T,
                          cmap=cm_cont, levels=levels)
    elif cont_var_name == 'temperature':
        cb = plt.colorbar(cf, shrink=0.5, ticks=np.arange(np.floor(varmin), np.floor(varmax)+1, .5), aspect=17)
        cf = plt.contourf(x_arr_, y_arr_, var_[:, :, k0].T,
                          cmap=cm_cont, levels=levels, extend='min')
    else:
        cb = plt.colorbar(cf, shrink=0.6, ticks=np.arange(np.floor(varmin), np.floor(varmax)+1, 1), aspect=12)
        cf = plt.contourf(x_arr_, y_arr_, var_[:, :, k0].T,
                          cmap=cm_cont, levels=levels)
    cb.ax.tick_params(width=1, size=4) #, labelsize=6)

    plt.plot(x_arr[ic],y_arr[jc],'go', markersize=10)
    plt.streamplot(x_arr_, y_arr_,
                   vel[0, ic-di:ic+di+1, jc-dj:jc+dj+1, k0].T,
                   vel[1, ic-di:ic+di+1, jc-dj:jc+dj+1, k0].T,
                   color='k', density=1., linewidth=lw[:,:].T)

    plt.plot([ x_arr[ic_arr[0]], x_arr[ic_arr[1]], x_arr[ic_arr[2]] ],
             [y_arr[jc_arr[0]], y_arr[jc_arr[1]], y_arr[jc_arr[2]] ], 'o', markersize=12, color='k')
    plt.plot([x_arr[ic_arr[0]], x_arr[ic_arr[1]] ], [y_arr[jc_arr[0]], y_arr[jc_arr[1]] ], 'k')
    plt.plot([x_arr[ic_arr[1]], x_arr[ic_arr[2]] ], [y_arr[jc_arr[1]], y_arr[jc_arr[2]] ], 'k')
    plt.plot([x_arr[ic_arr[0]], x_arr[ic_arr[2]] ], [y_arr[jc_arr[0]], y_arr[jc_arr[2]] ], 'k')

    x_ticks = ax.get_xticks() * 1e-3
    y_ticks = ax.get_yticks() * 1e-3
    ax.set_xticklabels(x_ticks)
    ax.set_yticklabels(y_ticks)
    plt.xlim([x_arr[ic - di], x_arr[ic + di]])
    plt.ylim(y_arr[jc - dj],y_arr[jc + dj])
    plt.xlabel('x [km]   (dx=' + str(np.round(dx,0)) + 'm)')
    plt.ylabel('y [km]   (dy=' + str(np.round(dy,0)) + 'm)')
    plt.title(cont_var_name + ' with horizontal velocity (t=' + str(t0) + ', z=' + str(dz * k0) + ')', fontsize=21)
    plt.savefig(os.path.join(path_out, cont_var_name+'_streamlines_xy_collision_lw_t' + str(t0) + '_k' + str(k0) + '.png'))
    plt.close()



    # Varying color along a streamline
    fig, ax = plt.subplots(figsize=(10.5, 8))
    ax.set_aspect('equal')  # ax.set_aspect(1.0)

    if cont_var_name == 's':
        cb = plt.colorbar(cf, shrink=0.5, ticks=np.arange(np.floor(varmin), np.floor(varmax)+1, 2), aspect=17)
        cf = plt.contourf(x_arr_, y_arr_, var_[:, :, k0].T,
                          cmap=cm_cont, levels=levels, alpha=1.)

    elif cont_var_name == 'temperature':
        cb = plt.colorbar(cf, shrink=0.5, ticks=np.arange(np.floor(varmin), np.floor(varmax)+1, .5), aspect=17)
        cf = plt.contourf(x_arr_, y_arr_, var_[:, :, k0].T,
                          cmap=cm_cont, levels=levels, alpha=1., extend='min')
    else:
        cb = plt.colorbar(cf, shrink=0.5, ticks=np.arange(np.floor(varmin), np.floor(varmax)+1, 1), aspect=17)
        cf = plt.contourf(x_arr_, y_arr_, var_[:, :, k0].T,
                          cmap=cm_cont, levels=levels, alpha=1.)
    cb.ax.tick_params(width=1, size=4)  # , labelsize=6)

    u = vel[0,ic-di:ic+di+1,jc-dj:jc+dj+1,k0].T
    v = vel[1,ic-di:ic+di+1,jc-dj:jc+dj+1,k0].T
    strm = plt.streamplot(x_arr_, y_arr_, u, v,
                          color=speed[ic - di:ic + di + 1, jc - dj:jc + dj + 1, k0].T, cmap=cm_lines,
                          linewidth=lw[:,:].T)
    cb = fig.colorbar(strm.lines, shrink=0.5, ticks=np.arange(0, np.amax(speed),1), aspect=17)
    cb.ax.tick_params(width=1, size=4)  # , labelsize=6)

    x_ticks = ax.get_xticks() * 1e-3
    y_ticks = ax.get_yticks() * 1e-3
    ax.set_xticklabels(x_ticks)
    ax.set_yticklabels(y_ticks)
    plt.xlim([x_arr[ic - di], x_arr[ic + di]])
    plt.ylim(y_arr[jc - dj], y_arr[jc + dj])
    plt.xlabel('x [km]   (dx=' + str(np.round(dx, 0)) + 'm)')
    plt.ylabel('y [km]   (dy=' + str(np.round(dy, 0)) + 'm)')
    plt.title(cont_var_name + ' with horizontal velocity (t=' + str(t0) + ', z=' + str(dz * k0) +')', fontsize=21)
    plt.savefig(os.path.join(path_out, cont_var_name + '_streamlines_xy_collision_col_t' + str(t0) + '_k' + str(k0) + '.png'))
    plt.close()
    return





def plot_streamplot_xy_varythickness(cont_var_name, var, vel, x_arr, y_arr, speed, k0, t0, path_out):
    cm = plt.cm.get_cmap('bwr')
    cm_lines = plt.cm.get_cmap('winter')

    cont_var_ = var[:, :, k0]
    if cont_var_name == 'w':
        varmax = np.amax(cont_var_)
        varmin = np.amin(cont_var_)
        varmax = np.maximum(np.abs(varmin), np.abs(varmax))
        varmin = -varmax
    elif cont_var_name == 'temperature':
        varmin = np.amin(cont_var_)
        varmax = np.amax(cont_var_)
    else:
        varmax = np.amax(cont_var_)
        varmin = np.amin(cont_var_)
    levels = np.linspace(varmin, varmax, 1e3)

    fig, ax = plt.subplots(figsize=(14,10))
    ax.set_aspect('equal')    # ax.set_aspect(1.0)
    if np.abs(speed[:,:,k0].max()) > 0.0:
        lw = 5 * speed[:,:,k0] / speed[:,:,k0].max()
    else:
        lw = 2 * np.ones(shape=speed[:,:,k0].shape)
    cf = plt.contourf(x_arr, y_arr, cont_var_.T, levels=levels, alpha=1., cmap = cm)
    if cont_var_name == 's':
        cb = plt.colorbar(cf, shrink=0.5, ticks=np.arange(np.floor(varmin), np.floor(varmax) + 1, 2), aspect=17)
    else:
        cb = plt.colorbar(cf, shrink=0.5, ticks=np.arange(np.floor(varmin), np.floor(varmax) + 1, 1), aspect=17)
    cb.ax.tick_params(width=1, size=4)  # , labelsize=6)
    strm = plt.streamplot(x_arr, y_arr, vel[0,:,:,k0].T ,vel[1,:,:,k0].T,
                          color=speed[:,:,k0].T, cmap=cm_lines, density=1.5, linewidth=lw[:,:].T)
    cb = plt.colorbar(strm.lines, shrink=0.5, ticks=np.arange(0, np.amax(speed),1), aspect=17)
    cb.ax.tick_params(width=1, size=4)  # , labelsize=6)

    x_ticks = ax.get_xticks() * 1e-3
    y_ticks = ax.get_yticks() * 1e-3
    ax.set_xticklabels(x_ticks)
    ax.set_yticklabels(y_ticks)
    plt.xlabel('x [m]   (dx='+str(dx)+')')
    plt.ylabel('y [m]   (dy='+str(dy)+')')
    plt.title('t='+str(t0) + ', z='+str(dz*k0), 'k='+str(k0), fontsize=label_size)
    plt.savefig(os.path.join(path_out, cont_var_name + '_streamlines_xy_lw_t'+str(t0)+'_k'+str(k0)+'.png'))
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

    global krange, kmax
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
        times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
                 and time_min <= np.int(name[:-3]) <= time_max]
        times.sort()
        files = [str(t) + '.nc' for t in times]
    else:  # 100_k120.nc
        times = [np.int(name[:-8]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
                 and time_min <= np.int(name[:-8]) <= time_max]
        times = times.sort()
        files = [str(t) + '.nc' for t in times]
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
    print('--- define geometry ---' + case_name)
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
    elif case_name[:21] == 'ColdPoolDry_single_3D':
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
    elif case_name[:21] == 'ColdPoolDry_double_3D':
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
    elif case_name[:21] == 'ColdPoolDry_triple_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        # zstar = nml['init']['h']
        # kstar = np.int(np.round(zstar / dx[2]))
        try:
            d = nml['init']['d']
        except:
            d = 10000
        id = np.round(d / dx)
        idhalf = np.int(np.round(id / 2))
        a = np.int(np.round(id * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        r_int = np.int(np.sqrt(3.) / 6 * id)  # radius of inscribed circle
        # point of 3-CP collision (ic, jc)
        ic = np.int(np.round(nx / 2))
        jc = np.int(np.round(ny / 2))
        ic1 = ic - r_int
        ic2 = ic1
        ic3 = ic + (a - r_int)
        jc1 = jc - idhalf
        jc2 = jc + idhalf
        jc3 = jc
        ic_arr = [ic1, ic2, ic3]
        jc_arr = [jc1, jc2, jc3]
        i0_coll = ic
        i0_center = ic_arr[0]
        j0_coll = jc
        j0_center = jc_arr[0]

        ic_arr = [ic1, ic2, ic3]
        jc_arr = [jc1, jc2, jc3]


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

    return i0_coll, j0_coll, ic_arr, jc_arr, x_half, y_half, z_half

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