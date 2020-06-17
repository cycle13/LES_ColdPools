import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import netCDF4 as nc
import argparse
import json as simplejson
import os
import time

from thermodynamic_functions import thetas_c


global tick_size, label_size
tick_size = 24
label_size = 32
plt.rcParams['xtick.labelsize'] = tick_size
plt.rcParams['ytick.labelsize'] = tick_size
plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['legend.fontsize'] = 10
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
    parser.add_argument("--times", nargs='+', type=int)
    parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    parser.add_argument("--imin")
    parser.add_argument("--imax")
    args = parser.parse_args()

    global cm_bwr
    cm_bwr = plt.cm.get_cmap('bwr')


    files, times, nml = set_input_parameters(args)
    # configuration (defined in define_geometry)
    # ic1 = ic - r_int; ic2 = ic1; ic3 = ic + (a - r_int)
    # jc1 = jc - idhalf; jc2 = jc + idhalf; jc3 = jc
    print('')
    ID = os.path.split(path_in)[1]
    if ID == '':
        ID = os.path.basename(path_in[:-1])
    print ('id: ', ID)
    i0_coll, j0_coll, ic_arr, jc_arr, d, x_half, y_half, z_half = define_geometry(case_name, nml, files)

    if args.times:
        time_range = [np.int(t) for t in args.times]
    else:
        if d == 12000:
            time_range = [600, 1200, 1500, 1800, 2000]
        elif d == 15000:
            time_range = [600, 1100, 1200, 1300, 1500]
        # time_range = [600, 1100]
    dj = (jc_arr[1] - jc_arr[0])  # side length of triangle
    h = (ic_arr[2] - ic_arr[0])  # height of triangle
    dplot = 1.7 * dj
    di = (dplot - h) / 2
    print('distance: ' + str(d))
    print('geometry: x=' + str(ic_arr) + ', y=' + str(jc_arr))
    print('plotting: dplot', dplot)
    k0 = 0
    print('k0=' + str(k0))

    imin = np.int(i0_coll - 1./3*h - di)
    imax = np.int(i0_coll + 2./3*h + di)
    jmin = np.int(j0_coll - .5*dplot)
    jmax = np.int(j0_coll + .5*dplot)
    # x_arr = x_half[i0_coll - di:i0_coll + di + 1]
    # y_arr = y_half[j0_coll - dj:j0_coll + dj + 1]
    x_arr = x_half[:(imax-imin)]
    y_arr = y_half[:(jmax-jmin)]
    print('')



    print('')
    ''' --- call plotting functions --- '''
    print('')
    print('path Figures: ' + pathout_figs)

    cm_cont = plt.cm.get_cmap('gnuplot')
    cm_lines = plt.cm.get_cmap('winter')
    cbar_pos_x = .94
    cbar_width = .012
    cbar_height = .3

    fig_name = 'triple_casestudy.png'
    nrow = 2
    ncol = len(time_range)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.9 * ncol, 5 * nrow), sharex='all', sharey='all')
    for it, t0 in enumerate(time_range):
        ax0 = axes[0, it]
        ax1 = axes[1, it]
        file = str(t0)+'.nc'
        print('--- time: ', t0, '('+str(it)+') ---')
        vel_h = np.ndarray((2, imax-imin, jmax-jmin, kmax+1))
        print(os.path.join(path_fields, file))

        print('')
        print('imin, etc', imin, imax, jmin, jmax)
        print('')

        root = nc.Dataset(os.path.join(path_fields, file), 'r')
        vel_h[0, :, :, :] = root.groups['fields'].variables['u'][imin:imax, jmin:jmax, :kmax+1]
        vel_h[1, :, :, :] = root.groups['fields'].variables['v'][imin:imax, jmin:jmax, :kmax+1]
        w = root.groups['fields'].variables['w'][imin:imax, jmin:jmax, k0]
        w_up = root.groups['fields'].variables['w'][imin:imax, jmin:jmax, 5]
        s = root.groups['fields'].variables['s'][imin:imax, jmin:jmax, :kmax+1]
        theta = thetas_c(s, 0.0)
        temp = root.groups['fields'].variables['temperature'][imin:imax, jmin:jmax, :kmax+1]
        root.close()

        speed_h = np.sqrt(vel_h[0, :] * vel_h[0, :] + vel_h[1, :] * vel_h[1, :])
        # t_mean = np.average(np.average(temp[:, :, :], axis=0), axis=0)
        # t_anomaly = temp[:, :, :] - t_mean
        if np.abs(speed_h[:,:, k0].max()) > 0.0:
            lw = 5 * speed_h[:,:,k0] / speed_h[:,:, k0].max()
        else:
            lw = 5 * np.ones(shape=speed_h[:,:,k0].shape)

        cm_cont_w = plt.cm.get_cmap('bwr')
        wmax = np.maximum(np.abs(np.amin(w)), np.abs(np.amax(w)))
        wmin = -wmax
        levels_w = np.linspace(wmin, wmax, 1e2)

        cm_cont_temp = plt.cm.get_cmap('gnuplot')
        tempmin = 298.
        tempmax = 300.
        levels_temp = np.linspace(tempmin, tempmax, 1e2)

        cf0 = ax0.contourf(x_arr, y_arr, w[:, :].T, cmap=cm_cont_w, levels=levels_w)#, extend='min')
        # cf1 = ax1.contourf(x_arr, y_arr, temp[:, :, k0].T, cmap=cm_cont_temp, levels=levels_temp, extend='min')
        cf1 = ax1.contourf(x_arr, y_arr, theta[:, :, k0].T, cmap=cm_cont_temp, levels=levels_temp, extend='min')

        # ax0.plot([x_arr[ic_arr[0]-imin], x_arr[ic_arr[1]-imin], x_arr[ic_arr[2]-imin]],
        #          [y_arr[jc_arr[0]-jmin], y_arr[jc_arr[1]-jmin], y_arr[jc_arr[2]-jmin]], 'o', markersize=12, color='k')
        # ax0.plot(x_arr[i0_coll-imin], y_arr[j0_coll-jmin], 'o', markersize=12, color='r')

        if it == ncol-1:
            cax = plt.axes([cbar_pos_x, 0.55, cbar_width, cbar_height])
            cb = plt.colorbar(cf0, cax=cax, ticks=np.arange(np.floor(wmin), np.floor(wmax)+1, 1))
            cb.set_label(r'w [m/s]', rotation=90)
            cb.ax.tick_params(width=1, size=4)  # , labelsize=6)
            cax = plt.axes([cbar_pos_x, 0.1, cbar_width, cbar_height])
            cb = plt.colorbar(cf1, cax=cax, ticks=np.arange(np.floor(tempmin), np.floor(tempmax)+.5, 1.))
            # text = cax.yaxis.label
            # font = plt.font_manager.FontProperties(family='times new roman', style='italic', size=16)
            # text.set_font_properties(font)
            cb.set_label(r'pot. temp. [K]', rotation=90)
            cb.ax.tick_params(width=1, size=4)  # , labelsize=6)

        # ax0.streamplot(x_arr, y_arr, vel_h[0,:,:,k0].T, vel_h[1,:,:,k0].T,
        #                 color='k', density=1., linewidth=lw[:, :].T)
        # ax0.contour(x_arr, y_arr, w_up[:, :].T, levels=np.arange(.5,5.1,.5))#, colors='k')
        ax0.contour(x_arr, y_arr, w_up[:, :].T, levels=np.arange(.5,.9,.5))#, colors='k')
        ax1.streamplot(x_arr, y_arr, vel_h[0, :, :, k0].T, vel_h[1, :, :, k0].T,
                       color='k', density=1., linewidth=lw[:, :].T)
        # ax1.streamplot(x_arr, y_arr, vel_h[0,:,:,k0].T, vel_h[1,:,:,k0].T,
        #                 color=speed_h[:,:,k0].T, cmap=cm_lines,
        #                 linewidth=lw[:, :].T)
        # cb = fig.colorbar(strm.lines, shrink=0.5, ticks=np.arange(0, np.amax(speed), 1), aspect=17)
        # cb.ax.tick_params(width=1, size=4)  # , labelsize=6)

        textstr = 't=' + str(t0) + 's'
        props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='white')
        ax0.text(0.1, 1.13, textstr, transform=ax0.transAxes, fontsize=24,
                verticalalignment='top', bbox=props)

    for ax in axes.flat:
        ax.set_aspect('equal')
        ax.set_xlim(x_arr[0], x_arr[imax-imin-1])
        ax.set_ylim(y_arr[0], y_arr[jmax-jmin-1])
    for ax in axes[-1,:].flat:
        x_ticks = ax.get_xticks() * 1e-3
        ax.set_xticklabels(x_ticks)
    for ax in axes[:,0].flat:
        y_ticks = ax.get_yticks() * 1e-3
        ax.set_yticklabels(y_ticks)
        ax.set_ylabel('')
    # lbls = ['a)']
    # for lbl,ax in enumerate(axes.flat):
    #     textstr = 't=' + str(t0) + 's'
    #     props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='white')
    #     ax0.text(0.1, 1.13, textstr, transform=ax0.transAxes, fontsize=24,
    #              verticalalignment='top', bbox=props)

    plt.subplots_adjust(bottom=0.05, right=.93, left=0.03, top=0.92, wspace=0.08, hspace=0.08)
    print('saving: ', os.path.join(pathout_figs, fig_name))
    fig.savefig(os.path.join(pathout_figs, fig_name))
    plt.close(fig)


    # ''' (a) xy-plane '''
    # for k0 in krange:
    #     print('-- k0=' + str(k0) + ' --')
    #     plot_streamplot_xy_collision(cont_var_name, cont_var[:,:,:kmax+1], vel_h, speed_h,
    #                                      x_half, y_half, ic_arr, jc_arr,
    #                                      i0_coll, di, j0_coll, dj, k0, t0, path_out)
    #     # plot_streamplot_xy_varythickness(cont_var_name, cont_var, vel_h,
    #     #                                  x_half, y_half, speed_h, k0, t0, path_out)



    return






def plot_streamplot_xy_collision(cont_var_name, cont_var, vel, speed,
                                 x_arr, y_arr, ic_arr, jc_arr, ic, di, jc, dj, k0, t0, path_out_figs):
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
    plt.savefig(os.path.join(path_out_figs, cont_var_name + '_streamlines_xy_collision_lw_t' + str(t0) + '_k' + str(k0) + '.png'))
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
    plt.savefig(os.path.join(path_out_figs, cont_var_name + '_streamlines_xy_collision_col_t' + str(t0) + '_k' + str(k0) + '.png'))
    plt.close()
    return





def plot_streamplot_xy_varythickness(cont_var_name, var, vel, x_arr, y_arr, speed, k0, t0, pathout_figs):
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
    plt.savefig(os.path.join(pathout_figs, cont_var_name + '_streamlines_xy_lw_t' + str(t0) + '_k' + str(k0) + '.png'))
    plt.close()
    return




# ----------------------------------
def set_input_parameters(args):
    global path_in, pathout_figs, path_fields
    path_in = args.path
    if os.path.exists(os.path.join(path_in, 'fields')):
        path_fields = os.path.join(path_in, 'fields')
    elif os.path.exists(os.path.join(path_in, 'fields_k120')):
        path_fields = os.path.join(path_in, 'fields_k120')
    pathout_figs = os.path.join(path_in, 'streamlines')
    if not os.path.exists(pathout_figs):
        os.mkdir(pathout_figs)
    print('path in: ' + path_in)
    print('path figures: ' + pathout_figs)
    print('')

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
    # if args.tmin:
    #     time_min = np.int(args.tmin)
    # else:
    #     try:
    #         time_min = np.int(files[0][:-3])
    #     except:
    #         time_min = 100
    # if args.tmax:
    #     time_max = np.int(args.tmax)
    # else:
    #     try:
    #         time_max = np.int(files[-1][:-3])
    #     except:
    #         time_max = 3000

    global krange, kmax, imin, imax
    if args.kmin:
        kmin = np.int(args.kmin)
    else:
        kmin = 1
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = 1
    krange = np.arange(kmin, kmax+1)
    if args.imin:
        imin = np.int(args.imin)
    else:
        imin = 50
    if args.imax:
        imax = np.int(args.imax)
    else:
        imax = nx - imin

    ''' determine file range '''
    if len(files[0]) < 7:  # 100.nc
        times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
                 #and time_min <= np.int(name[:-3]) <= time_max
                ]
        times.sort()
        files = [str(t) + '.nc' for t in times]
    else:  # 100_k120.nc
        times = [np.int(name[:-8]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
                 #and time_min <= np.int(name[:-8]) <= time_max
                ]
        times = times.sort()
        files = [str(t) + '.nc' for t in times]
        for it,t0 in enumerate(times):
            files[it] = str(t0)+files[0][3:]

    print('')
    print('files', files)
    # print('len', len(files[0]))
    # print('times', times)
    print('krange', krange)
    print('')

    return files, times, nml



def define_geometry(case_name, nml, files):
    print('--- define geometry: '  + case_name + '---')
    global rstar, irstar
    if case_name[:21] == 'ColdPoolDry_triple_3D':
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
            d = 1e5
            print('d not defined')
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
        print('ic_arr: ' + str(ic_arr))
        print('jc_arr: ' + str(jc_arr))
        print('collision: ' + str(ic) + ', ' + str(jc))
        print('d: ' + str(d))
        print('')


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

    return i0_coll, j0_coll, ic_arr, jc_arr, d, x_half, y_half, z_half

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
