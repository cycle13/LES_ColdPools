import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import netCDF4 as nc
import argparse
import json as simplejson
import os

label_size = 10
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['text.usetex'] = 'true'
plt.rcParams['legend.numpoints'] = 1
plt.rcParams["legend.handlelength"] = 2.0


def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--t_shift")
    parser.add_argument("--x_shift")
    args = parser.parse_args()
    # set_input_parameters(args)

    path_olga = '/nbi/ac/conv1/henneb/results/coldpool/lindp2K_13/output/cp/'
    path_out_figs = os.path.join('/nbi/home/meyerbe/paper_olga/figs_2D_real/')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    file_olga = 'coldpool_tracer_out_all.txt'
    CP_id_olga = 3
    col_id = 4  # column in textfile of CP-ID
    # tracer file:
    # - col=0:      timestep of simulation
    # - col=1:      age of CP (timestep since beginning of first CP)
    # - col=2:      tracer id
    # - col=3:      CP id
    # - col=4,5:
    # - col=8:      tracer radius
    # - col=14,15:  CP center (xc,yc)


    path_tracer_file = os.path.join(path_olga, file_olga)

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')

    dt_fields = 300
    dx = 200
    dy = dx
    print('???? dt_fields, dx', dt_fields, dx)

    if args.t_shift:
        t_shift = np.int(args.t_shift)
    else:
        t_shift = 0
    if args.x_shift:
        x_shift = np.int(args.x_shift)
    else:
        x_shift = 0
    print('t_shift: '+str(t_shift))
    print('x_shift: '+str(x_shift))


    print('')
    print('path figs:    ' + path_out_figs)
    print('path tracers: ' + path_tracer_file)
    print('')



    ''' ----- get CP properties ----- '''
    cp_id = CP_id_olga
    print('CP ID: ' + str(cp_id))
    # n_cps = get_number_cps(path_tracer_file)
    n_cps = 2265
    print('number of CPs: ' + str(n_cps))
    # n_tracers = get_number_tracers(path_tracer_file)
    n_tracers = 999
    print('number of tracers: ' + str(n_tracers))

    lifetime, n_lines = get_cp_lifetime(cp_id, n_tracers, path_tracer_file)
    t_ini = lifetime[0]
    t_end = lifetime[-1]
    print('CP lifetime computed: ' + str(tau), lifetime)
    print('line #: ', n_lines)
    tau = 5
    t_ini = 47
    t_end = 51
    lifetime = np.arange(t_ini, t_end+1)
    print('CP lifetime: ' + str(tau), lifetime)
    xc, yc = get_cp_center(cp_id, t_ini, tau, n_tracers, n_lines, path_tracer_file)
    # yc, xc = get_cp_center(cp_id, tau, n_tracers, path_tracer_file)
    # print('CP center: (!!! switching coordinates)')
    ic = np.asarray([np.int(i) for i in xc])
    jc = np.asarray([np.int(i) for i in yc])
    print('CP center: ')
    print('xc: ', xc)
    print('yc: ', yc)
    # print('           '+ str(ic)+', '+str(jc))
    print('')


    ''' ----- tracer coordinates ----- '''
    coordinates = get_tracer_coords(cp_id, n_cps, n_tracers, lifetime, path_tracer_file)
    # var_list = ['s', 'w']
    # shift = 0
    # fig_name = 's_w'
    # plot_tracers_field(coordinates, n_tracers, shift, var_list, lifetime,
    #                    path_fields, path_out_figs, fig_name)
    print('')


    ''' averaging tracer radius'''
    aux = np.zeros(shape=coordinates.shape)
    for it in range(tau):
        aux[it,:,0] = coordinates[it,:,0] - xc[it]
        aux[it,:,1] = coordinates[it,:,1] - yc[it]
    rad_ = np.sqrt(aux[:,:,0]**2 + aux[:,:,1]**2)
    tracer_dist = np.average(rad_, axis=1)
    # del aux
    print('averaging: ', rad_.shape, tracer_dist.shape, coordinates.shape)
    print(np.amax(rad_))
    print('av tracer dist: ', tracer_dist)
    # # fig_name = 'dist_tracers.png'
    # # fig, axis = plt.subplots(1, 2, figsize=(11, 6))
    # # ax0 = axis[0]
    # # ax0.plot(rad[:])
    # # plt.tight_layout()
    # # fig.savefig(os.path.join(path_out_figs, fig_name))
    # # plt.close(fig)


    ''' ----- read in input fields ----- '''
    k0 = 0
    path_fields = path_out_figs
    rootgrp = nc.Dataset(os.path.join(path_fields, 'input_u.nc'))
    u_in = rootgrp.variables['var1'][:, k0, :, :]
    rootgrp.close()
    rootgrp = nc.Dataset(os.path.join(path_fields, 'input_v.nc'))
    v_in = rootgrp.variables['var1'][:, k0, :, :]
    rootgrp.close()
    print('input fields: shapes', u_in.shape, v_in.shape)


    ''' ----- define subdomain for which to compute radial veloctiy wrt CP center (xc,yc) ----- '''
    [nx,ny] = u_in.shape[1:3]
    lx = 10e3
    irange = np.int(np.minimum(np.int(lx / dx), np.minimum(nx - ic[0], ic[0])))
    jrange = np.int(np.minimum(np.int(lx / dx), np.minimum(ny - jc[0], jc[0])))
    nx_ = 2 * irange
    ny_ = 2 * jrange
    rmax = np.minimum(irange, jrange)
    print('subdomain: ')
    print('nx: ' + str(nx) + ', nx_: ' + str(nx_))
    print('ny: ' + str(ny) + ', ny_: ' + str(ny_))
    print('')



    '''' plot input fields & tracer coordinates '''
    # plot_input_fields(coordinates, u_in, v_in, xc, yc, irange, jrange,
    #                   cp_id, n_tracers, x_shift, t_shift, lifetime,
    #                   path_out_figs)
    # print('')



    ''' ----- compute & average radial velocity field ----- '''
    print('compute theta, r')
    # th_field, r_field = compute_radius(ic[0], jc[0], irange, jrange, nx, nx, path_out_figs)
    th_field, r_field = compute_radius(irange, jrange, irange, jrange, nx_, nx_, path_out_figs)
    print('shapes: r, th', r_field.shape, th_field.shape)

    # interpolate horizontal velocity field to cell centre in subdomain
    # v_hor_int = compute_v_hor_int(u_in, v_in, ic, jc, irange, jrange, nx, ny, t_ini, t_end, tau, path_out_figs)
    v_hor_int = compute_v_hor_int(u_in, v_in, ic[0], jc[0], irange, jrange, nx_, ny_,
                                  t_ini, t_end, tau, path_out_figs)
    # compute and average radial and tangential velocity
    print('')
    v_rad_int = np.zeros((tau, nx_, ny_))
    v_tan_int = np.zeros((tau, nx_, ny_))
    v_rad_av = np.zeros((tau, rmax))
    v_tan_av = np.zeros((tau, rmax))

    # print(r_field.shape, ic[0], ic[0]-irange, ic[0] + irange, jc[0], jc-jrange, jc[0] + jrange)
    # fig = plt.figure()
    # plt.contourf(r_field[ic[0]-irange:ic[0]+irange, jc[0]-jrange:jc[0]+jrange])
    # fig.show()
    #
    for it, t0 in enumerate(lifetime):
        v_rad_int[it,:,:], v_tan_int[it,:,:] = compute_radial_vel(v_hor_int[:,it,:,:], th_field, coordinates, n_tracers,
                                                                  irange, jrange, nx_, ny_, ic[it], jc[it],
                                                                  t0, path_out_figs)
        # average radial velocity
        v_rad_av[it, :] = compute_average_var(v_rad_int[it,:,:], rmax, r_field, nx_, ny_)
        v_tan_av[it, :] = compute_average_var(v_tan_int[it,:,:], rmax, r_field, nx_, ny_)

        print('')
        fig_name = 'test_vrad_vtan_radial_av_t'+str(t0)+'.png'
        fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(11, 4), sharey='none')
        ax0.set_title('radial velocity')
        ax0.plot(np.arange(rmax)*dx, v_rad_av[it,:rmax])
        ax0.plot([tracer_dist[it]*dx, tracer_dist[it]*dx], [np.amin(v_rad_av[it,:]), np.amax(v_rad_av[it,:])], 'k', linewidth=1)
        ax0.set_xlabel('radius / m')
        ax0.set_ylabel(r'v$_r$ / ms$^{-1}$')

        ax1.set_title('tangential velocity')
        ax1.plot(np.arange(rmax)*dx, v_tan_av[it,:rmax])
        ax1.plot([tracer_dist[it]*dx, tracer_dist[it]*dx], [np.amin(v_rad_av[it,:]), np.amax(v_rad_av[it,:])], 'k', linewidth=1)
        ax1.set_ylabel(r'v$_{tan}$ / ms$^{-1}$')
        ax2.set_title('tracer distance  / m')
        ax2.plot(lifetime, tracer_dist*dx, 'o-')
        ax2.set_xlabel('time / [-]')
        plt.tight_layout()
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)

    ''' compute average spreading velocity '''
    print('')
    print('average spreading velocity')
    # based on tracer position
    v_spread1 = dx*np.asarray([(tracer_dist[i+1]-tracer_dist[i])/dt_fields for i in range(len(tracer_dist)-1)])
    # based on position of peak in v_r
    r_vmax = np.argmax(v_rad_av[:,:rmax], axis=1)*dx
    print(r_vmax.shape, tau)
    v_spread2 = np.asarray([(r_vmax[i+1] - r_vmax[i])/dt_fields for i in range(len(tracer_dist)-1)])
    print('CP velocity')
    # print(tracer_dist)
    print('v spread1', v_spread1)
    # print(v_rad_av)
    # print(v_spread2)
    print('')



    ''' PLOTTING '''
    colmap = plt.cm.winter

    fig_name = 'vrad_vtan_radial_av.png'
    fig, ax0 = plt.subplots(1, 1, figsize=(8, 4), sharey='none')
    min = np.amin(v_rad_av)
    max = np.amax(v_rad_av)+0.2
    for it, t0 in enumerate(lifetime):
        print('plotting ', it, t0)
        count_color = np.double(it) / (len(lifetime)-1)
        ax0.plot(np.arange(rmax)*dx, v_rad_av[it, :rmax], '-x', color=colmap(count_color), label='t=' + str((t0-t_ini)*dt_fields) + 's')
        # ax0.plot([tracer_dist[it] * dx, tracer_dist[it] * dx], [min, max], 'k', linewidth=1)
        ax0.plot([tracer_dist[it] * dx, tracer_dist[it] * dx], [min, max], '-',
                 linewidth=1, color=colmap(count_color))
        ax0.plot(tracer_dist[it]*dx, v_rad_av[it, tracer_dist[it]], 'ko', markersize=6)
        if it<4:
            i = 0
            while (v_rad_av[it,i] < v_spread1[it]) and (i<rmax-1):
                i+=1
            a = i
            while (v_rad_av[it,i] > v_spread1[it]) and (i<rmax-1):
                i+=1
            b = i
            ax0.plot([0,6e3], [v_spread1[it], v_spread1[it]], '-', linewidth=1, color=colmap(count_color))
            ax0.plot([a*dx,b*dx], [v_spread1[it], v_spread1[it]], '-', linewidth=2, color=colmap(count_color))
            print('-------a, b', a, b, aux.shape, v_spread1.shape)
        # ax0.plot([r_vmax[it], r_vmax[it]],[min, max], 'k')
        ax0.plot([a*dx, a*dx], [min, max], '--', linewidth=1, color=colmap(count_color), label='a')
        ax0.plot([b*dx, b*dx], [min, max], '-.', linewidth=1, color=colmap(count_color), label='b')
        ax0.plot([6*dx, 6*dx], [min, max], 'r-', linewidth=2)
        if it == tau-1:
            ax0.plot(tracer_dist[it]*dx, v_rad_av[it, tracer_dist[it]], 'ko', markersize=6, label='tracer')
        else:
            ax0.plot(tracer_dist[it]*dx, v_rad_av[it, tracer_dist[it]], 'ko', markersize=6)

    ax0.set_xlabel('radius r / km')
    ax0.set_ylabel('radial velocity / ms' + r'$^{-1}$')

    x_ticks = [np.int(ti * 1e-3) for ti in ax0.get_xticks()]
    ax0.set_xticklabels(x_ticks)
    y_ticks = [ti for ti in ax0.get_yticks()]
    ax0.set_yticklabels(y_ticks)
    # for label in ax0.xaxis.get_ticklabels()[1::2]:
    #     label.set_visible(False)

    # rect = mpatches.Rectangle((4.8e3, 1.4), 1.e3, 1., fill=True, linewidth=1, edgecolor='k', facecolor='white', zorder=10)
    # ax0.add_patch(rect)
    ax0.legend(loc='center left', bbox_to_anchor=(.8, 0.81), frameon=True)
    textprops = dict(facecolor='white', alpha=0.9, linewidth=0.)
    ax0.text(3e2, 2.1, 'c)', fontsize=18, bbox=textprops)
    ax0.set_xlim(0, 6e3)
    ax0.set_ylim(min, max)
    fig.subplots_adjust(top=0.97, bottom=0.08, left=0.11, right=0.95, hspace=0.2, wspace=0.25)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)


    # fig_name = 'vrad_vtan_radial_av.png'
    # fig, axis = plt.subplots(2, 1, figsize=(8, 8), sharey='none')
    # ax0 = axis[0]
    # ax1 = axis[1]
    # min = np.amin(v_rad_av)
    # max = np.amax(v_rad_av)+0.2
    # for it, t0 in enumerate(lifetime):
    #     count_color = np.double(it) / (len(lifetime)-1)
    #     ax0.plot(np.arange(rmax)*dx, v_rad_av[it, :rmax], '-x', color=colmap(count_color), label='t=' + str((t0-t_ini)*dt_fields) + 's')
    #     # ax0.plot([tracer_dist[it] * dx, tracer_dist[it] * dx], [min, max], 'k', linewidth=1)
    #     ax0.plot([tracer_dist[it] * dx, tracer_dist[it] * dx], [min, max], '-',
    #              linewidth=1, color=colmap(count_color))
    #     ax0.plot(tracer_dist[it]*dx, v_rad_av[it, tracer_dist[it]], 'ko', markersize=6)
    #     if it<4:
    #         aux = np.asarray([np.int(v_rad_av[it,i]) for i in range(rmax)])
    #         # a = np.where(np.int(v_rad_av[it,:15])==np.int(v_spread1[it]))
    #         a = np.where(aux[:6]==np.int(v_spread1[it]))[0][0]
    #         b = np.where(aux[6:]==np.int(v_spread1[it]))[0][0]
    #         ax0.plot([0,6e3], [v_spread1[it], v_spread1[it]], '-', linewidth=1, color=colmap(count_color))
    #         # ax0.plot([a,b], [v_spread1[it], v_spread1[it]], '-', linewidth=2, color=colmap(count_color))
    #         print('-------a, b', a, b, aux.shape, v_spread1.shape)
    #     # ax0.plot([r_vmax[it], r_vmax[it]],[min, max], 'k')
    #     ax0.plot([a*dx, a*dx], [min, max], '--', linewidth=1, color=colmap(count_color))
    #     ax0.plot([b*dx, b*dx], [min, max], '-.', linewidth=1, color=colmap(count_color))
    #     if it == tau-1:
    #         ax0.plot(tracer_dist[it]*dx, v_rad_av[it, tracer_dist[it]], 'ko', markersize=6, label='tracer')
    #     else:
    #         ax0.plot(tracer_dist[it]*dx, v_rad_av[it, tracer_dist[it]], 'ko', markersize=6)
    # ax1.plot(lifetime, tracer_dist * dx, 'o-')
    #
    # ax0.set_xlabel('radius r / km')
    # ax0.set_ylabel('radial velocity / ms' + r'$^{-1}$')
    # ax1.set_xlabel('time / [-]')
    # ax1.set_ylabel('radius of tracer position / km')
    #
    # x_ticks = [np.int(ti * 1e-3) for ti in ax0.get_xticks()]
    # ax0.set_xticklabels(x_ticks)
    # y_ticks = [ti for ti in ax0.get_yticks()]
    # ax0.set_yticklabels(y_ticks)
    # # for label in ax0.xaxis.get_ticklabels()[1::2]:
    # #     label.set_visible(False)
    # x_ticks = [np.round(ti, 1) for ti in ax1.get_xticks()]
    # ax1.set_xticklabels(x_ticks)
    # y_ticks = [np.round(ti*1e-3,1) for ti in ax1.get_yticks()]
    # ax1.set_yticklabels(y_ticks)
    #
    # # rect = mpatches.Rectangle((4.8e3, 1.4), 1.e3, 1., fill=True, linewidth=1, edgecolor='k', facecolor='white', zorder=10)
    # # ax0.add_patch(rect)
    # ax0.legend(loc='center left', bbox_to_anchor=(.8, 0.81), frameon=True)
    # textprops = dict(facecolor='white', alpha=0.9, linewidth=0.)
    # ax0.text(3e2, 2, 'a)', fontsize=18, bbox=textprops)
    # # ax1.text(0, 1, 'b)', fontsize=18, bbox=textprops)
    # # ax1.text(1, 50, 'b)', fontsize=18, bbox=textprops)
    # # ax1.text(1, 1, 'b)', fontsize=18, bbox=textprops)
    # # ax1.text(8e2, 1, 'b)', fontsize=18, bbox=textprops)
    # # ax1.text(50, 1, 'b)', fontsize=18, bbox=textprops)
    # ax1.text(47.2, 1.95e3, 'b)', fontsize=18, bbox=textprops)
    # ax0.set_xlim(0, 6e3)
    # ax0.set_ylim(min, max)
    # fig.subplots_adjust(top=0.97, bottom=0.07, left=0.11, right=0.95, hspace=0.2, wspace=0.25)
    # fig.savefig(os.path.join(path_out_figs, fig_name))
    # plt.close(fig)





    fig_name = 'vrad_vtan_radial_av2.png'
    fig, axis = plt.subplots(2, 1, figsize=(8, 8), sharey='none')
    ax0 = axis[0]
    ax1 = axis[1]
    min_av = np.amin(v_rad_av)
    max_av = np.amax(v_rad_av)+0.2
    for it, t0 in enumerate(lifetime):
        count_color = np.double(it) / (len(lifetime)-1)
        ax0.plot(np.arange(rmax) * dx, v_rad_av[it, :rmax], color=colmap(count_color),
                 label='t=' + str((t0-t_ini)*dt_fields) + 's')
        ax0.plot([tracer_dist[it] * dx, tracer_dist[it] * dx], [min_av, max_av],
                 linewidth=1, color=colmap(count_color))
        if it == tau-1:
            ax0.plot(tracer_dist[it]*dx, v_rad_av[it, tracer_dist[it]], 'ko', markersize=6, label='tracer')
        else:
            ax0.plot(tracer_dist[it]*dx, v_rad_av[it, tracer_dist[it]], 'ko', markersize=6)

    max = np.amax(np.abs(v_hor_int[:,-1,:,:]))
    lvls = np.linspace(-max, max, 20)
    cf = ax1.contourf(v_rad_int[-1,:,:], levels=lvls, alpha=0.5, cmap=cm_bwr)
    cbar_ax1 = fig.add_axes([0.75, 0.14, 0.015, 0.3])
    ax1.quiver(v_hor_int[0,-1,:,:], v_hor_int[1,-1,:,:], pivot='tip',
                  width=0.22, units='x')
    for i in range(n_tracers):
        ax1.plot(coordinates[-1,i,0]-xc[-1]+irange+x_shift, coordinates[-1, i, 1]-yc[-1]+irange + x_shift, 'ok', markersize=2)
    cb1 = fig.colorbar(cf, cax=cbar_ax1, ticks=np.arange(np.floor(-max),np.ceil(max)+1, 1))
    cb1.set_label('radial velocity / ms' + r'$^{-1}$')

    ax0.set_xlabel('radius r / km')
    ax0.set_ylabel('radial velocity / ms' + r'$^{-1}$')
    ax1.set_xlabel('distance / km')
    ax1.set_ylabel('distance / km')

    x_ticks = [np.int(ti * 1e-3) for ti in ax0.get_xticks()]
    ax0.set_xticklabels(x_ticks)
    y_ticks = [ti for ti in ax0.get_yticks()]
    ax0.set_yticklabels(y_ticks)
    # for label in ax0.xaxis.get_ticklabels()[1::2]:
    #     label.set_visible(False)
    x_ticks = [np.round(ti*dx*1e-3, 1) for ti in ax1.get_xticks()]
    ax1.set_xticklabels(x_ticks)
    y_ticks = [np.round(ti*dx*1e-3) for ti in ax1.get_yticks()]
    ax1.set_yticklabels(y_ticks)

    ax0.legend(loc='center left', bbox_to_anchor=(.77, 0.78), frameon=False)
    rect = mpatches.Rectangle((4.5e3, 1.25), 1.2e3, 1.15, fill=True, linewidth=1, edgecolor='k', facecolor='white')
    ax0.add_patch(rect)

    textprops = dict(facecolor='white', alpha=0.9, linewidth=0.)
    # ax0.text(1e3, 7.2, 'a)', fontsize=18, bbox=textprops)
    # ax1.text(1e3, 1.13, 'b)', fontsize=18)
    ax0.text(3e2,2, 'a)', fontsize=18, bbox=textprops)
    # ax1.text(10,40, 'b)', fontsize=18, bbox=textprops)      # ok
    ax1.text(25,72, 'b)', fontsize=18, bbox=textprops)

    ax0.set_xlim(0, 6e3)
    ax0.set_ylim(min_av, max_av)
    ax1.set_xlim(20,80)
    ax1.set_ylim(20,80)
    ax1.set_aspect('equal')
    fig.subplots_adjust(top=0.98, bottom=0.07, left=0.11, right=0.95, hspace=0.2, wspace=0.1)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    #
    # #
    # #
    # # var_list = ['v_rad', 'v_tan']
    # # rootgrp = nc.Dataset(os.path.join(path, 'fields_v_rad', 'v_rad.nc'))
    # # for it, t0 in enumerate(times):
    # #     print('-plot time: ' + str(t0))
    # #     fig_name = 'v_rad_tan' + '_t' + str(t0) + '_tracers.png'
    # #     fig, axis = plt.subplots(1, 2, figsize=(11, 6), sharey='all')
    # #     for j, var_name in enumerate(var_list):
    # #         print(var_name, j)
    # #         var = rootgrp.variables[var_name][it, :, :, k0]
    # #         max = np.amax(var)
    # #         if var_name in ['w', 'v_rad', 'v_tan']:
    # #             min = -max
    # #         else:
    # #             min = np.amin(var)
    # #         axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_bwr)
    # #         axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_bwr)
    # #         axis[j].set_title(var_name)
    # #         axis[j].set_aspect('equal')
    # #     for i in range(n_tracers):
    # #         for j in range(len(var_list)):
    # #             axis[j].plot(coordinates[it,i,0]+shift, coordinates[it,i,1]+shift, 'ok', markersize=2)
    # #     plt.tight_layout()
    # #     fig.savefig(os.path.join(path_out_figs, fig_name))
    # #     plt.close(fig)
    # # rootgrp.close()
    # #
    # #
    # #
    # # var_list = ['u', 'v']
    # # rootgrp = nc.Dataset(os.path.join(path_tracers, 'input', 'uv_alltimes.nc'))
    # # for it, t0 in enumerate(times):
    # #     print('-plot time: ' + str(t0))
    # #     fig_name = 'uv_alltimes' + '_t' + str(t0) + '_tracers.png'
    # #     fig, axis = plt.subplots(1, 2, figsize=(11, 6), sharey='all')
    # #     for j, var_name in enumerate(var_list):
    # #         print(var_name, j)
    # #         var = rootgrp.variables[var_name][it, :, :]
    # #         max = np.amax(var)
    # #         if var_name in ['w', 'v_rad', 'v_tan']:
    # #             min = -max
    # #         else:
    # #             min = np.amin(var)
    # #         axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_bwr)
    # #         axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_bwr)
    # #         axis[j].set_title(var_name)
    # #         axis[j].set_aspect('equal')
    # #     for i in range(n_tracers):
    # #         for j in range(len(var_list)):
    # #             axis[j].plot(coordinates[it,i,0]+shift, coordinates[it,i,1]+shift, 'ok', markersize=2)
    # #     plt.tight_layout()
    # #     fig.savefig(os.path.join(path_out_figs, fig_name))
    # #     plt.close(fig)
    # # rootgrp.close()


    return
# ----------------------------------------------------------------------





def plot_input_fields(coordinates, var1, var2, xc, yc, irange, jrange,
                      cp_id, n_tracers, x_shift, t_shift, lifetime,
                      path_out_figs):
    # k0 = 0
    # rootgrp = nc.Dataset(os.path.join(path_fields, 'input_u.nc'))
    # var1 = rootgrp.variables['var1'][:, k0, :, :]
    # rootgrp.close()
    # rootgrp = nc.Dataset(os.path.join(path_fields, 'input_v.nc'))
    # var2 = rootgrp.variables['var1'][:, k0, :, :]
    # rootgrp.close()
    print('plot input fields: ')

    imin = xc[0] - irange
    imax = xc[0] + irange
    jmin = yc[0] - jrange
    jmax = yc[0] + jrange
    imin = 580
    imax = 720
    jmin = 530
    jmax = 670


    # max = np.maximum(np.amax(var1[:,:,:]), np.amax(var2[:,:,:]))
    for it, t0 in enumerate(lifetime):
        print('-plot time: ' + str(t0), it)
        fig_name = 'uv_input_fields' + '_t' + str(t0) + '_tracers.png'
        fig, axis = plt.subplots(3, 2, figsize=(10, 10))
        max = 4.5
        lvls = np.linspace(-max, max, 20)
        for ax0 in axis[:2,0].flat:
            cf = ax0.contourf(var1[t0,:,:], levels=lvls, cmap=cm_bwr)
            ax0.set_title('u')
            plt.colorbar(cf, ax=ax0)
        ax0 = axis[2,0]
        cf = ax0.contourf(var1[t0,:,:].T, levels=lvls, cmap=cm_bwr)
        plt.colorbar(cf, ax=ax0)
        for ax1 in axis[:2,1].flat:
            cf = ax1.contourf(var2[t0,:,:], levels=lvls, cmap=cm_bwr)
            ax1.set_title('v')
            plt.colorbar(cf, ax=ax1)
        ax1 = axis[2,1]
        cf = ax1.contourf(var2[t0,:,:].T, levels=lvls, cmap=cm_bwr)
        plt.colorbar(cf, ax=ax1)

        for i in range(n_tracers):
            for ax in axis.flat:
                ax.plot(coordinates[it, i, 0] + x_shift, coordinates[it, i, 1] + x_shift, 'ok', markersize=2)
        for ax in axis[0,:].flat:
            # ax.set_aspect('equal')
            # ax.set_xlim(np.amin(coordinates[it, :, 0] - 1), np.amax(coordinates[it, :, 0] + 1))
            # ax.set_ylim(np.amin(coordinates[it, :, 1] - 1), np.amax(coordinates[it, :, 1] + 1))
            ax.set_xlim(imin, imax)
            ax.set_ylim(jmin, jmax)
            ax.plot(xc[it],yc[it], 'xk', markersize=20)
            ax.plot([xc[it]-irange,xc[it]-irange],[jmin,jmax],'-k')
            ax.plot([xc[it]+irange,xc[it]+irange],[jmin,jmax],'-k')
            ax.plot([imin,imax],[yc[it]-jrange,yc[it]-jrange],'-k')
            ax.plot([imin,imax],[yc[it]+jrange,yc[it]+jrange],'-k')
        for ax in axis[1:,:].flat:
            # ax.set_aspect('equal')
            # ax.set_xlim(np.amin(coordinates[it, :, 0] - 1), np.amax(coordinates[it, :, 0] + 1))
            # ax.set_ylim(np.amin(coordinates[it, :, 1] - 1), np.amax(coordinates[it, :, 1] + 1))
            ax.set_xlim(630, 670)
            ax.set_ylim(580, 610)
        plt.suptitle('t='+str(t0)+'s, t-shift '+str(t_shift)+', x-shift '+str(x_shift))
        #     plt.tight_layout()
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)

    return


def plot_tracers_field(coordinates, n_tracers, shift, var_list, times,
                       path_fields, path_out_figs, fig_name_prefix):
    # var_list = ['s', 'w']
    for it,t0 in enumerate(times):
        print('-plot time: '+str(t0))
        # print(os.path.join(path_fields, str(t0)+'.nc'))
        fig_name = fig_name_prefix + '_t' + str(t0) + '_tracers.png'
        fig, axis = plt.subplots(1, 2, figsize=(11, 6), sharey='all')
        rootgrp = nc.Dataset(os.path.join(path_fields, str(t0)+'.nc'))
        for j, var_name in enumerate(var_list):
            print(var_name, j)
            var = rootgrp.groups['fields'].variables[var_name][:,:,k0]
            max = np.amax(var)
            if var_name in ['w', 'v_rad', 'v_tan']:
                min = -max
                cm_ = cm_bwr
            else:
                min = np.amin(var)
                cm_ = cm_hsv
            # axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_)
            axis[j].imshow(var.T, vmin=min, vmax=max, cmap=cm_, origin='lower')
            axis[j].set_title(var_name)
            axis[j].set_aspect('equal')
        rootgrp.close()
        for i in range(n_tracers):
            for j in range(len(var_list)):
                axis[j].plot(coordinates[it,i,0]+shift, coordinates[it,i,1]+shift, 'ok', markersize=2)
        plt.tight_layout()
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)
    return

# ------------------------- CP STATISTICS ---------------------------------------------


def get_tracer_coords(cp_id, n_cps, n_tracers, times, fullpath_in):
    print('get tracer coordinates')
    f = open(fullpath_in, 'r')
    lines = f.readlines()

    nt = len(times)
    coords = np.zeros((nt, n_tracers, 2))

    count_start = 0
    for it,t0 in enumerate(times):
        print('----t0='+str(t0), it, '----', count_start)
        timestep = int(lines[count_start].split()[0])
        age = int(lines[count_start].split()[1])
        cp_ID = int(lines[count_start].split()[3])
        while (cp_ID < cp_id):
            count_start += 1
            cp_ID = int(lines[count_start].split()[3])
        # print('t0 - cp_ID:', t0, timestep, cp_ID, count_start)
        while (timestep < t0):
            count_start += 1
            timestep = int(lines[count_start].split()[0])
        age = int(lines[count_start].split()[1])
        cp_ID = int(lines[count_start].split()[3])
        # print('t_end - timestep: ', t_ini, timestep, cp_ID, count_start)

        count = count_start
        i = 0
        while (timestep == t0 and age == it and cp_ID == cp_id):
            columns = lines[count].split()
            coords[it,i,0] = float(columns[4])
            coords[it,i,1] = float(columns[5])
            i += 1
            count += 1
            timestep = int(lines[count].split()[0])
            age = int(lines[count].split()[1])
            cp_ID = int(lines[count].split()[3])

    f.close()
    print('')
    return coords


def get_number_cps(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    cp_number = int(lines[-1].split()[3])
    f.close()

    return cp_number


def get_number_tracers(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    count = 0
    cp_age = int(lines[count].split()[0])
    cp_ID = int(lines[count].split()[3])
    # print('cp_age', cp_age)
    # print('cp_ID', cp_ID)
    age = cp_age
    ID = cp_ID
    while (age == cp_age and ID == cp_ID):
        count += 1
        age = int(lines[count].split()[0])
        ID = int(lines[count].split()[3])
    n_tracers = count
    f.close()

    return n_tracers


def get_cp_lifetime(cp_ID, n_tracers, fullpath_in):
    print('get cp lifetime')
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    n_lines = []
    times = []
    count = 0
    ID = int(lines[count].split()[3])
    while (ID < cp_ID):
        count += 1
        ID = int(lines[count].split()[3])
    t0 = int(lines[count].split()[0])
    ID = int(lines[count].split()[3])
    # print(t0, ID, count)

    t1 = t0
    while (ID == cp_ID) and (count<len(lines)):
        id_tr = int(lines[count].split()[2])
        t1 = int(lines[count].split()[0])
        # print('loop', count, ID, id_tr, t1)
        if id_tr == 1:
            # print('append')
            n_lines.append(count)
            times.append(t1)
        elif id_tr > 1 and t1 == times[-1]:
            while id_tr > 1 and t1 == times[-1]:
                count += 1
                id_tr = int(lines[count].split()[2])
                t1 = int(lines[count].split()[0])
            if int(lines[count].split()[3]) == cp_ID:
                n_lines.append(count)
                times.append(t1)
                # print('diff append', count)
        # print(times)
        if t1 > times[-1]:
            break
        count += n_tracers
        ID = int(lines[count].split()[3])
        t2 = t1
        while (ID != cp_ID) and (t2 <= t1 + 1):
            count += n_tracers
            ID = int(lines[count].split()[3])
            t2 = int(lines[count].split()[0])
    # tau = t1-t0+1
    print('lifetime: ', len(times), t0, t1)
    print(times)

    f.close()
    return times, n_lines


def get_cp_center(cp_ID, t_ini, tau, n_tracers, n_lines, fullpath_in):
    xc = np.zeros(tau)
    yc = np.zeros(tau)

    f = open(fullpath_in, 'r')
    lines = f.readlines()
    # count = 0
    # ID = int(lines[count].split()[3])
    # while (ID < cp_ID):
    #     count += 1
    #     ID = int(lines[count].split()[3])
    # count_start = count
    #
    # for it in range(tau):
    #     xc[it] = float(lines[count_start+it*n_tracers].split()[14])
    #     yc[it] = float(lines[count_start+it*n_tracers].split()[15])

    for it in range(tau):
        count = n_lines[it]
        ID = int(lines[count].split()[3])
        t_ = int(lines[count].split()[0])
        print('cp center: ', it, count, ID, cp_ID, t_, it+t_ini)
        if ID == cp_ID and t_ == (it+t_ini):
            xc[it] = float(lines[count].split()[14])
            yc[it] = float(lines[count].split()[15])
        else:
            print('ERROR ERROR ERROR in reading out CP centre')
            print(ID, cp_ID, t_, it, count)

    f.close()
    # print('xc: ', xc)
    # print('yc: ', yc)
    return xc, yc


# def get_radius_vel(fullpath_in, t0, cp_id, n_tracers, n_cps):
#     f = open(fullpath_in, 'r')
#     # f = open(DIR+EXPID+'/'+child+'/output/irt_tracks_output_pure_sort.txt', 'r')
#     lines = f.readlines()
#     count = 0
#     dist = []
#     vel = []
#
#     count = t0 * n_cps * n_tracers + (cp_id - 1)*n_tracers
#     # while CP age is 0 and CP ID is cp_id
#     timestep = int(lines[count].split()[0])
#     cp_ID = int(lines[count].split()[3])
#     while (timestep-1 == t0 and int(lines[count].split()[3])==cp_id):
#         columns = lines[count].split()
#         dist.append(float(columns[8]))
#         # vel.append(np.sqrt(float(columns[10])**2 + float(columns[11])**2))
#         vel.append(float(columns[12]))
#         count += 1
#         timestep = int(lines[count].split()[0])
#     f.close()
#     r_av = np.average(dist)
#     vel_av = np.average(vel)
#
#     return r_av, vel_av


# ----------------------------------------------------------------------
def compute_radius(ic, jc, irange, jrange, nx, ny, path_out_figs):
    r_field = np.zeros((nx, ny), dtype=np.int)          # radius
    th_field = np.zeros((nx, ny), dtype=np.double)      # angle
    for i in range(irange):
        for j in range(jrange):
            r_field[ic+i, jc+j] = np.round(np.sqrt(i**2+j**2))
            r_field[ic-i, jc+j] = r_field[ic+i,jc+j]
            r_field[ic-i, jc-j] = r_field[ic+i,jc+j]
            r_field[ic+i, jc-j] = r_field[ic+i,jc+j]
            if i == 0:
                # i = 1e-9
                aux = np.arctan(np.double(j)/1e-9)
            else:
                aux = np.arctan(np.double(j)/i)
            # print('i,j', i, j, aux, np.pi-aux, np.pi + aux, 2*np.pi - aux)
            th_field[ic+i, jc+j] = aux
            th_field[ic-i, jc+j] = np.pi - aux
            th_field[ic-i, jc-j] = np.pi + aux
            th_field[ic+i, jc-j] = 2*np.pi - aux

    # fig, axis = plt.subplots(3, 2, figsize=(10,10))
    # # cf = ax1.imshow(r_field[ic - irange:ic + irange + 1, jc - jrange:jc + jrange + 1].T, origin='lower')
    # for ax in axis[:,0].flat:
    #     cf = ax.imshow(r_field.T, origin='lower')
    #     plt.colorbar(cf, ax=ax)
    #     ax.set_title('r(x,y)')
    # for ax in axis[:,1].flat:
    #     cf = ax.imshow(th_field.T, origin='lower')
    #     # cf = ax2.imshow(th_field[ic - irange:ic + irange + 1, jc - jrange:jc + jrange + 1].T, origin='lower')
    #     plt.colorbar(cf, ax=ax)
    #     ax.set_title('th(x,y)')
    #
    # for ax in axis[1:,:].flat:
    #     ax.plot([ic,ic], [0, ny], 'k', linewidth=1)
    #     ax.plot([0, nx], [jc,jc], 'k', linewidth=1)
    # # ax4.plot([ic,ic], [0, ny], 'k')
    # # ax4.plot([0, nx], [jc,jc], 'k')
    #
    # for ax in axis[:2,:].flat:
    #     ax.set_xlim(0, nx)
    #     ax.set_ylim(0, ny)
    # for ax in axis[2,:].flat:
    #     ax.set_xlim(ic-irange-5, ic+irange+5)
    #     ax.set_ylim(jc-jrange-5, jc+jrange+5)
    #
    # plt.suptitle('ic,jc='+str(ic) + ', '+str(jc) + ', (nxh = '+str(irange)+')')
    # plt.savefig(os.path.join(path_out_figs, 'r_th_field.png'))
    # plt.close()
    return th_field, r_field


def compute_v_hor_int(u_in, v_in, xc, yc, irange, jrange, nx_, ny_,
                      t_ini, t_end, tau, path_out_figs):
    print('')
    print('compute v hor')
    v_hor_int = np.zeros((2, tau, nx_, ny_))
    print('v_hor_int: ', v_hor_int.shape)
    # imin = xc[0] - irange
    # imax = xc[0] + irange
    # jmin = yc[0] - jrange
    # jmax = yc[0] + jrange
    #
    # for the setting, where v_hor_int only of size of subdomain & neglecting
    #        time variation in centre coordinates (xc, yc)
    imin = xc-irange
    imax = xc+irange
    jmin = yc-jrange
    jmax = yc+jrange
    [nt, nx, ny] = np.shape(u_in)
    print('imin: ', imin, xc, imax, irange)
    print('jmin: ', jmin, yc, jmax, jrange)

    # for i in range(irange):
    #     v_hor_int[0, :, xc[0] + i, jmin:jmax] = 0.5 * (
    #                 u_in[t_ini:t_end + 1, xc[0] + i, jmin:jmax] + u_in[t_ini:t_end + 1, xc[0] + i - 1, jmin:jmax])
    #     v_hor_int[0, :, xc[0] - i, jmin:jmax] = 0.5 * (
    #                 u_in[t_ini:t_end + 1, xc[0] - i, jmin:jmax] + u_in[t_ini:t_end + 1, xc[0] - i - 1, jmin:jmax])
    # for j in range(jrange):
    #     v_hor_int[1, :, imin:imax, yc[0] + j] = 0.5 * (
    #                 v_in[t_ini:t_end + 1, imin:imax, yc[0] + j] + v_in[t_ini:t_end + 1, imin:imax, yc[0] + j - 1])
    #     v_hor_int[1, :, imin:imax, yc[0] - j] = 0.5 * (
    #                 v_in[t_ini:t_end + 1, imin:imax, yc[0] - j] + v_in[t_ini:t_end + 1, imin:imax, yc[0] - j - 1])

    # for i in range(1, nx_ - 1):
    #     v_hor_int[0,:,i,:] = 0.5 * (u_in[t_ini:t_end+1,i,:] + u_in[t_ini:t_end+1,i-1,:])
    # for j in range(1, ny_ - 1):
    #     v_hor_int[1,:,:,j] = 0.5 * (v_in[t_ini:t_end+1,:,j] + v_in[t_ini:t_end+1,:,j-1])

    for it,t0 in enumerate(np.arange(t_ini, t_end+1)):
        for j in range(ny_):
            for i in range(nx_):
                v_hor_int[0, it, j, i] = 0.5 * ( u_in[t0, jmin+j, imin+i] + u_in[t0, jmin+j-1, imin+i] )
                v_hor_int[1, it, j, i] = 0.5 * ( v_in[t0, jmin+j, imin+i] + v_in[t0, jmin+j, imin+i-1] )

        print('it, t, i, imin+i', it, t0, i, imin, imin+i, j, jmin, jmin+j)


    # # plotting
    # for it, t0 in enumerate(np.arange(t_ini, t_end+1)):
    #     print('plot t='+str(t0), it)
    #     max_in = np.amax(u_in[t0, :, :])
    #     lvls_in = np.linspace(-max_in, max_in, 10)
    #
    #     fig, axis = plt.subplots(4, 2, figsize=(8, 12))
    #     ax = axis[0, 0]
    #     ax.set_title('u in')
    #     cf = ax.contourf(u_in[t0, :, :], cmap=cm_bwr, levels=lvls_in)
    #     plt.colorbar(cf, ax=ax, shrink=0.8)
    #     ax.plot(xc, yc, 'xk', markersize=20)
    #     ax.plot([xc+irange, xc+irange], [0, ny], 'k')
    #     ax.plot([xc-irange, xc-irange], [0, ny], 'k')
    #     ax = axis[0, 1]
    #     ax.set_title('v in')
    #     cf = ax.contourf(v_in[t0, :, :], cmap=cm_bwr, levels=lvls_in)
    #     plt.colorbar(cf, ax=ax, shrink=0.8)
    #     ax.plot(xc, yc, 'xk', markersize=20)
    #     ax = axis[1, 0]
    #     ax.set_title('u hor int')
    #     cf = ax.contourf(v_hor_int[0, it, :, :], cmap=cm_bwr, levels=lvls_in)
    #     plt.colorbar(cf, ax=ax, shrink=0.8)
    #     ax.plot(irange, jrange, 'xk', markersize=20)
    #     ax = axis[1, 1]
    #     ax.set_title('v hor int')
    #     cf = ax.contourf(v_hor_int[1, it, :, :], cmap=cm_bwr, levels=lvls_in)
    #     plt.colorbar(cf, ax=ax, shrink=0.8)
    #     ax.plot(irange, jrange, 'xk', markersize=20)
    #
    #     ax = axis[2, 0]
    #     ax.set_title('u hor int')
    #     cf = ax.contourf(v_hor_int[0, it, :, :], cmap=cm_bwr, levels=lvls_in)
    #     plt.colorbar(cf, ax=ax, shrink=0.8)
    #     ax.plot(irange, jrange, 'xk', markersize=20)
    #     ax = axis[2, 1]
    #     ax.set_title('v hor int')
    #     cf = ax.contourf(v_hor_int[1, it, :, :], cmap=cm_bwr, levels=lvls_in)
    #     plt.colorbar(cf, ax=ax, shrink=0.8)
    #     ax.plot(irange, jrange, 'xk', markersize=20)
    #     ax = axis[3, 0]
    #     ax.set_title('u in')
    #     cf = ax.contourf(u_in[t0, :, :], cmap=cm_bwr, levels=lvls_in)
    #     plt.colorbar(cf, ax=ax, shrink=0.8)
    #     ax.plot(xc, yc, 'xk', markersize=20)
    #     ax = axis[3, 1]
    #     ax.set_title('v in')
    #     cf = ax.contourf(v_in[t0, :, :], cmap=cm_bwr, levels=lvls_in)
    #     plt.colorbar(cf, ax=ax, shrink=0.8)
    #     ax.plot(xc, yc, 'xk', markersize=20)
    #
    #     for ax in axis[0,:].flat:
    #         ax.set_xlim(imin-5, imin+i+5)
    #         ax.set_ylim(jmin-5, jmin+j+5)
    #     for ax in axis[1,:].flat:
    #         ax.set_xlim(-5, nx_+5)
    #         ax.set_ylim(-5, ny_+5)
    #     for ax in axis[2,:].flat:
    #         ax.set_xlim(irange-(xc-625), irange+(665-xc))
    #         ax.set_ylim(jrange-(yc-570), jrange+(610-yc))
    #     for ax in axis[3,:].flat:
    #         ax.set_xlim(625, 665)
    #         ax.set_ylim(570, 610)
    #     # plt.tight_layout()
    #     plt.suptitle('t='+str(t0))
    #     plt.savefig(os.path.join(path_out_figs, 'test_field_vhor_int_t' + str(t0) + '.png'))
    #     plt.close()
    return v_hor_int


def compute_radial_vel(uv, th_field,
                       coordinates, n_tracers,
                       irange, jrange, nx, ny, ic, jc, t0, path_out_figs):
    x_shift_coords = -2
    ur = np.zeros((nx,ny), dtype=np.double)
    utan = np.zeros((nx,ny), dtype=np.double)
    print('compute radial velocity: ')
    print(uv.shape, th_field.shape)

    # for th of size of subdomain & adjusting for switched indices in velocity field
    for j in range(ny):
        for i in range(nx):
            th = th_field[i, j]
            # counter-clockwise rotation
            ur[j,i] = uv[0,j,i] * np.cos(th) + uv[1,j,i] * np.sin(th)
            utan[j,i] = -uv[0,j,i] * np.sin(th) + uv[1,j,i] * np.cos(th)
    # for th of size of large domain
    # for i in range(nx):
    #     for j in range(ny):
    #         ii = ic-np.int(nx/2)+i
    #         jj = jc-np.int(ny/2)+j
    #         th = th_field[ii,jj]
    #         # # clockwise rotation
    #         # ur[i,j] = uv[0,i,j]*np.cos(th) + uv[1,i,j]*np.sin(th)
    #         # counter-clockwise rotation
    #         ur[i,j] = uv[0,ii,jj]*np.cos(th) + uv[1,ii,jj]*np.sin(th)
    #         utan[i,j] = -uv[0,ii,jj]*np.sin(th) + uv[1,ii,jj]*np.cos(th)



    # fig, axis = plt.subplots(1, 3, figsize=(12, 5))
    # ax = axis[0]
    # ax.set_title('th(x,y)')
    # cf = ax.imshow(th_field.T, origin='lower')
    # plt.colorbar(cf, ax=ax, shrink=0.5)
    # # ax.plot([ic, ic], [0, ny], 'k', linewidth=1)
    # # ax.plot([0, nx], [jc, jc], 'k', linewidth=1)
    # ax = axis[1]
    # ax.set_title('cos')
    # cf = ax.imshow(np.cos(th_field), origin='lower')
    # plt.colorbar(cf, ax=ax, shrink=0.5)
    # ax = axis[2]
    # ax.set_title('sin')
    # cf = ax.imshow(np.sin(th_field), origin='lower')
    # plt.colorbar(cf, ax=ax, shrink=0.5)
    # plt.tight_layout()
    # plt.suptitle('ic,jc=' + str(ic) + ', ' + str(jc) + ', (nxh = ' + str(irange) + ')')
    # plt.savefig(os.path.join(path_out_figs, 'th_cos_sin_field.png'))
    # plt.close()


    max = np.maximum(np.amax(np.abs(uv)), np.maximum(np.amax(ur), np.amax(utan)))
    lvls = np.linspace(-max, max, 20)

    fig, axis = plt.subplots(2, 2, figsize=(11, 11))
    ax = axis[0,0]
    ax.set_title('u')
    cf = ax.contourf(uv[0, :, :], levels=lvls, alpha=0.5, cmap=cm_bwr)
    plt.colorbar(cf, ax=ax, shrink=0.8)
    ax = axis[1,0]
    ax.set_title('v')
    cf = ax.contourf(uv[1, :, :], levels=lvls, alpha=0.5, cmap=cm_bwr)
    plt.colorbar(cf, ax=ax, shrink=0.8)
    for ax in axis.flat:
        ax.quiver(uv[0, :, :], uv[1, :, :], pivot='tip',
                  width=0.22, units='x')

    ax = axis[0,1]
    ax.set_title('radial velocity')
    cf = ax.contourf(ur[:, :], levels=lvls, alpha=0.5, cmap=cm_bwr)
    plt.colorbar(cf, ax=ax, shrink=0.8)
    ax = axis[1,1]
    ax.set_title('tangential velocity')
    cf = ax.contourf(utan[:, :], levels=lvls, alpha=0.5, cmap=cm_bwr)
    plt.colorbar(cf, ax=ax, shrink=0.8)

    for ax in axis.flat:
        for i in range(n_tracers):
            ax.plot(coordinates[t0-47, i, 0]-ic+irange+x_shift_coords,
                     coordinates[t0-47, i, 1]-jc+irange+x_shift_coords, 'ok', markersize=2)

    for ax in axis.flat:
        ax.plot(irange, jrange, 'xk', markersize=20)
        ax.set_aspect('equal')
        ax.set_xlim(25, 70)
        ax.set_ylim(30, 75)
    plt.tight_layout()
    plt.savefig(os.path.join(path_out_figs, 'test_field_vrad_vtan_t' + str(t0) + '.png'))
    plt.close()


    # fig, axis = plt.subplots(2, 4, figsize=(30, 11))
    # ax01 = axis[0, 0]
    # ax02 = axis[1, 0]
    # ax11 = axis[0, 1]
    # ax12 = axis[1, 1]
    # ax21 = axis[0, 2]
    # ax22 = axis[1, 2]
    # ax31 = axis[0, 3]
    # ax32 = axis[1, 3]
    #
    # cf = ax01.contourf(uv[0, :, :], levels=lvls, alpha=0.3, cmap=cm_bwr)
    # plt.colorbar(cf, ax=ax01, shrink=0.8)
    # cf = ax02.contourf(uv[1, :, :], levels=lvls, alpha=0.3, cmap=cm_bwr)
    # plt.colorbar(cf, ax=ax02, shrink=0.8)
    # for ax in axis[:,0].flat:
    #     ax.quiver(uv[0,:,:], uv[1, :,:], pivot='tip',
    #               width=0.22, units='x')
    #
    # cf = ax11.contourf(uv[0, :, :], levels=lvls, cmap=cm_bwr)
    # plt.colorbar(cf, ax=ax11, shrink=0.8)
    # ax11.set_title('u')
    # cf = ax12.contourf(uv[1, :, :], levels=lvls, cmap=cm_bwr)
    # plt.colorbar(cf, ax=ax12, shrink=0.8)
    # ax12.set_title('v')
    #
    # circle1 = plt.Circle((ic, jc), rmax / 2, fill=False, color='k', linewidth=2)
    # ax21.set_title('radial velocity')
    # cf = ax21.contourf(ur[:, :], levels=lvls, cmap=cm_bwr)
    # plt.colorbar(cf, ax=ax21, shrink=0.8)
    # ax21.add_artist(circle1)
    #
    # ax22.set_title('tangential velocity')
    # # ax11.plot(xc[0] - 0.5, jc - 0.5, 'ow', markersize=7)
    # cf = ax22.contourf(utan[:, :], levels=lvls, cmap=cm_bwr)
    # plt.colorbar(cf, ax=ax22, shrink=0.8)
    #
    # ax31.set_title('radial velocity')
    # cf = ax31.contourf(ur[:, :], levels=lvls, alpha=0.5, cmap=cm_bwr)
    # plt.colorbar(cf, ax=ax31, shrink=0.8)
    # ax32.set_title('tangential velocity')
    # cf = ax32.contourf(utan[:, :], levels=lvls, alpha=0.5, cmap=cm_bwr)
    # plt.colorbar(cf, ax=ax32, shrink=0.8)
    # for ax in axis[:,3].flat:
    #     ax.quiver(uv[0, :, :], uv[1, :, :], pivot='tip',
    #               width=0.22, units='x')
    #
    # for ax in axis.flat:
    #     ax.plot(irange, jrange, 'xk', markersize=20)
    #     ax.set_aspect('equal')
    #     ax.set_xlim(25, 70)
    #     ax.set_ylim(30, 75)
    # plt.tight_layout()
    # plt.savefig(os.path.join(path_out_figs, 'test_field_uv_vrad_vtan_t' + str(t0) + '.png'))
    # plt.close()
    print('')
    return ur, utan


def compute_average_var(var, rmax, r_field, nx, ny):
    count = np.zeros(rmax, dtype=np.int)

    var_av = np.zeros((rmax), dtype=np.double)
    for i in range(nx):
        for j in range(ny):
            r = r_field[i, j]
            if r < rmax:
                count[r] += 1
                var_av[r] += var[i, j]
    # print('ij', i, j, nx, ny)
    # print('r', r, rmax)
    # print('r field', r_field[ic, :])

    for r in range(rmax):
        if count[r] > 0:
            var_av[r] /= count[r]

    return var_av
# ----------------------------------------------------------------------

def set_input_parameters(args):
    print('--- set input parameters ---')
    global path, path_fields, case_name

    return

# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()
