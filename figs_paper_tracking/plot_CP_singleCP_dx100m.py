import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# from matplotlib.widgets import TextBox
import netCDF4 as nc
import argparse
import json as simplejson
import os

label_size = 13
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 6
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['text.usetex'] = 'true'


# tracer file:
# - col=0:      timestep of simulation
# - col=1:      age of CP (timestep since beginning of first CP)
# - col=2:      tracer id
# - col=3:      CP id
# - col=4,5:    tracer position (coordinates)
# - col=6,7:    tracer position (coordinates), rounded
# - col=8,9:    tracer radius & angle
# - col=10,11:  u, v wind components
# - col=12,13:  radial & tangential wind components
# - col=14,15:  x-, y-distance of tracer to center
# - col=16,17:  CP center (xc,yc) (COM)

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--dx")
    parser.add_argument("--shift")
    parser.add_argument("--irmax")
    parser.add_argument("--nsub")
    args = parser.parse_args()

    global dx
    if args.dx:
        dx = np.int(args.dx)
    else:
        dx = 25
    if args.irmax:
        irmax = np.round(np.double(args.irmax) / dx, 1)
    else:
        irmax = np.round(10e3 / dx, 1)

    k0 = 0
    krange = [k0]
    nk = len(krange)

    if args.nsub:
        nsub = np.int(args.nsub)
    else:
        nsub = 0

    case_name = 'ColdPoolDry_single_3D'
    path_single_dx25m = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run4_dx25m/dTh3_z1000_r1000'
    path_single_dx50m = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run3_dx50m/dTh3_z1000_r1000'
    path_single_dx100m = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/dTh3_z1000_r1000'
    if dx == 25:
        path_data = path_single_dx25m
        path_tracer_file = os.path.join(path_data, 'tracer_k' + str(k0) + '/output/', 'coldpool_tracer_out.txt')
    elif dx == 50:
        path_data = path_single_dx50m
        # path_tracer_file = os.path.join(path_data, 'tracer_k' + str(k0) + '_nointerpol/output/', 'coldpool_tracer_out.txt')
        if nsub > 0:
            path_tracer_file = os.path.join(path_data, 'tracer_k' + str(k0) + '_interpol/output/', 'coldpool_tracer_out_nsub' + str(nsub) +'.txt')
            # path_tracer_file = os.path.join(path_data, 'tracer_k' + str(k0) + '/output/', 'coldpool_tracer_out_nsub' + str(nsub) +'.txt')
        else:
            path_tracer_file = os.path.join(path_data, 'tracer_k' + str(k0) + '/output/', 'coldpool_tracer_out.txt')
    elif dx == 100:
        path_data = path_single_dx100m
        path_tracer_file = os.path.join(path_data, 'tracer_k' + str(k0) + '/output/', 'coldpool_tracer_out.txt')
        if nsub > 0:
            # path_tracer_file = os.path.join(path_data, 'tracer_k' + str(k0) + '_interpol/output/', 'coldpool_tracer_out_nsub' + str(nsub) +'.txt')
            path_tracer_file = os.path.join(path_data, 'tracer_k' + str(k0) + '/output/', 'coldpool_tracer_out_nsub' + str(nsub) +'.txt')
        else:
            path_tracer_file = os.path.join(path_data, 'tracer_k' + str(k0) + '/output/', 'coldpool_tracer_out.txt')
    path_fields = os.path.join(path_data, 'fields')
    path_fields_merged = os.path.join(path_data, 'fields_merged')
    path_out_figs = os.path.join('/nbi/home/meyerbe/paper_olga/figs_profiles_1CP_dx'+str(dx)+'m/')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print('path figures: ' + path_out_figs)
    print('path tracers: ' + path_tracer_file)
    print('')


    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gray')
    cm_grey_r = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')
    cm_winter = plt.cm.get_cmap('winter')
    cm_summer = plt.cm.get_cmap('summer')
    # cm_fall = plt.cm.get_cmap('fall')

    global tmin, tmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 0
    if args.tmax:
        tmax = np.int(args.tmax) + 100
    else:
        tmax = 3600
    times = np.arange(tmin, tmax, 100)
    nt = len(times)
    print('tmin: ', tmin)
    print('tmax: ', tmax)
    print('times: ' + str(times))

    nml = simplejson.loads(open(os.path.join(path_data, case_name + '.in')).read())
    global dt_fields
    dt_fields = nml['fields_io']['frequency']
    global nx, ny, nz
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    # dx = np.zeros(3, dtype=np.int)
    # dx[0] = nml['grid']['dx']
    # dx[1] = nml['grid']['dy']
    # dx[2] = nml['grid']['dz']
    print('')




    ''' get tracer coordinates '''
    cp_id = 1
    n_cps = get_number_cps(path_tracer_file)
    n_tracers = get_number_tracers(path_tracer_file)
    print('CP ID: ' + str(cp_id))
    print('number of CPs: ' + str(n_cps))
    print('number of tracers: ' + str(n_tracers))
    print('')
    # coordinates = get_tracer_coords(cp_id, n_cps, n_tracers, times, dt_fields, path_tracer_file)
    coordinates = get_tracer_coords(cp_id, n_cps, n_tracers, times, dt_fields, path_tracer_file)

    # a) read in CP radius r_torus(t)=r_av and CP spreading velocity u_torus(t)=U_rad_av from tracer statistics
    print('path tracers: ', path_tracer_file)
    dist_av = np.zeros((nt, nk))
    U_rad_av = np.zeros((nt, nk))
    for it, t0 in enumerate(times):
        # print('---t0: ' + str(t0) + '---', it)
        dist_av[it, k0], U_rad_av[it, k0] = get_radius_vel(path_tracer_file, it, cp_id, n_tracers, n_cps)
    r_tracers_av = dist_av * dx
    del dist_av
    print ''

    ''' get azimuthally averaged fields '''
    rad_stats_file = nc.Dataset(os.path.join(path_data, 'data_analysis', 'stats_radial_averaged.nc'))
    v_rad_av = rad_stats_file.groups['stats'].variables['v_rad'][:,:,:]
    w_av = rad_stats_file.groups['stats'].variables['w'][:,:,:]                # nt, nr, nk
    r_av = rad_stats_file.groups['stats'].variables['r'][:]                    # nr
    time_av = rad_stats_file.groups['timeseries'].variables['time'][:]         # nt
    rad_stats_file.close()

    v_rad_av_at_rim = np.zeros((nt, nk), dtype=np.double)
    ir_tracer = np.zeros(nt)
    print('v_rad_av, v_rad_av_at_rim', v_rad_av.shape, v_rad_av_at_rim.shape, nt, times)
    for it, t0 in enumerate(times):
        for k0 in krange:
            ir_tracer[it] = np.where(r_av == np.int(np.round(r_tracers_av[it, k0], -2)))[0][0]
            v_rad_av_at_rim[it, k0] = v_rad_av[it, ir_tracer[it], k0]

    v_rad_gradient = np.zeros(shape=v_rad_av.shape)
    for ir, r in enumerate(r_av[1:-1]):
        v_rad_gradient[:, ir, :] = (v_rad_av[:, ir + 1, :] - v_rad_av[:, ir - 1, :]) / (r_av[ir + 1] - r_av[ir - 1])
    ir_vrad_max = np.argmax(v_rad_av, axis=1)  # nt, nk
    ir_vrad_grad_max = np.argmin(v_rad_gradient, axis=1)  # nt, nk
    ir_w_max = np.argmax(w_av[:,20:,:], axis=1)+20

    # test plot v_rad_av vs. r_av and U_rad_av
    # plot_vel_at_rim(r_tracers_av, U_rad_av, r_av, v_rad_av[:,:,k0], w_av[:,:,k0], v_rad_av_at_rim,
    #                 time_av, k0, krange, irmax,
    #                 times, path_out_figs)
    if dx == 50:
        timerange = [200, 600, 1200, 1800, 2400, 3000]
    elif dx == 100:
        timerange = [200, 600, 1200, 1800, 2400, 3000, 3600]
    plot_vrad_w(r_tracers_av, U_rad_av, r_av, v_rad_av[:,:,k0], w_av[:,:,k0], v_rad_av_at_rim, ir_vrad_grad_max,
                    time_av, k0, irmax,
                    timerange, path_out_figs)
    # plot_vrad_w_trange(r_tracers_av, U_rad_av, r_av, v_rad_av[:, :, k0], w_av[:,:,k0], v_rad_av_at_rim,
    #                    ir_w_max, ir_vrad_grad_max,
    #             time_av, k0, irmax,
    #             timerange, path_out_figs)
    print ''





    return


# ---------------------------- PLOTTING -----------------------

def plot_vrad_w_trange(r_av, U_rad_av, radius_rad_av,
                v_rad_av, w_av, vel_at_rim, ir_w_max, ir_vrad_grad_max,
                time_rad_av, k0_tracer, irmax, trange, path_out_figs):
    # r_av[t,k0]:              average radius of tracers
    # radius_rad_av[r]:        radius from azimuthally averaged profiles

    nt = len(trange)
    fig_name = 'w_v_rad_k' + str(k0_tracer) + '_test.png'
    fig, axis = plt.subplots(2, 2, figsize=(15, 8), sharex='col')
    ax0 = axis[0,0]
    ax1 = axis[1,0]
    ax2 = axis[0,1]
    ax3 = axis[1,1]
    for i,t0 in enumerate(trange):
        it = np.int(t0/dt_fields)
        print('it, t0', i, it, t0)
        count_color = np.double(i) / (len(trange)+1)
        ir_tracer_2 = np.where(radius_rad_av == np.int(np.round(r_av[it + 2, k0_tracer], -2)))[0][0]
        ir_tracer_1 = np.where(radius_rad_av == np.int(np.round(r_av[it + 1, k0_tracer], -2)))[0][0]
        ir_tracer_0 = np.where(radius_rad_av == np.int(np.round(r_av[it, k0_tracer], -2)))[0][0]
        ir_tracer_m1 = np.where(radius_rad_av == np.int(np.round(r_av[it -1, k0_tracer], -2)))[0][0]

    #     ax0.plot([r_av[dt * it, k0_tracer], r_av[dt * it, k0_tracer]], [-10, 10], ':', color='0.25',
    #              linewidth=1)  # line through tracer
        ax0.plot(radius_rad_av, v_rad_av[it, :], color=plt.cm.coolwarm(count_color), label='t=' + str(t0) + 's')
        if i == 0:
            ax0.plot(radius_rad_av[ir_w_max[it, k0_tracer]], v_rad_av[it, ir_w_max[it, k0_tracer]], 'wd', markersize=7,
                     label='max(w)')  # tracer
            ax0.plot(radius_rad_av[ir_vrad_grad_max[it, k0_tracer]], v_rad_av[it, ir_vrad_grad_max[it, k0_tracer]],
                     'r^', label='max(grad)', markersize=5)
            ax0.plot(radius_rad_av[ir_tracer_m1], v_rad_av[it, ir_tracer_m1], 'go', markersize=5,
                     label='tracer: t-1')  # tracer
            ax0.plot(radius_rad_av[ir_tracer_0], v_rad_av[it, ir_tracer_0], 'bo', markersize=5,
                     label='tracer: t')  # tracer
            ax0.plot(radius_rad_av[ir_tracer_1], v_rad_av[it, ir_tracer_1], 'ko', markersize=5,
                     label='tracer: t+1')  # tracer
            ax0.plot(radius_rad_av[ir_tracer_2], v_rad_av[it, ir_tracer_2], 'ro', markersize=5,
                     label='tracer: t+2')  # tracer

        else:
            ax0.plot(radius_rad_av[ir_w_max[it, k0_tracer]], v_rad_av[it, ir_w_max[it, k0_tracer]], 'wd', markersize=7)
            ax0.plot(radius_rad_av[ir_vrad_grad_max[it, k0_tracer]], v_rad_av[it, ir_vrad_grad_max[it, k0_tracer]],
                     'r^', markersize=5)
            ax0.plot(radius_rad_av[ir_tracer_m1], v_rad_av[it, ir_tracer_m1], 'bo', markersize=5)  # tracer
            ax0.plot(radius_rad_av[ir_tracer_0], v_rad_av[it, ir_tracer_0], 'go', markersize=5)  # tracer
            ax0.plot(radius_rad_av[ir_tracer_1], v_rad_av[it, ir_tracer_1], 'ko', markersize=5)  # tracer
            ax0.plot(radius_rad_av[ir_tracer_2], v_rad_av[it, ir_tracer_2], 'ro', markersize=5)  # tracer
        ax0.legend(loc='best', fontsize=8)

        ax2.plot(radius_rad_av, v_rad_av[it, :], color=plt.cm.coolwarm(count_color), label='t=' + str(t0) + 's')
        if i == 0:
            ax2.plot(radius_rad_av[ir_w_max[it,k0_tracer]], v_rad_av[it, ir_w_max[it,k0_tracer]], 'wd', markersize=7, label='max(w)')  # tracer
            ax2.plot(radius_rad_av[ir_vrad_grad_max[it,k0_tracer]], v_rad_av[it, ir_vrad_grad_max[it,k0_tracer]],
                     'r^', label='max(grad)', markersize=5)
            # ax2.plot(radius_rad_av[ir_tracer_m1], v_rad_av[it, ir_tracer_m1], 'bo', markersize=5, label='tracer: t-1')  # tracer
            ax2.plot(radius_rad_av[ir_tracer_0], v_rad_av[it, ir_tracer_0], 'go', markersize=5, label='tracer: t')  # tracer
            ax2.plot(radius_rad_av[ir_tracer_1], v_rad_av[it, ir_tracer_1], 'ko', markersize=5, label='tracer: t+1')  # tracer
            # ax2.plot(radius_rad_av[ir_tracer_2], v_rad_av[it, ir_tracer_2], 'ro', markersize=5, label='tracer: t+2')  # tracer

        else:
            ax2.plot(radius_rad_av[ir_w_max[it, k0_tracer]], v_rad_av[it, ir_w_max[it, k0_tracer]], 'wd', markersize=7)
            ax2.plot(radius_rad_av[ir_vrad_grad_max[it,k0_tracer]], v_rad_av[it, ir_vrad_grad_max[it,k0_tracer]], 'r^', markersize=5)
            # ax2.plot(radius_rad_av[ir_tracer_m1], v_rad_av[it, ir_tracer_m1], 'bo', markersize=5)  # tracer
            ax2.plot(radius_rad_av[ir_tracer_0], v_rad_av[it, ir_tracer_0], 'go', markersize=5)  # tracer
            ax2.plot(radius_rad_av[ir_tracer_1], v_rad_av[it, ir_tracer_1], 'ko', markersize=5)  # tracer
            # ax2.plot(radius_rad_av[ir_tracer_2], v_rad_av[it, ir_tracer_2], 'ro', markersize=5)  # tracer
        ax2.legend(loc='best', fontsize=8)


        ax1.plot(radius_rad_av, w_av[it, :], color=plt.cm.coolwarm(count_color), label='t=' + str(t0) + 's')
        # ax1.plot([r_av[it, k0_tracer], r_av[dt * it, k0_tracer]], [-10, 10], ':', color='0.25',
        #          linewidth=1)  # line through tracer
        ax1.plot(radius_rad_av[ir_w_max[it, k0_tracer]], w_av[it, ir_w_max[it, k0_tracer]], 'wd', markersize=7)
        ax1.plot(radius_rad_av[ir_vrad_grad_max[it, k0_tracer]], w_av[it, ir_vrad_grad_max[it, k0_tracer]], 'r^',
                 markersize=5)
        ax1.plot(radius_rad_av[ir_tracer_m1], w_av[it, ir_tracer_m1], 'bo', markersize=5)  # tracer
        ax1.plot(radius_rad_av[ir_tracer_0], w_av[it, ir_tracer_0], 'go', markersize=5)  # tracer
        ax1.plot(radius_rad_av[ir_tracer_1], w_av[it, ir_tracer_1], 'ko', markersize=5)  # tracer
        ax1.plot(radius_rad_av[ir_tracer_2], w_av[it, ir_tracer_2], 'ro', markersize=5)  # tracer


        # ax3.plot([r_av[it, k0_tracer], r_av[dt * it, k0_tracer]], [-10, 10], ':', color='0.25',
        #          linewidth=1)  # line through tracer
        ax3.plot(radius_rad_av, w_av[it, :], color=plt.cm.coolwarm(count_color), label='t=' + str(t0) + 's')
        ax3.plot(radius_rad_av[ir_w_max[it, k0_tracer]], w_av[it, ir_w_max[it, k0_tracer]], 'wd', markersize=7)
        ax3.plot(radius_rad_av[ir_vrad_grad_max[it, k0_tracer]], w_av[it, ir_vrad_grad_max[it, k0_tracer]], 'r^',
                 markersize=5)
        # ax3.plot(radius_rad_av[ir_tracer_m1], w_av[it, ir_tracer_m1], 'bo', markersize=5)  # tracer
        ax3.plot(radius_rad_av[ir_tracer_0], w_av[it, ir_tracer_0], 'go', markersize=5)  # tracer
        ax3.plot(radius_rad_av[ir_tracer_1], w_av[it, ir_tracer_1], 'ko', markersize=5)  # tracer
        # ax3.plot(radius_rad_av[ir_tracer_2], w_av[it, ir_tracer_2], 'ro', markersize=5)  # tracer

        for ax in [ax1, ax3]:
            ax.plot([])

    # ax1.legend(loc='center left', bbox_to_anchor=(1., 0.5))
    # ax0.set_ylim(np.amin(v_rad_av), np.amax(v_rad_av))
    # ax1.set_ylim(-1.2, 1.5)
    # # ax1.set_ylim(np.amin(w_av), np.amax(w_av))
    for ax in axis.flat:
        ax.set_xlim(0, irmax * dx)
        ax.grid()
    #     x_ticks = [np.int(ti * 1e-3) for ti in ax.get_xticks()]
    #     ax.set_xticklabels(x_ticks)
    #     # y_ticks = [np.int(ti * 1e-3) for ti in ax.get_yticks()]
    #     print('yticks', ax.get_yticks())
    #     y_ticks = ax.get_yticks()
    #     print('      ', y_ticks)
    #     ax.set_yticklabels(y_ticks)
    # #     for label in ax.xaxis.get_ticklabels()[1::2]:
    # #         label.set_visible(False)
    # #     for label in ax.yaxis.get_ticklabels()[1::2]:
    # #         label.set_visible(False)
    # # ax0.set_ylim(np.amin(v_rad_av), np.amax(v_rad_av))
    # # ax1.set_ylim(np.amin(w_av), np.amax(w_av))
    for i in range(2):
        axis[1, i].set_xlabel('radius r / km')
    ax0.set_ylabel('radial velocity / ms' + r'$^{-1}$')
    ax2.set_ylabel('radial velocity / ms' + r'$^{-1}$')
    ax1.set_ylabel('vertical velocity / ms' + r'$^{-1}$')
    fig.subplots_adjust(top=0.95, bottom=0.15, left=0.07, right=0.87, hspace=0.1, wspace=0.25)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return




def plot_vrad_w(r_av, U_rad_av, radius_rad_av,
                    v_rad_av, w_av, vel_at_rim, ir_vrad_grad_max,
                time_rad_av, k0_tracer, irmax, trange, path_out_figs):
    # r_av[t,k0]:              average radius of tracers
    # radius_rad_av[r]:        radius from azimuthally averaged profiles

    colmap = plt.cm.coolwarm
    colmap = plt.cm.winter

    fig_name = 'w_v_rad_k' + str(k0_tracer) + '_twotimesteps.png'
    fig, axis = plt.subplots(2, 1, figsize=(12, 8), sharex='col')
    ax0 = axis[0]
    ax1 = axis[1]
    for i, t0 in enumerate(trange):
        it = np.int(t0 / dt_fields)
        print('it, t0', it, t0, tmax)
        if t0 < tmax-2*dt_fields:
            print('go', t0, tmax-2)
            count_color = np.double(i) / (len(trange))
            ir_tracer_0 = np.where(radius_rad_av == np.int(np.round(r_av[it, k0_tracer], -2)))[0][0]
            ir_tracer_1 = np.where(radius_rad_av == np.int(np.round(r_av[it + 1, k0_tracer], -2)))[0][0]
            ir_tracer_2 = np.where(radius_rad_av == np.int(np.round(r_av[it + 2, k0_tracer], -2)))[0][0]
            ir_max_grad = ir_vrad_grad_max[it, k0_tracer]

            for ax in axis:
                ax.plot([radius_rad_av[ir_tracer_0], radius_rad_av[ir_tracer_0]], [-10, 10], ':', color='0.25',
                        linewidth=1, zorder=0)  # line through tracer
                ax.plot([radius_rad_av[ir_tracer_1], radius_rad_av[ir_tracer_1]], [-10, 10], ':', color='0.25',
                        linewidth=1, zorder=0)  # line through tracer
                ax.plot([radius_rad_av[ir_tracer_2], radius_rad_av[ir_tracer_2]], [-10, 10], ':', color='0.25',
                        linewidth=1, zorder=0)  # line through tracer
                ax.plot([radius_rad_av[ir_max_grad], radius_rad_av[ir_max_grad]], [-10, 10], '--', color='0.25',
                        linewidth=1, zorder=0)  # line through max gradient

            ax0.plot(radius_rad_av, v_rad_av[it, :], color=colmap(count_color))  # , label='t='+str(t0)+'s')
            if i == 0:
                ax0.plot(radius_rad_av[ir_tracer_0], U_rad_av[it, k0_tracer], 'gd', label='U tracer: t', markersize=8)  # tracer
                ax0.plot(radius_rad_av[ir_tracer_1], U_rad_av[it+1, k0_tracer], 'kd', label='U tracer: t+1', markersize=8)  # tracer
                ax0.plot(radius_rad_av[ir_tracer_2], U_rad_av[it+2, k0_tracer], 'bd', label='U tracer: t+2', markersize=8)  # tracer
                ax0.plot(radius_rad_av[ir_max_grad], v_rad_av[it, ir_max_grad], 'rv', label='max(grad)', markersize=8)
                ax0.plot(radius_rad_av[ir_tracer_0], v_rad_av[it, ir_tracer_0], 'go', label='tracer: t')
                ax0.plot(radius_rad_av[ir_tracer_1], v_rad_av[it, ir_tracer_1], 'ko', label='tracer: t+1')
                ax0.plot(radius_rad_av[ir_tracer_2], v_rad_av[it, ir_tracer_2], 'bo', label='tracer: t+2')
            else:
                ax0.plot(radius_rad_av[ir_tracer_0], U_rad_av[it, k0_tracer], 'gd', markersize=8)  # tracer
                ax0.plot(radius_rad_av[ir_tracer_1], U_rad_av[it + 1, k0_tracer], 'kd', markersize=8)  # tracer
                ax0.plot(radius_rad_av[ir_tracer_2], U_rad_av[it + 2, k0_tracer], 'bd', markersize=8)
                ax0.plot(radius_rad_av[ir_max_grad], v_rad_av[it, ir_max_grad],'rv', markersize=8)
                ax0.plot(radius_rad_av[ir_tracer_0], v_rad_av[it, ir_tracer_0], 'go')  # tracer
                ax0.plot(radius_rad_av[ir_tracer_1], v_rad_av[it, ir_tracer_1], 'ko')  # tracer
                ax0.plot(radius_rad_av[ir_tracer_2], v_rad_av[it, ir_tracer_2], 'bo')  # tracer

            ax1.plot(radius_rad_av, w_av[it, :], color=colmap(count_color), label='t=' + str(t0) + 's')
            ax1.plot(radius_rad_av[ir_max_grad], w_av[it, ir_max_grad], 'bv', markersize=7)
            ax1.plot(radius_rad_av[ir_tracer_0], w_av[it, ir_tracer_0], 'go')  # tracer
            ax1.plot(radius_rad_av[ir_tracer_1], w_av[it, ir_tracer_1], 'ko')  # tracer
            ax1.plot(radius_rad_av[ir_tracer_2], w_av[it, ir_tracer_2], 'bo')  # tracer

    ax0.legend(loc='center left', bbox_to_anchor=(.76, 0.6))
    ax1.legend(loc='center left', bbox_to_anchor=(.76, 0.2))#, frameon=False)
    # rect = mpatches.Rectangle(((irmax*dx-2e3), 2.1), 2.e3, 6.2, fill=True, linewidth=1, edgecolor='k', facecolor='white')
    # ax0.add_patch(rect)
    ax0.set_ylim(np.amin(v_rad_av), np.amax(v_rad_av))
    ax1.set_ylim(-1.2, 1.5)
    for ax in axis.flat:
        ax.set_xlim(0, irmax * dx)
        x_ticks = [np.int(ti * 1e-3) for ti in ax.get_xticks()]
        ax.set_xticklabels(x_ticks)
        # y_ticks = [np.int(ti * 1e-3) for ti in ax.get_yticks()]
        y_ticks = ax.get_yticks()
        ax.set_yticklabels(y_ticks)
        # for label in ax.xaxis.get_ticklabels()[1::2]:
        #     label.set_visible(False)
        # for label in ax.yaxis.get_ticklabels()[1::2]:
        #     label.set_visible(False)
    ax1.set_xlabel('radius r / km')
    ax0.set_ylabel('radial velocity / ms' + r'$^{-1}$')
    ax1.set_ylabel('vertical velocity / ms' + r'$^{-1}$')
    fig.subplots_adjust(top=0.95, bottom=0.1, left=0.1, right=0.9, hspace=0.1, wspace=0.25)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)





    fig_name = 'w_v_rad_k' + str(k0_tracer) + '.png'
    fig, axis = plt.subplots(2, 1, figsize=(8, 8), sharex='col')
    ax0 = axis[0]
    ax1 = axis[1]
    for i, t0 in enumerate(trange):
        it = np.int(t0 / dt_fields)
        print('it, t0', it, t0)
        count_color = np.double(i) / (len(trange)-1)
        # ir_tracer_1 = np.where(radius_rad_av == np.int(np.round(r_av[it + 1, k0_tracer], -2)))[0][0]
        ir_tracer_0 = np.where(radius_rad_av == np.int(np.round(r_av[it, k0_tracer], -2)))[0][0]
        ir_max_grad = ir_vrad_grad_max[it, k0_tracer]

        for ax in axis:
            ax.plot([radius_rad_av[ir_tracer_0], radius_rad_av[ir_tracer_0]], [-10, 10], ':', color='0.25',
                 linewidth=1, zorder=0)  # line through tracer
            # ax.plot([radius_rad_av[ir_tracer_1], radius_rad_av[ir_tracer_1]], [-10, 10], ':', color='0.25',
            #      linewidth=1, zorder=0)  # line through tracer
            ax.plot([radius_rad_av[ir_max_grad], radius_rad_av[ir_max_grad]], [-10, 10], '--', color='0.25',
                    linewidth=1, zorder=0)  # line through max gradient

        ax0.plot(radius_rad_av, v_rad_av[it, :], color=colmap(count_color))  # , label='t='+str(t0)+'s')
        ax1.plot(radius_rad_av, w_av[it, :], color=colmap(count_color), label='t=' + str(t0) + 's')

        if i == 0:
            ax0.plot(radius_rad_av[ir_max_grad], v_rad_av[it, ir_max_grad], 'bv', label='max(grad)', markersize=7)
            ax0.plot(radius_rad_av[ir_tracer_0], v_rad_av[it, ir_tracer_0], 'ko', label='tracer')  # tracer
            # ax0.plot(radius_rad_av[ir_tracer_1], v_rad_av[it, ir_tracer_1], 'ko', label='tracer')  # tracer
            # ax0.plot(radius_rad_av[ir_tracer_1], U_rad_av[it, k0_tracer], 'ko', label='tracer')  # tracer
        else:
            ax0.plot(radius_rad_av[ir_max_grad], v_rad_av[it, ir_max_grad], 'bv', markersize=7)
            ax0.plot(radius_rad_av[ir_tracer_0], v_rad_av[it, ir_tracer_0], 'ko')  # tracer
            # ax0.plot(radius_rad_av[ir_tracer_1], v_rad_av[it, ir_tracer_1], 'ko')  # tracer

        ax1.plot(radius_rad_av[ir_max_grad], w_av[it, ir_max_grad], 'bv', markersize=7)
        ax1.plot(radius_rad_av[ir_tracer_0], w_av[it, ir_tracer_0], 'ko')  # tracer
        # ax1.plot(radius_rad_av[ir_tracer_1], w_av[it, ir_tracer_1], 'ko')  # tracer

    ax0.legend(loc='center left', bbox_to_anchor=(.78, 0.90), frameon=False)
    ax1.legend(loc='center left', bbox_to_anchor=(.78, 1.67), frameon=False)
    rect = mpatches.Rectangle((7.e3, 1.8), 1.9e3, 4.5, fill=True, linewidth=1, edgecolor='k', facecolor='white')
    ax0.add_patch(rect)

    textprops = dict(facecolor='white', alpha=0.9, linewidth=0.)
    ax0.text(4.9e2, 5.7, 'a)', fontsize=18, bbox=textprops)
    ax1.text(4.9e2, 1.13, 'b)', fontsize=18)

    ax0.set_ylim(np.amin(v_rad_av), np.amax(v_rad_av))
    # ax0.set_ylim(-0.5, 8.)
    ax1.set_ylim(-1., 1.5)
    for ax in axis.flat:
        ax.set_xlim(0, irmax * dx)
        x_ticks = [np.int(ti * 1e-3) for ti in ax.get_xticks()]
        ax.set_xticklabels(x_ticks)
        # y_ticks = [np.int(ti * 1e-3) for ti in ax.get_yticks()]
        y_ticks = ax.get_yticks()
        ax.set_yticklabels(y_ticks)
        # for label in ax.xaxis.get_ticklabels()[1::2]:
        #     label.set_visible(False)
        # for label in ax.yaxis.get_ticklabels()[1::2]:
        #     label.set_visible(False)
    ax1.set_xlabel('radius r / km')
    ax0.set_ylabel('radial velocity / ms' + r'$^{-1}$')
    ax1.set_ylabel('vertical velocity / ms' + r'$^{-1}$')
    fig.subplots_adjust(top=0.97, bottom=0.07, left=0.11, right=0.95, hspace=0.1, wspace=0.25)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return






def plot_vel_at_rim(r_av, U_rad_av, radius_rad_av,
                    v_rad_av, w_av, vel_at_rim, time_rad_av,
                    k0_tracer, krange, irmax, times, path_out_figs):
    nt = len(times)
    for k0 in krange:

        fig_name = 'v_rad_vel_k'+str(k0)+'.png'
        fig, axis = plt.subplots(1, 3, figsize=(18, 5))
        ax0 = axis[0]
        ax1 = axis[1]
        ax2 = axis[2]
        for it, t0 in enumerate(times[1::2]):
            count_color = 2 * np.double(it) / len(time_rad_av)
            ir_tracer = np.where(radius_rad_av == np.int(np.round(r_av[2 * it + 1, k0_tracer], -2)))[0][0]
            ax0.plot(radius_rad_av[:irmax], v_rad_av[2 * it + 1, :irmax], color=plt.cm.jet(count_color), label='t='+str(t0)+'s')
            ax0.plot([r_av[2 * it + 1, k0_tracer], r_av[2 * it + 1, k0_tracer]], [-10, 10],':', color='0.25', linewidth=1)  # line through tracer
            ax0.plot(radius_rad_av[ir_tracer], v_rad_av[2 * it + 1, ir_tracer], 'ko', markersize=5)     # tracer

            ax1.plot(radius_rad_av[:irmax], w_av[2 * it + 1, :irmax], color=plt.cm.jet(count_color), label='t='+str(t0)+'s')
        ax2.plot(times, vel_at_rim[:,k0], label='radial velocity at tracer position')
        ax2.plot(times, np.amax(v_rad_av[:nt,:irmax], axis=1), label='max(radial velocity)')
        ax2.plot(times, U_rad_av[:nt,k0_tracer], label='rim vel from tracer')
        # ax2.fill_between(times, 0.6*np.ones(nt), 0.4*np.ones(nt), alpha=0.2, color='0.3')
        ax0.legend(loc=1, ncol=2, fontsize=6)
        ax2.legend()
        ax0.set_ylim(np.amin(v_rad_av[:,:irmax]), np.amax(v_rad_av[:,:irmax]))
        ax2.set_xlim(times[0], times[-1])
        ax1.grid()
        ax2.grid()
        # ax0.set_title('radial velocity: v_rad')
        # ax1.set_title('radial velocity at rim')
        # ax2.set_title('(radial velocity at rim)/max(v_rad)')
        ax0.set_xlabel('radius r  [m]')
        ax1.set_xlabel('time  [s]')
        ax2.set_xlabel('time  [s]')
        # ax0.set_ylabel('radial velocity  [m/s]')
        ax0.set_xlabel('radial velocity / ms' + r'$^{-1}$')
        # ax1.set_ylabel('w    [m/s]')
        ax1.set_xlabel('vertical velocity / ms' + r'$^{-1}$')
        ax2.set_ylabel(r'v$_{rad}$    [m/s]')
        # plt.tight_layout()
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)

    return

# ---------------------------- TRACER STATISTICS -----------------------

def get_tracer_coords(cp_id, n_cps, n_tracers, times, dt_fields, fullpath_in):
    print('get tracer coordinates')
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    column = lines[0].split()

    nt = len(times)
    coords = np.zeros((nt, n_tracers, 2))
    coords_pol = np.zeros((nt, n_tracers, 2))

    for it, t0 in enumerate(times):
        print('----t0='+str(t0), it, '----')
        i = 0
        # count = t0 * n_cps * n_tracers + (cp_id - 1) * n_tracers
        # count = it * n_cps * n_tracers + (cp_id - 1) * n_tracers
        count = np.int(t0/dt_fields) * n_cps * n_tracers + (cp_id - 1) * n_tracers
        # while CP age is 0 and CP ID is cp_id
        timestep = int(lines[count].split()[0])
        cp_ID = int(lines[count].split()[3])
        while (timestep - 1 == t0/dt_fields and cp_ID == cp_id):
            columns = lines[count].split()
            coords[it,i,0] = float(columns[4])
            coords[it,i,1] = float(columns[5])
            coords_pol[it,i,0] = float(columns[8])
            coords_pol[it,i,1] = float(columns[9])
            i += 1
            count += 1
            cp_ID = int(lines[count].split()[3])
            timestep = int(lines[count].split()[0])

    f.close()
    # print ''
    return coords, coords_pol


def get_number_tracers(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    count = 0
    # while CP age is 0 and CP ID is 1
    cp_age = int(lines[count].split()[0])
    cp_ID = int(lines[count].split()[3])
    print('cp_age', cp_age)
    while (cp_age == 1 and cp_ID == 1):
        count += 1
        cp_age = int(lines[count].split()[0])
        cp_ID = int(lines[count].split()[3])
    n_tracers = count
    f.close()

    return n_tracers



def get_number_cps(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    # count = 0
    # # while CP age is 0 and CP ID is 1
    # while (int(lines[count].split()[0]) == 1):
    #     count += 1
    # cp_number = int(lines[count-1].split()[3])
    cp_number = int(lines[-1].split()[3])
    f.close()

    return cp_number



def get_radius_vel(fullpath_in, t0, cp_id, n_tracers, n_cps):
    # print('in', fullpath_in)
    f = open(fullpath_in, 'r')
    # f = open(DIR+EXPID+'/'+child+'/output/irt_tracks_output_pure_sort.txt', 'r')
    lines = f.readlines()
    count = 0
    dist = []
    vel = []

    count = t0 * n_cps * n_tracers + (cp_id - 1)*n_tracers
    # while CP age is 0 and CP ID is cp_id
    timestep = int(lines[count].split()[0])
    cp_ID = int(lines[count].split()[3])
    # print(timestep, cp_ID)
    while (timestep-1 == t0 and int(lines[count].split()[3])==cp_id):
        columns = lines[count].split()
        dist.append(float(columns[8]))
        # vel.append(np.sqrt(float(columns[10])**2 + float(columns[11])**2))
        vel.append(float(columns[12]))
        count += 1
        timestep = int(lines[count].split()[0])
    f.close()
    r_av = np.average(dist)
    vel_av = np.average(vel)

    return r_av, vel_av

# ----------------------------------------------------------------------

# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()