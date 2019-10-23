import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os
import scipy
from scipy import optimize

label_size = 10
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.labelsize'] = 15
# plt.rcParams['text.usetex'] = 'true'


def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path_root")
    parser.add_argument("dTh", type=int)
    parser.add_argument("--zparams", nargs='+', type=int)
    parser.add_argument('--rparams', nargs='+', type=int)
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--k0")
    args = parser.parse_args()

    # level of tracers
    if args.k0:
        k0 = np.int(args.k0)
    else:
        k0 = 0
    dTh, z_params, r_params = set_input_parameters(args)
    n_params = len(r_params)

    global cm_bwr, cm_grey, cm_vir, cm_hsv, cm_grey2
    cm_bwr = plt.cm.get_cmap('seismic')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')
    cm_grey2 = plt.cm.get_cmap('bone_r')
    cm_fall = plt.cm.get_cmap('winter')
    cm_summer = plt.cm.get_cmap('spring')

    # # test colorlist
    colorlist_all = ['darkred', 'maroon', 'r', 'tomato', 'indianred', 'orange', 'gold',
                 'limegreen', 'forestgreen', 'g', 'darkgreen', 'seagreen', 'lightseagreen', 'darkcyan',
                 'mediumblue', 'darkblue', 'midnightblue', 'navy']
    colorlist = ['maroon', 'indianred', 'orange',
                     'limegreen', 'darkgreen', 'darkcyan', 'lightseagreen',
                     'mediumblue', 'navy']
    colorlist5 = ['maroon', 'indianred', 'orange', 'darkcyan', 'navy']
    # colorlist5 = ['orange', 'indianred', 'maroon', 'navy', 'lightseagreen']
    colorlist4 = ['indianred', 'orange', 'darkcyan', 'navy']
    # not in list: siena,
    fig, (ax0, ax1, ax2) = plt.subplots(1,3)
    x_range = np.arange(10)
    for i in range(len(colorlist_all)):
        ax0.plot(x_range, i + x_range, color=colorlist_all[i])
    for i in range(len(colorlist)):
        ax1.plot(x_range, i + x_range, color=colorlist[i])
    for i in range(5):
        ax2.plot(x_range, i + x_range, color=colorlist5[i])
    plt.savefig(os.path.join(path_out_figs, 'test_colors.png'))
    plt.close()
    # --------------------------------------



    # reference case: dTh3_z1000_r1000
    rstar_ref = 1000
    zstar_ref = 1000
    id_ref = 'dTh3_z'+str(zstar_ref) + '_r' + str(rstar_ref)
    # path_ref = os.path.join(path_root, id_ref)
    path_ref = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/' + id_ref
    # path_ref = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run3_dx50m/' + id_ref
    dt_fields = 100
    cp_id = 2  # circle ID that is used for statistics

    id0 = 'dTh' + str(dTh) + '_z' + str(z_params[0]) + '_r' + str(r_params[0])
    fullpath_in = os.path.join(path_root, id0, 'tracer_k' + str(k0), 'output')
    n_tracers = get_number_tracers(fullpath_in)
    n_cps = get_number_cps(fullpath_in)
    print('number of CPs: ', n_cps)
    print('number of tracers per CP: ', n_tracers)
    print ''
    nt = len(times)
    krange = [k0]
    nk = len(krange)

    # --------------------------------------
    ''' (a) read in data from tracer output (text-file)'''
    dist_av = np.zeros((n_params, nt, nk))
    r_av = np.zeros((n_params, nt, nk))         # absolute radius
    r_av_abs = np.zeros((n_params, nt, nk))     # radius minus initial radius (r(t) - r(t=0))
    drdt_av = np.zeros((n_params, nt, nk))
    U_rad_av = np.zeros((n_params, nt, nk))

    print('--- reading in tracer data')
    for istar in range(n_params):
        zstar = z_params[0]
        rstar = r_params[istar]
        if rstar == 1000:
            id = id_ref
            fullpath_in = os.path.join(path_ref, 'tracer_k' + str(k0), 'output')
        else:
            id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
            fullpath_in = os.path.join(path_root, id, 'tracer_k'+str(k0), 'output')
        print('id', id)
        # read_in_txtfile(fullpath_in)
        for it, t0 in enumerate(times):
            print('---t0: '+str(t0)+'---', it)
            dist_av[istar, it, k0], U_rad_av[istar, it, k0] = get_radius_vel(fullpath_in, it, cp_id, n_tracers, n_cps)
        r_av = dist_av * dx[0]
        r_av_abs[istar,:,:] = r_av[istar,:,:] - rstar
        r_av_abs[istar,:,:] = r_av[istar,:,:] - r_av[istar,0,:]
        print('rstar: ', dx[0], rstar, r_av[istar,0,k0], r_av_abs[istar,0,k0])
        for it, t0 in enumerate(times[1:]):
            drdt_av[:,it,:] = 1./dt_fields * (r_av[:,it,:] - r_av[:,it-1,:])
        print ''
    print ''

    # --------------------------------------
    ''' (b) plot r_av, dtdt_av, U_rad_av'''
    figname = 'CP_rim_dTh' + str(dTh) + '.png'
    title = 'CP rim (dTh=' + str(dTh) + 'K, dx='+str(dx[0]) +'m)'
    plot_dist_vel(r_av, drdt_av, U_rad_av,
                  [dTh], z_params, r_params, n_params, k0, id_ref, title, colorlist5, figname)
    ''' r_av_abs'''
    figname = 'CP_rim_dTh' + str(dTh) + '_r_abs.png'
    title = 'CP rim (dTh=' + str(dTh) + 'K, dx=' + str(dx[0]) + 'm)'
    plot_dist_vel(r_av_abs, drdt_av, U_rad_av,
                  [dTh], z_params, r_params, n_params, k0, id_ref, title, colorlist5, figname)


    ''' (c) fitting '''
    # func(x, a, b, c) = a + b * x ** c
    # func2(x, a, b, c) = a * (x / b) ** c  # a = y(x=b); a=R0=R(t=t0), t0: start time for fitting
    # func_log1(x, a, b, c) = a + b * np.log(x)
    # func_log2(x, a, b, c) = a + b * np.log(c*x)
    # func_log_romps_r(x, a, b, c) = a + b * np.log(1 + x * c)
    # # ''' fit function to r_av '''
    # # tmin = 4        # fitting start
    # # figname = 'CP_rim.png'
    # # p1_r, p_log1_r, p_log2_r, p_log_romps_r = plot_dist_fitting(r_av,
    # #                 [dTh], z_params, r_params, n_params, tmin, k0, id_ref, colorlist5, figname)
    # # plot_parameters(p1_r, p_log1_r, p_log2_r, p_log_romps_r, 'r')
    # ''' fit function to r_av_abs '''
    # figname = 'CP_rim_abs.png'
    tmin_r = 4
    # tmin = tmin_r
    # p1_r, p_log1_r, p_log2_r, p_log_romps_r = plot_dist_fitting(r_av_abs,
    #                   [dTh], z_params, r_params, n_params, tmin, k0, id_ref, colorlist5, figname)
    # plot_parameters(p1_r, p_log1_r, p_log2_r, p_log_romps_r, 'r_abs')
    # # ''' fit function to U_rad '''
    # # tmin = 5
    # # figname = 'CP_rim_fit_vrad'
    # # p0 = [0, 1, -1]
    # # p1_vel, p_log1_vel, p_log2_vel, p_log_romps_vel = plot_vel_fitting(r_av, U_rad_av, p0,
    # #                 [dTh], z_params, r_params, n_params, tmin, k0, id_ref, colorlist5, figname)
    # # plot_parameters(p1_vel, p_log1_vel, p_log2_vel, p_log_romps_vel, 'U_rad')
    # ''' fit function to dr/dt '''
    tmin_drdt = 7
    # tmin = tmin_drdt
    # figname = 'CP_rim_fit_dRdt'
    # p1_vel, p_log1_vel, p_log2_vel, p_log_romps_vel = plot_vel_fitting(r_av, drdt_av, [0,1,-1],
    #                  [dTh], z_params, r_params, n_params, tmin, k0, id_ref, colorlist5, figname)
    # plot_parameters(p1_vel, p_log1_vel, p_log2_vel, p_log_romps_vel, 'dRdt')


    # ''' (d) error '''
    # fig_name = 'Fit_comparison_r.png'
    # test_fit('R', r_av[:,:,k0], p1_r, p_log2_r, p_log_romps_r,
    #          n_params, [dTh], r_params, z_params, id_ref, tmin, colorlist5, fig_name)
    # fig_name = 'Fit_comparison_drdt.png'
    # test_fit('U', drdt_av[:, :, k0], p1_vel, p_log2_vel, p_log_romps_vel,
    #          n_params, [dTh], r_params, z_params, id_ref, tmin, colorlist5, fig_name)




    # # ''' (d) plot normalized radius / velocity'''
    # # # (i) for vertical velocity from crosssection in 3D field
    # # # (ii) for azimuthally averaged vertical vleocity
    # # print('plotting normalized')
    # # for istar in range(n_params):
    # #     zstar = z_params[0]
    # #     rstar = r_params[istar]
    # #     id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
    # #     print('id', id)
    # #     nml = simplejson.loads(open(os.path.join(path_root, id, case_name + '.in')).read())
    # #     nx = nml['grid']['nx']
    # #     ny = nml['grid']['ny']
    # #     ic = np.int(nx / 2)
    # #     jc = np.int(ny / 2)
    # #     print(ic, jc)
    # #     figname_norm = 'CP_rim_normalized_' + id + '.png'
    # #     plot_vel_normalized(r_av, times, istar, k0, nx, ic, jc, id, figname_norm)
    # #     # figname_norm = 'CP_rim_normalized_' + id + '_av.png'
    # #     # plot_vel_normalized_w_av(r_av, times, istar, k0, id, figname_norm)
    #
    #
    # ''' (e) scaling with PE'''
    # trange = [600, 1200, 1800, 2400, 3000, 3500 ]
    # scaling = [-1, 0, 1, 2, 3]
    # fig_name = 'PE_scaling_dTh' + str(dTh) + '.png'
    # plot_PE_scaling(r_av, U_rad_av, scaling[:n_params+1], k0, trange, nt, n_params, fig_name, colorlist)
    # # fig_name = 'PE_scaling_dTh' + str(dTh) + '_log2.png'
    # # plot_PE_scaling_U(r_av, U_rad_av, drdt_av, scaling[:n_params+1], k0, trange, nt, n_params, fig_name)

    print('')
    print('!!!!!!!!!!!!!!!!!!!!!')
    print('')
    ''' compute slope m '''
    m_r = np.ndarray(len(times))
    m_r_ = np.ndarray((n_params, len(times)))
    m_U = np.ndarray(len(times))
    m_U_ = np.ndarray(len(times))
    # for it, t0 in enumerate(times):
    #     m_r_[2:, it] = np.log(r_av[2:, it, k0] / R0[it]) / scaling[2:]
    #     m_U_[2:, it] = np.log(drdt_av[2:, it, k0] / r0[it]) / scaling[2:]
    #     # m10[it] = np.mean(np.log(r_av[:,it,k0]/r0[it])/scaling[:])
    # m_r[:] = np.mean(m_r_[2:, :], axis=0)
    # m_r = m_r * np.log(2) / np.log(10)
    #
    # fig, (ax0, ax1, ax2) = plt.subplots(1, 3, sharex='none', figsize=(18, 5))
    # ax0.plot(times, m10_r, '-o')
    # ax1.plot(times, m_r, '-o')
    # ax0.set_title('m10')
    # ax1.set_title('m')
    # for ax in (ax0, ax1):
    #     ax.grid(True)
    # fig.tight_layout()
    # fig.savefig(os.path.join(path_out_figs, 'test_slope.png'))
    # plt.close(fig)



    ''' scaling of R'''
    tmin = tmin_r
    fig_name = 'test_R.png'
    R0 = r_av[0, tmin, k0]
    dR_log = np.log(r_av[:, tmin, k0]) - np.log(R0)
    dR = R0/r_av[:, tmin, k0]
    r_av_ens_mean = np.average(r_av[:,:,k0], axis=0)
    dR_ens_mean = R0/r_av_ens_mean[tmin]
    t0 = np.double(tmin*dt_fields)

    fig, axes = plt.subplots(6, 4, sharex='none', figsize=(18, 20))
    [ax0, ax1, ax2, ax3] = [axes[0, i] for i in range(4)]
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar])
        ax0.plot([t0, t0], [0, np.amax(r_av)], 'k', linewidth=1)
        ax2.plot([np.log(times[1]), np.log(times[1])], [np.log(np.amin(r_av)), np.log(np.amax(r_av))], 'k',
                 linewidth=1)
        ax0.plot(times[1:-1], r_av[istar, 1:-1, k0], 'o-', color=colorlist5[istar], label=lbl)
        ax1.loglog(times[1:-1], r_av[istar, 1:-1, k0], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[1:-1]), np.log(r_av[istar, 1:-1, k0]), 'o-', color=colorlist5[istar], label=lbl)
    ax3.plot(np.arange(-1, 4), 1. / dR[:], '-ok', markersize=4)
    ax0.plot(times[1:-1], r_av_ens_mean[1:-1], 'b')
    ax0.set_xlabel('time')
    ax1.set_xlabel('time')
    ax2.set_xlabel('log(t)')
    ax0.set_ylabel('time')
    ax0.set_ylabel('R')
    ax1.set_ylabel('R')
    ax2.set_ylabel('log(R)')
    ax3.set_xlabel('log2(PE/PE_ref)')
    ax3.set_ylabel('R0')

    [ax0, ax1, ax2, ax3] = [axes[1, i] for i in range(4)]
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar])
        ax0.plot([0, 0], [0, np.amax(r_av)], 'k', linewidth=1)
        ax0.plot(times[4:-1] - t0, r_av[istar, 4:-1, k0], 'o-', color=colorlist5[istar], label=lbl)
        ax1.loglog(times[4:-1] - t0, r_av[istar, 4:-1, k0], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[4:-1] - t0), np.log(r_av[istar, 4:-1, k0]), 'o-',
                 color=colorlist5[istar], label=lbl)
    ax0.plot(times[4:-1]-t0, r_av_ens_mean[4:-1], 'b')
    ax0.set_xlabel('t-t0')
    ax1.set_xlabel('t-t0')
    ax2.set_xlabel('log(t-t0)')
    ax0.set_ylabel('R')
    ax1.set_ylabel('R')
    ax2.set_ylabel('log(R)')
    ax3.set_ylabel('log(R)')

    [ax0, ax1, ax2, ax3] = [axes[2, i] for i in range(4)]
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar])
        ax0.plot([0, 0], [0, np.amax(r_av)], 'k', linewidth=1)
        ax0.plot(times[4:-1]/t0, r_av[istar, 4:-1, k0], 'o-', color=colorlist5[istar], label=lbl)
        ax1.loglog(times[4:-1]/t0, r_av[istar, 4:-1, k0], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[4:-1])-np.log(t0), np.log(r_av[istar, 4:-1, k0]), 'o-',
                 color=colorlist5[istar], label=lbl)
        ax3.plot(np.log(times[4:-1]/t0), np.log(r_av[istar, 4:-1, k0]), 'o-',
                 color=colorlist5[istar], label=lbl)
    ax0.plot(times[4:-1]/t0, r_av_ens_mean[4:-1], 'b')
    ax0.set_xlabel('t/t0')
    ax1.set_xlabel('t/t0')
    ax2.set_xlabel('log(t)-log(t0)')
    ax3.set_xlabel('log(t/t0)')
    ax0.set_ylabel('R')
    ax1.set_ylabel('R')
    ax2.set_ylabel('log(R)')
    ax3.set_ylabel('log(R)')

    [ax0, ax1, ax2, ax3] = [axes[3, i] for i in range(4)]
    m = 0.55
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar])
        ax0.plot([0, 0], [0, np.amax(r_av)], 'k', linewidth=1)
        ax0.plot(times[4:-1] - t0, r_av[istar, 4:-1, k0]*dR[istar], 'o-',
                 color=colorlist5[istar], label=lbl)
        ax1.loglog(times[4:-1] - t0, r_av[istar, 4:-1, k0]*dR[istar], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[4:-1]-t0), np.log(r_av[istar, 4:-1, k0]) - dR_log[istar], 'o-',
                 color=colorlist5[istar], label=lbl)
        ax3.plot(np.log(times[4:-1]-t0), m*(np.log(r_av[istar, 4:-1, k0]) - dR_log[istar]), 'o-',
                 color=colorlist5[istar], label=lbl)

        ax2.plot(np.log(times[4:-1]-t0),
                 np.log(R0) + m* (np.log(times[4:-1]) - np.log(t0)), '-r', label='m=0.54')

        ax0.plot(times[4:tmin] - tmin * dt_fields, r_av[istar, 4:tmin, k0] - dR_log[istar], 'ow')
        ax1.loglog(times[4:tmin] - tmin * dt_fields, r_av[istar, 4:tmin, k0], 'ow')
        ax2.plot(np.log(times[4:tmin]) - np.log(tmin * dt_fields), np.log(r_av[istar, 4:tmin, k0]) - dR_log[istar],
                 'ow')
    ax0.plot(times[4:-1]-t0, r_av_ens_mean[4:-1]*dR_ens_mean, 'b')
    ax0.set_xlabel('t-t0')
    ax1.set_xlabel('t-t0')
    ax2.set_xlabel('log(t-t0)')
    ax0.set_ylabel('R * R0(t0)/R(t0)')
    ax1.set_ylabel('R * R0(t0)/R(t0)')
    ax2.set_ylabel(r'log(R)-dlog(R)')
    ax3.set_ylabel(r'log(R)-dlog(R)')
    # ax3.set_ylabel('log(R)-[log(R(t0)-log(R0(t0)]')
    # ax3.set_xlim(np.log(times[tmin]), np.log(times[-2]))
    # ax2.set_xlim(0, 2.3)
    ax3.set_xlim(0, 2.3)
    ax0.set_ylim(0, 8e3)

    [ax0, ax1, ax2, ax3] = [axes[4, i] for i in range(4)]
    m = 0.55
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar])
        ax0.plot([0, 0], [0, np.amax(r_av)], 'k', linewidth=1)
        ax0.plot(times[4:-1]/t0, r_av[istar, 4:-1, k0] * dR[istar], 'o-',
                 color=colorlist5[istar], label=lbl)
        ax1.loglog(times[4:-1]/t0, r_av[istar, 4:-1, k0] * dR[istar], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[4:-1]) - np.log(t0), np.log(r_av[istar, 4:-1, k0]) - dR_log[istar], 'o-',
                 color=colorlist5[istar], label=lbl)
        ax3.plot(np.log(times[4:-1]) - np.log(t0), m * (np.log(r_av[istar, 4:-1, k0]) - dR_log[istar]), 'o-',
                 color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[4:-1]) - np.log(t0),
                 np.log(R0) + m * (np.log(times[4:-1]) - np.log(t0)), '-r', label='m=0.54')

        ax0.plot(times[4:tmin]-t0, r_av[istar, 4:tmin, k0] - dR_log[istar], 'ow')
        ax1.loglog(times[4:tmin]-t0, r_av[istar, 4:tmin, k0], 'ow')
        ax2.plot(np.log(times[4:tmin]-t0), np.log(r_av[istar, 4:tmin, k0]) - dR_log[istar],
                 'ow')
    ax0.plot(times[4:-1] / t0, r_av_ens_mean[4:-1] * dR_ens_mean, 'b')
    ax0.set_xlabel('t/t0')
    ax1.set_xlabel('t/t0')
    ax2.set_xlabel('log(t)-log(t0)')
    ax0.set_ylabel('R * R0(t0)/R(t0)')
    ax1.set_ylabel('R * R0(t0)/R(t0)')
    ax2.set_ylabel(r'log(R)-dlog(R)')
    ax3.set_ylabel(r'log(R)-dlog(R)')
    # ax3.set_ylabel('log(R)-[log(R(t0)-log(R0(t0)]')
    # ax3.set_xlim(np.log(times[tmin]), np.log(times[-2]))
    ax2.set_xlim(0, 2.3)
    ax3.set_xlim(0, 2.3)
    ax0.set_ylim(0, 8e3)


    [ax0, ax1, ax2, ax3] = [axes[5, i] for i in range(4)]
    fit_r = np.log(R0) + m * (np.log(times/t0))
    for istar in range(n_params):
        ax0.plot(times[1:-1]/t0, r_av[istar, 1:-1, k0]*dR[istar], '-', color='0.5')
        ax1.loglog(times[1:-1]/t0, r_av[istar, 1:-1, k0] * dR[istar], '-', color='0.5')
        ax2.plot(np.log(times[1:-1]/t0), np.log(r_av[istar, 1:-1, k0]) - dR_log[istar], '-', color='0.5')
    lbl = 'mean(r_av)'
    ax0.plot(times[1:-1]/t0, r_av_ens_mean[1:-1]*dR_ens_mean, '-', color='b', linewidth=2, label=lbl)
    ax1.loglog(times[1:-1]/t0, r_av_ens_mean[1:-1]*dR_ens_mean, '-', color='b', linewidth=2, label=lbl)
    ax2.plot(np.log(times[1:-1]/t0), np.log(r_av_ens_mean[1:-1]*dR_ens_mean), '-', color='b', linewidth=2, label=lbl)
    ax3.plot(np.log(times[4:-1]/t0), np.log(r_av_ens_mean[4:-1]*dR_ens_mean) - np.log(R0)+ m*(np.log(times[4:-1]/t0)),
             '-', color='k', linewidth=2, label='diff')

    ax2.plot(np.log(times[1:-1]/t0), fit_r[1:-1], '-r', label='m='+str(m))
    ax2.legend(loc='best')

    axes[0, 0].legend(loc='best')
    for ax in axes.flatten():
        ax.grid(True)

    plt.suptitle('scaling R')
    plt.subplots_adjust(bottom=0.05, right=.95, left=0.06, top=0.95, wspace=0.25, hspace=0.2)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)





    ''' scaling of U '''
    tmin = tmin_drdt
    fig_name = 'test_U.png'
    U0 = drdt_av[0,tmin, k0]
    dU_log = np.log(drdt_av[:,tmin,k0]) - np.log(U0)
    dU = U0/drdt_av[:,tmin,k0]
    t0 = np.double(tmin*dt_fields)
    U_ens_mean = np.average(drdt_av[:,:,k0], axis=0)
    dU_ens_mean = U0/U_ens_mean[tmin]

    fig, axes = plt.subplots(4, 4, sharex='none', figsize=(18, 20))
    [ax0, ax1, ax2, ax3] = [axes[0, i] for i in range(4)]
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar])
        ax0.plot([t0,t0], [0, np.amax(drdt_av)], 'k', linewidth=1)
        ax1.plot([np.log(t0), np.log(t0)], [0, np.log(np.amax(drdt_av))], 'k', linewidth=1)
        ax0.plot(times[1:-1], drdt_av[istar,1:-1,k0], 'o-', color=colorlist5[istar], label=lbl)
        ax1.loglog(times[1:-1], drdt_av[istar, 1:-1,k0], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[1:-1]), np.log(drdt_av[istar,1:-1,k0]), 'o-', color=colorlist5[istar], label=lbl)

    ax3.plot(np.arange(-1,4),1./dU[:],'-ok', markersize=4)
    ax0.plot(times[1:-1], U_ens_mean[1:-1], 'b')
    ax0.set_xlabel('time')
    ax1.set_xlabel('time')
    ax2.set_xlabel('log(t)')
    ax0.set_ylabel('time')
    ax0.set_ylabel('U')
    ax1.set_ylabel('U')
    ax2.set_ylabel('log(U)')
    ax3.set_xlabel('log2(PE/PE_ref)')
    ax3.set_ylabel('U0')
    ax1.set_xlim(1e2, 5e3)

    [ax0, ax1, ax2, ax3] = [axes[1, i] for i in range(4)]
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar])
        ax0.plot([0,0], [0,np.amax(drdt_av)], 'k', linewidth=1)
        ax0.plot(times[4:-1]/t0, drdt_av[istar,4:-1,k0], 'o-', color=colorlist5[istar], label=lbl)
        ax1.loglog(times[4:-1]/t0, drdt_av[istar, 4:-1,k0], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[4:-1]/t0), np.log(drdt_av[istar,4:-1,k0]), 'o-', color=colorlist5[istar], label=lbl)
    # ax0.plot(times[4:-1]/t0, r_av_ens_mean[4:-1], 'b')
    ax0.set_xlabel('t/t0')
    ax1.set_xlabel('t/t0')
    ax2.set_xlabel('log(t/t0)')
    ax0.set_ylabel('U')
    ax1.set_ylabel('U')
    ax2.set_ylabel('log(U)')
    ax1.set_xlim(5e-1,6)


    [ax0, ax1, ax2, ax3] = [axes[2, i] for i in range(4)]
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar])
        # ax0.plot([0, 0], [0, np.amax(drdt_av)], 'k', linewidth=1)
        ax0.plot(times[4:-1]/t0, drdt_av[istar, 4:-1, k0]*dU[istar], 'o-', color=colorlist5[istar], label=lbl)
        ax1.loglog(times[4:-1]/t0, drdt_av[istar, 4:-1, k0]*dU[istar], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[4:-1]/t0), np.log(drdt_av[istar, 4:-1, k0])-dU_log[istar], 'o-', color=colorlist5[istar], label=lbl)
        ax3.plot(np.log(times[4:-1]/t0), np.log(drdt_av[istar, 4:-1, k0])-dU_log[istar], 'o-', color=colorlist5[istar], label=lbl)

        ax0.plot(times[4:tmin]/t0, drdt_av[istar, 4:tmin, k0]*dU[istar], 'ow')
        ax1.loglog(times[4:tmin]/t0, drdt_av[istar, 4:tmin, k0]*dU[istar], 'ow')
        ax2.plot(np.log(times[4:tmin]/t0), np.log(drdt_av[istar, 4:tmin, k0]) - dU_log[istar], 'ow')
    ax3.plot(np.log(times[4:-1] / t0), np.log(U0) - 0.8 * (np.log(times[4:-1] / t0)), '-r', label='m=-0.8')
    # ax3.set_xlim(np.log(times[tmin]), np.log(times[-2]))
    ax0.set_xlabel('t/t0')
    ax1.set_xlabel('t/t0')
    ax2.set_xlabel('log(t/t0)')
    ax3.set_xlabel('log(t/t0)')
    ax0.set_ylabel('U * U0(t0)/U(t0)')
    ax1.set_ylabel('U * U0(t0)/U(t0)')
    ax2.set_ylabel('log(U) - dlog(U)')
    ax3.set_ylabel('log(U) - dlog(U)')
    ax1.set_xlim(5e-1,6)
    ax3.set_xlim(0, 1.6)
    ax3.legend(loc='best')

    [ax0, ax1, ax2, ax3] = [axes[3, i] for i in range(4)]
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar])
        ax0.plot(times[4:-1] / t0, drdt_av[istar, 4:-1, k0] * dU[istar], '-', color='0.5')
        ax1.loglog(times[4:-1] / t0, drdt_av[istar, 4:-1, k0] * dU[istar], '-', color='0.5')
        ax2.plot(np.log(times[4:-1] / t0), np.log(drdt_av[istar, 4:-1, k0]) - dU_log[istar], '-', color='0.5')
        ax3.plot(np.log(times[4:-1] / t0), np.log(drdt_av[istar, 4:-1, k0]) - dU_log[istar], '-', color='0.5')
    lbl = 'mean(r_av)'
    ax0.plot(times[4:-1] / t0, U_ens_mean[4:-1] * dU_ens_mean, '-', color='b', linewidth=2, label=lbl)
    ax1.loglog(times[4:-1]/t0, U_ens_mean[4:-1]*dU_ens_mean, '-', color='b', linewidth=2, label=lbl)
    ax2.plot(np.log(times[4:-1]/t0), np.log(U_ens_mean[4:-1]*dU_ens_mean), '-', color='b', linewidth=3, label=lbl)
    ax3.plot(np.log(times[4:-1]/t0), np.log(U_ens_mean[4:-1]*dU_ens_mean), '-', color='b', linewidth=3, label=lbl)
    for m in np.arange(0.5,1,0.1):
        ax3.plot(np.log(times[4:-1] / t0), np.log(U0) - m* (np.log(times[4:-1] / t0)), '-r', linewidth=1, label='m=-'+str(m))
    ax0.set_xlabel('t/t0')
    ax1.set_xlabel('t/t0')
    ax2.set_xlabel('log(t/t0)')
    ax3.set_xlabel('log(t/t0)')
    ax0.set_ylabel('U * U0(t0)/U(t0)')
    ax1.set_ylabel('U * U0(t0)/U(t0)')
    ax2.set_ylabel('log(U) - dlog(U)')
    ax3.set_ylabel('log(U) - dlog(U)')
    ax1.set_xlim(5e-1, 6)
    ax3.set_xlim(0, 1.6)
    ax3.legend(loc='best')

    axes[0,0].legend(loc='best')
    for ax in axes.flatten():
        ax.grid(True)

    plt.suptitle('scaling U')
    plt.subplots_adjust(bottom=0.075, right=.95, left=0.06, top=0.9, wspace=0.25, hspace=0.2)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)




    return


# ----------------------------------------------------------------------
def plot_PE_scaling_U(r_av, U_rad_av, drdt_av, scaling, k0, trange, nt, n_params, fig_name):
    fig, axes = plt.subplots(1, 3, sharex='none', figsize=(15, 5))
    ax0 = axes[0]
    ax1 = axes[1]
    ax2 = axes[2]
    # # scaling = [-1, 0, 1, 2, 3]
    for t0 in trange:
        it = np.where(times == t0)[0][0]
        # ax0.plot(scaling, r_av[1,it,k0] + 9.5e2 * np.asarray(scaling), 'k', linewidth=1)
        # ax0.plot(scaling, r_av[3,it,k0] + 12e2 * (np.asarray(scaling) - 2), 'k--', linewidth=1)
        # ax0.plot(scaling, 2e3 + 0.1 * 1e3 * (np.asarray(scaling) + 2) ** 2, 'r', linewidth=1)
        ax0.plot(scaling, r_av[:, it, k0], 'o-', label='t=' + str(t0))
        ax1.plot(scaling, U_rad_av[:, it, k0], 'o-', label='t=' + str(t0))
        ax2.plot(scaling, drdt_av[:, it, k0], 'o-', label='t=' + str(t0))
    ax0.set_xlabel('log2(PE/PE0)')
    ax1.set_xlabel('log2(PE/PE0)')
    ax2.set_xlabel('log2(PE/PE0)')
    ax0.set_ylabel('r_av  [m]')
    ax1.set_ylabel('U_rad_av [m/s]')
    ax2.set_ylabel('d(r_av)/dt [m/s]')
    ax1.legend(loc='best')
    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    return




def plot_PE_scaling(r_av, U_rad_av, scaling, k0, trange, nt, n_params, fig_name, colorlist):
    print('plot PE scaling')

    PEPE0 = [2**s for s in scaling]
    it0 = np.where(times == trange[0])[0][0]
    it1 = np.where(times == trange[-1])[0][0]
    rmax = 1.1*np.amax(r_av[:,it1,k0])
    r0 = r_av[1, :, k0]     # radius of PE=PE0
    dr0 = r_av[1, :, k0]-r_av[1, it0, k0]
    dr0_log = np.log10(r_av[1,:,k0])-np.log10(r_av[1,it0,k0])

    fig, axes = plt.subplots(3, 4, sharex='none', figsize=(18, 10))
    [ax0, ax1, ax2, ax3] = [axes[0,i] for i in range(4)]
    for ic,t0 in enumerate(trange):
        it = np.where(times == t0)[0][0]
        ax0.plot(PEPE0, r_av[:, it, k0], 'o-', color=colorlist[ic], label='t=' + str(t0))
        ax1.plot(scaling, r_av[:, it, k0], 'o-', color=colorlist[ic], label='t=' + str(t0))
        ax2.plot(scaling, np.log10(r_av[:, it, k0]), 'o-', color=colorlist[ic], )
        ax3.loglog(PEPE0, r_av[:, it, k0], 'o-', color=colorlist[ic], )
        ax0.plot(PEPE0[1], (r0[it]), 'wo')
        ax1.plot(scaling[1], (r0[it]), 'wo')
        ax2.plot(scaling[1], np.log10(r0[it]), 'wo')
        ax3.loglog(PEPE0[1], (r0[it]), 'wo')
    ax0.legend(loc='best')

    [ax0, ax1, ax2, ax3] = [axes[1,i] for i in range(4)]
    ax2.fill_between(np.arange(-1,4), 5e4*np.ones(5))
    for ic,t0 in enumerate(trange):
        it = np.where(times == t0)[0][0]
        ax0.plot(PEPE0, r_av[:, it, k0]-dr0[it], 'o-', color=colorlist[ic], label='t=' + str(t0))
        ax1.plot(scaling, r_av[:, it, k0]-dr0[it], 'o-', color=colorlist[ic], label='t=' + str(t0))
        ax2.semilogy(scaling, r_av[:, it, k0]-dr0[it], 'o-', color=colorlist[ic], )
        ax3.loglog(PEPE0, r_av[:, it, k0]-dr0[it], 'o-', color=colorlist[ic], )

    [ax0, ax1, ax2, ax3] = [axes[2,i] for i in range(4)]
    for ic,t0 in enumerate(trange):
        it = np.where(times == t0)[0][0]
        ax0.plot(PEPE0, r_av[:,it,k0]*r0[it0]/r0[it], 'o-', color=colorlist[ic], label='t=' + str(t0))
        ax1.plot(scaling, r_av[:,it,k0]*r0[it0]/r0[it], 'o-', color=colorlist[ic], label='t=' + str(t0))
        ax2.plot(scaling, np.log10(r_av[:, it, k0])-dr0_log[it], 'o-', color=colorlist[ic], )
        ax3.semilogx(PEPE0, np.log10(r_av[:, it, k0])-dr0_log[it], 'o-', color=colorlist[ic], )

    r_fit = np.ndarray((5))
    # for i,m in enumerate(np.arange(0.0, 0.12, 0.02)):
    #     r_fit[i] = r0[it0]+m*scaling
    for i,m in enumerate(np.arange(0.04, 0.12, 0.02)):
        ax2.plot(scaling, np.log10(r0[it0])+m*np.asarray(scaling), color=str(np.double(i)/5), linewidth=1)

    ax0.set_xlabel('PE/PE0')
    ax1.set_xlabel('log2(PE/PE0)')
    ax2.set_xlabel('log2(PE/PE0)')
    ax3.set_xlabel('log(PE/PE0)')
    for i in range(3):
        axes[i, 0].set_ylabel('r_av  [m]')
        axes[i, 1].set_ylabel('r_av  [m]')
        axes[i, 2].set_ylabel('log(r_av)  [m]')
        axes[i, 3].set_ylabel('log(r_av)  [m]')
    axes[0, 2].set_ylabel('log10(r_av)  [m]')
    axes[2, 2].set_ylabel('log10(r_av)  [m]')
    axes[2, 3].set_ylabel('log10(r_av)  [m]')
    # axes[0,2].set_ylim(3,5)
    for ax in axes[:,3].flatten():
        ax.set_xlim(4e-1,1e1)
    for ax in axes[:, 2].flatten():
        ax.set_xlim(-1,3)
    for ax in axes[:2,:2].flatten():
        ax.set_ylim(1e3,rmax)
    for ax in axes.flatten():
        ax.grid(True, which='major', axis='x')
        ax.grid(True, which='both', axis='y')
    for ax in axes[:,0].flatten():
        ax.legend(loc='best')
    plt.suptitle('scaling of radius with PE', fontsize=18)
    # fig.tight_layout()
    plt.subplots_adjust(bottom=0.075, right=.95, left=0.06, top=0.9, wspace=0.25, hspace=0.2)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)



    ''' compute slope m '''
    m10_r = np.ndarray(len(times))
    m10_r_ = np.ndarray((n_params, len(times)))
    for it,t0 in enumerate(times):
        m10_r_[2:,it] = np.log(r_av[2:,it,k0]/r0[it]) / scaling[2:]
        # m10[it] = np.mean(np.log(r_av[:,it,k0]/r0[it])/scaling[:])
    m10_r[:] = np.mean(m10_r_[2:,:], axis=0)
    m_r = m10_r * np.log(2)/np.log(10)

    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, sharex='none', figsize=(18, 5))
    ax0.plot(times, m10_r, '-o')
    ax1.plot(times, m_r, '-o')
    ax0.set_title('m10')
    ax1.set_title('m')
    for ax in (ax0, ax1):
        ax.grid(True)
    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, fig_name[:-4] + '_slope.png'))
    plt.close(fig)






    return

# ----------------------------------------------------------------------
def plot_vel_normalized_w_av(r_av, times, istar, k0, id, fig_name):
    # read in azimuthally averaged vertical velocity
    rootgrp = nc.Dataset(os.path.join(path_root, id, 'data_analysis', 'stats_radial_averaged.nc'))
    ts_grp = rootgrp.groups['timeseries']
    t_out = ts_grp.variables['time'][:]
    stats_grp = rootgrp.groups['stats']
    w = stats_grp.variables['w'][:,:,:] # dimensions (nt, nr, nz)
    nr = stats_grp.dimensions['nr'].size
    n_times = ts_grp.dimensions['nt'].size
    rootgrp.close()

    # defining geometry for single CP
    x_array = np.arange(0, nx) * dx[0]
    ic = nx / 2
    # r_array = x_array - dx[0] * ic
    r_array = np.arange(0,nr)*dx[0]

    # plotting parameters
    imin = 0
    imax = 80
    k0 = 0

    fig, axes = plt.subplots(2, 2, sharex='none', figsize=(20, 12))
    ax0 = axes[0, 0]
    ax1 = axes[1, 0]
    ax2 = axes[0, 1]
    ax3 = axes[1, 1]
    ax0.plot(r_av[istar, :, k0], times, '-o', label=id)
    ax0.set_xlim(r_array[imin], r_array[imax])
    ax0.grid()
    it = 0
    for t0 in t_out[1::2]:
        if t0 >= np.int(times[0]) and t0 <= times[-1]:
            count = np.double(t0) / times[-1]
            print '>> it', it, t0, count
            iR = r_av[istar, it, k0] / dx[0]
            ax1.plot(r_array[imin:imax], w[it, imin:imax, k0], '-', color=cm_grey2(count),
                     label='t=' + str(t0) + 's')
            ax1.plot(r_array[iR], w[it, iR, k0], 'd', color=cm_grey2(count), markersize=12)
            R_array = r_array / r_av[istar, it, k0]
            ax2.plot(R_array[imin:imax], w[it, imin:imax, k0], '-', color=cm_grey2(count), linewidth=3,
                     label='t=' + str(t0) + 's')
            ax3.plot(R_array[imin:imax], w[it, imin:imax, k0], '-', color=cm_grey2(count), linewidth=3,
                     label='t=' + str(t0) + 's')
            ax3.plot(R_array[iR], w[it, iR, k0], 'd', color=cm_grey(count), markersize=12)
            it += 1
    fig.suptitle(id)
    ax0.legend(loc=3)
    ax1.legend(loc='lower center', bbox_to_anchor=(0.5, 0.1),
               fancybox=True, shadow=True, ncol=6, fontsize=11)
    ax1.set_xlim(r_array[imin], r_array[imax])
    ax2.set_xlim([0., 1.5])
    ax3.set_xlim([0.75, 1.25])
    ax3.set_ylim([-1.1, 2.])
    ax0.set_xlabel('r_av  [m]')
    ax1.set_xlabel('r [m]')
    ax2.set_xlabel('r/r_av [m]')
    ax3.set_xlabel('r/r_av [m]')
    ax0.set_ylabel('times [s]')
    ax1.set_ylabel('w  [m/s]')
    ax2.set_ylabel('w  [m/s]')
    ax3.set_ylabel('w  [m/s]')
    ax2.set_title('w normalized by tracer radius')
    ax3.set_title('w normalized by tracer radius')
    ax1.grid()
    ax2.grid()
    ax3.grid()
    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    return


def plot_vel_normalized(r_av, times, istar, k0, nx, ic, jc, id, fig_name):
    # defining geometry for single CP
    x_array = np.arange(0, nx) * dx[0]
    r_array = x_array - dx[0] * ic

    # read in vertical velocity for all times
    print(os.path.join(path_root, id, 'fields_merged', 'fields_allt_xz_j'+str(jc)+'.nc'))
    rootgrp = nc.Dataset(os.path.join(path_root, id, 'fields_merged', 'fields_allt_xz_j'+str(jc)+'.nc'))
    w = rootgrp.variables['w'][:, :, :]
    t_out = rootgrp.variables['time'][:]
    rootgrp.close()

    #plotting parameters
    imin = ic
    imax = ic + 80
    print('plotting: ', ic, jc, imin, imax, nx, ny)
    print(r_array)

    fig, axes = plt.subplots(2, 2, sharex='none', figsize=(20, 12))
    ax0 = axes[0,0]
    ax1 = axes[1,0]
    ax2 = axes[0,1]
    ax3 = axes[1,1]
    ax0.plot(r_av[istar, :, k0], times, '-o', label=id)
    ax0.set_xlim(r_array[imin], r_array[imax])
    ax0.grid()
    it = 0
    for t0 in t_out[0::2]:
        if t0 >= np.int(times[0]) and t0 <= times[-1]:
            count = np.double(t0) / times[-1]
            print '>> it', it, t0, count
            iR = r_av[istar, it, k0] / dx[0] + ic
            ax1.plot(r_array[imin:imax], w[it, imin:imax, k0], '-', color=cm_grey2(count),
                     label='t=' + str(t0) + 's')
            ax1.plot(r_array[iR], w[it, iR, k0], 'd', color=cm_grey2(count), markersize=12)
            R_array = r_array / r_av[istar, it, k0]
            ax2.plot(R_array[imin:imax], w[it, imin:imax, k0], '-', color=cm_grey2(count), linewidth=3,
                     label='t=' + str(t0) + 's')
            ax3.plot(R_array[imin:imax], w[it, imin:imax, k0], '-', color=cm_grey2(count), linewidth=3,
                     label='t=' + str(t0) + 's')
            ax3.plot(R_array[iR], w[it, iR, k0], 'd', color=cm_grey(count), markersize=12)
            it += 1
    fig.suptitle(id)
    ax0.legend(loc=3)
    ax1.legend(loc='lower center', bbox_to_anchor=(0.5, 0.1),
               fancybox=True, shadow=True, ncol=6, fontsize=11)
    # ax2.legend()
    ax1.set_xlim(r_array[imin], r_array[imax])
    ax2.set_xlim([-0.2, 1.2])
    ax3.set_xlim([0.75, 1.2])
    ax3.set_ylim([-1.1, 2.])
    ax0.set_xlabel('r_av  [m]')
    ax1.set_xlabel('r [m]')
    ax2.set_xlabel('r/r_av [m]')
    ax3.set_xlabel('r/r_av [m]')
    ax0.set_ylabel('times [s]')
    ax1.set_ylabel('w  [m/s]')
    ax2.set_ylabel('w  [m/s]')
    ax3.set_ylabel('w  [m/s]')
    ax2.set_title('w normalized by tracer radius')
    ax3.set_title('w normalized by tracer radius')
    ax1.grid()
    ax2.grid()
    ax3.grid()
    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    print('')
    return


# ----------------------------------------------------------------------
def plot_dist_vel(r_av, drdt_av, U_rad_av,
                  dTh_params, z_params, r_params, n_params, k0,
                  id_ref, title, colorlist, fig_name):

    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, sharex='all', figsize=(18, 5))
    for istar in range(n_params):
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[0]
        rstar = r_params[istar]
        if rstar == 1000:
            id = id_ref
        else:
            id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        ax0.plot(times, r_av[istar, :, k0], 'o-', color=colorlist[istar], label=id)
        ax1.plot(times[1:], drdt_av[istar, 1:, k0], 'o-', color=colorlist[istar], label=id)
        ax2.plot(times, U_rad_av[istar, :, k0], 'o-', color=colorlist[istar], label=id)
    # ax0.plot(times, r_av_ref[:, k0], 'ko-', label=id_ref)
    # ax1.plot(times[1:], drdt_av_ref[1:, k0], 'ko-', label=id_ref)
    # ax2.plot(times, U_rad_av_ref[:, k0], 'ko-', label=id_ref)

    ax0.set_xlabel('times [s]')
    ax1.set_xlabel('times [s]')
    ax2.set_xlabel('times [s]')
    ax0.set_ylabel('r_av  [m]')
    ax1.set_ylabel('drdt_av')
    ax2.set_ylabel('U_rad_av  [m/s]')
    ax2.legend()
    # fig.suptitle('CP rim (dTh=' + str(dTh) + ')')
    fig.suptitle(title, fontsize=18)
    fig.tight_layout()
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.06, top=0.9, wspace=0.15)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    ''' plot logarithmically '''
    fig, axis = plt.subplots(3, 3, sharex='none', figsize=(18, 15))
    for istar in range(n_params):
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[0]
        rstar = r_params[istar]
        if rstar == 1000:
            id = id_ref
        else:
            id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        axis[0,0].semilogx(times, r_av[istar, :, k0], 'o-', color=colorlist[istar], label=id)
        axis[0,1].semilogx(times, drdt_av[istar, :, k0], 'o-', color=colorlist[istar], label=id)
        axis[0,2].semilogx(times, U_rad_av[istar, :, k0], 'o-', color=colorlist[istar], label=id)
        axis[1,0].loglog(times, r_av[istar, :, k0], 'o-', color=colorlist[istar], label=id)
        axis[1,1].loglog(times, drdt_av[istar, :, k0], 'o-', color=colorlist[istar], label=id)
        axis[1,2].loglog(times, U_rad_av[istar, :, k0], 'o-', color=colorlist[istar], label=id)
        axis[2,0].semilogy(times, r_av[istar, :, k0], 'o-', color=colorlist[istar], label=id)
        axis[2,1].semilogy(times, drdt_av[istar, :, k0], 'o-', color=colorlist[istar], label=id)
        axis[2,2].semilogy(times, U_rad_av[istar, :, k0], 'o-', color=colorlist[istar], label=id)
    # axis[0, 0].semilogx(times, r_av_ref[:, k0], 'ko-', label=id)
    # axis[0, 1].semilogx(times, drdt_av_ref[:, k0], 'ko-', label=id)
    # axis[0, 2].semilogx(times, U_rad_av_ref[:, k0], 'ko-', label=id)
    # axis[1, 0].loglog(times, r_av_ref[:, k0], 'ko-', label=id)
    # axis[1, 1].loglog(times, drdt_av_ref[:, k0], 'ko-', label=id)
    # axis[1, 2].loglog(times, U_rad_av_ref[:, k0], 'ko-', label=id)
    # axis[2, 0].semilogy(times, r_av_ref[:, k0], 'ko-', label=id)
    # axis[2, 1].semilogy(times, drdt_av_ref[:, k0], 'ko-', label=id)
    # axis[2, 2].semilogy(times, U_rad_av_ref[:, k0], 'ko-', label=id)
    axis[0, 0].set_title('CP radius (r_av)', fontsize=18)
    axis[0, 1].set_title('d(r_av)/dt', fontsize=18)
    axis[0, 2].set_title('radial spreading velocity (U_av)', fontsize=18)
    # axis[0,0].set_xlabel('log(times) [s]')
    for i in range(3):
        axis[0, i].set_xlabel('log(times) [s]')
        axis[1, i].set_xlabel('log(times) [s]')
        axis[2, i].set_xlabel('times [s]')
    axis[0, 0].set_ylabel('r_av  [m]')
    axis[1, 0].set_ylabel('log(r_av)  [m]')
    axis[2, 0].set_ylabel('log(r_av)  [m]')
    axis[0, 1].set_ylabel('dr/dt  [m/s]')
    axis[1, 1].set_ylabel('log(dr/dt)  [m/s]')
    axis[2, 1].set_ylabel('log(dr/dt)  [m/s]')
    axis[0, 2].set_ylabel('U_rad_av  [m/s]')
    axis[1, 2].set_ylabel('log(U_rad_av)  [m/s]')
    axis[2, 2].set_ylabel('log(U_rad_av)  [m/s]')
    axis[0,1].set_ylim(-0.1,6)
    axis[0,2].set_ylim(-0.1,6)
    fig.tight_layout()
    plt.subplots_adjust(bottom=0.075, right=.95, left=0.07, top=0.9, wspace=0.25)
    axis[1,1].legend(loc='best')
    fig.suptitle(title, fontsize=21)
    fig.savefig(os.path.join(path_out_figs, fig_name[:-4]+'_log.png'))
    plt.close(fig)
    return



def plot_dist_fitting(r_av,
                      dTh_params, z_params, r_params, n_params, tmin, k0,
                      id_ref, colorlist, fig_name):
    # Fit the first set
    # fitfunc = lambda p, x: p[0] * np.cos(2 * np.pi / p[1] * x + p[2]) + p[3] * x  # Target function
    fitfunc1 = lambda p, x: p[0] + p[1] * x ** p[2]  # Target function
    errfunc = lambda p, x, y: fitfunc1(p, x) - y  # Distance to the target function
    fitfunc2 = lambda p, x: p[0] + p[1] * np.log(x)  # Target function


    # matrix of optimal parameters
    f1 = np.zeros((n_params, 3))
    f_log1 = np.zeros((n_params, 3))
    f_log2 = np.zeros((n_params, 3))
    f_log_romps = np.zeros((n_params, 3))

    # initial guess for the parameters of power-law relation
    n_init = 4
    p0 = np.zeros((n_init, 3))
    p0[:,0] = [0., -1e4, -1e3, -1e2]
    p0[:,2] = [0., 1e-2, 1e-1, 1e0]
    p1 = np.zeros((n_init**2, 3))
    # initial guess for the parameters of log-relation
    p0_log = np.zeros((3))
    p0_log[0] = 0.
    p0_log[1] = 1.
    p0_log[2] = .5

    rmin = np.amin(r_av[:, :, k0])
    for istar in range(n_params):
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[0]
        rstar = r_params[istar]
        if rstar == 1000:
            id = id_ref
        else:
            id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        print('fitting radius R: '+id)
        rmax = np.amax(r_av[istar, tmin:, k0])


        ''' fitting '''
        popt, pcov = optimize.curve_fit(func, times[tmin:], r_av[istar, tmin:, k0], p0[0, :])
        f1[istar, :] = popt
        popt_log1, pcov = optimize.curve_fit(func_log1, times[tmin:], r_av[istar, tmin:, k0], p0_log[:])
        f_log1[istar, :] = popt_log1
        popt_log2, pcov = optimize.curve_fit(func_log2, times[tmin:], r_av[istar, tmin:, k0], p0_log[:])
        f_log2[istar, :] = popt_log2
        popt_log_romps, pcov = optimize.curve_fit(func_log_romps_r, times[tmin:], r_av[istar, tmin:, k0], p0_log[:])
        f_log_romps[istar, :] = popt_log_romps


        ''' plotting separately for each case '''
        fig, axis = plt.subplots(2, 4, sharex='none', figsize=(18, 12))
        for i in range(2):
            axis[i,0].plot(times, r_av[istar, :, k0], 'ko-', label=id)
            axis[i,1].semilogx(times, r_av[istar, :, k0], 'ko-', label=id)
            axis[i,2].loglog(times, r_av[istar, :, k0], 'ko-', label=id)
            axis[i,3].semilogy(times, r_av[istar, :, k0], 'ko-', label=id)
            # mark the ones used for fitting
            axis[i,0].plot(times[:tmin], r_av[istar, :tmin, k0], 'ow', alpha=0.75)
            axis[i,1].semilogx(times[:tmin], r_av[istar, :tmin, k0], 'ow', alpha=0.75)
            axis[i,2].loglog(times[:tmin], r_av[istar, :tmin, k0], 'ow', alpha=0.75)
            axis[i,3].semilogy(times[:tmin], r_av[istar, :tmin, k0], 'ow', alpha=0.75)
        # reference CP
        axis[0,0].plot(times, r_av[1, :, k0], 'd-', color='0.7', label=id_ref)
        axis[0,1].semilogx(times, r_av[1, :, k0], 'd-', color='0.7', label=id_ref)
        axis[0,2].loglog(times, r_av[1, :, k0], 'd-', color='0.7', label=id_ref)
        axis[0,3].semilogy(times, r_av[1, :, k0], 'd-', color='0.7', label=id_ref)
        # fitted functions
        axis[1,0].plot(times, func(times, *popt), 'r-', label='a+b*x^c ('+str(popt)+')')
        axis[1,1].semilogx(times, func(times, *popt), 'r-', label='a+b*x^c ('+str(popt)+')')
        axis[1,2].loglog(times, func(times, *popt), 'r-')
        axis[1,3].semilogy(times, func(times, *popt), 'r-')
        axis[1,0].plot(times, func_log1(times, *popt_log1), '-', linewidth=3, color='limegreen', label='b*log(x) ('+str(popt_log1)+')')
        axis[1,1].semilogx(times, func_log1(times, *popt_log1), '-', color='limegreen', label='b*log(x) ('+str(popt_log1)+')')
        axis[1,2].loglog(times, func_log1(times, *popt_log1), '-', color='limegreen', )
        axis[1,3].semilogy(times, func_log1(times, *popt_log1), '-', color='limegreen', )
        axis[1,0].plot(times, func_log2(times, *popt_log2), '-', color='darkgreen', label='a+b*log(x) ('+str(popt_log2)+')')
        axis[1,1].semilogx(times, func_log2(times, *popt_log2), '-', color='darkgreen', label='a+b*log(x) ('+str(popt_log2)+')')
        axis[1,2].loglog(times, func_log2(times, *popt_log2), '-', color='darkgreen', )
        axis[1,3].semilogy(times, func_log2(times, *popt_log2), '-', color='darkgreen', )
        axis[1,0].plot(times, func_log_romps_r(times, *popt_log_romps), '-', color='navy', label='a+b*log(1+c*x) (' + str(popt_log_romps) + ')')
        axis[1,1].semilogx(times, func_log_romps_r(times, *popt_log_romps), '-', color='navy', label='a+b*log(1+c*x) (' + str(popt_log_romps) + ')')
        axis[1,2].loglog(times, func_log_romps_r(times, *popt_log_romps), '-', color='navy')
        axis[1,3].semilogy(times, func_log_romps_r(times, *popt_log_romps), '-', color='navy')
        for i in range(n_init):
            for j in range(n_init):
                ij = i*n_init+j
                p1[ij, :], success = optimize.leastsq(errfunc, [p0[i,0], p0[i,1], p0[j,2]], args=(times[tmin:], r_av[istar, tmin:, k0]))
                # print(i, 'success: ', success)
                axis[0,0].plot(times[tmin:], fitfunc1(p1[i,:], times[tmin:]), "-", label='p='+str(p1[ij,:]))  # Plot of the data and the fit
                axis[0,1].semilogx(times, fitfunc1(p1[i,:], times), "-", label='p='+str(p1[ij,:]))  # Plot of the data and the fit
                axis[0,2].loglog(times, fitfunc1(p1[i,:], times), "-", label='p='+str(p1[ij,:]))  # Plot of the data and the fit
                axis[0,3].semilogy(times, fitfunc1(p1[i,:], times), "-", label='p='+str(p1[ij,:]))  # Plot of the data and the fit
        axis[0,0].set_title('scipy.optizmize.leastsq: a + b*x^c', fontsize=18)
        axis[1,0].set_title('scipy.optizmize.curve_fit: a + b*x^c, a+b*log(c), a+b*log(1+c)', fontsize=18)
        for i in range(2):
            axis[i,0].set_xlabel('time [s]')
            axis[i,1].set_xlabel('log(time) [s]')
            axis[i,2].set_xlabel('log(time) [s]')
            axis[i,3].set_xlabel('time [s]')
            axis[i,0].set_ylabel('r_av  [m/s]')
            axis[i,2].set_ylabel('log(r_av)  [m/s]')
            axis[i,1].set_xlim([1e2,4e3])
            axis[i,2].set_xlim([1e2,4e3])
            axis[i,0].set_ylim([0,rmax])
            axis[i,1].set_ylim([0,rmax])
            axis[i,2].set_ylim([1.5e3,rmax])
            axis[i,3].set_ylim([1.5e3,rmax])
            # axis[i,3].set_ylim([1e3,1e4])
        plt.subplots_adjust(bottom=0.15, right=.95, left=0.07, top=0.93, wspace=0.25, hspace=0.7)
        axis[0,0].legend(loc='upper left', bbox_to_anchor=(0., -0.2),
                   fancybox=True, shadow=False, ncol=3)
        axis[1,0].legend(loc='upper left', bbox_to_anchor=(0., -0.2),
                          fancybox=True, shadow=False, ncol=2)
        fig.suptitle('CP radius (' + id +', dx='+str(dx[0]) +'m)', fontsize=21)
        fig.savefig(os.path.join(path_out_figs, fig_name[:-4] + '_fit_r_' + id + '.png'))
        plt.close(fig)




    ''' all '''
    fig, axis = plt.subplots(3, 4, sharex='none', figsize=(18, 16))
    for istar in range(n_params):
        col = colorlist[istar]
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[0]
        rstar = r_params[istar]
        if rstar == 1000:
            id = id_ref
        else:
            id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        for i in range(3):
            axis[i,0].plot(times, r_av[istar, :, k0], 'o-', linewidth=1, color=col, label=id)#, color=cmap(i_color), label=id)
            axis[i,1].semilogx(times, r_av[istar, :, k0], 'o-', linewidth=1, color=col, label=id)#, color=cmap(i_color), label=id)
            axis[i,2].loglog(times, r_av[istar, :, k0], 'o-', linewidth=1, color=col, label=id)#, color=cmap(i_color), label=id)
            axis[i,3].semilogy(times, r_av[istar, :, k0], 'o-', linewidth=1, color=col, label=id)#, color=cmap(i_color), label=id)
            # mark points thatt are not used for fitting
            axis[i,0].plot(times[:tmin], r_av[istar, :tmin, k0], 'ow', alpha=0.7)
            axis[i,1].semilogx(times[:tmin], r_av[istar, :tmin, k0], 'ow', alpha=0.7)
            axis[i,2].loglog(times[:tmin], r_av[istar, :tmin, k0], 'ow', alpha=0.7)
            axis[i,3].semilogy(times[:tmin], r_av[istar, :tmin, k0], 'ow', alpha=0.7)
        axis[0, 0].plot(times, func(times, *f1[istar, :]), '-', linewidth=3, color=col,
                       label='a+b*x^c (' + str(f1[istar,:]) + ')')
        axis[0, 1].semilogx(times, func(times, *f1[istar, :]), '-', linewidth=3, color=col,
                           label='a+b*x^c (' + str(f1[istar,:]) + ')')
        axis[0, 2].loglog(times, func(times, *f1[istar, :]), '-', linewidth=3, color=col,
                         label='a+b*x^c (' + str(f1[istar,:]) + ')')
        axis[0, 3].semilogy(times, func(times, *f1[istar, :]), '-', linewidth=3, color=col,
                            label='a+b*x^c (' + str(f1[istar, :]) + ')')
        axis[1, 0].plot(times, func_log_romps_r(times, *f_log_romps[istar, :]), '-', linewidth=3, color=col,
                            label='a+b*log(1+c*x) (' + str(np.round(f_log_romps[istar, :], 5)) + ')')
        axis[1, 1].semilogx(times, func_log_romps_r(times, *f_log_romps[istar, :]), '-', linewidth=3, color=col,
                            label='a+b*log(1+c*x) (' + str(np.round(f_log_romps[istar, :], 5)) + ')')
        axis[1, 2].loglog(times, func_log_romps_r(times, *f_log_romps[istar, :]), '-', linewidth=3, color=col,
                            label='a+b*log(1+c*x) (' + str(np.round(f_log_romps[istar, :], 5)) + ')')
        axis[1, 3].semilogy(times, func_log_romps_r(times, *f_log_romps[istar, :]), '-', linewidth=3, color=col,
                            label='a+b*log(1+c*x) (' + str(np.round(f_log_romps[istar, :], 5)) + ')')
        axis[2, 0].plot(times, func_log2(times, *f_log2[istar, :]), '-', linewidth=3, color=col,
                            label='a+b*log(x) (' + str(f_log2[istar, :]) + ')')
        axis[2, 1].semilogx(times, func_log2(times, *f_log2[istar, :]), '-', linewidth=3, color=col,
                            label='a+b*log(x) (' + str(f_log2[istar, :]) + ')')
        axis[2, 2].loglog(times, func_log2(times, *f_log2[istar, :]), '-', linewidth=3, color=col,
                            label='a+b*log(x) (' + str(f_log2[istar, :]) + ')')
        axis[2, 3].semilogy(times, func_log2(times, *f_log2[istar, :]), '-', linewidth=3, color=col,
                            label='a+b*log(x) (' + str(f_log2[istar, :]) + ')')
    for i in range(2,3):
        axis[i,0].set_xlabel('time [s]')
        axis[i,1].set_xlabel('log(time) [s]')
        axis[i,2].set_xlabel('log(time) [s]')
        axis[i,3].set_xlabel('time [s]')
    for i in range(3):
        axis[i,0].set_ylabel('r_av  [m]')
        axis[i,2].set_ylabel('log(r_av)  [m]')
        axis[i,1].set_xlim(0,4e3)
        axis[i,2].set_xlim(0,4e3)
        axis[i,2].set_ylim(rmin,2e4)
        axis[i,3].set_ylim(rmin,2e4)

    textprops = dict(facecolor='white', alpha=0.5, linewidth=0.)
    axis[0, 0].text(0.05, 0.95, 'Power Law: R=a+b*t^c', transform=axis[0,0].transAxes, fontsize=16,
                    verticalalignment='top', bbox=textprops)
    axis[1, 0].text(0.05, 0.95, 'Romps: R=a+b*log(1+c*t)', transform=axis[1,0].transAxes, fontsize=16,
                    verticalalignment='top', bbox=textprops)
    axis[2, 0].text(0.05, 0.95, 'Logarithm: R=a+b*log(t)', transform=axis[2,0].transAxes, fontsize=16,
                    verticalalignment='top', bbox=textprops)

    plt.subplots_adjust(bottom=0.1, right=.97, left=0.07, top=0.95, wspace=0.25, hspace=0.35)
    axis[0, 0].legend(loc='upper left', bbox_to_anchor=(-0.1, -0.07),
                     fancybox=True, shadow=False, ncol=3)
    axis[1,0].legend(loc='upper left', bbox_to_anchor=(-0.1, -0.07),
                     fancybox=True, shadow=False, ncol=3)
    axis[2, 0].legend(loc='upper left', bbox_to_anchor=(-0.1, -0.15),
                      fancybox=True, shadow=False, ncol=3)
    fig.suptitle('CP radius (r_av) (dx=' + str(dx[0]) + 'm)', fontsize=21)
    fig.savefig(os.path.join(path_out_figs, fig_name[:-4] + '_fit_r_all.png'))
    plt.close(fig)







    return f1, f_log1, f_log2, f_log_romps





def plot_vel_fitting(r_av, vel, init_params,
                      dTh_params, z_params, r_params, n_params, tmin, k0,
                      id_ref, colorlist, fig_name):
    print('fitting velocity')

    # Fit the first set
    # fitfunc = lambda p, x: p[0] * np.cos(2 * np.pi / p[1] * x + p[2]) + p[3] * x  # Target function
    fitfunc1 = lambda p, x: p[0] + p[1] * x ** p[2]  # Target function
    errfunc = lambda p, x, y: fitfunc1(p, x) - y  # Distance to the target function

    # matrix of optimal parameters
    f1 = np.zeros((n_params, 3))
    f_log1 = np.zeros((n_params, 3))
    f_log2 = np.zeros((n_params, 3))
    f_log_romps = np.zeros((n_params, 3))

    # initial guess for the parameters of power-law relation
    n_init = 3
    p0 = np.zeros((n_init, 3))
    p1 = np.zeros((n_init**3, 3))
    a = [1e0, 1e2, 1e3]
    p0[:,0] = a
    b = [-1e1, 1e0, 1e1]
    p0[:,1] = b
    c = [-1e-1, -1e0, 1e-1]
    p0[:, 2] = c
    # initial guess for the parameters of log-relation
    p0_log = np.zeros((3))
    p0_log[0] = 0.
    p0_log[1] = 1.
    p0_log[2] = .5

    vmax = np.amax(vel[:, :, k0])

    for istar in range(n_params):
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[0]
        rstar = r_params[istar]
        if rstar == 1000:
            id = id_ref
        else:
            id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        print('fitting vel: ' + id, 'tmin: '+str(tmin))
        tmax = -1


        ''' fitting '''
        print('init', init_params, tmin)
        try:
            popt, pcov = optimize.curve_fit(func, times[tmin:tmax], vel[istar, tmin:tmax, k0], init_params)
            f1[istar, :] = popt
        except:
            print('no fitting possible')
            f1[istar,:] = [1,0,1]
        popt_log1, pcov = optimize.curve_fit(func_log1, times[tmin:], vel[istar, tmin:, k0], p0_log[:])
        f_log1[istar, :] = popt_log1
        popt_log2, pcov = optimize.curve_fit(func_log2, times[tmin:], vel[istar, tmin:, k0], p0_log[:])
        f_log2[istar, :] = popt_log2
        popt_log_romps, pcov = optimize.curve_fit(func_log_romps_U, times[tmin:], vel[istar, tmin:, k0], p0_log[:])
        f_log_romps[istar, :] = popt_log_romps


        ''' plotting separately for each case '''
        fig, axis = plt.subplots(2, 4, sharex='none', figsize=(18, 12))
        for i in range(2):
            axis[i,0].plot(times, vel[istar, :, k0], 'ko-', label=id)
            axis[i,1].semilogx(times, vel[istar, :, k0], 'ko-', label=id)
            axis[i,2].loglog(times, vel[istar, :, k0], 'ko-', label=id)
            axis[i,3].semilogy(times, vel[istar, :, k0], 'ko-', label=id)
            # mark the ones used for fitting
            axis[i, 0].plot(times[:tmin], vel[istar, :tmin, k0], 'ow', alpha=0.7)
            axis[i, 1].semilogx(times[:tmin], vel[istar, :tmin, k0], 'ow', alpha=0.7)
            axis[i, 2].loglog(times[:tmin], vel[istar, :tmin, k0], 'ow', alpha=0.7)
            axis[i, 3].semilogy(times[:tmin], vel[istar, :tmin, k0], 'ow', alpha=0.7)
            axis[i, 0].plot(times[tmax:], vel[istar, tmax:, k0], 'ow', alpha=0.7)
            axis[i, 1].semilogx(times[tmax:], vel[istar, tmax:, k0], 'ow', alpha=0.7)
            axis[i, 2].loglog(times[tmax:], vel[istar, tmax:, k0], 'ow', alpha=0.7)
            axis[i, 3].semilogy(times[tmax:], vel[istar, tmax:, k0], 'ow', alpha=0.7)
        # reference CP
        axis[0,0].plot(times, vel[1, :, k0], 'd-', color='0.7', label=id_ref)
        axis[0,1].semilogx(times, vel[1, :, k0], 'd-', color='0.7', label=id_ref)
        axis[0,2].loglog(times, vel[1, :, k0], 'd-', color='0.7', label=id_ref)
        axis[0,3].semilogy(times, vel[1, :, k0], 'd-', color='0.7', label=id_ref)
        # fitted functions
        axis[1, 0].plot(times, func(times, *popt), 'r-', label='a+b*x^c (' + str(popt) + ')')
        axis[1, 1].semilogx(times, func(times, *popt), 'r-')
        axis[1, 2].loglog(times, func(times, *popt), 'r-')
        axis[1, 3].semilogy(times, func(times, *popt), 'r-')
        axis[1, 0].plot(times, func_log1(times, *popt_log1), color='limegreen', linewidth=3, label='a+b*log(x) (' + str(popt_log1) + ')')
        axis[1, 1].semilogx(times, func_log1(times, *popt_log1), color='limegreen', linewidth=3)
        axis[1, 2].loglog(times, func_log1(times, *popt_log1), color='limegreen', linewidth=3)
        axis[1, 3].semilogy(times, func_log1(times, *popt_log1), color='limegreen', linewidth=3)
        axis[1, 0].plot(times, func_log2(times, *popt_log2), color='darkgreen', label='a+b*log(c*x) (' + str(popt_log2) + ')')
        axis[1, 1].semilogx(times, func_log2(times, *popt_log2), color='darkgreen')
        axis[1, 2].loglog(times, func_log2(times, *popt_log2), color='darkgreen')
        axis[1, 3].semilogy(times, func_log2(times, *popt_log2), color='darkgreen')
        # axis[1, 0].plot(times[1:], func_log_romps_U(times[1:], *popt_log_romps), color='navy',
        #                 label='a*(1+bx)^(-1) (' + str(popt_log_romps) + ')')
        # axis[1, 1].semilogx(times[1:], func_log_romps_U(times[1:], *popt_log_romps), color='navy')
        # axis[1, 2].loglog(times[1:], func_log_romps_U(times[1:], *popt_log_romps), color='navy')
        # axis[1, 3].semilogy(times[1:], func_log_romps_U(times[1:], *popt_log_romps), color='navy')
        # for i in range(n_init):
        #     for j in range(n_init):
        #         for k in range(n_init):
        #             ij = i*n_init**2+j*n_init+k
        #             # print(i, j, k, p0[i,0], p0[j,1], p0[k,2])
        #             p1[ij, :], success = optimize.leastsq(errfunc, [p0[i,0], p0[j,1], p0[k,2]], args=(times[tmin:], r_av[istar, tmin:, k0]))
        #             axis[0,0].plot(times[tmin:], fitfunc1(p1[ij,:], times[tmin:]), "-", label='p='+str(p1[i,:]))  # Plot of the data and the fit
        #             axis[0,1].semilogx(times, fitfunc1(p1[ij,:], times), "-", label='p='+str(p1[i,:]))  # Plot of the data and the fit
        #             axis[0,2].loglog(times, fitfunc1(p1[ij,:], times), "-", label='p='+str(p1[i,:]))  # Plot of the data and the fit
        # # # axis[0, 1].set_title('radial spreading velocity (U_av)', fontsize=18)
        axis[0,0].set_title('scipy.optizmize.leastsq: a + b*x^c', fontsize=18)
        axis[1,0].set_title('scipy.optizmize.curve_fit: a + b*x^c, a+b*log(c), a+b*log(1+c)', fontsize=18)
        for i in range(2):
            axis[i, 0].set_xlabel('time [s]')
            axis[i, 1].set_xlabel('log(time) [s]')
            axis[i, 2].set_xlabel('log(time) [s]')
            axis[i, 3].set_xlabel('time [s]')
            axis[i, 0].set_ylabel('U  [m/s]')
            axis[i, 2].set_ylabel('U  [m/s]')
            axis[i, 0].set_ylim(0.,1e1)
            axis[i, 1].set_ylim(0.,1e1)
        axis[1,3].set_ylim(1.1e-1,1e1)
        axis[1,2].set_ylim(1.1e-1,1e1)
        axis[1,3].set_ylim(1.1e-1,1e1)
        plt.subplots_adjust(bottom=0.1, right=.95, left=0.07, top=0.9, wspace=0.25, hspace=0.9)
        axis[0, 0].legend(loc='upper left', bbox_to_anchor=(0.5, -0.2),
                          fancybox=True, shadow=True, ncol=3)
        axis[1, 0].legend(loc='upper left', bbox_to_anchor=(0.5, -0.2),
                          fancybox=True, shadow=True, ncol=3)
        fig.suptitle('CP rim velocity (' + id +', dx='+str(dx[0]) +'m)', fontsize=21)
        fig.savefig(os.path.join(path_out_figs, fig_name + '_' + id + '.png'))
        plt.close(fig)






    ''' all '''
    print('plotting all')
    cmap = cm_hsv
    fig, axis = plt.subplots(3, 4, sharex='none', figsize=(18, 16))
    for istar in range(n_params):
        col = colorlist[istar]
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[0]
        rstar = r_params[istar]
        tmin_plot = 3
        if rstar == 1000:
            id = id_ref
        else:
            id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        for i in range(3):
            axis[i, 0].plot(times, vel[istar, :, k0], 'o-', linewidth=1, color=col, label=id)
            axis[i, 1].semilogx(times, vel[istar, :, k0], 'o-', linewidth=1, color=col, label=id)
            axis[i, 2].loglog(times, vel[istar, :, k0], 'o-', linewidth=1, color=col, label=id)
            axis[i, 3].semilogy(times, vel[istar, :, k0], 'o-', linewidth=1, color=col, label=id)
            axis[i, 0].plot(times[:tmin], vel[istar, :tmin, k0], 'ow', alpha=0.7)
            axis[i, 1].semilogx(times[:tmin], vel[istar, :tmin, k0], 'ow', alpha=0.7)
            axis[i, 2].loglog(times[:tmin], vel[istar, :tmin, k0], 'ow', alpha=0.7)
            axis[i, 3].semilogy(times[:tmin], vel[istar, :tmin, k0], 'ow', alpha=0.7)
        axis[0, 0].plot(times[tmin_plot:], func(times[tmin_plot:], *f1[istar, :]), '-', linewidth=3, color=col,
                         label='a+b*x^c (' + str(f1[istar, :]) + ')')
        axis[0, 1].semilogx(times[tmin_plot:], func(times[tmin_plot:], *f1[istar, :]), '-', linewidth=3, color=col)
        axis[0, 2].loglog(times[tmin_plot:], func(times[tmin_plot:], *f1[istar, :]), '-', linewidth=3, color=col)
        axis[0, 3].semilogy(times[tmin_plot:], func(times[tmin_plot:], *f1[istar, :]), '-', linewidth=3, color=col)
        axis[1, 0].plot(times[tmin_plot:], func_log_romps_U(times[tmin_plot:], *f_log_romps[istar, :]), '-', linewidth=3, color=col,
                        label='b/(1+c*x) (' + str(f_log_romps[istar, :]) + ')')
        axis[1, 1].semilogx(times[tmin_plot:], func_log_romps_U(times[tmin_plot:], *f_log_romps[istar, :]), '-', linewidth=3, color=col)
        axis[1, 2].loglog(times[tmin_plot:], func_log_romps_U(times[tmin_plot:], *f_log_romps[istar, :]), '-', linewidth=3, color=col)
        axis[1, 3].semilogy(times[tmin_plot:], func_log_romps_U(times[tmin_plot:], *f_log_romps[istar, :]), '-', linewidth=3, color=col)
        axis[2, 0].plot(times[tmin_plot:], func_log2(times[tmin_plot:], *f_log2[istar, :]), '-', linewidth=3, color=col,
                        label='a+b*log(c*x) (' + str(f_log2[istar, :]) + ')')
        axis[2, 1].semilogx(times[tmin_plot:], func_log2(times[tmin_plot:], *f_log2[istar, :]), '-', linewidth=3, color=col)
        axis[2, 2].loglog(times[tmin_plot:], func_log2(times[tmin_plot:], *f_log2[istar, :]), '-', linewidth=3, color=col)
        axis[2, 3].semilogy(times[tmin_plot:], func_log2(times[tmin_plot:], *f_log2[istar, :]), '-', linewidth=3, color=col)
    for i in range(2,3):
        axis[i,0].set_xlabel('time [s]')
        axis[i,1].set_xlabel('log(time) [s]')
        axis[i,2].set_xlabel('log(time) [s]')
        axis[i,3].set_xlabel('time [s]')
    for i in range(3):
        axis[i,0].set_ylabel('U_av  [m/s]')
        axis[i,2].set_ylabel('log(U_av)  [m/s]')
        axis[i,1].set_xlim(0,3.6e3)
        axis[i,2].set_xlim(0,3.6e3)
        axis[i,0].set_ylim(0., vmax)
        axis[i,1].set_ylim(0., vmax)
        axis[i,2].set_ylim(2e-1, 1e1)
        axis[i,3].set_ylim(2e-1, 1e1)
    textprops = dict(facecolor='white', alpha=0.5, linewidth=0.)
    axis[0, 0].text(0.05, 0.95, 'Power Law: U=a+b*t^c', transform=axis[0, 0].transAxes, fontsize=16,
                    verticalalignment='top', bbox=textprops)
    axis[1, 0].text(0.05, 0.95, 'Romps: U=b/(1+c*t)', transform=axis[1, 0].transAxes, fontsize=16,
                    verticalalignment='top', bbox=textprops)
    axis[2, 0].text(0.05, 0.95, 'Logarithm: U=a+b*log(c*t)', transform=axis[2, 0].transAxes, fontsize=16,
                    verticalalignment='top', bbox=textprops)

    # # plt.subplots_adjust(bottom=0.3, right=.95, left=0.07, top=0.9, wspace=0.25)
    plt.subplots_adjust(bottom=0.1, right=.97, left=0.07, top=0.95, wspace=0.25, hspace=0.35)
    axis[0, 0].legend(loc='upper left', bbox_to_anchor=(-0.1, -0.07),
                      fancybox=True, shadow=False, ncol=3)
    axis[1, 0].legend(loc='upper left', bbox_to_anchor=(-0.1, -0.07),
                      fancybox=True, shadow=False, ncol=3)
    axis[2, 0].legend(loc='upper left', bbox_to_anchor=(-0.1, -0.15),
                      fancybox=True, shadow=False, ncol=3)
    fig.suptitle('radial spreading velocity (' + fig_name[-4:] + ') (dx=' + str(dx[0]) + 'm)', fontsize=21)
    fig.savefig(os.path.join(path_out_figs, fig_name + '_all.png'))
    plt.close(fig)




    ''' comparison fitting '''
    fig, axis = plt.subplots(2, 3, sharex='all', figsize=(18, 7))
    for istar in range(n_params):
        if rstar == 1000:
            id = id_ref
        else:
            id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        popt = f1[istar, :]
        popt_log2 = f_log2[istar, :]
        popt_log_romps = f_log_romps[istar, :]
        diff1 = func(times, *popt) - func_log2(times, *popt_log2)
        diff2 = func(times, *popt) - func_log_romps_U(times, *popt_log_romps)
        axis[0, 0].plot(times, diff1, '-', linewidth=2, color=colorlist[istar], label=id)
        axis[0, 1].plot(times, diff2, '-', linewidth=2, color=colorlist[istar], label=id)
        axis[0, 2].plot(times, diff2, '-', linewidth=2, color=colorlist[istar], label=id)
        axis[1, 0].plot(times[tmin:], diff1[tmin:], '-', linewidth=2, color=colorlist[istar], label=id)
        axis[1, 1].plot(times[tmin:], diff2[tmin:], '-', linewidth=2, color=colorlist[istar], label=id)
    axis[0, 0].set_title('a+b*x^c  -  a+b*log(x)')
    axis[0, 1].set_title('a+b*x^c  -  b/(1+c*x)')
    axis[0, 2].set_title('a+b*x^c  -  b/(1+c*x)')
    # axis[0, 2].set_ylim(-1e2, 1e2)
    axis[0,0].legend(loc='best')

    for i in range(3):
        axis[-1, i].set_xlabel('time')
    for i in range(2):
        axis[i, 0].set_ylabel('diff in radius [m]')
    fig.suptitle('comparison of fittings', fontsize=21)
    fig.savefig(os.path.join(path_out_figs, fig_name + '_comparison.png'))
    plt.close(fig)

    return f1, f_log1, f_log2, f_log_romps




def plot_parameters(f1, f_log1, f_log2, f_log_romps, var):
    ''' plot parameters '''
    print('plot parameters')
    PE_range_log = [-1, 0, 1, 2, 3]
    PE_range = [0.5, 1, 2, 4, 8]

    fig, axis = plt.subplots(3, 3, sharex='all', figsize=(18, 7))
    axis[0, 0].plot(PE_range, f1[:, 0], '-o')
    axis[0, 1].plot(PE_range, f1[:, 1], '-o')
    axis[0, 2].plot(PE_range, f1[:, 2], '-o')
    axis[1, 0].plot(PE_range, f_log2[:, 0], '-o')
    axis[1, 1].plot(PE_range, f_log2[:, 1], '-o')
    # axis[1, 2].plot(PE_range, f_log2[:, 2], '-o')
    axis[2, 0].plot(PE_range, f_log_romps[:, 0], '-o')
    axis[2, 1].plot(PE_range, f_log_romps[:, 1], '-o')
    axis[2, 2].plot(PE_range, f_log_romps[:, 2], '-o')
    # axis[1,2].plot(f_log2[:,2], '-o')
    axis[0, 0].set_title('a+b*x^c: a')
    axis[0, 1].set_title('a+b*x^c: b')
    axis[0, 2].set_title('a+b*x^c: c')
    axis[1, 0].set_title('a+b*log(c*x): a')
    axis[1, 1].set_title('a+b*log(c*x): b')
    axis[1, 2].set_title('a+b*log(c*x): c')
    if var[0] in ['r','R']:
        axis[2, 0].set_title('a+b*log(1+c*x): a')
        axis[2, 1].set_title('a+b*log(1+c*x): b')
        axis[2, 2].set_title('a+b*log(1+c*x): c')
    elif var[0] == 'U' or var == 'dRdt':
        axis[2, 0].set_title('b/(1+c*x): a')
        axis[2, 1].set_title('b/(1+c*x): b')
        axis[2, 2].set_title('b/(1+c*x): c')
    for i in range(3):
        axis[2, i].set_xlabel('PE/PE_ref')
    fig.suptitle('fitting parameters ' + var, fontsize=21)
    fig.savefig(os.path.join(path_out_figs, 'parameters_' + var + '.png'))
    plt.close(fig)

    fig, axis = plt.subplots(3, 3, sharex='all', figsize=(18, 7))
    axis[0, 0].plot(PE_range_log, f1[:, 0], '-o')
    axis[0, 1].plot(PE_range_log, f1[:, 1], '-o')
    axis[0, 2].plot(PE_range_log, f1[:, 2], '-o')
    axis[1, 0].plot(PE_range_log, f_log2[:, 0], '-o')
    axis[1, 1].plot(PE_range_log, f_log2[:, 1], '-o')
    axis[1,2].plot(PE_range_log, f_log2[:,2], '-o')
    axis[2, 0].plot(PE_range_log, f_log_romps[:, 0], '-o')
    axis[2, 1].plot(PE_range_log, f_log_romps[:, 1], '-o')
    axis[2, 2].plot(PE_range_log, f_log_romps[:, 2], '-o')
    axis[0, 0].set_title('a+b*x^c: a')
    axis[0, 1].set_title('a+b*x^c: b')
    axis[0, 2].set_title('a+b*x^c: c')
    axis[1, 0].set_title('a+b*log(x): a')
    axis[1, 1].set_title('a+b*log(x): b')
    axis[1, 1].set_title('a+b*log(c*x): c')
    # axis[1,2].set_title('a+b*x^c: c')
    if var[0] in ['r','R']:
        axis[2, 0].set_title('a+b*log(1+c*x): a')
        axis[2, 1].set_title('a+b*log(1+c*x): b')
        axis[2, 2].set_title('a+b*log(1+c*x): c')
    elif var[0] == 'U' or var == 'dRdt':
        axis[2, 0].set_title('b/(1+c*x): a')
        axis[2, 1].set_title('b/(1+c*x): b')
        axis[2, 2].set_title('b/(1+c*x): c')
    for i in range(3):
        axis[2, i].set_xlabel('log2(PE/PE_ref)')
    fig.suptitle('fitting parameters ' + var, fontsize=21)
    fig.savefig(os.path.join(path_out_figs, 'parameters_' + var + '_log.png'))
    plt.close(fig)
    return


# ----------------------------------------------------------------------
def test_fit(var_name, data, f1, f_log2, f_log_romps,
             n_params, dTh_params, r_params, z_params, id_ref, tmin, colorlist, figname):
    print('')
    print('test fitting')
    print(times, times[-1])

    err1 = np.ndarray((n_params, len(times)))
    err2 = np.ndarray((n_params, len(times)))
    err3 = np.ndarray((n_params, len(times)))
    diff1 = np.ndarray((n_params, len(times)))
    diff2 = np.ndarray((n_params, len(times)))
    for istar in range(n_params):
        popt = f1[istar,:]
        popt_log2 = f_log2[istar,:]
        popt_log_romps = f_log_romps[istar,:]
        diff1[istar,:] = func(times, *popt) - func_log2(times, *popt_log2)
        err1[istar, :] = np.sqrt((data[istar, :] - func(times, *popt)) ** 2)
        err2[istar, :] = np.sqrt((data[istar, :] - func_log2(times, *popt_log2)) ** 2)
        if var_name == 'R':
            diff2[istar,:] = func(times, *popt) - func_log_romps_r(times, *popt_log_romps)
            err3[istar, :] = np.sqrt((data[istar, :] - func_log_romps_r(times, *popt_log_romps)) ** 2)
        elif var_name == 'U':
            diff2[istar,:] = func(times, *popt) - func_log_romps_U(times, *popt_log_romps)
            err3[istar, :] = np.sqrt((data[istar, :] - func_log_romps_U(times, *popt_log_romps)) ** 2)
    err1_mean = np.mean(err1[:,tmin:-2], axis=1)
    err2_mean = np.mean(err2[:,tmin:-2], axis=1)
    err3_mean = np.mean(err3[:,tmin:-2], axis=1)




    fig, axis = plt.subplots(3, 3, sharex='all', figsize=(18, 12))
    if len(dTh_params) == 1:
        dTh = dTh_params[0]
    else:
        dTh = dTh_params[istar]
    for istar in range(n_params):
        rstar = r_params[istar]
        zstar = z_params[0]
        if rstar == 1000:
            id = id_ref
        else:
            id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        axis[2, 0].plot([0, times[-1]], [err1_mean, err1_mean], '-', color=colorlist[istar], linewidth=0.5)
        axis[2, 1].plot([0, times[-1]], [err2_mean, err2_mean], '-', color=colorlist[istar], linewidth=0.5)
        axis[2, 2].plot([0, times[-1]], [err3_mean, err3_mean], '-', color=colorlist[istar], linewidth=0.5)
        axis[0, 1].plot(times, diff1[istar,:], '-', linewidth=2, color=colorlist[istar], label=id)
        axis[0, 2].plot(times, diff2[istar,:], '-', linewidth=2, color=colorlist[istar], label=id)
        axis[1, 1].plot(times[tmin:], diff1[istar,tmin:], '-', linewidth=2, color=colorlist[istar], label=id)
        axis[1, 2].plot(times[tmin:], diff2[istar,tmin:], '-', linewidth=2, color=colorlist[istar], label=id)
        axis[2, 0].plot(times[tmin:-2], err1[istar,tmin:-2], '-', linewidth=2, color=colorlist[istar], label=id)
        axis[2, 1].plot(times[tmin:-2], err2[istar,tmin:-2], '-', linewidth=2, color=colorlist[istar], label=id)
        axis[2, 2].plot(times[tmin:-2], err3[istar,tmin:-2], '-', linewidth=2, color=colorlist[istar], label=id)
        # axis[2, 0].fill_between(times, err1_mean*np.ones(len(times)),
        #                  facecolor='k', alpha=0.2, linewidth=0)
    axis[0, 1].set_title('a+b*x^c  -  a+b*log(c*x)')
    axis[0, 2].set_title('a+b*x^c  -  a+b*log(1+c*x)')
    axis[2, 0].set_title('err(a+b*x^c)')
    axis[2, 1].set_title('err(a+b*log(c*x)')
    if var_name == 'R':
        axis[2, 2].set_title('err(a+b*log(1+c*x))')
    elif var_name == 'U':
        axis[2, 2].set_title('err(b/(1+c*x)')
    axis[0, 1].legend(loc='best')

    for i in range(2):
        axis[-1, i].set_xlabel('time')
    for i in range(2):
        axis[i, 0].set_ylabel('diff in dRdt [m/s]')
    fig.suptitle('Comparison & errors of fittings ('+var_name+')', fontsize=21)
    fig.savefig(os.path.join(path_out_figs, figname))
    plt.close(fig)

    return


# ----------------------------------------------------------------------
def plot_vel(r_av, U_rad_av, dTh_params, z_params, r_params, n_params, k0, fig_name):

    fig, (ax0, ax1) = plt.subplots(1, 2, sharex='all', figsize=(11, 5))
    for istar in range(n_params):
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[0]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        ax0.plot(times, r_av[istar, :, k0], 'o-', label=id)
        ax1.plot(times, U_rad_av[istar, :, k0], 'o-', label=id)
    ax0.set_xlabel('times [s]')
    ax1.set_xlabel('times [s]')
    ax0.set_ylabel('r_av')
    ax1.set_ylabel('U_rad_av')
    ax1.legend()
    fig.suptitle('CPs (r=1km)')
    fig.tight_layout()
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, wspace=0.25)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return

# ----------------------------------------------------------------------
def read_in_txtfile(fullpath_in):
    f = open(fullpath_in+'/coldpool_tracer_out.txt', 'r')
    # f = open(DIR+EXPID+'/'+child+'/output/irt_tracks_output_pure_sort.txt', 'r')
    lines = f.readlines()
    column =lines[0].split()
    # print column
    start_c = int(column[1])   #get the first timestep when tracking starts
    print 'precip starts at timestep', start_c
    column = lines[-1].split()
    end_c = int(column[1])
    print 'last timestep in precip tracking', end_c
    for line in lines:
        columns = line.split()
        tist = (int(columns[0]))    # timestep
        # age = (int(columns[1]))    # age = timestep - tstart
        tracerID = (int(columns[2]))
        cpID = (int(columns[3]))
        # position tracer: columns 4,5
        # rounded position tracer: columns 6, 7
        dist = (float(columns[8])) # distance of tracer from COG
        # angle = (float(columns[9])) # angle of vector (position tracer)-COG
        # u =(float(columns[10]))  # u- wind component
        # v =(float(columns[11]))  # v- wind component
        vel_rad =(float(columns[12]))  # radial wind component
        # vel_tang =(float(columns[13]))  # tangential wind component
        # dist_x = (float(columns[14])) # x-component of distance from COG
        # dist_y = (float(columns[15])) # x-component of distance from COG
        COGx =(float(columns[16]))  # center of CP (and initialized circle)
        COGy = (float(columns[17]))

        if tracerID == 1:
            print('t', tist)
            print('vel_rad', vel_rad, dist)

    f.close()

    return


def get_radius_vel(fullpath_in, t0, cp_id, n_tracers, n_cps):
    # print('in', fullpath_in)
    f = open(fullpath_in + '/coldpool_tracer_out.txt', 'r')
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


def get_number_tracers(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in + '/coldpool_tracer_out.txt', 'r')
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
    f = open(fullpath_in + '/coldpool_tracer_out.txt', 'r')
    lines = f.readlines()
    # count = 0
    # # while CP age is 0 and CP ID is 1
    # while (int(lines[count].split()[0]) == 1):
    #     count += 1
    # cp_number = int(lines[count-1].split()[3])
    cp_number = int(lines[-1].split()[3])

    f.close()

    return cp_number


# ----------------------------------------------------------------------

def set_input_parameters(args):
    print('--- set input parameters ---')
    global case_name
    global path_root, path_out_figs
    global times

    path_root = args.path_root
    path_out_figs = os.path.join(path_root, 'figs_tracers')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print 'path_root: ', path_root
    print('path figs: ', path_out_figs)
    print ''

    dTh = args.dTh
    if args.zparams:
        z_params = args.zparams
    else:
        z_params = [1000]  # run5, 6
    if args.rparams:
        r_params = args.rparams
    else:
        r_params = [500, 1000, 1100, 1600, 2300]
    n_params = len(r_params)


    print('dTh: ', dTh)
    print('z*: ', z_params)
    print('r*: ', r_params)
    print('n_params: ', n_params)

    case_name = args.casename
    id0 = 'dTh' + str(dTh) + '_z' + str(z_params[0]) + '_r' + str(r_params[0])
    nml = simplejson.loads(open(os.path.join(path_root, id0, case_name + '.in')).read())
    global nx, ny, nz, dx, dV, gw
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = np.zeros(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx[0] * dx[1] * dx[2]

    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = np.int(100)
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = np.int(3600)
    times = np.arange(tmin, tmax + 100, 100)
    # times = [np.int(name[:-3]) for name in files]
    times.sort()
    print('times', times)
    print ''

    return dTh, z_params, r_params



# ----------------------------------------------------------------------
# use optimize.curve_fit
def func(x, a, b, c):
    return a + b * x ** c
def func2(x, a, b, c):
    return a * (x / b) ** c  # a = y(x=b); a=R0=R(t=t0), t0: start time for fitting
def func_log1(x, a, b, c):
    return a + b * np.log(x)
def func_log2(x, a, b, c):
    return a + b * np.log(c*x)
def func_log_romps_r(x, a, b, c):
    return a + b * np.log(1 + x * c)
def func_log_romps_U(x, a, b, c):
    return b * (1 + x * c)**(-1)
# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()
