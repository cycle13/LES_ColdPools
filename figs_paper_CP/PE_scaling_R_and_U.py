import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os
import scipy
from scipy import optimize

execfile('settings.py')


label_size = 18
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['axes.labelsize'] = 24
# plt.rcParams['text.usetex'] = 'true'
# plt.rcParams["legend.facecolor"] = 'w'

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    case_name = 'ColdPoolDry_single_3D'
    path_root = '/nbi/ac/conv4/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/single_3D_noise/run5_PE_scaling_dx100m/'
    # path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single'
    path_out_figs = '/nbi/home/meyerbe/paper_CP'
    print('path figures: '  + path_out_figs)
    dTh, z_params, r_params = set_input_parameters(args, path_root, case_name)

    # reference case: dTh3_z1000_r1000
    rstar_ref = 1000
    zstar_ref = 1000
    id_ref = 'dTh3_z' + str(zstar_ref) + '_r' + str(rstar_ref)
    path_ref = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/' + id_ref
    dt_fields = 100
    cp_id = 1  # circle ID that is used for statistics

    id0 = 'dTh' + str(dTh) + '_z' + str(z_params[0]) + '_r' + str(r_params[0])
    fullpath_in = os.path.join(path_root, id0, 'tracer_k' + str(k0), 'output')
    n_tracers = get_number_tracers(fullpath_in)
    n_cps = get_number_cps(fullpath_in)
    print('number of CPs: ', n_cps)
    print('number of tracers per CP: ', n_tracers)
    print ''


    colorlist5 = ['maroon', 'indianred', 'orange', 'darkcyan', 'navy']

    # --------------------------------------
    ''' (a) read in data from tracer output (text-file)'''
    dist_av = np.zeros((n_params, nt))
    r_av = np.zeros((n_params, nt))  # absolute radius
    r_av_abs = np.zeros((n_params, nt))  # radius minus initial radius (r(t) - r(t=0))
    drdt_av = np.zeros((n_params, nt))
    U_rad_av = np.zeros((n_params, nt))

    print('--- reading in tracer data')
    for istar in range(n_params):
        zstar = z_params[0]
        rstar = r_params[istar]
        if rstar == 1000:
            id = id_ref
            fullpath_in = os.path.join(path_ref, 'tracer_k' + str(k0), 'output')
        else:
            id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
            fullpath_in = os.path.join(path_root, id, 'tracer_k' + str(k0), 'output')
        print('id', id)
        # read_in_txtfile(fullpath_in)
        for it, t0 in enumerate(times):
            print('---t0: ' + str(t0) + '---', it)
            dist_av[istar, it], U_rad_av[istar, it] = get_radius_vel(fullpath_in, it, cp_id, n_tracers, n_cps)
        r_av = dist_av * dx[0]
        r_av_abs[istar, :] = r_av[istar, :] - rstar
        r_av_abs[istar, :] = r_av[istar, :] - r_av[istar, 0]
        print('rstar: ', dx[0], rstar, r_av[istar, 0], r_av_abs[istar, 0])
        for it, t0 in enumerate(times[1:]):
            drdt_av[:, it] = 1. / dt_fields * (r_av[:, it] - r_av[:, it - 1])
        print('')
    print('')

    # --------------------------------------
    tmin_r = 6
    tmin_u = 7
    # --------------------------------------
    fig_name = 'test_R.png'
    test_R_figure(r_av, tmin_r, r_params, dt_fields, path_out_figs, fig_name)


    def fit_r_powerlaw(t, m):
        return m * (np.log(t / t0_r))
    def fit_r_test_istar0(t, m):
        istar0 = 0
        R0_log = np.log(r_av[istar0, tmin_r])  # log(R[t0])
        # print('fitting 0:', r_av.shape, R0_log.shape, t.shape, t0_r)
        return R0_log + m * (np.log(t / t0_r))
    def fit_r_test_istar1(t, m):
        istar0 = 1
        R0_log = np.log(r_av[istar0, tmin_r])  # log(R[t0])
        # print('fitting 1:', r_av.shape, R0_log.shape, t.shape, t0_r)
        return R0_log + m * (np.log(t / t0_r))
    def fit_r_test_istar2(t, m):
        istar0 = 2
        R0_log = np.log(r_av[istar0, tmin_r])  # log(R[t0])
        # print('fitting 2:', r_av.shape, R0_log.shape, t.shape, t0_r)
        return R0_log + m * (np.log(t / t0_r))
    def fit_r_test_istar3(t, m):
        istar0 = 3
        R0_log = np.log(r_av[istar0, tmin_r])  # log(R[t0])
        # print('fitting 3:', r_av.shape, R0_log.shape, t.shape, t0_r)
        return R0_log + m * (np.log(t / t0_r))
    def fit_r_test_istar4(t, m):
        istar0 = 4
        R0_log = np.log(r_av[istar0, tmin_r])  # log(R[t0])
        # print('fitting 4:', r_av.shape, R0_log.shape, t.shape, t0_r)
        return R0_log + m * (np.log(t / t0_r))


    # logarithmic fit by Romps
    def fit_log_romps(t, a, b, c):
        return a + b*np.log(1+c*t)





    # --------------------------------------
    print('Fitting R:')
    fig_name = 'R_scaling.png'

    t0_r = np.double(tmin_r * dt_fields)
    print('       t0:' + str(t0_r), tmin_r)

    # R00 = r_av[0, tmin_r]
    R0_log = np.log(r_av[:, tmin_r])  # log(R[t0])
    # dR = 1. / r_av[:, tmin_r]
    r_av_ens_mean = np.average(r_av[:, :], axis=0)  # ensemble mean
    dR_ens_mean = 1. / r_av_ens_mean[tmin_r]

    m0 = -0.5
    istar = 0
    # popt, pcov = curve_fit(func, xdata, ydata)
    m_fitted = np.zeros(5, dtype=np.double)
    for istar in range(n_params):
        m_fitted[istar], pcov = optimize.curve_fit(fit_r_test, times[1:-1], np.log(r_av[istar, 1:-1]) - R0_log[istar])
    # m_fitted[1], pcov = optimize.curve_fit(fit_r_test, times[1:-1], np.log(r_av[istar, 1:-1]) - R0_log[istar])
    # m_fitted[2], pcov = optimize.curve_fit(fit_r_test, times[1:-1], np.log(r_av[istar, 1:-1]) - R0_log[istar])
    # m_fitted[3], pcov = optimize.curve_fit(fit_r_test, times[1:-1], np.log(r_av[istar, 1:-1]) - R0_log[istar])
    # m_fitted[4], pcov = optimize.curve_fit(fit_r_test, times[1:-1], np.log(r_av[istar, 1:-1]) - R0_log[istar])
    print('Parameters: ', m_fitted)
    print('')

    collist = ['b', 'k', 'g', 'r', 'cyan']
    plot_R_1x2(r_av, r_av_ens_mean, dR_ens_mean, R0_log, t0_r, m_fitted,
               r_params, dt_fields, collist, path_out_figs, fig_name)
>>>>>>> db4da1ab62cb28f3964774546f8423898ee4d545



    # --------------------------------------
    print('Fitting R:')
    fig_name = 'R_scaling.png'
    #
    t0_r = np.double(tmin_r * dt_fields)
    # print('       t0:' + str(t0_r), tmin_r)
    #
    R0 = r_av[:, tmin_r]
    R0_log = np.log(r_av[:, tmin_r])  # log(R[t0])
    # dR = 1. / r_av[:, tmin_r]
    r_av_ens_mean = np.average(r_av[:, :], axis=0)  # ensemble mean
    dR_ens_mean = 1. / r_av_ens_mean[tmin_r]

    r_av_rel_mean = 0.
    for istar in range(n_params):
        r_av_rel_mean += r_av[istar, :] / R0[istar]
    r_av_rel_mean /= n_params


    m0 = -0.5
    istar = 0
    # popt, pcov = curve_fit(func, xdata, ydata)
    m_fitted = np.zeros(5, dtype=np.double)
    for istar in range(n_params):
        m_fitted[istar], pcov = optimize.curve_fit(fit_r_powerlaw, times[1:-1], np.log(r_av[istar, 1:-1]) - R0_log[istar])
    print('Parameters: ', m_fitted)
    print('')

    # collist = ['b', 'k', 'g', 'r', 'cyan']
    # plot_R_1x2(r_av, r_av_ens_mean, dR_ens_mean, R0_log, t0_r, m_fitted,
    #            r_params, dt_fields, collist, path_out_figs, fig_name)

    # --------------------------------------

    params_Romps = []
    # params_Romps = np.zeros(5, dtype=np.double)
    # tmin_r = 4
    p0 = np.zeros((3))
    p0[0] = 0.
    p0[1] = 1.
    p0[2] = .5
    p0_log = np.copy(p0)
    # p0_log = np.zeros((3))
    p0_log[0] = np.log(p0[0])
    # p0_log[1] = 1.
    # p0_log[2] = .5
    for istar in range(n_params):
        popt_log_romps, pcov = optimize.curve_fit(fit_log_romps, times[tmin_r:], r_av[istar, tmin_r:], p0[:])
        params_Romps.append(popt_log_romps)
    params_Romps_mean, pcov = optimize.curve_fit(fit_log_romps, times[tmin_r:], r_av_rel_mean[tmin_r:], p0[:])

    collist = colorlist5
    fig_name = 'R_scaling_Romps.png'
    plot_Romps(r_av, r_av_ens_mean, dR_ens_mean, R0_log, tmin_r, params_Romps, params_Romps_mean,
               r_params, dt_fields, collist, path_out_figs, fig_name)

    fig_name = 'R_scaling_pow_log.png'
    plot_R_1x3(r_av, r_av_ens_mean, dR_ens_mean, R0_log, tmin_r, m_fitted,
               params_Romps, params_Romps_mean,
               r_params, dt_fields, collist, path_out_figs, fig_name)


    # --------------------------------------
    # fig_name = 'R_U_scaling.png'
    #
    # t0_r = np.double(tmin_r * dt_fields)
    # t0_u = np.double(tmin_u * dt_fields)
    #
    # R00 = r_av[0, tmin_r]
    # R0_log = np.log(r_av[:, tmin_r])                 # log(R[t0])
    # dR = 1. / r_av[:, tmin_r]
    # r_av_ens_mean = np.average(r_av[:, :], axis=0)   # ensemble mean
    # dR_ens_mean = 1. / r_av_ens_mean[tmin_r]
    #
    #
    # U00 = drdt_av[0,tmin_r]
    # dU_log = np.log(drdt_av[:, tmin_r])
    # dU = 1. / drdt_av[:, tmin_r]
    # U_ens_mean = np.average(drdt_av[:, :], axis=0)
    # dU_ens_mean = 1. / U_ens_mean[tmin_r]
    #
    # plot_R_U_2x2(r_av, r_av_ens_mean, dR_ens_mean, R0_log, t0_r,
    #              drdt_av, U_ens_mean, dU_ens_mean, dU_log,
    #              r_params, dt_fields, path_out_figs, fig_name)

    return


# --------------------------------------
def fit_r(R0, m, times, t0):
    return np.log(R0) + m * (np.log(times / t0))
    # return m * (np.log(times / t0))

# def fit_r_test(m, times, t0):
#     R0_log = np.log(r_av[:, tmin_r])  # log(R[t0])
#     return np.log(R0) + m * (np.log(times / t0))
#     # return m * (np.log(times / t0))
# --------------------------------------

def plot_Romps(r_av, r_av_ens_mean, dR_ens_mean, R0_log, tmin_r, params_fitted, params_fitted_mean,
               r_params, dt_fields, collist, path_out_figs, fig_name):
    t0_r = np.double(tmin_r * dt_fields)
    R0 = r_av[:, tmin_r]

    fig, axes = plt.subplots(3, 3, sharex='none', figsize=(18, 16))
    ax0 = axes[0,0]
    ax1 = axes[0,1]
    ax2 = axes[0,2]

    ax4 = axes[1,0]

    ax7 = axes[2,0]
    ax8 = axes[2,1]
    ax9 = axes[2,2]

    ax10 = axes[1,2]

    ax0.plot([t0_r, t0_r], [0, 13e3], '.5', linewidth=.5)
    ax1.semilogy(t0_r*np.ones(5), np.logspace(1e2, 5e4, 5), '.5', linewidth=.5)
    # ax2.loglog(t0_r*np.ones(5), np.logspace(1e2, 5e4, 5), '.5', linewidth=.5)
    ax4.plot([t0_r, t0_r], [0, 13e3], '.5', linewidth=.5)
    # ax7.plot([t0_r, t0_r], [-2e3, 8e3], '.5', linewidth=.5)
    # ax8.semilogy([t0_r, t0_r], [0, 1e2], '.5', linewidth=.5)
    # # ax9.plot([0, 0], [-1, 1.], '.5', linewidth=.5)
    # # ax9.plot([-.5, 2.], [0, 0], '.5', linewidth=.5)
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar]) + 'm'
        ax0.plot(times[:-1], r_av[istar, :-1], 'o-', color=collist[istar], label=lbl)
        ax1.semilogy(times[1:-1], r_av[istar, 1:-1], 'o-', color=collist[istar])
        ax2.loglog(times[1:-1], r_av[istar, 1:-1], 'o-', color=collist[istar])
        ax4.plot(times[1:-1], r_av[istar, 1:-1] - R0[istar], 'o-', color=collist[istar])
        ax7.plot(times[1:-1], r_av[istar, 1:-1]/R0[istar], 'o-', color=collist[istar])
        ax8.semilogy(times[1:-1], r_av[istar, 1:-1]/R0[istar], 'o-', color=collist[istar])
        # ax9.loglog(times[1:-1] / t0_r, r_av[istar, 1:-1]/R0[istar], '-', color=collist[istar])
        ax9.loglog(times[1:-1], r_av[istar, 1:-1]/R0[istar], 'o-', color=collist[istar])

    ''' mean values'''
    ax10.plot(np.log(times[3:-1] / t0_r), np.log(r_av_ens_mean[3:-1] * dR_ens_mean), '-', color='k', linewidth=3,
             label='ensemble mean')
    [a, b, c] = params_fitted_mean
    r_av_rel_mean = 0.
    for istar in range(n_params):
        r_av_rel_mean += r_av[istar, :] / R0[istar]
    r_av_rel_mean /= n_params
    ax8.semilogy(times[1:-1], r_av_rel_mean[1:-1], 'k', linewidth=3, label='ensemble mean <R_i/R_i0>')
    lbl = 'a='+str(np.round(a,1))+',b='+str(np.round(b,1))+',c='+str(np.round(c, 4))
    ax8.semilogy(times[1:-1], a+b*np.log(1+c*times[1:-1]), 'g-', linewidth=3,
                 label='a+b*log(1+c*t) (fit on ensemble mean), ' + lbl)

    a = np.zeros(n_params)
    b = np.zeros(n_params)
    c = np.zeros(n_params)
    for istar,rstar in enumerate(r_params):
        [a[istar],b[istar],c[istar]] = params_fitted[istar]

        lbl = 'a='+str(np.round(a[istar],1))+',b='+str(np.round(b[istar],1))+',c='+str(np.round(c[istar], 4))
        if istar == 0:
            ax0.plot(times[1:-1], a[istar]+b[istar]*np.log(1+c[istar]*times[1:-1]), 'r-', label='R=a+b*log(1+cx)')
        else:
            ax0.plot(times[1:-1], a[istar]+b[istar]*np.log(1+c[istar]*times[1:-1]), 'r-')
        # ax1.plot(times[1:-1], a[istar]+b[istar]*np.log(1+c[istar]*times[1:-1]), 'r-')
        # ax3_.semilogy(times[1:-1], a[istar]+b[istar]*np.log(1+c[istar]*times[1:-1]), 'g-', label=lbl)
        ax1.semilogy(times[1:-1], (a[istar]+b[istar]*np.log(1+c[istar]*times[1:-1])), 'r-', label=lbl)
        ax2.loglog(times[1:-1], (a[istar]+b[istar]*np.log(1+c[istar]*times[1:-1])), 'r-', label=lbl)
        ax4.plot(times[1:-1], (a[istar]+b[istar]*np.log(1+c[istar]*times[1:-1]))-R0[istar], 'r-', label=lbl+'-R0')
    #     ax7.plot(times[1:-1], (a[istar]+b[istar]*np.log(1+c[istar]*times[1:-1]))/ R0[istar], 'r-', label=lbl)
    #
    #     ax7.plot(times[1:-1], (b[istar]*np.log(1+c[istar]*times[1:-1])) / R0[istar], 'g-')
    #     ax8.semilogy(times[1:-1], (b[istar]*np.log(1+c[istar]*times[1:-1])) / R0[istar], 'g-', label='b='+str(np.round(b[istar],1))+',c='+str(np.round(c[istar], 4)))
    # # print('a=' + str(np.round(a, 1)) + ', a/R0=' + str(np.round(a / R0, 1)))
    # # # b = np.mean(b)
    # # # c = np.mean(c)
    # # # for istar,rstar in enumerate(r_params):
    # # #     [a,aux1,aux2] = params_fitted[istar]
    # # #     ax7.plot(times[1:-1], (a + b * np.log(1 + c * times[1:-1])) / R0[istar], 'g-')
    # #
    # #
    # # # m = 0.6
    # # print('---- m: ', params_fitted)
    # m = np.mean(params_fitted)
    # ax9.semilogx(times[3:-1] / t0_r, m * (np.log(times[3:-1] / t0_r)), '-r', label='y=m*x, m=' + str(m))
    # # # for istar, m_ in enumerate(m_fitted):
    # # #     ax9.plot(np.log(times[1:-1] / t0_r), m_ * np.log(times[1:-1] / t0_r), color=collist[istar], linewidth=1,
    # # #              label='y=m*x, m=' + str(np.round(m_, 2)))
    # # a = 0
    # # b = np.mean(b)
    # # c = np.mean(c)
    # # lbl = 'a=0, b=' + str(np.round(b, 1)) + ',c=' + str(np.round(c,4))
    # # ax7.plot(times[1:-1], (b*np.log(1+c*times[1:-1])) / R0[istar], 'g--', label=lbl)
    # # ax8.semilogy(times[1:-1], (b*np.log(1+c*times[1:-1])) / R0[istar], 'g--', label=lbl)
    #
    # # ax9.plot(np.log(times[3:-1] / (t0_r - 3)), m * (np.log(times[3:-1] / (t0_r - 3))), '-b', label='m=' + str(m))
    # # ax9.plot(np.log(times[3:-1]/(t0_r-3)), m * (np.log(times[3:-1] / (t0_r-3))), '-b', label='m='+str(m))
    #
    #
    ax0.set_xlabel('t  [min]')
    ax0.set_ylabel('R  [km]')
    ax1.set_xlabel('t  [min]')
    ax1.set_ylabel('log(R)  [-]')
    ax2.set_xlabel('log(t)  [min]')
    ax2.set_ylabel('log(R)  [-]')
    ax7.set_xlabel('t  [min]')
    ax7.set_ylabel('R/R0  [-]')
    ax8.set_xlabel('t  [min]')
    ax8.set_ylabel('log(R/R0)  [-]')
    ax4.set_xlabel('t  [min]')
    ax4.set_ylabel('R - R0  [m]')
    ax9.set_xlabel('log(t/t0)  [-]')
    ax9.set_ylabel('log(R/R0)  [-]')
    leg = ax0.legend(bbox_to_anchor=(.012, .99), loc='upper left', fontsize=15,
                     borderpad=.1)  # , edgecolor='w')#, facecolor='w')
    leg.get_frame().set_edgecolor('w')
    leg = ax4.legend(bbox_to_anchor=(.012, .99), loc='upper left', fontsize=10, fancybox=False, shadow=False, ncol=1,
                     borderpad=.1)  # , frameon=False)
    leg.get_frame().set_edgecolor('w')
    # leg = ax7.legend(bbox_to_anchor=(.01, .99), loc='upper left', fontsize=12, fancybox=False, shadow=False, ncol=1,
    #                  borderpad=.1)  # , frameon=False)
    # leg.get_frame().set_edgecolor('w')
    leg = ax8.legend(bbox_to_anchor=(.01, .99), loc='upper left', fontsize=12, fancybox=False, shadow=False, ncol=1,
                     borderpad=.1)  # , frameon=False)
    leg.get_frame().set_edgecolor('w')
    # # leg = ax9.legend(bbox_to_anchor=(.012, .99), loc='upper left', fontsize=15, fancybox=False, shadow=False, ncol=1,
    # #                  borderpad=.1)  # , frameon=False)
    # # leg.get_frame().set_edgecolor('w')
    ax1.set_ylim(6e2, 2e4)
    ax2.set_xlim(2e2, 5e3)
    ax2.set_ylim(6e2, 2e4)
    ax7.set_ylim(0, 3)
    ax8.set_ylim(0, 3)
    ax9.set_xlim(2e2, 5e3)
    # # ax9.set_xlim(-.5, 1.6)
    # # ax9.set_ylim(-.6, 1.)
    # #
    ax0.set_xticks(np.arange(0, 3600, step=600))
    ax0.set_xticklabels([np.int(n / 60) for n in ax0.get_xticks()])
    ax0.set_yticklabels([np.int(n * 1e-3) for n in ax0.get_yticks()])
    ax1.set_xticks(np.arange(0, 3600, step=600))
    ax1.set_xticklabels([np.int(n / 60) for n in ax1.get_xticks()])
    ax4.set_xticks(np.arange(0, 3600, step=600))
    ax4.set_xticklabels([np.int(n / 60) for n in ax4.get_xticks()])
    ax7.set_xticks(np.arange(0, 3600, step=600))
    ax7.set_xticklabels([np.int(n / 60) for n in ax7.get_xticks()])
    ax8.set_xticks(np.arange(0, 3600, step=600))
    ax8.set_xticklabels([np.int(n / 60) for n in ax8.get_xticks()])
    ax9.set_xticklabels([n for n in ax9.get_xticks()])
    ax9.set_yticklabels([n for n in ax9.get_yticks()])
    for label in ax0.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    # for label in ax4.yaxis.get_ticklabels()[0::2]:
    #     label.set_visible(False)
    # for label in ax4.yaxis.get_ticklabels()[0::2]:
    #     label.set_visible(False)
    plt.subplots_adjust(bottom=0.05, right=.95, left=0.06, top=0.95, wspace=0.3, hspace=0.3)
    print('saving: '+str(os.path.join(path_out_figs, fig_name)))
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.show(fig)
    return

# --------------------------------------
def plot_R_1x3(r_av, r_av_ens_mean, dR_ens_mean, R0_log, tmin_r, m_fitted,
            log_params_fitted, log_params_fitted_mean,
               r_params, dt_fields, collist, path_out_figs, fig_name):
    print('')
    print('plotting R (1x3)')
    R0 = r_av[:, tmin_r]
    t0_r = np.double(tmin_r * dt_fields)

    fig, axes = plt.subplots(1, 3, sharex='none', figsize=(18, 5.5))

    ax0 = axes[0]
    ax1 = axes[1]
    ax2 = axes[2]
    ax0.plot([t0_r,t0_r],[0,14e3], '.5', linewidth=.5)
    ax1.plot([t0_r,t0_r],[0,3], '.5', linewidth=.5)
    ax2.plot([0,0],[-1,1.], '.5', linewidth=.5)
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar]) + 'm'
        ax0.plot(times[:-1], r_av[istar, :-1], 'o-', color=collist[istar], label=lbl)
        ax1.plot(times[:-1], r_av[istar, :-1] / R0[istar], 'o-', color=collist[istar])
        ax2.plot(np.log(times[1:-1] / t0_r), np.log(r_av[istar, 1:-1]) - R0_log[istar], 'o-', color=collist[istar])

    ''' mean R '''
    r_av_rel_mean = 0.
    for istar in range(n_params):
        r_av_rel_mean += r_av[istar, :] / R0[istar]
    r_av_rel_mean /= n_params
    ax1.plot(times[1:-1], r_av_rel_mean[1:-1], '-', color='k', linewidth=4, label='ensemble mean')

    ''' fits '''
    for istar,m_ in enumerate(m_fitted):
        ax2.plot(np.log(times[1:-1]/t0_r), m_*np.log(times[1:-1]/t0_r), color=collist[istar], linewidth=1,
                 # label='y=m*x, m='+str(np.round(m_,2)))
                 label=r'y=m$\cdot$ x, m='+str(np.round(m_,2)))
                 # label=r'$y=m\cdot x, m=$'+str(np.round(m_,2)))
    a = np.zeros(n_params)
    b = np.zeros(n_params)
    c = np.zeros(n_params)
    for istar, rstar in enumerate(r_params):

        [a[istar], b[istar], c[istar]] = log_params_fitted[istar]
        lbl = 'a=' + str(np.round(a[istar], 1)) + ',b=' + str(np.round(b[istar], 1)) + ',c=' + str(np.round(c[istar], 4))
        print('r*='+str(rstar)+': ' + lbl)
        if istar == 0:
            ax0.plot(times[1:-1], a[istar] + b[istar] * np.log(1 + c[istar] * times[1:-1]), 'r-',
                     #label='R=a+b*log(1+cx)')
                    # label = r'$R=a+b\,\log(1+c\,t)$')
                    label = r'R=a+b$\,$log(1+c$\,$t)')
        else:
            ax0.plot(times[1:-1], a[istar] + b[istar] * np.log(1 + c[istar] * times[1:-1]), 'r-')
    [a, b, c] = log_params_fitted_mean
    lbl = 'a=' + str(np.round(a, 1)) + ',b=' + str(np.round(b, 1)) + ',c=' + str(np.round(c, 4))
    print('ensemble mean: ' + lbl)
    ax1.plot(times[1:-1], a + b * np.log(1 + c * times[1:-1]), 'r-',
                 # label=r'$R=a+b\,\log(1+c\,t)$ (fit ens. mean)')
                 label=r'R=a+b$\,$log(1+c$\,$t) (fit ensemble mean)')

    ax0.set_xlabel('t  [min]')
    ax0.set_ylabel('R  [km]')
    ax1.set_xlabel('t  [min]')
    ax1.set_ylabel('R/R0  [-]')
    ax2.set_xlabel('log(t/t0)  [-]')
    ax2.set_ylabel('log(R/R0)  [-]')
    leg = ax0.legend(bbox_to_anchor=(.012, .99), loc='upper left', fontsize=15, borderpad=.1)#, edgecolor='w')#, facecolor='w')
    leg.get_frame().set_edgecolor('w')
    leg = ax1.legend(bbox_to_anchor=(.012, .99), loc='upper left', fontsize=15, fancybox=False, shadow=False, ncol=1, borderpad=.1)#, frameon=False)
    leg.get_frame().set_edgecolor('w')
    leg = ax2.legend(bbox_to_anchor=(.012, .99), loc='upper left', fontsize=15, fancybox=False, shadow=False, ncol=1, borderpad=.1)#, frameon=False)
    leg.get_frame().set_edgecolor('w')
    ax2.set_xlim(-.5, 1.6)
    ax2.set_ylim(-.6, 1.)

    ax0.set_xticks(np.arange(0, 3600, step=600))
    ax0.set_xticklabels([np.int(n / 60) for n in ax0.get_xticks()])
    ax0.set_yticklabels([np.int(n * 1e-3) for n in ax0.get_yticks()])
    ax1.set_xticks(np.arange(0, 3600, step=600))
    ax1.set_xticklabels([np.int(n / 60) for n in ax1.get_xticks()])
    ax1.set_yticklabels([np.int(n) for n in ax1.get_yticks()])
    ax2.set_xticklabels([n for n in ax2.get_xticks()])
    ax2.set_yticklabels([n for n in ax2.get_yticks()])
    for label in ax0.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    for label in ax1.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    for label in ax2.yaxis.get_ticklabels()[0::2]:
        label.set_visible(False)

    textprops = dict(facecolor='white', alpha=0.9, linewidth=0.)
    t_pos_x = [30, .5, -.5]
    t_pos_y = [14.9e3, 3.2, 1.1]
    labels = ['a)', 'b)', 'c)', 'd)']
    for i, ax in enumerate(axes.flat):
        ax.text(t_pos_x[i], t_pos_y[i], labels[i], fontsize=24, horizontalalignment='left', bbox=textprops)

    plt.subplots_adjust(bottom=0.12, right=.99, left=0.05, top=0.87, wspace=0.25, hspace=0.2)
    print('saving: ' + str(os.path.join(path_out_figs, fig_name)))

    fig.savefig(os.path.join(path_out_figs, fig_name))
    fig.savefig(os.path.join(path_out_figs, fig_name[:-4]+'.pdf'))
    plt.close(fig)
    return

# --------------------------------------
def plot_R_1x2(r_av, r_av_ens_mean, dR_ens_mean, R0_log, t0_r, m_fitted,
               r_params, dt_fields, collist, path_out_figs, fig_name):

    fig, axes = plt.subplots(1, 2, sharex='none', figsize=(13, 5))

    ax0 = axes[0]
    ax1 = axes[1]
    ax0.plot([t0_r,t0_r],[0,13e3], '.5', linewidth=.5)
    ax1.plot([0,0],[-1,1.], '.5', linewidth=.5)
    ax1.plot([-.5,2.],[0,0],'.5', linewidth=.5)
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar]) + 'm'
        ax0.plot(times[:-1], r_av[istar, :-1], 'o-', color=collist[istar], label=lbl)
        # ax0.plot(times[1:-1], r_av[istar, 1:-1], 'o-', color=colorlist5[istar], label=lbl)
        # ax1.plot(np.log(times[1:-1] / t0_r), np.log(r_av[istar, 1:-1]) - R0_log[istar], '-', color='0.5')
        ax1.plot(np.log(times[1:-1] / t0_r), np.log(r_av[istar, 1:-1]) - R0_log[istar], 'o-', color=collist[istar])

    ''' mean values '''
    ax1.plot(np.log(times[3:-1] / t0_r), np.log(r_av_ens_mean[3:-1] * dR_ens_mean), '-', color='k', linewidth=3,
             label='ensemble mean')
    # m = 0.6
    # ax1.plot(np.log(times[3:-1] / t0_r), m * (np.log(times[3:-1] / t0_r)), '-r', label='y=a*x, a=' + str(m))
    # m = 0.54
    # ax1.plot(np.log(times[4:-1]) - np.log(t0_r), m * (np.log(times[4:-1]) - np.log(t0_r)), '-g', label='m=0.54')
    for istar,m_ in enumerate(m_fitted):
        ax1.plot(np.log(times[1:-1]/t0_r), m_*np.log(times[1:-1]/t0_r), color=collist[istar], linewidth=1, label='y=m*x, m='+str(np.round(m_,2)))

    # ax1.plot(np.log(times[3:-1] / (t0_r - 3)), m * (np.log(times[3:-1] / (t0_r - 3))), '-b', label='m=' + str(m))
    # ax1.plot(np.log(times[3:-1]/(t0_r-3)), m * (np.log(times[3:-1] / (t0_r-3))), '-b', label='m='+str(m))

    ax0.set_xlabel('t  [min]')
    ax0.set_ylabel('R  [km]')
    ax1.set_xlabel('log(t/t0)  [-]')
    ax1.set_ylabel('log(R/R0)  [-]')
    leg = ax0.legend(bbox_to_anchor=(.012, .99), loc='upper left', fontsize=15, borderpad=.1)#, edgecolor='w')#, facecolor='w')
    leg.get_frame().set_edgecolor('w')
    leg = ax1.legend(bbox_to_anchor=(.012, .99), loc='upper left', fontsize=15, fancybox=False, shadow=False, ncol=1, borderpad=.1)#, frameon=False)
    leg.get_frame().set_edgecolor('w')
    ax1.set_xlim(-.5, 1.6)
    ax1.set_ylim(-.6, 1.)

    ax0.set_xticks(np.arange(0, 3600, step=600))
    ax0.set_xticklabels([np.int(n / 60) for n in ax0.get_xticks()])
    ax0.set_yticklabels([np.int(n * 1e-3) for n in ax0.get_yticks()])
    ax1.set_xticklabels([n for n in ax1.get_xticks()])
    ax1.set_yticklabels([n for n in ax1.get_yticks()])
    for label in ax0.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    for label in ax1.yaxis.get_ticklabels()[0::2]:
        label.set_visible(False)

    textprops = dict(facecolor='white', alpha=0.9, linewidth=0.)
    t_pos_x = [110, -.5]
    # t_pos_y = [12.5e3, .85, 5.3, .67]
    t_pos_y = [14.5e3, 1.08]
    labels = ['a)', 'b)', 'c)', 'd)']
    for i, ax in enumerate(axes.flat):
        ax.text(t_pos_x[i], t_pos_y[i], labels[i], fontsize=21, horizontalalignment='left', bbox=textprops)

    plt.subplots_adjust(bottom=0.12, right=.95, left=0.06, top=0.9, wspace=0.25, hspace=0.2)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    return
# --------------------------------------


def plot_R_U_2x2(r_av, r_av_ens_mean, dR_ens_mean, R0_log, t0_r,
                 drdt_av, U_ens_mean, dU_ens_mean, dU_log,
                 r_params, dt_fields, path_out_figs, fig_name):
    fig, axes = plt.subplots(2, 2, sharex='none', figsize=(15, 10))

    ax0 = axes[0,0]
    ax2 = axes[0,1]
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar]) + 'm'
        ax0.plot(times[1:-1], r_av[istar, 1:-1], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[1:-1] / t0_r), np.log(r_av[istar, 1:-1]) - R0_log[istar], '-', color='0.5')
    ax2.plot(np.log(times[3:-1]/t0_r), np.log(r_av_ens_mean[3:-1] * dR_ens_mean), '-', color='k', linewidth=3, label='ensemble mean')
    m = 0.6
    ax2.plot(np.log(times[3:-1]/t0_r), m * (np.log(times[3:-1] / t0_r)), '-r', label='y=a*x, a='+str(m))
    ax2.plot(np.log(times[4:-1]) - np.log(t0_r), m * (np.log(times[4:-1]) - np.log(t0_r)), '-g', label='m=0.54')
    # ax2.plot(np.log(times[3:-1]/(t0_r-3)), m * (np.log(times[3:-1] / (t0_r-3))), '-b', label='m='+str(m))
    # ax2.plot(np.log(times[3:-1]/(t0_r-3)), m * (np.log(times[3:-1] / (t0_r-3))), '-b', label='m='+str(m))
    ax0.set_xlabel('t  [min]')
    ax0.set_ylabel('R  [km]')
    ax2.set_xlabel('log(t/t0)  [-]')
    ax2.set_ylabel('log(R/R0)  [-]')
    ax0.legend(loc=4)
    ax2.legend(loc=4)


    ax0 = axes[1, 0]
    ax2 = axes[1, 1]
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar]) + 'm'
        ax0.plot(times[1:-1], drdt_av[istar, 1:-1], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[3:-1] / t0_r), np.log(drdt_av[istar, 3:-1]) - dU_log[istar], '-', color='0.5')
    ax2.plot(np.log(times[3:-1] / t0_r), np.log(U_ens_mean[3:-1] * dU_ens_mean), '-', color='k', linewidth=3, label=lbl)
    times_ = np.append(times, np.arange(times[-1]+dt_fields, times[-1]+20*dt_fields, dt_fields))
    for m in np.arange(0.5,1,0.1):
        #ax2.plot(np.log(times_[7:] / t0_u), (np.log(U00)) - m*(np.log((times_[7:])/t0_u)-0.5), '-r', linewidth=1, label='m=-'+str(m))
        ax2.plot(np.log(times[3:-1] / t0_r), -m*(np.log((times[3:-1])/t0_r)), '-r', linewidth=1, label='m=-'+str(m))
    ax0.set_xlabel('t  [min]')
    ax0.set_ylabel('U  [m/s]')
    ax2.set_xlabel('log(t/t0)  [-]')
    ax2.set_ylabel('log(U/U0)  [-]')
    # ax0.plot([t0_u, t0_u], [0,6], 'k')
    ax2.legend(bbox_to_anchor=(1., 1.0), loc='upper right', fontsize=15, frameon=False)


    for ax in axes[:,1].flat:
        ax.set_xlim(-.5, 1.6)

    for ax in axes[:,0].flat:
        ax.set_xticks(np.arange(0, 3600, step=600))
        x_ticks = [np.int(n/60) for n in ax.get_xticks()]
        ax.set_xticklabels(x_ticks)
    y_ticks = [np.int(n * 1e-3) for n in axes[0,0].get_yticks()]
    axes[0,0].set_yticklabels(y_ticks)
    y_ticks = [np.int(n) for n in axes[1,0].get_yticks()]
    axes[1,0].set_yticklabels(y_ticks)
    for ax in axes[:,1].flat:
        ax.set_xticklabels([n for n in ax.get_xticks()])
        ax.set_yticklabels([n for n in ax.get_yticks()])

    textprops = dict(facecolor='white', alpha=0.9, linewidth=0.)
    t_pos_x = [120, -.43, 120, -.43]
    t_pos_y = [12.5e3, .85, 5.3, .67]
    labels = ['a)', 'b)', 'c)', 'd)']
    for i,ax in enumerate(axes.flat):
        ax.text(t_pos_x[i], t_pos_y[i], labels[i], fontsize=21, horizontalalignment='left', bbox=textprops)

    plt.subplots_adjust(bottom=0.05, right=.95, left=0.06, top=0.95, wspace=0.25, hspace=0.2)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)


    return



# ----------------------------------------------------------------------
def test_R_figure(r_av, tmin_r, r_params, dt_fields, path_out_figs, fig_name):
    ''' scaling of R'''
    tmin = tmin_r

    R0 = r_av[0, tmin]
    dR_log = np.log(r_av[:, tmin]) - np.log(R0)
    dR = R0 / r_av[:, tmin]
    r_av_ens_mean = np.average(r_av[:, :], axis=0)
    dR_ens_mean = R0 / r_av_ens_mean[tmin]
    t0 = np.double(tmin * dt_fields)

    fig, axes = plt.subplots(6, 4, sharex='none', figsize=(18, 20))
    [ax0, ax1, ax2, ax3] = [axes[0, i] for i in range(4)]
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar]) + 'm'
        ax0.plot([t0, t0], [0, np.amax(r_av)], 'k', linewidth=1)
        ax2.plot([np.log(times[1]), np.log(times[1])], [np.log(np.amin(r_av)), np.log(np.amax(r_av))], 'k', linewidth=1)
        ax0.plot(times[1:-1], r_av[istar, 1:-1], 'o-', color=colorlist5[istar], label=lbl)
        ax1.loglog(times[1:-1], r_av[istar, 1:-1], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[1:-1]), np.log(r_av[istar, 1:-1]), 'o-', color=colorlist5[istar], label=lbl)
    # ax3.plot(np.arange(-1, 4), 1. / dR[:], '-ok', markersize=4)
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
        ax0.plot(times[4:-1] - t0, r_av[istar, 4:-1], 'o-', color=colorlist5[istar], label=lbl)
        ax1.loglog(times[4:-1] - t0, r_av[istar, 4:-1], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[4:-1] - t0), np.log(r_av[istar, 4:-1]), 'o-',
                 color=colorlist5[istar], label=lbl)
    ax0.plot(times[4:-1] - t0, r_av_ens_mean[4:-1], 'b')
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
        ax0.plot(times[4:-1] / t0, r_av[istar, 4:-1], 'o-', color=colorlist5[istar], label=lbl)
        ax1.loglog(times[4:-1] / t0, r_av[istar, 4:-1], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[4:-1]) - np.log(t0), np.log(r_av[istar, 4:-1]), 'o-',
                 color=colorlist5[istar], label=lbl)
        ax3.plot(np.log(times[4:-1] / t0), np.log(r_av[istar, 4:-1]), 'o-',
                 color=colorlist5[istar], label=lbl)
    ax0.plot(times[4:-1] / t0, r_av_ens_mean[4:-1], 'b')
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
        ax0.plot(times[4:-1] / t0, r_av[istar, 4:-1] / r_av[istar, tmin], 'o-',
                 color=colorlist5[istar], label=lbl)
        ax1.loglog(times[4:-1] / t0, r_av[istar, 4:-1] / r_av[istar, tmin], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[4:-1]) - np.log(t0), np.log(r_av[istar, 4:-1]) - np.log(r_av[istar, tmin]), 'o-',
                 color=colorlist5[istar], label=lbl)
    #     ax3.plot(np.log(times[4:-1]) - np.log(t0), m * (np.log(r_av[istar, 4:-1]) - dR_log[istar]), 'o-',
    #              color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[4:-1]) - np.log(t0), m * (np.log(times[4:-1]) - np.log(t0)), '-r', label='m=0.54')
    #
    #     ax0.plot(times[4:tmin] - t0, r_av[istar, 4:tmin] - dR_log[istar], 'ow')
    #     ax1.loglog(times[4:tmin] - t0, r_av[istar, 4:tmin], 'ow')
    #     ax2.plot(np.log(times[4:tmin] - t0), np.log(r_av[istar, 4:tmin]) - dR_log[istar],
    #              'ow')
    # ax0.plot(times[4:-1] / t0, r_av_ens_mean[4:-1] * dR_ens_mean, 'b')
    ax0.set_xlabel('t/t0')
    ax1.set_xlabel('t/t0')
    ax2.set_xlabel('log(t)-log(t0)')
    ax0.set_ylabel('R/R(t0)')
    ax1.set_ylabel('R/R(t0)')
    ax2.set_ylabel(r'log(R)-log(R(t0))')
    ax3.set_ylabel(r'log(R)-log(R(t0))')
    # ax3.set_ylabel('log(R)-[log(R(t0)-log(R0(t0)]')
    # ax3.set_xlim(np.log(times[tmin]), np.log(times[-2]))
    ax2.set_xlim(0, 2.3)
    ax3.set_xlim(0, 2.3)
    ax0.set_ylim(0, 3)

    [ax0, ax1, ax2, ax3] = [axes[4, i] for i in range(4)]
    m = 0.55
    for istar in range(n_params):
        lbl = 'r*=' + str(r_params[istar])
        ax0.plot([0, 0], [0, np.amax(r_av)], 'k', linewidth=1)
        ax0.plot(times[4:-1] / t0, r_av[istar, 4:-1] * dR[istar], 'o-',
                 color=colorlist5[istar], label=lbl)
        ax1.loglog(times[4:-1] / t0, r_av[istar, 4:-1] * dR[istar], 'o-', color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[4:-1]) - np.log(t0), np.log(r_av[istar, 4:-1]) - dR_log[istar], 'o-',
                 color=colorlist5[istar], label=lbl)
        ax3.plot(np.log(times[4:-1]) - np.log(t0), m * (np.log(r_av[istar, 4:-1]) - dR_log[istar]), 'o-',
                 color=colorlist5[istar], label=lbl)
        ax2.plot(np.log(times[4:-1]) - np.log(t0),
                 np.log(R0) + m * (np.log(times[4:-1]) - np.log(t0)), '-r', label='m=0.54')

        # ax0.plot(times[4:tmin] - t0, r_av[istar, 4:tmin] - dR_log[istar], 'ow')
        # ax1.loglog(times[4:tmin] - t0, r_av[istar, 4:tmin], 'ow')
        # ax2.plot(np.log(times[4:tmin] - t0), np.log(r_av[istar, 4:tmin]) - dR_log[istar], 'ow')
    ax0.plot(times[4:-1] / t0, r_av_ens_mean[4:-1] * dR_ens_mean, 'b')
    ax0.set_xlabel('t/t0')
    ax1.set_xlabel('t/t0')
    ax2.set_xlabel('log(t)-log(t0)')
    ax0.set_ylabel('R * R0(t0)/R(t0)')
    ax1.set_ylabel('R * R0(t0)/R(t0)')
    ax2.set_ylabel(r'log(R)-log(R0(t0)/R(t0))')
    ax3.set_ylabel(r'log(R)-log(R0(t0)/R(t0))')
    # ax3.set_ylabel('log(R)-[log(R(t0)-log(R0(t0)]')
    # ax3.set_xlim(np.log(times[tmin]), np.log(times[-2]))
    ax2.set_xlim(0, 2.3)
    ax3.set_xlim(0, 2.3)
    ax0.set_ylim(0, 8e3)
    #
    # [ax0, ax1, ax2, ax3] = [axes[5, i] for i in range(4)]
    # fit_r = np.log(R0) + m * (np.log(times / t0))
    # for istar in range(n_params):
    #     ax0.plot(times[1:-1] / t0, r_av[istar, 1:-1] * dR[istar], '-', color='0.5')
    #     ax1.loglog(times[1:-1] / t0, r_av[istar, 1:-1] * dR[istar], '-', color='0.5')
    #     ax2.plot(np.log(times[1:-1] / t0), np.log(r_av[istar, 1:-1]) - dR_log[istar], '-', color='0.5')
    # lbl = 'mean(r_av)'
    # ax0.plot(times[1:-1] / t0, r_av_ens_mean[1:-1] * dR_ens_mean, '-', color='b', linewidth=2, label=lbl)
    # ax1.loglog(times[1:-1] / t0, r_av_ens_mean[1:-1] * dR_ens_mean, '-', color='b', linewidth=2, label=lbl)
    # ax2.plot(np.log(times[1:-1] / t0), np.log(r_av_ens_mean[1:-1] * dR_ens_mean), '-', color='b', linewidth=2,
    #          label=lbl)
    # ax3.plot(np.log(times[4:-1] / t0),
    #          np.log(r_av_ens_mean[4:-1] * dR_ens_mean) - np.log(R0) + m * (np.log(times[4:-1] / t0)),
    #          '-', color='k', linewidth=2, label='diff')
    #
    # ax2.plot(np.log(times[1:-1] / t0), fit_r[1:-1], '-r', label='m=' + str(m))
    # ax2.legend(loc='best')
    #
    # axes[0, 0].legend(loc='best')
    # for ax in axes.flatten():
    #     ax.grid(True)
    #
    # plt.suptitle('scaling R')
    # plt.subplots_adjust(bottom=0.05, right=.95, left=0.06, top=0.95, wspace=0.25, hspace=0.2)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    return

# ----------------------------------------------------------------------

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
    cp_number = int(lines[-1].split()[3])

    f.close()

    return cp_number

# ----------------------------------------------------------------------

def set_input_parameters(args, path_root, case_name):
    print('--- set input parameters ---')

    global n_params
    dTh = 5
    z_params = [1000]  # run5, 6
    r_params = [500, 1000, 1100, 1600, 2300]
    n_params = len(r_params)
    print('dTh: ', dTh)
    print('z*: ', z_params)
    print('r*: ', r_params)
    print('n_params: ', n_params)

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

    global times, nt, k0
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = np.int(100)
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = np.int(3600)
    times = np.arange(tmin, tmax + 100, 100)
    times.sort()

    nt = len(times)
    k0 = 0
    # krange = [0]
    # nk = len(krange)

    print('times', times)
    print('nt: '+str(nt))
    print('k0: ' + str(k0))
    print ''

    return dTh, z_params, r_params

# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()
