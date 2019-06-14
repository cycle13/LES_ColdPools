import numpy as np
import matplotlib
import matplotlib.pyplot as plt
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

    global cm_bwr, cm_grey, cm_vir, cm_hsv, cm_grey2
    cm_bwr = plt.cm.get_cmap('seismic')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')
    cm_grey2 = plt.cm.get_cmap('bone_r')
    cm_fall = plt.cm.get_cmap('winter')
    cm_summer = plt.cm.get_cmap('spring')

    dTh, z_params, r_params = set_input_parameters(args)

    # reference case: dTh3_z1000_r1000
    id_ref = 'dTh3_z1000_r1000'
    # path_ref = os.path.join(path_root, id_ref)
    # path_ref = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/' + id_ref
    # path_ref = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run5_PE_scaling_dx100m/' + id_ref
    path_ref = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run3_dx50m/' + id_ref
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
    krange = [0]
    nk = len(krange)
    k0 = 0

    # --------------------------------------
    ''' (a) read in data from tracer output (text-file)'''
    n_params = len(r_params)
    dist_av = np.zeros((n_params, nt, nk))
    r_av = np.zeros((n_params, nt, nk))
    drdt_av = np.zeros((n_params, nt, nk))
    U_rad_av = np.zeros((n_params, nt, nk))
    dU_rad_av = np.zeros((n_params, nt, nk))
    dist_av_ref = np.zeros((nt, nk))
    U_rad_av_ref = np.zeros((nt, nk))

    for istar in range(n_params):
        zstar = z_params[0]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        print('id', id)
        fullpath_in = os.path.join(path_root, id, 'tracer_k'+str(k0), 'output')
        fullpath_in_ref = os.path.join(path_ref, 'tracer_k'+str(k0), 'output')
        # read_in_txtfile(fullpath_in)
        for it, t0 in enumerate(times):
            print('---t0: '+str(t0)+'---', it)
            dist_av[istar, it, k0], U_rad_av[istar, it, k0] = get_radius_vel(fullpath_in, it, cp_id, n_tracers, n_cps)
            dist_av_ref[it, k0], U_rad_av_ref[it, k0] = get_radius_vel(fullpath_in_ref, it, cp_id, n_tracers, n_cps)
        r_av = dist_av * dx[0]
        r_av_ref = dist_av_ref * dx[0]
        for it, t0 in enumerate(times[1:]):
            drdt_av[:,it,:] = 1./dt_fields * (r_av[:,it,:] - r_av[:,it-1,:])
        print ''
    print ''


    ''' (b) plot r_av, dtdt_av, U_rad_av'''
    figname = 'CP_rim_dTh' + str(dTh) + '.png'
    title = 'CP rim (dTh=' + str(dTh) + 'K, dx='+str(dx[0]) +'m)'
    plot_dist_vel(r_av, drdt_av, U_rad_av, r_av_ref, U_rad_av_ref,
                  [dTh], z_params, r_params, n_params, k0, id_ref, title, figname)
    ''' fit function to U_rad '''
    plot_vel_fitting(r_av, drdt_av, U_rad_av, r_av_ref, U_rad_av_ref,
                  [dTh], z_params, r_params, n_params, k0, id_ref, figname)



    ''' (c) plot normalized radius / velocity'''
    # (i) for vertical velocity from crosssection in 3D field
    # (ii) for azimuthally averaged vertical vleocity
    print('plotting normalized')
    for istar in range(n_params):
        zstar = z_params[0]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        print('id', id)
        nml = simplejson.loads(open(os.path.join(path_root, id, case_name + '.in')).read())
        nx = nml['grid']['nx']
        ny = nml['grid']['ny']
        ic = np.int(nx / 2)
        jc = np.int(ny / 2)
        print(ic, jc)
        figname_norm = 'CP_rim_normalized_' + id + '.png'
        plot_vel_normalized(r_av, times, istar, k0, nx, ic, jc, id, figname_norm)
        # figname_norm = 'CP_rim_normalized_' + id + '_av.png'
        # plot_vel_normalized_w_av(r_av, times, istar, k0, id, figname_norm)


    trange = [600, 1200, 1800, 2400, 3000, 3600 ]
    fig_name = 'PE_scaling_dTh' + str(dTh) + '.png'
    fig, axes = plt.subplots(1, 2, sharex='none', figsize=(12, 5))
    ax0 = axes[0]
    ax1 = axes[1]
    scaling = [-1, 0, 1, 2, 3]
    r_av_ = np.zeros(shape=len(scaling))
    U_rad_av_ = np.zeros(shape=len(scaling))
    for t0 in trange:
        it = np.where(times == t0)[0][0]
        for i,s in enumerate(scaling):
            if s < 0:
                r_av_[i] = r_av[i, it, k0]
                U_rad_av_[i] = U_rad_av[i, it, k0]
            elif s == 0:
                r_av_[i] = r_av_ref[it, k0]
                U_rad_av_[i] = U_rad_av_ref[it, k0]
            else:
                r_av_[i] = r_av[i-1, it, k0]
                U_rad_av_[i] = U_rad_av[i-1, it, k0]
        ax0.plot(scaling, r_av_[1] + 9.5e2*np.asarray(scaling), 'k', linewidth=1)
        ax0.plot(scaling, r_av_[3]+ 12e2*(np.asarray(scaling)-2), 'k--', linewidth=1)
        ax0.plot(scaling, 2e3 + 0.1*1e3*(np.asarray(scaling)+2)**2, 'r', linewidth=1)
        ax0.plot(scaling, r_av_[:], 'o-', label='t='+str(t0))
        ax1.plot(scaling, U_rad_av_[:], 'o-', label='t='+str(t0))
    ax0.set_xlabel('log2(PE/PE0)')
    ax1.set_xlabel('log2(PE/PE0)')
    ax0.set_ylabel('r_av  [m]')
    ax1.set_ylabel('U_rad_av [m/s]')
    ax1.legend(loc='best')
    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, fig_name))
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
def plot_dist_vel(r_av, drdt_av, U_rad_av, r_av_ref, U_rad_av_ref,
                  dTh_params, z_params, r_params, n_params, k0,
                  id_ref, title, fig_name):

    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, sharex='all', figsize=(18, 5))
    for istar in range(n_params):
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[0]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        ax0.plot(times, r_av[istar, :, k0], 'o-', label=id)
        ax1.plot(times[1:], drdt_av[istar, 1:, k0], 'o-', label=id)
        ax2.plot(times, U_rad_av[istar, :, k0], 'o-', label=id)
    ax0.plot(times, r_av_ref[:, k0], 'ko-', label=id_ref)
    ax2.plot(times, U_rad_av_ref[:, k0], 'ko-', label=id_ref)

    # ax0.set_title('r_av')
    # ax1.set_title('drdt_av')
    # ax2.set_title('U_av')
    ax0.set_xlabel('times [s]')
    ax1.set_xlabel('times [s]')
    ax2.set_xlabel('times [s]')
    ax0.set_ylabel('r_av  [m]')
    ax1.set_ylabel('drdt_av')
    ax2.set_ylabel('U_rad_av  [m/s]')
    ax2.legend()
    # fig.suptitle('CP rim (dTh=' + str(dTh) + ')')
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, fig_name))
    # fig.savefig(os.path.join(path_out_figs, 'CP_rim_dTh' + str(dTh) + '.png'))
    plt.close(fig)


    fig, axis = plt.subplots(2, 2, sharex='none', figsize=(18, 10))
    for istar in range(n_params):
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[0]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        axis[0,0].semilogx(times, r_av[istar, :, k0], 'o-', label=id)
        axis[0,1].semilogx(times, U_rad_av[istar, :, k0], 'o-', label=id)
        axis[1,0].loglog(times, r_av[istar, :, k0], 'o-', label=id)
        axis[1,1].loglog(times, U_rad_av[istar, :, k0], 'o-', label=id)
    axis[0, 0].semilogx(times, r_av_ref[:, k0], 'ko-', label=id)
    axis[0, 1].semilogx(times, U_rad_av_ref[:, k0], 'ko-', label=id)
    axis[1, 0].loglog(times, r_av_ref[:, k0], 'ko-', label=id)
    axis[1, 1].loglog(times, U_rad_av_ref[:, k0], 'ko-', label=id)
    axis[0, 0].set_title('CP radius (r_av)', fontsize=18)
    axis[0, 1].set_title('radial spreading velocity (U_av)', fontsize=18)
    axis[0,0].set_xlabel('times [s]')
    axis[1, 0].set_xlabel('time [s]')
    axis[1, 1].set_xlabel('time [s]')
    axis[0, 0].set_ylabel('r_av  [m]')
    axis[1, 0].set_ylabel('r_av  [m]')
    axis[0, 1].set_ylabel('U_rad_av  [m/s]')
    axis[1, 1].set_ylabel('U_rad_av  [m/s]')
    fig.tight_layout()
    plt.subplots_adjust(bottom=0.075, right=.95, left=0.07, top=0.9, wspace=0.25)
    axis[1,1].legend(loc='best')
    # fig.suptitle('CP rim (dTh=' + str(dTh) + ')')
    fig.suptitle(title, fontsize=21)
    fig.savefig(os.path.join(path_out_figs, fig_name[:-4]+'_log.png'))
    # fig.savefig(os.path.join(path_out_figs, 'CP_rim_dTh' + str(dTh) + '.png'))
    plt.close(fig)
    return



def plot_vel_fitting(r_av, drdt_av, U_rad_av, r_av_ref, U_rad_av_ref,
                      dTh_params, z_params, r_params, n_params, k0,
                      id_ref, fig_name):
    import scipy
    from scipy import optimize

    # Fit the first set
    # fitfunc = lambda p, x: p[0] * np.cos(2 * np.pi / p[1] * x + p[2]) + p[3] * x  # Target function
    fitfunc1 = lambda p, x: p[0] + p[1] * x ** p[2]  # Target function
    errfunc = lambda p, x, y: fitfunc1(p, x) - y  # Distance to the target function
    n_init = 5
    p0 = np.zeros((n_init, 3))
    p1 = np.zeros((n_init, 3))
    p0[0, :] = [0., 0., -3.]  # Initial guess for the parameters
    p0[1, :] = [0., 0., 2.]  # Initial guess for the parameters
    p0[2, :] = [10., 0., 2.]  # Initial guess for the parameters
    p0[3, :] = [0., 0., -2.]  # Initial guess for the parameters
    p0[4, :] = [10., 0., -2.]  # Initial guess for the parameters


    for istar in range(n_params):
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[0]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)

        fig, axis = plt.subplots(1, 3, sharex='none', figsize=(18, 7))
        axis[0].plot(times, U_rad_av[istar, :, k0], 'o-', label=id)
        axis[1].semilogx(times, U_rad_av[istar, :, k0], 'o-', label=id)
        axis[2].loglog(times, U_rad_av[istar, :, k0], 'o-', label=id)
        tmin = 4
        axis[0].plot(times, U_rad_av_ref[:, k0], 'ko-', label=id)
        axis[1].semilogx(times, U_rad_av_ref[:, k0], 'ko-', label=id)
        axis[2].loglog(times, U_rad_av_ref[:, k0], 'ko-', label=id)
        for i in range(n_init):
            p1[i, :], success = optimize.leastsq(errfunc, p0[i, :], args=(times[tmin:], U_rad_av[istar, tmin:, k0]))
            axis[0].plot(times[tmin:], fitfunc1(p1[i,:], times[tmin:]), "-", label='p='+str(p1[i,:]))  # Plot of the data and the fit
            axis[1].semilogx(times, fitfunc1(p1[i,:], times), "-", label='p='+str(p1[i,:]))  # Plot of the data and the fit
            axis[2].loglog(times, fitfunc1(p1[i,:], times), "-", label='p='+str(p1[i,:]))  # Plot of the data and the fit
        # # axis[0, 1].set_title('radial spreading velocity (U_av)', fontsize=18)
        axis[0].set_xlabel('time [s]')
        axis[1].set_xlabel('time [s]')
        axis[2].set_xlabel('time [s]')
        axis[0].set_ylabel('U_rad_av  [m/s]')
        # fig.tight_layout()
        plt.subplots_adjust(bottom=0.3, right=.95, left=0.07, top=0.9, wspace=0.25)
        axis[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
                   fancybox=True, shadow=True, ncol=3)
        fig.suptitle('CP rim (' + id +', dx='+str(dx[0]) +'m)', fontsize=21)
        fig.savefig(os.path.join(path_out_figs, fig_name[:-4] + '_fit_' + id + '.png'))
        plt.close(fig)



    ''' reference '''
    fig, axis = plt.subplots(1, 3, sharex='none', figsize=(18, 7))
    tmin = 2
    axis[0].plot(times, U_rad_av_ref[:, k0], 'ko-', label=id)
    axis[1].semilogx(times, U_rad_av_ref[:, k0], 'ko-', label=id)
    axis[2].loglog(times, U_rad_av_ref[:, k0], 'ko-', label=id)
    for istar in range(n_params):
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[0]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        axis[0].plot(times, U_rad_av[istar, :, k0], 'o-', label=id)
        axis[1].semilogx(times, U_rad_av[istar, :, k0], 'o-', label=id)
        axis[2].loglog(times, U_rad_av[istar, :, k0], 'o-', label=id)
    for i in range(n_init):
        p1[i, :], success = optimize.leastsq(errfunc, p0[i, :], args=(times[tmin:], U_rad_av_ref[tmin:, k0]))
        axis[0].plot(times[tmin:], fitfunc1(p1[i, :], times[tmin:]), "-",
                     label='p=' + str(p1[i, :]))  # Plot of the data and the fit
        axis[1].semilogx(times, fitfunc1(p1[i, :], times), "-",
                         label='p=' + str(p1[i, :]))  # Plot of the data and the fit
        axis[2].loglog(times, fitfunc1(p1[i, :], times), "-",
                       label='p=' + str(p1[i, :]))  # Plot of the data and the fit
    # # axis[0, 1].set_title('radial spreading velocity (U_av)', fontsize=18)
    axis[0].set_xlabel('time [s]')
    axis[1].set_xlabel('time [s]')
    axis[2].set_xlabel('time [s]')
    axis[0].set_ylabel('U_rad_av  [m/s]')
    # axis[1, 1].set_ylabel('U_rad_av  [m/s]')
    # fig.tight_layout()
    plt.subplots_adjust(bottom=0.3, right=.95, left=0.07, top=0.9, wspace=0.25)
    axis[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
                   fancybox=True, shadow=True, ncol=3)
    fig.suptitle('CP rim (' + id + ', dx=' + str(dx[0]) + 'm)', fontsize=21)
    fig.savefig(os.path.join(path_out_figs, fig_name[:-4] + '_fit_' + id_ref + '.png'))
    plt.close(fig)




    ''' all '''
    cmap = cm_hsv
    fig, axis = plt.subplots(1, 3, sharex='none', figsize=(18, 7))
    axis[0].plot(times, U_rad_av_ref[:, k0], 'ko-', linewidth=1, label=id)
    axis[1].semilogx(times, U_rad_av_ref[:, k0], 'ko-', linewidth=1, label=id)
    axis[2].loglog(times, U_rad_av_ref[:, k0], 'ko-', linewidth=1, label=id)
    p1[i, :], success = optimize.leastsq(errfunc, p0[i, :], args=(times[tmin:], U_rad_av_ref[tmin:, k0]))
    axis[0].plot(times[tmin:], fitfunc1(p1[i, :], times[tmin:]), "k-", linewidth=3,
                 label='p=' + str(p1[i, :]))  # Plot of the data and the fit
    axis[1].semilogx(times[tmin:], fitfunc1(p1[i, :], times[tmin:]), "k-", linewidth=3,
                     label='p=' + str(p1[i, :]))  # Plot of the data and the fit
    axis[2].loglog(times[tmin:], fitfunc1(p1[i, :], times)[tmin:], "k-", linewidth=3,
                   label='p=' + str(p1[i, :]))  # Plot of the data and the fit
    for istar in range(n_params):
        i_color = np.double(istar)/n_params
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[0]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        tmin = 4
        i = 4   # parameters
        axis[0].plot(times, U_rad_av[istar, :, k0], 'o-', linewidth=1, color=cmap(i_color), label=id)
        axis[1].semilogx(times, U_rad_av[istar, :, k0], 'o-', linewidth=1, color=cmap(i_color), label=id)
        axis[2].loglog(times, U_rad_av[istar, :, k0], 'o-', linewidth=1, color=cmap(i_color), label=id)
        p1[i, :], success = optimize.leastsq(errfunc, p0[i, :], args=(times[tmin:], U_rad_av[istar, tmin:, k0]))
        axis[0].plot(times[tmin:], fitfunc1(p1[i, :], times[tmin:]), "-", linewidth=3, color=cmap(i_color), label='p=' + str(p1[i, :]))  # Plot of the data and the fit
        axis[1].semilogx(times[tmin:], fitfunc1(p1[i, :], times[tmin:]), "-", linewidth=3, color=cmap(i_color), label='p=' + str(p1[i, :]))  # Plot of the data and the fit
        axis[2].loglog(times[tmin:], fitfunc1(p1[i, :], times[tmin:]), "-", linewidth=3, color=cmap(i_color), label='p=' + str(p1[i, :]))  # Plot of the data and the fit
    axis[0].set_xlabel('time [s]')
    axis[1].set_xlabel('log(time) [s]')
    axis[2].set_xlabel('log(time) [s]')
    axis[0].set_ylabel('U_rad_av  [m/s]')
    axis[1].set_ylabel('U_rad_av  [m/s]')
    axis[2].set_ylabel('log(U_rad_av)  [m/s]')
    # axis[1, 1].set_ylabel('U_rad_av  [m/s]')
    # fig.tight_layout()
    plt.subplots_adjust(bottom=0.3, right=.95, left=0.07, top=0.9, wspace=0.25)
    axis[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
                   fancybox=True, shadow=True, ncol=3)
    fig.suptitle('radial spreading velocity (U_av) (' + id + ', dx=' + str(dx[0]) + 'm)', fontsize=21)
    fig.savefig(os.path.join(path_out_figs, fig_name[:-4] + '_fit_all.png'))
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
            z_params = [1000]  # run1
    if args.rparams:
        r_params = args.rparams
    else:
        r_params = [500, 1100, 1600, 2300]


    print('dTh: ', dTh)
    print('z*: ', z_params)
    print('r*: ', r_params)

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

if __name__ == '__main__':
    main()
