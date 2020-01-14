import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import scipy
from scipy import stats


execfile('settings.py')

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    case_name = 'ColdPoolDry_single_3D'
    dx = 100
    if dx == 100:
        path_root = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/'
        path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single/figs_run2_dx100m/'
    elif dx == 50:
        path_root = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run3_dx50m/'
        path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single/figs_run3_dx50m/'
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print('path in: '+ path_root)
    print('path figs: '+path_out_figs)
    print('')

    # loop through all cases
    # use azimuthally averaged statistics
    # statistics:
    # - radius
    # (- spreading velocity)
    # - max. radial velocity component ?
    # - max vertical velocity
    # - min vertical velocity ?
    # - min theta
    # - CP volume
    # - CP height
    # - max vorticity
    # - mean over CP volume of
    #     - vorticity
    #     - theta

    # in stats_radial_averaged.nc:
    # - v_rad
    # - v_tan
    # - w
    # - s
    # - phi
    # - temperature

    # - CP height: in data_analysis/CP_height_dTh3_z1000_r1000_sth0.5.nc
    # - CP volume: in data_analysis/CP_volume_dTh3_z1000_r1000.nc:

    ''' for run2 (dx=100m) '''
    dTh_range_A = [2, 3, 4]
    rstar_range_A = [1300, 1000, 900]
    zstar_range_A = [900, 1000, 900]

    dTh_range_B = [3, 3, 3, 3, 3]
    rstar_range_B = [500, 600, 700, 1000, 1500]
    zstar_range_B = [2500, 2000, 1600, 1000, 500]

    filename_stats = 'stats_radial_averaged.nc'
    filename_vort = 'Stats_vorticity_phi.nc'
    lvl_w = 1
    lvl = 0


    # plot_sensitivity_plots_temperature(dTh_range_A, rstar_range_A, zstar_range_A, path_root,
    #                                   lvl, lvl_w,
    #                                    filename_stats, filename_vort,
    #                                    colorlist3, path_out_figs)
    #
    # plot_sensitivity_plots_geometry(dTh_range_B, rstar_range_B, zstar_range_B, path_root,
    #                                 lvl, lvl_w,
    #                                 filename_stats, filename_vort,
    #                                 colorlist3, path_out_figs)

    # ''' test entropy / theta '''
    # rootname = 'dTh' + str(dTh_range_A[0]) + '_z' + str(zstar_range_A[0])\
    #            + '_r' + str(rstar_range_A[0])
    # path_in = os.path.join(path_root, rootname)
    # stats_root = nc.Dataset(os.path.join(path_in, 'data_analysis', filename_stats))
    # time_stats = stats_root.groups['timeseries'].variables['time'][:]
    # time_stats = stats_root.groups['stats'].variables['r'][:]
    # s_ = stats_root.groups['stats'].variables['s'][:, :, :]
    # s_min = np.amin(np.amin(s_, axis=-1), axis=-1)
    # theta_ = thetas_c(s_, 0)
    # theta_min = np.amin(np.amin(theta_, axis=-1), axis=-1)
    # theta_min_k0 = np.amin(theta_[:, :, lvl], axis=1)
    # # del s_, theta_
    # stats_root.close()
    # print('TEST: ', s_.shape, time_stats.shape)
    # ncols=10
    # fig_name = 'testfigure_theta.png'
    # fig, axis = plt.subplots(1, ncols, sharex='all', sharey='all', figsize=(ncols * 3, 5))
    # for i,ax in enumerate(axis):
    #     cf = ax.contourf(s_[i,:-1,:])
    #     plt.colorbar(cf, ax=ax)
    # plt.subplots_adjust(bottom=0.075, right=.99, left=0.01, top=0.95, wspace=0.2, hspace=0.1)
    # fig.savefig(os.path.join(path_out_figs, fig_name))
    # ''''''


    plot_sensitivity_plots_all(dTh_range_A, rstar_range_A, zstar_range_A,
                               dTh_range_B, rstar_range_B, zstar_range_B,
                               lvl, lvl_w, dx,
                               filename_stats, filename_vort,
                               colorlist3, path_root, path_out_figs)
    return





def plot_sensitivity_plots_temperature(dTh_range_A, rstar_range_A, zstar_range_A, path_root,
                                    lvl, lvl_w,
                                       filename_stats, filename_vort,
                                       colorlist3, path_out_figs):

    fig_name = 'sensitivity_plots_temperature.png'
    ncols = 6
    fig, axis = plt.subplots(1, ncols, sharex='all', figsize=(ncols*3, 5))
    ax1 = axis[1]
    ax2 = axis[2]
    ax3 = axis[3]

    for i, dTh in enumerate(dTh_range_A):
        rstar = rstar_range_A[i]
        zstar = zstar_range_A[i]
        rootname = 'dTh'+str(dTh) + '_z'+str(zstar) + '_r'+str(rstar)
        filename_CPheight = 'CP_height_'+rootname+'_sth0.5.nc'
        filename_CPvol = 'CP_volume_'+rootname+'.nc'
        path_in = os.path.join(path_root, rootname)

        print(os.path.join(path_in, 'data_analysis', filename_stats))
        stats_root = nc.Dataset(os.path.join(path_in, 'data_analysis', filename_stats))
        time_stats = stats_root.groups['timeseries'].variables['time'][:]
        w_ = stats_root.groups['stats'].variables['w'][:,:,:]
        w_max = np.amax(np.amax(w_, axis=-1), axis=-1)
        w_max_k0 = np.amax(w_[:,:,lvl_w], axis=1)
        w_min = np.amin(np.amin(w_, axis=-1), axis=-1)
        w_min_k0 = np.amin(w_[:,:,lvl_w], axis=1)
        del w_
        v_rad_ = stats_root.groups['stats'].variables['v_rad'][:,:,:]
        v_rad_max = np.amax(np.amax(v_rad_, axis=-1), axis=-1)
        v_rad_max_k0 = np.amax(v_rad_[:,:,lvl], axis=1)
        del v_rad_
        s_ = stats_root.groups['stats'].variables['s'][:,:,:]
        s_min = np.amin(np.amin(s_, axis=-1), axis=-1)
        theta_ = thetas_c(s_, 0)
        theta_min = np.amin(np.amin(theta_, axis=-1), axis=-1)
        theta_min_k0 = np.amin(theta_[:,:,lvl], axis=1)
        del s_, theta_
        stats_root.close()

        # # vorticity from azimuthally averaged velocity fields (v_rad, v_tan, w)
        # print(os.path.join(path_in, 'fields_vorticity', filename_vort))
        # vort_root = nc.Dataset(os.path.join(path_in, 'fields_vorticity', filename_vort))
        # time_vort = vort_root.variables['time'][:]
        # vort_phi_max = vort_root.variables['vort_phi_max'][:]
        # vort_phi_min = vort_root.variables['vort_phi_min'][:]
        # vort_root.close()

        print(os.path.join(path_in, 'data_analysis', filename_CPheight))
        root = nc.Dataset(os.path.join(path_in, 'data_analysis', filename_CPheight))
        time_geom = root.groups['timeseries'].variables['time'][:]
        CP_height_2D = root.groups['fields_2D'].variables['CP_height_2d'][:,:,:]
        CP_height_max = root.groups['timeseries'].variables['CP_height_max'][:]
        root.close()
        del CP_height_
        root = nc.Dataset(os.path.join(path_in, 'data_analysis', filename_CPvol))
        time_geom_ = root.groups['timeseries'].variables['time'][:]
        if time_geom != time_geom_:
            print('!!!!! different timeseries')
        CP_vol = root.grous['timeseries'].variables['CP_vol_sth0.5'][:]
        root.close()

        # ax0.plot(time_stats, s_min, '-o', color=colorlist3[i], label=rootname)
        ax1.plot(time_stats, theta_min, '-o', color=colorlist3[i], label=rootname)
        ax2.plot(time_stats, w_max, '-o', color=colorlist3[i], label=rootname)
        # ax3.plot(time_vort, vort_phi_max, '-o', color=colorlist3[i], label=rootname)
        ax4.plot(time_geom, CP_height_max, '-o', color=colorlist3[i], label=rootname)
        ax5.plot(time_geom, CP_vol, '-o', color=colorlist3[i], label=rootname)

    ax1.legend(loc='best')
    # fig.suptitle(title, fontsize=18)
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.06, top=0.9, wspace=0.15)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return



def plot_sensitivity_plots_geometry(dTh_range_B, rstar_range_B, zstar_range_B, path_root,
                                lvl, lvl_w,
                                           filename_stats, filename_vort,
                                           colorlist3, path_out_figs):
    fig_name = 'sensitivity_plots_geometry.png'
    ncols = 6
    fig, axis = plt.subplots(1, ncols, sharex='all', figsize=(ncols*3, 5))
    ax1 = axis[1]
    ax2 = axis[2]
    ax3 = axis[3]
    ax4 = axis[4]
    ax5 = axis[5]

    for i, dTh in enumerate(dTh_range_B):
        rstar = rstar_range_B[i]
        zstar = zstar_range_B[i]
        rootname = 'dTh'+str(dTh) + '_z'+str(zstar) + '_r'+str(rstar)
        filename_CPheight = 'CP_height_'+rootname+'_sth0.5.nc'
        filename_CPvol = 'CP_volume_'+rootname+'.nc'
        path_in = os.path.join(path_root, rootname)

        print(os.path.join(path_in, 'data_analysis', filename_stats))
        stats_root = nc.Dataset(os.path.join(path_in, 'data_analysis', filename_stats))
        time_stats = stats_root.groups['timeseries'].variables['time'][:]
        w_ = stats_root.groups['stats'].variables['w'][:,:,:]
        w_max = np.amax(np.amax(w_, axis=-1), axis=-1)
        w_max_k0 = np.amax(w_[:,:,lvl_w], axis=1)
        w_min = np.amin(np.amin(w_, axis=-1), axis=-1)
        w_min_k0 = np.amin(w_[:,:,lvl_w], axis=1)
        del w_
        v_rad_ = stats_root.groups['stats'].variables['v_rad'][:,:,:]
        v_rad_max = np.amax(np.amax(v_rad_, axis=-1), axis=-1)
        v_rad_max_k0 = np.amax(v_rad_[:,:,lvl], axis=1)
        del v_rad_
        s_ = stats_root.groups['stats'].variables['s'][:,:,:]
        # s_min = np.amin(np.amin(s_, axis=-1), axis=-1)
        theta_ = thetas_c(s_, 0)
        theta_min = np.amin(np.amin(theta_, axis=-1), axis=-1)
        theta_min_k0 = np.amin(theta_[:,:,lvl], axis=1)
        del s_, theta_
        stats_root.close()

        # # vorticity from azimuthally averaged velocity fields (v_rad, v_tan, w)
        # print(os.path.join(path_in, 'fields_vorticity', filename_vort))
        # vort_root = nc.Dataset(os.path.join(path_in, 'fields_vorticity', filename_vort))
        # time_vort = vort_root.variables['time'][:]
        # vort_phi_max = vort_root.variables['vort_phi_max'][:]
        # vort_phi_min = vort_root.variables['vort_phi_min'][:]
        # vort_root.close()

        print(os.path.join(path_in, 'data_analysis', filename_CPheight))
        root = nc.Dataset(os.path.join(path_in, 'data_analysis', filename_CPheight))
        time_geom = root.groups['timeseries'].variables['time'][:]
        CP_height_2D = root.groups['fields_2D'].variables['CP_height_2d'][:,:,:]
        CP_height_max = root.groups['timeseries'].variables['CP_height_max'][:]
        root.close()
        del CP_height_
        root = nc.Dataset(os.path.join(path_in, 'data_analysis', filename_CPvol))
        time_geom_ = root.groups['timeseries'].variables['time'][:]
        if time_geom != time_geom_:
            print('!!!!! different timeseries')
        CP_vol = root.grous['timeseries'].variables['CP_vol_sth0.5'][:]
        root.close()

        # ax0.plot(time_stats, s_min, '-o', color=colorlist3[i], label=rootname)
        ax1.plot(time_stats, theta_min, '-o', color=colorlist3[i], label=rootname)
        ax2.plot(time_stats, w_max, '-o', color=colorlist3[i], label=rootname)
        # ax3.plot(time_vort, vort_phi_max, '-o', color=colorlist3[i], label=rootname)
        ax4.plot(time_geom, CP_height_max, '-o', color=colorlist3[i], label=rootname)
        ax5.plot(time_geom, CP_vol, '-o', color=colorlist3[i], label=rootname)

    ax1.legend(loc='best')
    # fig.suptitle(title, fontsize=18)
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.06, top=0.9, wspace=0.15)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return





def plot_sensitivity_plots_all(dTh_range_A, rstar_range_A, zstar_range_A,
                               dTh_range_B, rstar_range_B, zstar_range_B,
                               lvl, lvl_w, dx,
                               filename_stats, filename_vort,
                               colorlist3, path_root, path_out_figs):


    ncols = 6
    nrows = 2
    max_range = np.zeros(ncols, dtype=np.double)
    min_range = 9999.9*np.ones(ncols, dtype=np.double)
    fig_name = 'sensitivity_plots_all_dx'+str(dx)+'m.png'
    fig, axis = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*5), sharex='col')
    ax0 = axis[0,0]
    ax1 = axis[0,1]
    ax2 = axis[0,2]
    ax3 = axis[0,3]
    ax4 = axis[0,4]
    ax5 = axis[0,5]
    ax1.set_title('max. pot. temp.')
    ax2.set_title('max. w')
    ax3.set_title('max. vorticity')
    ax4.set_title('max. CP height')
    ax5.set_title('CP volume')
    for i, dTh in enumerate(dTh_range_A):
        rstar = rstar_range_A[i]
        zstar = zstar_range_A[i]
        rootname = 'dTh'+str(dTh) + '_z'+str(zstar) + '_r'+str(rstar)
        lbl = 'dTh='+str(dTh) + 'K, z*='+str(zstar) + ', r*='+str(rstar)
        filename_CPheight = 'CP_height_'+rootname+'_sth0.5.nc'
        filename_CPvol = 'CP_volume_'+rootname+'.nc'
        path_in = os.path.join(path_root, rootname)
        print('')
        print(rootname)
        print(os.path.join(path_in, 'data_analysis', filename_stats))
        stats_root = nc.Dataset(os.path.join(path_in, 'data_analysis', filename_stats))
        time_stats = stats_root.groups['timeseries'].variables['time'][:]
        w_ = stats_root.groups['stats'].variables['w'][:,:-1,:]
        w_max = np.amax(np.amax(w_, axis=-1), axis=-1)
        w_max_k0 = np.amax(w_[:,:,lvl_w], axis=1)
        w_min = np.amin(np.amin(w_, axis=-1), axis=-1)
        w_min_k0 = np.amin(w_[:,:,lvl_w], axis=1)
        del w_
        v_rad_ = stats_root.groups['stats'].variables['v_rad'][:,:-1,:]
        v_rad_max = np.amax(np.amax(v_rad_, axis=-1), axis=-1)
        v_rad_max_k0 = np.amax(v_rad_[:,:,lvl], axis=1)
        del v_rad_
        s_ = stats_root.groups['stats'].variables['s'][:,:-1,:]
        # s_min = np.amin(np.amin(s_, axis=-1), axis=-1)
        theta_ = thetas_c(s_, 0)
        theta_min = np.amin(np.amin(theta_, axis=-1), axis=-1)
        theta_min_k0 = np.amin(theta_[:,:,lvl], axis=1)
        del s_, theta_
        stats_root.close()

        # vorticity from azimuthally averaged velocity fields (v_rad, v_tan, w)
        print(os.path.join(path_in, 'fields_vorticity', filename_vort))
        vort_root = nc.Dataset(os.path.join(path_in, 'fields_vorticity', filename_vort))
        time_vort = vort_root.groups['timeseries'].variables['time'][:]
        vort_phi_max = vort_root.groups['timeseries'].variables['vort_phi_max'][:]
        vort_phi_min = vort_root.groups['timeseries'].variables['vort_phi_min'][:]
        vort_root.close()

        print(os.path.join(path_in, 'data_analysis', filename_CPheight))
        root = nc.Dataset(os.path.join(path_in, 'data_analysis', filename_CPheight))
        time_geom = root.groups['timeseries'].variables['time'][:]
        CP_height_2D = root.groups['fields_2D'].variables['CP_height_2d'][:,:,:]
        CP_height_max = root.groups['timeseries'].variables['CP_height_max'][:]
        root.close()
    #     del CP_height_
        root = nc.Dataset(os.path.join(path_in, 'data_analysis', filename_CPvol))
        time_geom_ = root.groups['timeseries'].variables['time'][:]
        if time_geom.any() != time_geom_.any():
            print('!!!!! different timeseries')
        CP_vol = root.groups['timeseries'].variables['CP_vol_sth0.5'][:]
        root.close()

        max_range[1] = np.maximum(np.amax(theta_min_k0), max_range[1])
        max_range[2] = np.maximum(np.amax(w_max), max_range[2])
        max_range[3] = np.maximum(np.amax(vort_phi_max), max_range[3])
        min_range[1] = np.minimum(np.amin(theta_min), min_range[1])
        # ax0.plot(time_stats, s_min, '-o', color=colorlist3[i], label=lbl)
        ax1.plot(time_stats, theta_min, '-o', color=colorlist3[i], label=lbl)
        ax2.plot(time_stats, w_max, '-o', color=colorlist3[i], label=lbl)
        ax3.plot(time_vort, vort_phi_max, '-o', color=colorlist3[i], label=lbl)
        ax4.plot(time_geom, CP_height_max, '-o', color=colorlist3[i], label=lbl)
        ax5.plot(time_geom, 1e-9*CP_vol, '-o', color=colorlist3[i], label=lbl)

    # ax0.set_ylabel('entropy [J/K]')
    ax1.set_ylabel('pot. temperature [K]')
    ax2.set_ylabel('max(w) [m/s]')
    ax3.set_ylabel('max(vorticity) [1/s]')
    ax4.set_ylabel('max(CP height) [m]')
    ax5.set_ylabel('CP volume [m^3]')

    # ax1.legend(loc='best')
    ax1.legend(loc='upper center', bbox_to_anchor=(-2., .75),
               fancybox=True, shadow=False, ncol=1, fontsize=10)

    ax0 = axis[1, 0]
    ax1 = axis[1, 1]
    ax2 = axis[1, 2]
    ax3 = axis[1, 3]
    ax4 = axis[1, 4]
    ax5 = axis[1, 5]
    for i, dTh in enumerate(dTh_range_B):
        rstar = rstar_range_B[i]
        zstar = zstar_range_B[i]
        rootname = 'dTh'+str(dTh) + '_z'+str(zstar) + '_r'+str(rstar)
        lbl = 'dTh='+str(dTh) + 'K, z*='+str(zstar) + ', r*='+str(rstar)
        filename_CPheight = 'CP_height_'+rootname+'_sth0.5.nc'
        filename_CPvol = 'CP_volume_'+rootname+'.nc'
        path_in = os.path.join(path_root, rootname)
        print('')
        print(rootname)
        print(os.path.join(path_in, 'data_analysis', filename_stats))
        # filename stats: variables = [nt,nr,nz]
        stats_root = nc.Dataset(os.path.join(path_in, 'data_analysis', filename_stats))
        time_stats = stats_root.groups['timeseries'].variables['time'][:]
        w_ = stats_root.groups['stats'].variables['w'][:,:-1,:]
        w_max = np.amax(np.amax(w_, axis=-1), axis=-1)
        w_max_k0 = np.amax(w_[:,:,lvl_w], axis=1)
        w_min = np.amin(np.amin(w_, axis=-1), axis=-1)
        w_min_k0 = np.amin(w_[:,:,lvl_w], axis=1)
        del w_
        v_rad_ = stats_root.groups['stats'].variables['v_rad'][:,:-1,:]
        v_rad_max = np.amax(np.amax(v_rad_, axis=-1), axis=-1)
        v_rad_max_k0 = np.amax(v_rad_[:,:,lvl], axis=1)
        del v_rad_
        s_ = stats_root.groups['stats'].variables['s'][:,:-1,:]
        # s_min = np.amin(np.amin(s_, axis=-1), axis=-1)
        theta_ = thetas_c(s_, 0)
        theta_min = np.amin(np.amin(theta_, axis=-1), axis=-1)
        theta_min_k0 = np.amin(theta_[:,:,lvl], axis=1)
        del s_, theta_
        stats_root.close()

        # vorticity from azimuthally averaged velocity fields (v_rad, v_tan, w)
        print(os.path.join(path_in, 'fields_vorticity', filename_vort))
        vort_root = nc.Dataset(os.path.join(path_in, 'fields_vorticity', filename_vort))
        time_vort = vort_root.groups['timeseries'].variables['time'][:]
        vort_phi_max = vort_root.groups['timeseries'].variables['vort_phi_max'][:]
        vort_phi_min = vort_root.groups['timeseries'].variables['vort_phi_min'][:]
        vort_root.close()

        root = nc.Dataset(os.path.join(path_in, 'data_analysis', filename_CPheight))
        time_geom = root.groups['timeseries'].variables['time'][:]
        CP_height_2D = root.groups['fields_2D'].variables['CP_height_2d'][:,:,:]
        CP_height_max = root.groups['timeseries'].variables['CP_height_max'][:]
        root.close()
        print(os.path.join(path_in, 'data_analysis', filename_CPvol))
        root = nc.Dataset(os.path.join(path_in, 'data_analysis', filename_CPvol))
        time_geom_ = root.groups['timeseries'].variables['time'][:]
        if time_geom.any() != time_geom_.any():
            print('!!!!! different timeseries')
        CP_vol = root.groups['timeseries'].variables['CP_vol_sth0.5'][:]
        root.close()

        max_range[1] = np.maximum(np.amax(theta_min_k0), max_range[1])
        max_range[2] = np.maximum(np.amax(w_max), max_range[2])
        max_range[3] = np.maximum(np.amax(vort_phi_max), max_range[3])
        max_range[4] = np.maximum(np.amax(CP_height_max), max_range[4])
        max_range[5] = np.maximum(1e-9*np.amax(CP_vol), max_range[5])
        min_range[1] = np.minimum(np.amin(theta_min), min_range[1])
        min_range[2] = 0.
        min_range[3] = 0.
        min_range[4] = 0.
        min_range[5] = 0.
        # ax0.plot(time_stats, s_min, '-o', color=colorlist5[i], label=lbl)
        ax1.plot(time_stats, theta_min, '-o', color=colorlist5[i], label=lbl)
        ax2.plot(time_stats, w_max, '-o', color=colorlist5[i], label=lbl)
        ax3.plot(time_vort, vort_phi_max, '-o', color=colorlist5[i], label=lbl)
        ax4.plot(time_geom, CP_height_max, '-o', color=colorlist5[i], label=lbl)
        ax5.plot(time_geom, 1e-9*CP_vol, '-o', color=colorlist5[i], label=lbl)
    # ax0.set_ylabel('entropy [J/K]')
    ax1.set_ylabel('pot. temperature [K]')
    ax2.set_ylabel('max(w) [m/s]')
    ax3.set_ylabel('max(vorticity) [1/s]')
    ax4.set_ylabel('max(CP height) [m]')
    ax5.set_ylabel('CP volume [km^3]')

    print('')
    ax1.legend(loc='upper center', bbox_to_anchor=(-2., .75),
               fancybox=True, shadow=False, ncol=1, fontsize=10)
    for ax in axis.flat:
        ax.set_xlim(0,3600)
        ax.set_xticks(np.arange(0, 3700, step=900.))
        x_ticks = [ti/3600 for ti in ax.get_xticks()]
        # print('ticks ', ax.get_xticks())
        # print('      ', x_ticks)
        ax.set_xticklabels(x_ticks)
        ax.set_ylim(min_range[i], np.amax(max_range[i]))
        # y_ticks = ax.get_yticks()
        # ax.set_yticklabels(y_ticks)
        # for label in ax.xaxis.get_ticklabels()[1::2]:
        #     label.set_visible(False)
    print('')
    for ax in axis[1,1:].flat:
        ax.set_xlabel('time [h]')
    max_range[1] += .5
    max_range[3] = .13
    max_range[5] += 0.1
    for i in range(1,ncols):
        axis[0,i].set_ylim(min_range[i],np.amax(max_range[i]))
        axis[1,i].set_ylim(min_range[i],np.amax(max_range[i]))
    for ax in axis.flat:
        y_ticks = ax.get_yticks()
        ax.set_yticklabels(y_ticks)
    plt.subplots_adjust(bottom=0.075, right=.99, left=0.14, top=0.95, wspace=0.3, hspace=0.1)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return


# ----------------------------------------------------------------------
def cpm_c(qt):
    cpd = 1004.0
    cpv = 1859.0
    return (1.0-qt) * cpd + qt * cpv

def thetas_c(s, qt):
    T_tilde = 298.15
    sd_tilde = 6864.8
    sv_tilde = 10513.6
    return T_tilde*np.exp((s-(1.0-qt)*sd_tilde - qt*sv_tilde)/cpm_c(qt))

# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()