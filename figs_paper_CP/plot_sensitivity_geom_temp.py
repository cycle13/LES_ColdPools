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
    path_root = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/'
    path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single/figs_run2_dx100m/'
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)

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
    zstar_range_B = [500, 1000, 1600, 2000, 2500]

    filename_stats = 'stats_radial_averaged.nc'
    filename_vort = 'Stats_vorticity_phi.nc'
    lvl_w = 1
    lvl = 0


    # plot_sensitivity_plots_temperature(dTh_range_A, rstar_range_A, zstar_range_A, path_root,
    #                                    filename_stats, filename_vort, filename_CP_height, filename_CP_vol,
    #                                    colorlist3, path_out_figs)
    #
    # plot_sensitivity_plots_geometry(dTh_range_B, rstar_range_B, zstar_range_B, path_root,
    #                                 filename_stats, filename_vort, filename_CP_height, filename_CP_vol,
    #                                 colorlist3, path_out_figs)

    plot_sensitivity_plots_all(dTh_range_A, rstar_range_A, zstar_range_A,
                               dTh_range_B, rstar_range_B, zstar_range_B,
                               filename_stats, filename_vort, filename_CP_height, filename_CP_vol,
                               colorlist3, path_root, path_out_figs)
    return





def plot_sensitivity_plots_temperature(dTh_range_A, rstar_range_A, zstar_range_A, path_root,
                                       filename_stats, filename_vort, filename_CP_height, filename_CP_vol,
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

        ax0.plot(time_stats, s_min, '-o', color=colorlist3[i], label=rootname)
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
                                           filename_stats, filename_vort, filename_CP_height, filename_CP_vol,
                                           colorlist3, path_out_figs):
    fig_name = 'sensitivity_plots_geometry.png'
    ncols = 6
    fig, axis = plt.subplots(1, ncols, sharex='all', figsize=(ncols*3, 5))
    ax1 = axis[1]
    ax2 = axis[2]
    ax3 = axis[3]

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

        ax0.plot(time_stats, s_min, '-o', color=colorlist3[i], label=rootname)
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
                               filename_stats, filename_vort, filename_CP_height, filename_CP_vol,
                               colorlist3, path_root, path_out_figs):

    fig_name = 'sensitivity_plots_all.png'
    ncols = 6
    nrows = 2
    fig, axis = plt.subplots(nrows, ncols, figsize=(ncols*3, nrows*5), sharex='col')
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

        ax0.plot(time_stats, s_min, '-o', color=colorlist3[i], label=rootname)
        ax1.plot(time_stats, theta_min, '-o', color=colorlist3[i], label=rootname)
        ax2.plot(time_stats, w_max, '-o', color=colorlist3[i], label=rootname)
        # ax3.plot(time_vort, vort_phi_max, '-o', color=colorlist3[i], label=rootname)
        ax4.plot(time_geom, CP_height_max, '-o', color=colorlist3[i], label=rootname)
        ax5.plot(time_geom, CP_vol, '-o', color=colorlist3[i], label=rootname)
        ax0.set_ylabel('entropy [J/K]')
        ax1.set_ylabel('pot. temperature [K]')
        ax2.set_ylabel('max(w) [m/s]')
        ax3.set_ylabel('max(vorticity) [1/s]')
        ax4.set_ylabel('max(CP height) [m]')
        ax5.set_ylabel('CP volume [m^3]')

    ax1.legend(loc='best')

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

        ax0.plot(time_stats, s_min, '-o', color=colorlist3[i], label=rootname)
        ax1.plot(time_stats, theta_min, '-o', color=colorlist3[i], label=rootname)
        ax2.plot(time_stats, w_max, '-o', color=colorlist3[i], label=rootname)
        # ax3.plot(time_vort, vort_phi_max, '-o', color=colorlist3[i], label=rootname)
        ax4.plot(time_geom, CP_height_max, '-o', color=colorlist3[i], label=rootname)
        ax5.plot(time_geom, CP_vol, '-o', color=colorlist3[i], label=rootname)
        ax0.set_ylabel('entropy [J/K]')
        ax1.set_ylabel('pot. temperature [K]')
        ax2.set_ylabel('max(w) [m/s]')
        ax3.set_ylabel('max(vorticity) [1/s]')
        ax4.set_ylabel('max(CP height) [m]')
        ax5.set_ylabel('CP volume [m^3]')

    # ax1.legend(loc='best')
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
               fancybox=True, shadow=True, ncol=1, fontsize=10)
    for ax in axis[1,:].flat:
        ax.set_xlabel('time [h]')
    # fig.suptitle(title, fontsize=18)
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.06, top=0.9, wspace=0.15)
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