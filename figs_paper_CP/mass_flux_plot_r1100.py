import netCDF4 as nc
import argparse
import json as simplejson
import os
import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import matplotlib.cbook as cbook
import time

execfile('settings.py')
label_size = 15
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['axes.labelsize'] = 18

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("--level")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    case_name_1CP = 'ColdPoolDry_single_3D'
    case_name_2CP = 'ColdPoolDry_double_3D'
    case_name_3CP = 'ColdPoolDry_triple_3D'

    dTh = 5
    zstar = 1000
    rstar = 1100
    rst = str(rstar)
    sep = d_range[rst][0]
    print('Parameters: ')
    print('d-range: ' + str(d_range[rst]))
    print('')

    case = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
    case_xCP = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar) + '_d' + str(sep) + 'km'
    path_1CP = '/nbi/ac/conv4/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/single_3D_noise/run5_PE_scaling_dx100m/'
    path_2CP = '/nbi/ac/coag/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/double_3D/'
    path_3CP = '/nbi/ac/cond2/meyerbe/ColdPools/3D_sfc_fluxes_off/triple_3D/'
    print('')
    print('Case: ' + case)
    print('path 1CP:   ' + path_1CP)
    print('path 2CP:   ' + path_2CP)
    print('path 3CP:   ' + path_3CP)
    path_data = os.path.join(path_3CP, 'data_analysis')
    # path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single'
    path_out_figs = '/nbi/home/meyerbe/paper_CP'
    print('Path Figures: ' + path_out_figs)
    print('')

    if args.level:
        zlevel = np.int(args.level)
    else:
        zlevel = 1000
    print('')

    # plotting parameters
    xmin_3CP = np.zeros(3, dtype=np.int)
    xmax_3CP = np.zeros(3, dtype=np.int)
    ymin_3CP = np.zeros(3, dtype=np.int)
    ymax_3CP = np.zeros(3, dtype=np.int)
    if rstar == 1100:
        xmin_3CP[0] = 100
        xmax_3CP[0] = 350
        ymin_3CP[0] = 50
        xmin_3CP[1] = 200
        xmax_3CP[1] = 500
        ymin_3CP[1] = 100
        xmin_3CP[2] = 200
        xmax_3CP[2] = 600
        ymin_3CP[2] = 50


    for i,sep in enumerate(d_range[rst][:]):
        print('--- d='+str(sep)+'km ---')
        case_xCP = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar) + '_d' + str(sep) + 'km'

        nx_1CP, nx_2CP, nx_3CP, dt_fields_1CP, dt_fields_2CP, dt_fields_3CP = set_input_output_parameters(args, case_name_1CP, case_name_2CP, case_name_3CP,
                                                                     path_1CP, path_2CP, path_3CP,
                                                                     case, case_xCP)

        ic_arr_1CP, jc_arr_1CP, ic_arr_2CP, jc_arr_2CP, \
                ic_arr_3CP, jc_arr_3CP, \
                ic_2CP, jc_2CP, ic_3CP, jc_3CP = define_geometry(nx_1CP, nx_2CP, nx_3CP, sep)

        times = np.arange(t_ini[rst][i], t_final[rst][i], 100)
        print('times: ', times)
        # print('CP geometry:')
        # print(ic_arr)
        # print(jc_arr)
        print('')


        ''' Mass Flux '''
        # read in flux (2D-field)
        filename_data = 'mass_flux_z' + str(zlevel) + '.nc'
        path_in = os.path.join(path_1CP, case, 'data_analysis', filename_data)
        print(path_in)
        root_in = nc.Dataset(path_in, 'r')
        mass_flux_1CP = root_in.groups['fields_2D'].variables['mass_flux_2D'][:, :, :]
        mass_flux_pos_1CP = root_in.groups['fields_2D'].variables['mass_flux_2D_positive'][:, :, :]
        time_mf_1CP = root_in.groups['timeseries'].variables['time'][:]
        root_in.close()
        path_in = os.path.join(path_2CP, case_xCP, 'data_analysis', filename_data)
        print(path_in)
        root_in = nc.Dataset(path_in, 'r')
        mass_flux_2CP = root_in.groups['fields_2D'].variables['mass_flux_2D'][:, :, :]
        mass_flux_pos_2CP = root_in.groups['fields_2D'].variables['mass_flux_2D_positive'][:, :, :]
        time_mf_2CP = root_in.groups['timeseries'].variables['time'][:]
        root_in.close()
        path_in = os.path.join(path_3CP, case_xCP, 'data_analysis', filename_data)
        print(path_in)
        root_in = nc.Dataset(path_in, 'r')
        mass_flux_3CP = root_in.groups['fields_2D'].variables['mass_flux_2D'][:, :, :]
        mass_flux_pos_3CP = root_in.groups['fields_2D'].variables['mass_flux_2D_positive'][:, :, :]
        time_mf_3CP = root_in.groups['timeseries'].variables['time'][:]
        root_in.close()
        ''' averaged mass flux '''
        tmin = t_ini[rst][i]
        tmax = t_final[rst][i]
        delta = 5       # averaging over band of withd delta_y=2*delta*dy=600m
        print('Mass flux: averaging from tmin='+str(tmin)+' to tmax='+str(tmax))
        it0 = np.where(time_mf_1CP == tmin)[0][0]#np.int(tmin / dt_fields_1CP)
        it1 = np.where(time_mf_1CP == tmax)[0][0]#np.int(tmax / dt_fields_1CP)
        print('MF times 1CP: ', it0, it1)
        MF_1CP = np.sum(mass_flux_1CP[it0:it1,:,:], axis=0)     # accumulate over time
        MF_pos_1CP = np.sum(mass_flux_pos_1CP[it0:it1,:,:], axis=0)     # accumulate over time
        MF_mean_1CP = np.mean(MF_1CP[:,jc_arr_1CP[0]-delta:jc_arr_1CP[0]+delta], axis=1)
        MF_mean_1CP_ = np.mean(MF_1CP[ic_arr_1CP[0]-delta:ic_arr_1CP[0]+delta,:], axis=0)
        MF_pos_mean_1CP = np.mean(MF_pos_1CP[:,jc_arr_1CP[0]-delta:jc_arr_1CP[0]+delta], axis=1)
        MF_pos_mean_1CP_ = np.mean(MF_pos_1CP[ic_arr_1CP[0]-delta:ic_arr_1CP[0]+delta,:], axis=0)
        it0 = np.where(time_mf_2CP == tmin)[0][0]  # np.int(tmin / dt_fields_2CP)
        it1 = np.where(time_mf_2CP == tmax)[0][0]  # np.int(tmax / dt_fields_2CP)
        print('MF times 2CP: ', it0, it1)
        MF_2CP = np.sum(mass_flux_2CP[it0:it1, :, :], axis=0)  # accumulate over time
        MF_pos_2CP = np.sum(mass_flux_pos_2CP[it0:it1, :, :], axis=0)  # accumulate over time
        MF_mean_2CP = np.mean(MF_2CP[ic_2CP-delta:ic_2CP+delta,:], axis=0)
        MF_mean_2CP_ = np.mean(MF_2CP[:,jc_2CP-delta:jc_2CP+delta], axis=1)
        MF_pos_mean_2CP = np.mean(MF_pos_2CP[ic_2CP-delta:ic_2CP+delta,:], axis=0)
        MF_pos_mean_2CP_ = np.mean(MF_pos_2CP[:,jc_2CP-delta:jc_2CP+delta], axis=1)
        it0 = np.where(time_mf_3CP == tmin)[0][0]  # np.int(tmin / dt_fields_3CP)
        it1 = np.where(time_mf_3CP == tmax)[0][0]  # np.int(tmax / dt_fields_3CP)
        print('MF times 3CP: ', it0, it1)
        MF_3CP = np.sum(mass_flux_3CP[it0:it1, :, :], axis=0)  # accumulate over time
        MF_pos_3CP = np.sum(mass_flux_pos_3CP[it0:it1, :, :], axis=0)  # accumulate over time
        MF_mean_3CP = np.mean(MF_3CP[:,jc_3CP-delta:jc_3CP+delta], axis=1)
        MF_mean_3CP_ = np.mean(MF_3CP[ic_3CP-delta:ic_3CP+delta,:], axis=0)
        MF_pos_mean_3CP = np.mean(MF_pos_3CP[:,jc_3CP-delta:jc_3CP+delta], axis=1)
        MF_pos_mean_3CP_ = np.mean(MF_pos_3CP[ic_3CP-delta:ic_3CP+delta,:], axis=0)

        print('')

        ''' CP height '''
        filename_CP_height = 'CP_height_' + case + '_sth0.5' + '.nc'
        print(os.path.join(path_1CP, case, 'data_analysis', filename_CP_height))
        root_in = nc.Dataset(os.path.join(path_1CP, case, 'data_analysis', filename_CP_height), 'r')
        time_CPheight_1CP = root_in.groups['timeseries'].variables['time'][:]
        CP_height_1CP = root_in.groups['fields_2D'].variables['CP_height_2d'][:, :, :]
        root_in.close()
        filename_CP_height = 'CP_height_' + case_xCP + '_sth0.5' + '.nc'
        print(os.path.join(path_2CP, case_xCP, 'data_analysis', filename_CP_height))
        root_in = nc.Dataset(os.path.join(path_2CP, case_xCP, 'data_analysis', filename_CP_height), 'r')
        time_CPheight_2CP = root_in.groups['timeseries'].variables['time'][:]
        CP_height_2CP = root_in.groups['fields_2D'].variables['CP_height_2d'][:, :, :]
        root_in.close()
        print(os.path.join(path_3CP, case_xCP, 'data_analysis', filename_CP_height))
        root_in = nc.Dataset(os.path.join(path_3CP, case_xCP, 'data_analysis', filename_CP_height), 'r')
        time_CPheight_3CP = root_in.groups['timeseries'].variables['time'][:]
        CP_height_3CP = root_in.groups['fields_2D'].variables['CP_height_2d'][:, :, :]
        root_in.close()

        global dt_CP_height_1CP, dt_CP_height_2CP, dt_CP_height_3CP
        dt_CP_height_1CP = time_CPheight_1CP[1] - time_CPheight_1CP[0]
        dt_CP_height_2CP = time_CPheight_2CP[1] - time_CPheight_2CP[0]
        dt_CP_height_3CP = time_CPheight_3CP[1] - time_CPheight_3CP[0]
        print('dt CP height: ', dt_CP_height_1CP, dt_CP_height_2CP, dt_CP_height_3CP)
        print('')



        plot_collision_massflux_CPheight_test(CP_height_1CP, CP_height_2CP, CP_height_3CP,
                                              time_CPheight_2CP, time_CPheight_3CP,
                                              t_ini[rst][i], t_2CP[rst][i], t_3CP[rst][i], t_final[rst][i],
                                              ic_2CP, jc_2CP, ic_3CP, jc_3CP,
                                              xmin_3CP[i], xmax_3CP[i],
                                              case, case_xCP,
                                              times, path_out_figs)


        fig_name = 'collisions_massflux_CPheight_' + case_xCP + '_123.png'
        plot_collision_massflux_CPheight(CP_height_1CP, CP_height_2CP, CP_height_3CP,
                                         time_CPheight_1CP, time_CPheight_2CP, time_CPheight_3CP,
                                         MF_1CP, MF_2CP, MF_3CP, MF_mean_1CP, MF_mean_2CP, MF_mean_3CP,
                                         t_ini[rst][i], t_2CP[rst][i], t_3CP[rst][i], t_final[rst][i],
                                         delta, ic_2CP, jc_2CP, ic_3CP, jc_3CP,
                                         ic_arr_1CP, jc_arr_1CP, ic_arr_2CP, jc_arr_2CP, ic_arr_3CP, jc_arr_3CP,
                                         xmin_3CP[i], xmax_3CP[i], ymin_3CP[i],
                                         times, nx_1CP, nx_2CP, nx_3CP, path_out_figs, fig_name)


        fig_name = 'collisions_massflux_pos_CPheight_' + case_xCP + '_123.png'
        plot_collision_massflux_CPheight(MF_pos_1CP, MF_pos_2CP, MF_pos_3CP, MF_pos_mean_1CP, MF_pos_mean_2CP, MF_pos_mean_3CP,
                                         t_ini[rst][i], t_2CP[rst][i], t_3CP[rst][i], t_final[rst][i],
                                         delta, ic_2CP, jc_2CP, ic_3CP, jc_3CP,
                                         ic_arr_1CP, jc_arr_1CP, ic_arr_2CP, jc_arr_2CP, ic_arr_3CP, jc_arr_3CP,
                                         xmin_3CP[i], xmax_3CP[i], ymin_3CP[i],
                                         times, nx_1CP, nx_2CP, nx_3CP, path_out_figs, fig_name)


        fig_name = 'collisions_massflux_CPheight_' + case_xCP + '_3panels.png'
        plot_collision_massflux_CPheight_3panels(CP_height_1CP, CP_height_2CP, CP_height_3CP,
                                         time_CPheight_1CP, time_CPheight_2CP, time_CPheight_3CP,
                                         MF_1CP, MF_2CP, MF_3CP, MF_mean_1CP, MF_mean_2CP, MF_mean_3CP,
                                         t_ini[rst][i], t_2CP[rst][i], t_3CP[rst][i], t_final[rst][i],
                                         delta, ic_2CP, jc_2CP, ic_3CP, jc_3CP,
                                         ic_arr_2CP, jc_arr_2CP, ic_arr_3CP, jc_arr_3CP,
                                         xmin_3CP[i], xmax_3CP[i], ymin_3CP[i],
                                         case, case_xCP,
                                         times, path_out_figs, fig_name)
    return


# ----------------------------------------------------------------------
def plot_collision_massflux_CPheight_test(CP_height_1CP, CP_height_2CP, CP_height_3CP,
                                          time_CPheight_2CP, time_CPheight_3CP,
                                          t_ini, t_2CP, t_3CP, t_final,
                                          ic_2CP, jc_2CP, ic_3CP, jc_3CP,
                                          xmin_3CP, xmax_3CP,
                                          case, case_xCP,
                                          times, path_out_figs):
    ncol = 4
    nrow = 5
    fig, axis = plt.subplots(nrow, ncol, figsize=(ncol*5, nrow*5))
    # ax = axis[0, 0]
    # # for it_CP,t0 in enumerate(time_CPheight_2CP):
    # #     ax.plot(CP_height_2CP[it_CP, ic_arr_2CP, :], color=cm_bwr(count_color), label='t=' + str(t0))
    ax = axis[1, 0]
    it = np.int(t_ini / dt_CP_height_1CP)
    cf = ax.contourf(CP_height_1CP[it, :, :], cmap=cm_gray)
    # cf = ax.pcolor(CP_height_2CP[it,:,:], cmap=cm_gray)
    plt.colorbar(cf, ax=ax, shrink=0.75)
    ax.set_title('t=' + str(t_ini))
    # ax.plot([0, nx_2CP[1]], [ic_2CP, ic_2CP], '-y')
    ax = axis[1, 1]
    it = np.int(t_2CP / dt_CP_height_1CP)
    cf = ax.contourf(CP_height_1CP[it, :, :], cmap=cm_gray)
    # cf = ax.pcolor(CP_height_2CP[it,:,:], cmap=cm_gray)
    plt.colorbar(cf, ax=ax, shrink=0.75)
    ax.set_title('t=' + str(t_2CP))
    ax = axis[1, 2]
    it = np.int(t_3CP / dt_CP_height_1CP)
    cf = ax.contourf(CP_height_1CP[it, :, :], cmap=cm_gray)
    # cf = ax.pcolor(CP_height_2CP[it,:,:], cmap=cm_gray)
    plt.colorbar(cf, ax=ax, shrink=0.75)
    ax.set_title('t=' + str(t_3CP))
    ax = axis[1, 3]
    it = np.int(t_final / dt_CP_height_1CP)
    cf = ax.contourf(CP_height_1CP[it, :, :], cmap=cm_gray)
    # cf = ax.pcolor(CP_height_2CP[it,:,:], cmap=cm_gray)
    plt.colorbar(cf, ax=ax, shrink=0.75)
    ax.set_title('t=' + str(t_final))
    for ax in axis[1, :].flat:
        #     ax.set_ylim(100, 400)
        ax.set_aspect('equal')

    ax = axis[2, 0]
    it = np.int(t_ini / dt_CP_height_2CP)
    cf = ax.contourf(CP_height_2CP[it, :, :], cmap=cm_gray)
    # cf = ax.pcolor(CP_height_2CP[it,:,:], cmap=cm_gray)
    plt.colorbar(cf, ax=ax, shrink=0.75)
    ax.set_title('t=' + str(t_ini))
    ax = axis[2, 1]
    it = np.int(t_2CP / dt_CP_height_2CP)
    cf = ax.contourf(CP_height_2CP[it, :, :], cmap=cm_gray)
    # cf = ax.pcolor(CP_height_2CP[it,:,:], cmap=cm_gray)
    plt.colorbar(cf, ax=ax, shrink=0.75)
    ax.set_title('t=' + str(t_2CP))
    ax = axis[2, 2]
    it = np.int(t_3CP / dt_CP_height_2CP)
    cf = ax.contourf(CP_height_2CP[it, :, :], cmap=cm_gray)
    plt.colorbar(cf, ax=ax, shrink=0.75)
    ax.set_title('t=' + str(t_3CP))
    ax = axis[2, 3]
    it = np.int(t_final / dt_CP_height_2CP)
    cf = ax.contourf(CP_height_2CP[it, :, :], cmap=cm_gray)
    plt.colorbar(cf, ax=ax, shrink=0.75)
    ax.set_title('t=' + str(t_final))
    for ax in axis[2,:].flat:
        ax.plot([0, nx_2CP[1]], [ic_2CP, ic_2CP], '-y')
        # ax.set_ylim(100,400)
        ax.set_aspect('equal')

    ax = axis[3, 0]
    it = np.int(t_ini / dt_CP_height_3CP)
    cf = ax.contourf(CP_height_3CP[it, :, :], cmap=cm_gray)
    # cf = ax.pcolor(CP_height_2CP[it,:,:], cmap=cm_gray)
    plt.colorbar(cf, ax=ax, shrink=0.75)
    ax.set_title('t=' + str(t_ini))
    ax = axis[3, 1]
    it = np.int(t_2CP / dt_CP_height_3CP)
    cf = ax.contourf(CP_height_3CP[it, :, :], cmap=cm_gray)
    plt.colorbar(cf, ax=ax, shrink=0.75)
    ax.set_title('t=' + str(t_2CP))
    ax = axis[3, 2]
    it = np.int(t_3CP / dt_CP_height_3CP)
    cf = ax.contourf(CP_height_3CP[it, :, :], cmap=cm_gray)
    plt.colorbar(cf, ax=ax, shrink=0.75)
    ax.set_title('t=' + str(t_3CP))
    ax = axis[3, 3]
    it = np.int(t_final / dt_CP_height_3CP)
    cf = ax.contourf(CP_height_3CP[it, :, :], cmap=cm_gray)
    plt.colorbar(cf, ax=ax, shrink=0.75)
    ax.set_title('t=' + str(t_final))
    for ax in axis[3,:].flat:
        ax.plot([0, nx_3CP[1]], [ic_3CP, ic_3CP], '-y')
        ax.set_xlim(xmin_3CP, xmax_3CP)
        ax.set_ylim(xmin_3CP, xmax_3CP)
        ax.set_aspect('equal')

    #     it = np.int(t0 / dt_CP_height_3CP)
    #     ax = axis[2, i]
    #     cf = ax.contourf(CP_height_3CP[it, :, :].T, cmap=cm_gray)
    #     # cf = ax.pcolor(CP_height_3CP[it, :, :].T, cmap=cm_gray)
    #     plt.colorbar(cf, ax=ax, shrink=0.75)
    #     ax.set_title('t=' + str(time_CPheight_3CP[it]))
    #     ax.plot([0, nx_3CP[0]], [jc_3CP, jc_3CP], '-y')
    #     ax.set_xlim(160, 500)
    #     ax.set_ylim(160, 500)
    #
    # for ax in axis[1:3, :].flatten():
    #     ax.set_aspect('equal')
    #
    print('times: ', times)
    for it, t0 in enumerate(times):
        it_2CP = np.where(time_CPheight_2CP == t0)[0][0]
        it_2CP_ = np.int(t0 / dt_CP_height_2CP)
        it_3CP = np.where(time_CPheight_3CP == t0)[0][0]
        #     it_3CP_ = np.int(t0 / dt_CP_height_3CP)
        #     print('comp: ', it_2CP, it_2CP_, it_3CP, it_3CP_)
        #     if t0 > time_bins[0] and t0 <= time_bins[-1]:
        if t0 > t_ini and t0 <= t_final:
            # var = read_in_netcdf_fields(var_name, os.path.join(path_fields, str(t0)+'.nc'))
            #         if t0 <= time_bins[1]:
            if t0 < t_2CP:
                cm = cm_gray
                count_color = (np.double(t0) - t_ini) / (t_2CP - t_ini) * 0.8
                for ax in axis[0,:].flat:
                    ax.plot(CP_height_2CP[it_2CP, ic_2CP, :], color=cm(count_color), label='t=' + str(t0))
                for ax in axis[4,:].flat:
                    ax.plot(CP_height_3CP[it_3CP, :, jc_3CP], color=cm(count_color), label='t=' + str(t0))
            elif t0 <= t_3CP:
                cm = cm_bwr_r
                count_color = (np.double(t0) - t_2CP) / (t_3CP - t_2CP) * 0.5 + 0.5
                for ax in axis[0,1:].flat:
                    ax.plot(CP_height_2CP[it_2CP, ic_2CP, :], color=cm(count_color), label='t=' + str(t0))
                for ax in axis[4,1:].flat:
                    ax.plot(CP_height_3CP[it_3CP, :, jc_3CP], color=cm(count_color), label='t=' + str(t0))
            elif t0 <= t_final:
                cm = cm_bwr
                count_color = (np.double(t0) - t_3CP) / (t_final - t_3CP) * 0.45 + 0.5
                for ax in axis[0,2:].flat:
                    ax.plot(CP_height_2CP[it_2CP, ic_2CP, :], color=cm(count_color), label='t=' + str(t0))
                for ax in axis[4,2:].flat:
                    ax.plot(CP_height_3CP[it_3CP, :, jc_3CP], color=cm(count_color), label='t=' + str(t0))

            # for i in range(ncol):
                # ax = axis[0, i]
                # # ax.plot([0, nx_2CP[0]], [0, 0], '0.5', linewidth=1)
                # ax.plot(CP_height_2CP[it_2CP, ic_2CP, :], color=cm(count_color), label='t=' + str(t0))
                # ax = axis[4, i]
                # ax.plot(CP_height_3CP[it_3CP, :, jc_3CP], color=cm(count_color), label='t=' + str(t0))
        for ax in axis[4, :].flat:
            ax.set_xlim(xmin_3CP, xmax_3CP)
    plt.subplots_adjust(bottom=0.05, right=.95, left=0.07, top=0.9, hspace=0.2, wspace=0.1)
    fig.savefig(os.path.join(path_out_figs, 'CP_height_test_' + case_xCP + '_123.png'))
    plt.close(fig)
    return




def plot_collision_massflux_CPheight(CP_height_1CP, CP_height_2CP, CP_height_3CP,
                                     time_CPheight_1CP, time_CPheight_2CP, time_CPheight_3CP,
                                     MF_1CP, MF_2CP, MF_3CP, MF_mean_1CP, MF_mean_2CP, MF_mean_3CP,
                                     t_ini, t_2CP, t_3CP, t_final,
                                     delta, ic_2CP, jc_2CP, ic_3CP, jc_3CP,
                                     ic_arr_1CP, jc_arr_1CP, ic_arr_2CP, jc_arr_2CP, ic_arr_3CP, jc_arr_3CP,
                                     xmin, xmax, ymin,
                                     times, nx_1CP, nx_2CP, nx_3CP, path_out_figs, fig_name):

    ncol = 5
    vmin = 1e-2
    vmin = np.amin(MF_3CP)
    vmax = 6.5
    vmax = np.amax(MF_3CP)
    lvls = np.linspace(vmin, vmax, 10)
    fig, axis = plt.subplots(1, ncol, figsize=(ncol * 5, 5))

    ''' CP Height '''
    print('time windows: ', t_ini, t_2CP, t_3CP, t_final)
    for it, t0 in enumerate(times[::2]):
        if t0 >= t_ini and t0 <= t_final:
            # it_3CP = np.where(time_CPheight_3CP == t0)[0][0]
            it_3CP = np.int(t0/dt_CP_height_3CP)
            print('-- it, t0: ', it, t0, '--', t_ini, t_2CP, t_3CP, t_final)
            count_color = 0
            if t0 < t_2CP:
                cm = cm_gray
                count_color = (np.double(t0) - t_ini) / (t_2CP - t_ini) * 0.8
            elif t0 < t_3CP:
                cm = cm_bwr_r
                count_color = (np.double(t0) - t_2CP) / (t_3CP - t_2CP) * 0.5 + 0.5
            elif t0 <= t_final:
                cm = cm_bwr
                count_color = (np.double(t0) - t_3CP) / (t_final - t_3CP) * 0.45 + 0.5
        ax = axis[0]
        if t0 == t_ini or t0 == t_ini+100:
            print('single')
            ax.plot(np.arange(nx_3CP[0]) - ic_3CP, CP_height_3CP[it_3CP, :, jc_3CP], color=cm(count_color), label='single CP')
        elif t0 == t_3CP-100 or t0 == t_3CP-200:
            print('double')
            ax.plot(np.arange(nx_3CP[0]) - ic_3CP, CP_height_3CP[it_3CP, :, jc_3CP], color=cm(count_color), label='double collison')
        elif t0 == t_final-100 or t0 == t_final-200:
            print('triple')
            ax.plot(np.arange(nx_3CP[0]) - ic_3CP, CP_height_3CP[it_3CP, :, jc_3CP], color=cm(count_color), label='triple collison')
        else:
            ax.plot(np.arange(nx_3CP[0]) - ic_3CP, CP_height_3CP[it_3CP, :, jc_3CP], color=cm(count_color))
        ax.set_xlim(xmin-ic_3CP, xmax-ic_3CP)
        ax.legend()
    axis[0].set_ylabel(r'CP Height  [m]')

    ''' Mass Flux '''
    ax = axis[1]
    jc_1CP = np.int(nx_1CP[1] / 2)
    ax.plot(np.arange(nx_1CP[0]) - jc_1CP, MF_mean_1CP, label='single CP', color=cm_gray(.1))  # color=colorlist2[0])
    ax.plot(np.arange(nx_2CP[1]) - jc_2CP, MF_mean_2CP, label='double collision', color=cm_bwr_r(.8))  # color=colorlist2[0])
    ax.plot(np.arange(nx_3CP[0]) - ic_arr_3CP[0], MF_mean_3CP, label='triple collision', color=cm_bwr(.8))  # color=colorlist2[1])
    ax.legend(loc=4)
    ax.set_xlim(xmin - ic_3CP, xmax - ic_3CP)
    ax.set_ylim(-2, np.ceil(np.amax(MF_mean_3CP)))
    axis[1].set_ylabel(r'Integrated Mass Flux  [kg/m$^2$]')


    ''' Mass Flux 2D '''
    ax = axis[2]
    b = ic_arr_3CP[2]-ic_arr_3CP[0]
    pcm = ax.contourf(MF_1CP.T, levels=lvls, cmap=cm_bw)  # cmap='RdBu_r')
    # pcm = ax.pcolormesh(np.arange(nx_1CP[0])-ic_arr_1CP[0]+b, np.arange(nx_1CP[1])-jc_arr_1CP[0], MF_1CP.T,
    #                     norm=colors.LogNorm(vmin=vmin, vmax=vmax), cmap=cm_bw)  # cmap='RdBu_r')
    rect1 = mpatches.Rectangle((-ic_arr_1CP[0], -delta), nx_1CP[0], 2*delta, fill=True,
                               linewidth=0, edgecolor='r', facecolor='lightyellow', alpha=0.3)
    ax.add_patch(rect1)
    ax.set_xlim(xmin-ic_3CP, ic_3CP-xmin)
    ax.set_ylim(ymin-ic_3CP, ic_3CP-ymin)
    ax = axis[3]
    pcm = ax.contourf(MF_2CP, levels=lvls, cmap=cm_bw)
    # pcm = ax.pcolormesh(np.arange(nx_2CP[1])-jc_2CP, np.arange(nx_2CP[0])-ic_2CP, MF_2CP,
    #                     # norm = colors.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-max,vmax=max),
    #                     # norm = colors.PowerNorm(),
    #                     norm=colors.LogNorm(vmin=vmin, vmax=vmax), cmap=cm_bw)  # cmap='RdBu_r')
    rect = mpatches.Rectangle((-jc_2CP, -delta), nx_2CP[1], 2*delta, fill=True,
                              linewidth=0, edgecolor='r', facecolor='lightyellow', alpha=0.3)
    ax.add_patch(rect)
    ax.set_xlim(xmin-ic_3CP, ic_3CP-xmin)
    ax.set_ylim(ymin-ic_3CP, ic_3CP-ymin)
    ax = axis[4]
    pcm = ax.contourf(MF_3CP.T, levels=lvls, cmap=cm_bw)  # cmap='RdBu_r')
    # pcm = ax.pcolormesh(np.arange(nx_3CP[0])-ic_3CP, np.arange(nx_3CP[1])-jc_3CP, MF_3CP.T,
    #                     norm=colors.LogNorm(vmin=vmin, vmax=vmax), cmap=cm_bw)  # cmap='RdBu_r')
    plt.colorbar(pcm, ax=ax, extend='both')
    # rect2 = mpatches.Rectangle((-100, jc_arr_3CP[2] - delta), 300, 2*delta, fill=True,
    rect2 = mpatches.Rectangle((-100, -delta), 300, 2*delta, fill=True,
                               linewidth=0, edgecolor='r', facecolor='lightyellow', alpha=0.3)
    ax.add_patch(rect2)
    ax.set_xlim(xmin-ic_3CP, xmax-ic_3CP)
    ax.set_ylim(ymin-ic_3CP, ic_3CP-ymin)
    for ax in axis[2:].flat:
        ax.set_aspect('equal')
    axis[2].set_ylabel('y [km]')

    for ax in axis:
        ax.set_xlabel('x [km]')
        ax.set_xticklabels([np.int(ti*1e-3*dx[0]) for ti in ax.get_xticks()])
    axis[0].set_yticklabels([np.round(ti*1e-3,1) for ti in axis[0].get_yticks()])
    axis[1].set_yticklabels([np.round(ti,1) for ti in axis[1].get_yticks()])
    for ax in axis[2:]:
        ax.set_yticklabels([np.int(ti*1e-3*dx[0]) for ti in ax.get_yticks()])

    textprops = dict(facecolor='white', alpha=0.9, linewidth=0.)
    axis[0].text(xmin-ic_3CP+15, 900, 'a)', fontsize=18, bbox=textprops)
    axis[1].text(xmin-ic_3CP+15, np.amax(MF_mean_3CP)-.3, 'b)', fontsize=18, bbox=textprops)
    axis[2].text(xmin-ic_3CP+15, ic_3CP-ymin-30, 'c)', fontsize=21, bbox=textprops)
    axis[3].text(xmin-ic_3CP+15, ic_3CP-ymin-30, 'd)', fontsize=21, bbox=textprops)
    axis[4].text(xmin-ic_3CP+15, ic_3CP-ymin-30, 'e)', fontsize=21, bbox=textprops)

    plt.subplots_adjust(bottom=0.1, right=.99, left=0.03, top=0.95, hspace=0.2, wspace=0.25)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    print('')
    print('')

    return



def plot_collision_massflux_CPheight_3panels(CP_height_1CP, CP_height_2CP, CP_height_3CP,
                                     time_CPheight_1CP, time_CPheight_2CP, time_CPheight_3CP,
                                     MF_1CP, MF_2CP, MF_3CP, MF_mean_1CP, MF_mean_2CP, MF_mean_3CP,
                                     t_ini, t_2CP, t_3CP, t_final,
                                     delta, ic_2CP, jc_2CP, ic_3CP, jc_3CP,
                                     ic_arr_2CP, jc_arr_2CP, ic_arr_3CP, jc_arr_3CP,
                                     xmin, xmax, ymin,
                                     case, case_xCP,
                                     times, path_out_figs, fig_name):


    ncol = 3
    vmin = 1e-2
    vmax = 6.5
    fig, axis = plt.subplots(1, ncol, figsize=(ncol * 5, 5))

    ''' CP Height '''
    print('time windows: ', t_ini, t_2CP, t_3CP, t_final)
    for it, t0 in enumerate(times[::2]):
        if t0 >= t_ini and t0 <= t_final:
            # it_3CP = np.where(time_CPheight_3CP == t0)[0][0]
            it_3CP = np.int(t0/dt_CP_height_3CP)
            print('-- it, t0: ', it, t0, '--', t_ini, t_2CP, t_3CP, t_final)
            count_color = 0
            if t0 < t_2CP:
                cm = cm_gray
                count_color = (np.double(t0) - t_ini) / (t_2CP - t_ini) * 0.8
            elif t0 < t_3CP:
                cm = cm_bwr_r
                count_color = (np.double(t0) - t_2CP) / (t_3CP - t_2CP) * 0.5 + 0.5
            elif t0 <= t_final:
                cm = cm_bwr
                count_color = (np.double(t0) - t_3CP) / (t_final - t_3CP) * 0.45 + 0.5
        ax = axis[0]
        if t0 == t_ini or t0 == t_ini+100:
            print('single')
            ax.plot(np.arange(nx_3CP[0]) - ic_3CP, CP_height_3CP[it_3CP, :, jc_3CP], color=cm(count_color), label='single CP')
        elif t0 == t_3CP-100 or t0 == t_3CP-200:
            print('double')
            ax.plot(np.arange(nx_3CP[0]) - ic_3CP, CP_height_3CP[it_3CP, :, jc_3CP], color=cm(count_color), label='double collison')
        elif t0 == t_final-100 or t0 == t_final-200:
            print('triple')
            ax.plot(np.arange(nx_3CP[0]) - ic_3CP, CP_height_3CP[it_3CP, :, jc_3CP], color=cm(count_color), label='triple collison')
        else:
            ax.plot(np.arange(nx_3CP[0]) - ic_3CP, CP_height_3CP[it_3CP, :, jc_3CP], color=cm(count_color))
        ax.set_xlim(xmin-ic_3CP, xmax-ic_3CP)
        ax.legend()
    axis[0].set_ylabel(r'CP Height  [m]')

    ''' Mass Flux '''
    ax = axis[1]
    jc_1CP = np.int(nx_1CP[1] / 2)
    ax.plot(np.arange(nx_1CP[1]) - jc_1CP, MF_mean_1CP, label='single CP', color=cm_gray(.1))  # color=colorlist2[0])
    ax.plot(np.arange(nx_2CP[1]) - jc_2CP, MF_mean_2CP, label='double collision', color=cm_bwr_r(.8))  # color=colorlist2[0])
    ax.plot(np.arange(nx_3CP[0]) - ic_arr_3CP[0], MF_mean_3CP, label='triple collision', color=cm_bwr(.8))  # color=colorlist2[1])
    ax.legend(loc=4)
    ax.set_xlim(xmin - ic_3CP, xmax - ic_3CP)
    ax.set_ylim(-2, np.ceil(np.amax(MF_mean_3CP)))
    axis[1].set_ylabel(r'Integrated Mass Flux  [kg/m$^2$]')


    ''' Mass Flux 2D '''
    ax = axis[2]
    pcm = ax.pcolormesh(np.arange(nx_3CP[0]) - ic_3CP, np.arange(nx_3CP[1]) - jc_3CP, MF_3CP.T,
                        norm=colors.LogNorm(vmin=vmin, vmax=vmax), cmap=cm_bw)  # cmap='RdBu_r')
    plt.colorbar(pcm, ax=ax, extend='both')
    # rect2 = mpatches.Rectangle((-100, jc_arr_3CP[2] - delta), 300, 2*delta, fill=True,
    rect2 = mpatches.Rectangle((-100, -delta), 300, 2 * delta, fill=True,
                               linewidth=0, edgecolor='r', facecolor='lightyellow', alpha=0.3)
    ax.add_patch(rect2)
    ax.set_xlim(xmin - ic_3CP, xmax - ic_3CP)
    ax.set_ylim(ymin - ic_3CP, ic_3CP - ymin)
    for ax in axis[2:].flat:
        ax.set_aspect('equal')
    axis[2].set_ylabel('y [km]')

    for ax in axis:
        ax.set_xlabel('x [km]')
        ax.set_xticklabels([np.int(ti*1e-3*dx[0]) for ti in ax.get_xticks()])
    axis[0].set_yticklabels([np.round(ti*1e-3,1) for ti in axis[0].get_yticks()])
    axis[1].set_yticklabels([np.round(ti,1) for ti in axis[1].get_yticks()])
    for ax in axis[2:]:
        ax.set_yticklabels([np.int(ti*1e-3*dx[0]) for ti in ax.get_yticks()])

    textprops = dict(facecolor='white', alpha=0.9, linewidth=0.)
    axis[0].text(xmin-ic_3CP+15, 900, 'a)', fontsize=18, bbox=textprops)
    axis[1].text(xmin-ic_3CP+15, np.amax(MF_mean_3CP)-.3, 'b)', fontsize=18, bbox=textprops)
    axis[2].text(xmin-ic_3CP+15, ic_3CP-ymin-30, 'c)', fontsize=21, bbox=textprops)

    plt.subplots_adjust(bottom=0.1, right=.99, left=0.03, top=0.95, hspace=0.2, wspace=0.25)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return
# ----------------------------------------------------------------------

def set_input_output_parameters(args, case_name_1CP, case_name_2CP, case_name_3CP,
                                path_1CP, path_2CP, path_3CP,
                                case, case_xCP):
    print('')
    print('--- set input parameters ---')
    # print(os.path.join(path_1CP, case_name_1CP + '.in'))
    # print(os.path.join(path_2CP, case_name_2CP + '.in'))
    # print(os.path.join(path_3CP, case_name_3CP + '.in'))
    # print('')
    nml_1CP = simplejson.loads(open(os.path.join(path_1CP, case, case_name_1CP + '.in')).read())
    nml_2CP = simplejson.loads(open(os.path.join(path_2CP, case_xCP, case_name_2CP + '.in')).read())
    nml_3CP = simplejson.loads(open(os.path.join(path_3CP, case_xCP, case_name_3CP + '.in')).read())
    nx_1CP = np.zeros(3, dtype=np.int)
    nx_2CP = np.zeros(3, dtype=np.int)
    nx_3CP = np.zeros(3, dtype=np.int)
    nx_1CP[0] = nml_1CP['grid']['nx']
    nx_1CP[1] = nml_1CP['grid']['ny']
    nx_1CP[2] = nml_1CP['grid']['nz']
    nx_2CP[0] = nml_2CP['grid']['nx']
    nx_2CP[1] = nml_2CP['grid']['ny']
    nx_2CP[2] = nml_2CP['grid']['nz']
    nx_3CP[0] = nml_3CP['grid']['nx']
    nx_3CP[1] = nml_3CP['grid']['ny']
    nx_3CP[2] = nml_3CP['grid']['nz']
    global dx
    dx = np.zeros(3, dtype=np.int)
    dx[0] = nml_3CP['grid']['dx']
    dx[1] = nml_3CP['grid']['dy']
    dx[2] = nml_3CP['grid']['dz']
    # gw = nml_3CP['grid']['gw']
    # dV = dx[0] * dx[1] * dx[2]
    print('grid 1CP: ', nx_1CP)
    print('grid 2CP: ', nx_2CP)
    print('grid 3CP: ', nx_3CP)


    # global dt_fields_1CP, dt_fields_2CP, dt_fields_3CP, dt_stats, dt
    dt_fields_1CP = np.int(nml_1CP['fields_io']['frequency'])
    dt_fields_2CP = np.int(nml_2CP['fields_io']['frequency'])
    dt_fields_3CP = np.int(nml_3CP['fields_io']['frequency'])
    print('Output timestep:')
    print('1CP: ', dt_fields_1CP)
    print('2CP: ', dt_fields_2CP)
    print('3CP: ', dt_fields_3CP)

    return nx_1CP, nx_2CP, nx_3CP, dt_fields_1CP, dt_fields_2CP, dt_fields_3CP


# _______________________________________________________

def define_geometry(nx_1CP, nx_2CP, nx_3CP, sep):
    # 1-CP simulation
    ic_arr_1CP = [np.int(nx_1CP[0]/2)]
    jc_arr_1CP = [np.int(nx_1CP[1]/2)]

    # 2-CP simulation
    isep = np.int(np.round(sep / dx[0]))
    jsep = 0
    ic_2CP = np.int(np.round(nx_2CP[0] / 2))
    jc_2CP = np.int(np.round(nx_2CP[1] / 2))
    ic1 = ic_2CP - np.int(np.round(isep / 2))
    jc1 = jc_2CP
    ic2 = ic1 + isep
    jc2 = jc1 + jsep
    ic_arr_2CP = [ic1, ic2]
    jc_arr_2CP = [jc1, jc2]

    # 3-CP simulation
    i_d = np.int(np.round(sep / dx[0]))
    idhalf = np.int(np.round(i_d / 2))
    a = np.int(np.round(i_d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
    r_int = np.int(np.sqrt(3.) / 6 * i_d)  # radius of inscribed circle
    # point of 3-CP collision (ic, jc)
    ic_3CP = np.int(np.round(nx_3CP[0] / 2))
    jc_3CP = np.int(np.round(nx_3CP[1] / 2))
    ic1 = ic_3CP - r_int
    ic2 = ic1
    ic3 = ic_3CP + (a - r_int)
    jc1 = jc_3CP - idhalf
    jc2 = jc_3CP + idhalf
    jc3 = jc_3CP
    ic_arr_3CP = [ic1, ic2, ic3]
    jc_arr_3CP = [jc1, jc2, jc3]
    print('')
    print('CP configuration: ')
    print('1CP: ', ic_arr_1CP, jc_arr_1CP)
    print('2CP: ', ic_arr_2CP, jc_arr_2CP)
    print('3CP: ', ic_arr_3CP, jc_arr_3CP)
    print('CP collision points: ')
    print('2CP: ', ic_2CP, jc_2CP)
    print('3CP: ', ic_3CP, jc_3CP)
    print('')

    return ic_arr_1CP, jc_arr_1CP, ic_arr_2CP, jc_arr_2CP, ic_arr_3CP, jc_arr_3CP, ic_2CP, jc_2CP, ic_3CP, jc_3CP

# _______________________________________________________

if __name__ == '__main__':
    main()