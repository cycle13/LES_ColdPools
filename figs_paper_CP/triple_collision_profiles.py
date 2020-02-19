import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35
plt.rcParams['text.usetex'] = 'true'

execfile('settings.py')


def main():
    path_root = '/nbi/ac/cond2/meyerbe/ColdPools/3D_sfc_fluxes_off/triple_3D/'
    # path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single'
    path_out_figs = '/nbi/home/meyerbe/paper_CP'
    case_name = 'ColdPoolDry_triple_3D'

    dTh = 5
    id_list = set_input_parameters(case_name, dTh)
    ncases = len(id_list)
    dt_fields = define_geometry(case_name, path_root, id_list)

    ''' determine time range '''
    global tmin, tmax
    tmin = np.int(0)
    tmax = np.int(7200)
    times = np.arange(tmin, tmax + dt_fields, dt_fields)
    print('tmin: ' + str(tmin))
    print('tmax: ' + str(tmax))
    print('times: ' + str(times))
    print('')


    ''' plot from crosssections stats-file '''
    # var_name = 'w'
    # filename_in = 'stats_crosssections_kmax20.nc'
    # test_file = nc.Dataset(os.path.join(path_root, id_list[0], 'data_analysis', filename_in), 'r', format='NETCDF4')
    # nz = test_file.groups[var_name].dimensions['nz'].size
    # test_file.close()
    # times_array_s = np.ndarray((ncases, 2), dtype=np.double)
    # times_array_w = np.ndarray((ncases, 2), dtype=np.double)
    # min_array_s = np.ndarray((ncases, 3, nz), dtype=np.double)
    # max_array_w = np.ndarray((ncases, 3, nz), dtype=np.double)
    #
    # for i_id, ID in enumerate(id_list[:]):
    #     print('--- '+ID)
    #     fullpath_in = os.path.join(path_root, ID, 'data_analysis', filename_in)
    #     print(fullpath_in)
    #     print('')
    #     file = nc.Dataset(fullpath_in, 'r', format='NETCDF4')
    #     min_array_s[i_id, 0, :] = file.groups['s'].variables['min_single'][:nz]
    #     min_array_s[i_id, 1, :] = file.groups['s'].variables['min_double'][:nz]
    #     min_array_s[i_id, 2, :] = file.groups['s'].variables['min_triple'][:nz]
    #     max_array_w[i_id, 0, :] = file.groups['w'].variables['max_single'][:nz]
    #     max_array_w[i_id, 1, :] = file.groups['w'].variables['max_double'][:nz]
    #     max_array_w[i_id, 2, :] = file.groups['w'].variables['max_triple'][:nz]
    #     times_array_s[i_id, 0] = file.groups['s'].variables['time_single'][:]
    #     times_array_s[i_id, 1] = file.groups['s'].variables['time_double'][:]
    #     times_array_w[i_id, 0] = file.groups['w'].variables['time_single'][:]
    #     times_array_w[i_id, 1] = file.groups['w'].variables['time_double'][:]
    #     file.close()
    #
    # fig_name = 'collisions_minmax_profiles_windows.png'
    # zrange = np.arange(nz) * dx[2]
    # fig, axis = plt.subplots(1, 3, figsize=(14, 10))
    # ax0 = axis[0]
    # ax1 = axis[1]
    # for i_id, ID in enumerate(id_list):
    #     al = np.double(i_id + 1) / (ncases + 1)
    #     ax0.plot(max_array_w[i_id, 0, :nz], zrange, color='b', alpha=al, label='single, ' + ID)
    #     ax0.plot(max_array_w[i_id, 1, :nz], zrange, '-', color='g', alpha=al, label='double')
    #     ax0.plot(max_array_w[i_id, 2, :nz], zrange, '-', color='r', alpha=al, label='triple')
    #     ax1.plot(min_array_s[i_id, 0, :nz], zrange, color='b', alpha=al, label='single, ' + ID)
    #     ax1.plot(min_array_s[i_id, 1, :nz], zrange, '-', color='g', alpha=al, label='double')
    #     ax1.plot(min_array_s[i_id, 2, :nz], zrange, '-', color='r', alpha=al, label='triple')
    # ax0.legend(loc='upper left', bbox_to_anchor=(-0.1, -0.1),
    #            fancybox=False, shadow=False, ncol=3, fontsize=9)
    # # ax1.legend(loc='upper left', bbox_to_anchor=(0.1, -0.1),
    # #            fancybox=False, shadow=False, ncol=1, fontsize=9)
    # ax0.set_ylabel('height [m]')
    # plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.9, hspace=0.4)
    # plt.savefig(os.path.join(path_out_figs, fig_name))
    # plt.close()

    # -------------------------------------------------------
    filename_in = 'minmax_'+id_list[0]+'.nc'
    test_file = nc.Dataset(os.path.join(path_root, id_list[0], 'data_analysis', filename_in), 'r', format='NETCDF4')
    nz = test_file.groups['profiles'].dimensions['nz'].size
    nt = test_file.groups['profiles'].dimensions['nt'].size
    nt = 72
    test_file.close()
    print('nt', nt, 'nz', nz)
    times_array = np.ndarray((ncases, nt), dtype=np.double)
    ts_max_w = np.ndarray((ncases, nz), dtype=np.double)
    ts_min_th = np.ndarray((ncases, nz), dtype=np.double)
    prof_max_w = np.ndarray((ncases, nt, nz), dtype=np.double)
    prof_min_th = np.ndarray((ncases, nt, nz), dtype=np.double)

    for i_id, ID in enumerate(id_list[:]):
        print('--- ' + ID)
        filename_in = 'minmax_' + ID + '.nc'
        fullpath_in = os.path.join(path_root, ID, 'data_analysis', filename_in)
        print(fullpath_in)
        print('')
        file = nc.Dataset(fullpath_in, 'r', format='NETCDF4')
        times_array[i_id, :] = file.groups['timeseries'].variables['time'][:nt]

        ts_max_w[i_id, :] = file.groups['timeseries'].variables['w_max'][:nz]
        ts_min_th[i_id, :] = file.groups['timeseries'].variables['theta_min'][:nz]

        prof_max_w[i_id, :,:] = file.groups['profiles'].variables['w_max'][:nt,:nz]
        # prof_min_th[i_id, :,:] = file.groups['profiles'].variables['theta_min'][:nt,nz]


        file.close()



    fig_name = 'collisions_minmax_profiles_domain.png'
    zrange = np.arange(nz) * dx[2]
    fig, axis = plt.subplots(1, 3, figsize=(14, 5))
    ax0 = axis[0]
    ax1 = axis[1]
    for i_id, ID in enumerate(id_list):
        al = np.double(i_id + 1) / (ncases + 1)
        ax0.plot(ts_max_w[i_id, :nz], zrange, color='b', alpha=al, label='single, ' + ID)
        ax1.plot(ts_min_th[i_id, :nz], zrange, color='b', alpha=al, label='single, ' + ID)
    ax0.legend(loc='upper left', bbox_to_anchor=(-0.1, -0.1),
               fancybox=False, shadow=False, ncol=3, fontsize=9)
    # ax1.legend(loc='upper left', bbox_to_anchor=(0.1, -0.1),
    #            fancybox=False, shadow=False, ncol=1, fontsize=9)
    ax0.set_ylabel('height [km]')
    ax0.set_xlabel('max(w) [m/s]')
    ax1.set_xlabel(r'min(potential temperature  ) [K]')
    for ax in axis.flatten():
        y_ticks = [np.int(ti * 1e-3) for ti in ax.get_yticks()]
        ax.set_yticklabels(y_ticks)
    plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.9, hspace=0.4)
    plt.savefig(os.path.join(path_out_figs, fig_name))
    plt.close()

    time_bins = np.ndarray((ncases,3), dtype=np.double)
    time_bins[0,:] = [1000, 1400, times[-1]]
    time_bins[1,:] = [4500, 5900, times[-1]]
    # time_bins[1,:] = [2300, 3200, times[-1]]
    # time_bins[2,:] = [4500, 5900, times[-1]]
    time_bins_ = {}
    time_bins_['d10km'] = [1000, 1400, times[-1]]
    time_bins_['d15km'] = [2300, 3200, times[-1]]
    time_bins_['d20km'] = [4500, 5900, times[-1]]

    fig_name = 'collisions_minmax_profiles.png'
    zrange = np.arange(nz) * dx[2]
    fig, axis = plt.subplots(1, 3, figsize=(14, 5))
    ax0 = axis[0]
    ax1 = axis[1]
    for i_id, ID in enumerate(id_list):
        print(ID[-5:])
        # it_s = np.int(times_array_w[i_id, 0]/dt_fields)
        # it_d = np.int(times_array_w[i_id, 1]/dt_fields)
        # it_s = np.int(time_bins[i_id, 0]/dt_fields)
        # it_d = np.int(time_bins[i_id, 1]/dt_fields)
        it_s = np.int(time_bins_[ID[-5:]][0]/dt_fields)
        it_d = np.int(time_bins_[ID[-5:]][1]/dt_fields)
        # it_t = times_array_w[i_id, 2]
        print('it single, double: ', it_s, it_d)
        al = np.double(i_id + 1) / (ncases + 1)
        aux = np.amax(prof_max_w[i_id, :it_s, :nz], axis=0)
        print('', prof_max_w.shape, aux.shape, nz, zrange.shape)
        ax0.plot(np.amax(prof_max_w[i_id, :it_s, :nz], axis=0), zrange, color='b', alpha=al, label='single, ' + ID)
        # ax0.plot(prof_max_w[i_id, :nz], zrange, '-', color='g', alpha=al, label='double')
        # ax0.plot(prof_max_w[i_id, :nz], zrange, '-', color='r', alpha=al, label='triple')
    #     ax1.plot(prof_min_th[i_id, :nz], zrange, color='b', alpha=al, label='single, ' + ID)
    #     # ax1.plot(prof_min_th[i_id, :nz], zrange, '-', color='g', alpha=al, label='double')
    #     # ax1.plot(prof_min_th[i_id, :nz], zrange, '-', color='r', alpha=al, label='triple')
    # ax0.legend(loc='upper left', bbox_to_anchor=(-0.1, -0.1),
    #            fancybox=False, shadow=False, ncol=3, fontsize=9)
    # # ax1.legend(loc='upper left', bbox_to_anchor=(0.1, -0.1),
    # #            fancybox=False, shadow=False, ncol=1, fontsize=9)
    ax0.set_ylabel('height [km]')
    ax0.set_xlabel('max(w) [m/s]')
    ax1.set_xlabel(r'min(potential temperature  ) [K]')
    # ax1.set_xlim(250,350)
    for ax in axis.flatten():
        y_ticks = [np.int(ti * 1e-3) for ti in ax.get_yticks()]
        ax.set_yticklabels(y_ticks)
    plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.9, hspace=0.4)
    plt.savefig(os.path.join(path_out_figs, fig_name))
    plt.close()

    return


# --------------------------------------------------------------------


# --------------------------------------------------------------------



def set_input_parameters(case_name, dTh):

    if case_name == 'ColdPoolDry_triple_3D':
        if dTh == 3:
            zstar = 2000
            rstar = 2000
        elif dTh == 5:
            zstar = 1000
            rstar = 1100
        else:
            print('dTh not defined !!!!')
        # d_range = [10, 15, 20]
        d_range = [10, 20]
        id_list = []
        for dstar in d_range:
            id_list.append('dTh' + str(dTh) +'_z' + str(zstar) + '_r' + str(rstar) + '_d' + str(dstar) + 'km')
    print('ID-list: ', id_list)
    print('')



    return id_list



def define_geometry(case_name, path_root, id_list):
    print 'define geometry'
    global nx, ny, nz, dx, gw
    nml = simplejson.loads(open(os.path.join(path_root, id_list[0], case_name + '.in')).read())
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = np.ndarray(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']
    global dt_fields
    dt_fields = np.int(nml['fields_io']['frequency'])

    return dt_fields


# ----------------------------------




if __name__ == '__main__':
    print('')
    main()

