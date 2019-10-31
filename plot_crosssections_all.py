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

from thermodynamic_functions import thetas_c

def main():
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path_root")
    parser.add_argument("dTh")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    kmax, krange, id_list = set_input_parameters(args)
    dt_fields, times = define_geometry(id_list)
    nt = len(times)
    ncases = len(id_list)
    print('')


    ''' (A) from stats_crosssections.nc '''
    filename_in = os.path.join('data_analysis', 'stats_crosssections.nc')
    var_list = ['w', 's']
    # var_list = ['w']
    for var_name in var_list:
        test_file = nc.Dataset(os.path.join(path_root, id_list[0], filename_in), 'r', format='NETCDF4')
        nz = test_file.groups[var_name].dimensions['nz'].size
        test_file.close()
        times_array = np.ndarray((ncases, 2), dtype=np.double)
        tw = np.ndarray((ncases, 2), dtype=np.double)
        min_array = np.ndarray((ncases, 3, nz), dtype=np.double)
        max_array = np.ndarray((ncases, 3, nz), dtype=np.double)
        min_array_tw = np.ndarray((ncases, 3, nz), dtype=np.double)
        max_array_tw = np.ndarray((ncases, 3, nz), dtype=np.double)
        for i_id, ID in enumerate(id_list[:]):
            print('--- ' + ID)
            fullpath_in = os.path.join(path_root, ID, filename_in)
            print(fullpath_in)
            print('')
            file = nc.Dataset(fullpath_in, 'r', format='NETCDF4')
            tw[i_id, 0] = (file.groups['w'].variables['time_single'][:] / dt_fields)
            tw[i_id, 1] = file.groups['w'].variables['time_double'][:]
            i_tw_s = np.int(file.groups['w'].variables['time_single'][:] / dt_fields)
            i_tw_d = np.int(file.groups['w'].variables['time_double'][:] / dt_fields)
            times_array[i_id, 0] = file.groups[var_name].variables['time_single'][:]
            times_array[i_id, 1] = file.groups[var_name].variables['time_double'][:]
            min_array[i_id, 0, :] = file.groups[var_name].variables['min_single'][:nz]
            min_array[i_id, 1, :] = file.groups[var_name].variables['min_double'][:nz]
            min_array[i_id, 2, :] = file.groups[var_name].variables['min_triple'][:nz]
            max_array[i_id, 0, :] = file.groups[var_name].variables['max_single'][:nz]
            max_array[i_id, 1, :] = file.groups[var_name].variables['max_double'][:nz]
            max_array[i_id, 2, :] = file.groups[var_name].variables['max_triple'][:nz]
            print(nz, tw[i_id, 0], i_tw_s)
            min_array_tw[i_id, 0, :] = np.amin(file.groups[var_name].variables['min_s'][:i_tw_s, :nz], axis=0)
            min_array_tw[i_id, 1, :] = np.amin(file.groups[var_name].variables['min_d'][i_tw_s:i_tw_d, :nz], axis=0)
            min_array_tw[i_id, 2, :] = np.amin(file.groups[var_name].variables['min_t'][i_tw_d:, :nz], axis=0)
            max_array_tw[i_id, 0, :] = np.amax(file.groups[var_name].variables['max_s'][:i_tw_s, :nz], axis=0)
            max_array_tw[i_id, 1, :] = np.amax(file.groups[var_name].variables['max_d'][i_tw_s:i_tw_d, :nz], axis=0)
            max_array_tw[i_id, 2, :] = np.amax(file.groups[var_name].variables['max_t'][i_tw_d:, :nz], axis=0)
            file.close()
        zrange = np.arange(nz) * dx[2]
        print('')
        if var_name == 's':
            max_array = thetas_c(max_array, 0.)
            print(max_array_tw[0, 0, :])
            max_array_tw = thetas_c(max_array, 0.)
            print(max_array_tw[0, 0, :])

        print('var: ', var_name)
        if nz < kmax:
            print('!!! nz < kmax !!!')
            print('')
        # print(nz)
        # print(np.amax(np.abs(min_array)))
        print(np.amax(np.abs(max_array)))
        # print(np.amax(np.abs(times_array)))
        print('')
        print('')

        fig_name_root = 'crosssections_profiles_test_' + var_name
        if var_name == 'w':
            plot_profiles_crosssections_test(var_name, max_array, 'max', id_list, zrange, kmax, filename_in, fig_name_root)
        elif var_name == 's':
            plot_profiles_crosssections_test(var_name, min_array_tw, 'min', id_list, zrange, kmax, filename_in, fig_name_root)

        fig_name_root = 'crosssections_profiles_all_' + var_name
        plot_profiles_crosssections(var_name, max_array, max_array_tw, min_array, min_array_tw,
                      times_array, id_list, ncases, zrange, kmax, fig_name_root)




    # # plot times
    # times_array = np.ndarray((ncases, len(var_list), 2), dtype=np.double)
    # for i_id, ID in enumerate(id_list[:]):
    #     print('--- ' + ID)
    #     fullpath_in = os.path.join(path_root, ID, filename_in)
    #     print(fullpath_in)
    #     print('')
    #     file = nc.Dataset(fullpath_in, 'r', format='NETCDF4')
    #     for i,var_name in enumerate(var_list):
    #         times_array[i_id, i, 0] = file.groups[var_name].variables['time_single'][:]
    #         times_array[i_id, i, 1] = file.groups[var_name].variables['time_double'][:]
    #     file.close()
    #
    # fig_name = 'times_' + case_name + '.png'
    # fig, axis = plt.subplots(1, 2, figsize=(9, 10), sharey='row')
    # for i_id, ID in enumerate(id_list):
    #     al = np.double(i_id + 1) / (ncases + 1)
    #     for i, var_name in enumerate(var_list):
    #         ax = axis[i]
    #         ax.plot([1,2],times_array[i_id, i, :], '-o', color='b', alpha=al, label='single, ' + ID)
    #         ax.set_title('from ' + var_list[0])
    # for ax in axis:
    #     ax.set_xlabel('collision (#CPs)')
    #     ax.set_ylabel('time [s]')
    #     ax.grid()
    # axis[0].legend(loc='upper left', bbox_to_anchor=(-0.1, -0.1),
    #            fancybox=False, shadow=False, ncol=3, fontsize=9)
    # plt.suptitle('times for 2-CP & 2-CP collision', fontsize=21)
    # plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.9, hspace=0.4)
    # plt.savefig(os.path.join(path_out_figs, fig_name))
    # plt.close()

    ''' (B) from minmax_ID.nc '''
    var_list = ['w', 's']
    for var_name in var_list:
        filename_in = os.path.join('data_analysis', 'minmax_'+id_list[0] + '.nc')
        test_file = nc.Dataset(os.path.join(path_root, id_list[0], filename_in), 'r', format='NETCDF4')
        nz = test_file.groups['profiles'].dimensions['nz'].size
        nt = test_file.groups['profiles'].dimensions['nt'].size
        nt = 72
        test_file.close()
        # times_array = np.ndarray((ncases, 2), dtype=np.double)
        # tw = np.ndarray((ncases, 2), dtype=np.double)
        min_array = np.ndarray((ncases, nt, nz), dtype=np.double)
        max_array = np.ndarray((ncases, nt, nz), dtype=np.double)
        # min_array_tw = np.ndarray((ncases, 3, nz), dtype=np.double)
        # max_array_tw = np.ndarray((ncases, 3, nz), dtype=np.double)
        for i_id, ID in enumerate(id_list[:]):
            print('--- ' + ID)
            filename_in = 'minmax_' + ID + '.nc'
            fullpath_in = os.path.join(path_root, ID, 'data_analysis', filename_in)
            print(fullpath_in)
            print('')
            file = nc.Dataset(fullpath_in, 'r', format='NETCDF4')
            prof_grp = file.groups['profiles']
            # tw[i_id, 0] = (file.groups['w'].variables['time_single'][:] / dt_fields)
            # tw[i_id, 1] = file.groups['w'].variables['time_double'][:]
            # i_tw_s = np.int(file.groups['w'].variables['time_single'][:] / dt_fields)
            # i_tw_d = np.int(file.groups['w'].variables['time_double'][:] / dt_fields)
            # times_array[i_id, 0] = file.groups[var_name].variables['time_single'][:]
            # times_array[i_id, 1] = file.groups[var_name].variables['time_double'][:]
            test =  prof_grp.variables[var_name + '_min'][:,:]
            print('test', test.shape)
            min_array[i_id, :nt, :] = prof_grp.variables[var_name + '_min'][:nt,:]
            max_array[i_id, :nt, :] = prof_grp.variables[var_name + '_max'][:nt,:]
            # print(nz, tw[i_id, 0], i_tw_s)
            # min_array_tw[i_id, 0, :] = np.amin(file.groups[var_name].variables['min_s'][:i_tw_s, :nz], axis=0)
            # min_array_tw[i_id, 1, :] = np.amin(file.groups[var_name].variables['min_d'][i_tw_s:i_tw_d, :nz], axis=0)
            # min_array_tw[i_id, 2, :] = np.amin(file.groups[var_name].variables['min_t'][i_tw_d:, :nz], axis=0)
            # max_array_tw[i_id, 0, :] = np.amax(file.groups[var_name].variables['max_s'][:i_tw_s, :nz], axis=0)
            # max_array_tw[i_id, 1, :] = np.amax(file.groups[var_name].variables['max_d'][i_tw_s:i_tw_d, :nz], axis=0)
            # max_array_tw[i_id, 2, :] = np.amax(file.groups[var_name].variables['max_t'][i_tw_d:, :nz], axis=0)
            file.close()
        zrange = np.arange(nz) * dx[2]
        # print('')
        # if var_name == 's':
        #     max_array = thetas_c(max_array, 0.)
        #     print(max_array_tw[0, 0, :])
        #     max_array_tw = thetas_c(max_array, 0.)
        #     print(max_array_tw[0, 0, :])


        fig_name_root = 'minmax_profiles_test_' + var_name
        plot_profiles_minmax_test(var_name, max_array, 'max', id_list, zrange, kmax, fig_name_root)

    return


#
# def plot_profiles_minmax_crosssections_test():
#     for i_id, ID in enumerate(id_list[:]):
#         fig, axis = plt.subplots(3, 4, figsize=(18, 15), sharex='all')
#         plt.savefig(os.path.join(path_out_figs, fig_name_ + '_' + ID + '.png'))
#         plt.close(fig)
#     return


def plot_profiles_minmax_test(var_name, minmax_array, type_min_max, id_list, zrange, kmax, fig_name_):
    delta_s = [7, 21, 36]  # 10: 7, 15: 21
    delta_d = [1, 4, 4]  # 10: 1, 15: 4
    it_s = [10, 25, 45]
    it_d = [15, 33, 57]
    it_t = [i+10 for i in it_d]
    for i_id, ID in enumerate(id_list[:]):
        fig, axis = plt.subplots(3, 4, figsize=(18, 15), sharex='all')
        ax00 = axis[0, 0]
        ax01 = axis[0, 1]
        ax02 = axis[0, 2]
        ax03 = axis[0, 3]
        ax10 = axis[1, 0]
        ax11 = axis[1, 1]
        ax12 = axis[1, 2]
        ax13 = axis[1, 3]
        ax20 = axis[2, 0]
        ax21 = axis[2, 1]
        ax22 = axis[2, 2]
        ax23 = axis[2, 3]
        print('')
        ax00.set_title('maxima over time windows')
        ax00.plot(np.amax(minmax_array[i_id, :it_s[i_id], :kmax], axis=0), zrange[:kmax], color='b', label='single, ' + ID)
        ax00.plot(np.amax(minmax_array[i_id, it_s[i_id]:it_d[i_id], :kmax], axis=0), zrange[:kmax], '-', color='g', label='double')
        ax00.plot(np.amax(minmax_array[i_id, it_d[i_id]:it_t[i_id], :kmax], axis=0), zrange[:kmax], '-', color='r', label='triple')
        for i, delta in enumerate(np.arange(0, 50, 1)):
            it_s_ = it_s[i_id] - delta
            if it_s_ > 0:
                al = 1. - np.double(delta) / (1 + it_s[i_id])
                max_s = np.amax(minmax_array[i_id, :it_s_, :], axis=0)
                ax01.plot(max_s[:kmax], zrange[:kmax], color='b', alpha=al, label='single, ' + ID)
                delta_max = delta
        ax01.plot(np.amax(minmax_array[i_id, :it_s[i_id], :kmax], axis=0), zrange[:kmax], '--k')
        ax01.set_title('it_s=' + str(it_s[i_id]) + ',...,' + str(it_s[i_id] - delta_max))
        for i, delta in enumerate(np.arange(0, 10, 1)):
            it_s_ = it_s[i_id] - delta
            it_d_ = it_d[i_id] - delta
            max_d = np.amax(minmax_array[i_id, it_s_:it_d_, :], axis=0)
            max_t = np.amax(minmax_array[i_id, it_d_:, :], axis=0)
            al = 1. - np.double(delta) / 15
            ax02.plot(max_d[:kmax], zrange[:kmax], color='g', alpha=al, label='double, ' + ID)
            ax03.plot(max_t[:kmax], zrange[:kmax], color='r', alpha=al, label='triple, ' + ID)
        ax02.plot(np.amax(minmax_array[i_id, it_s[i_id]:it_d[i_id], :kmax], axis=0), zrange[:kmax], '--k', label='double')
        ax03.plot(np.amax(minmax_array[i_id, it_d[i_id]:, :kmax], axis=0), zrange[:kmax], '--k', label='triple')
        for ax in axis[0, 2:].flatten():
            ax.set_title('delta=1,..,' + str(delta))

        ax10.set_title('instantaneous maxima')
        ax10.plot(minmax_array[i_id, it_s[i_id], :kmax], zrange[:kmax], color='b', label='single, ' + ID)
        ax10.plot(minmax_array[i_id, it_d[i_id], :kmax], zrange[:kmax], '-', color='g', label='double')
        ax10.plot(minmax_array[i_id, it_t[i_id], :kmax], zrange[:kmax], '-', color='r', label='triple')
        for i, delta in enumerate(np.arange(0, 50, 1)):
            it_s_ = it_s[i_id] - delta
            if it_s_ > 0:
                al = 1. - np.double(delta) / (1 + it_s[i_id])
                ax11.plot(minmax_array[i_id, it_s_, :kmax], zrange[:kmax], 'b', alpha=al, label='t='+str(it_s_))
                delta_max = delta
        ax11.plot(minmax_array[i_id, it_s[i_id], :kmax], zrange[:kmax], 'k--')
        for i, delta in enumerate(np.arange(0, 10, 1)):
            it_d_ = it_d[i_id] - delta
            it_t_ = it_t[i_id] - delta
            al = 1. - np.double(delta) / 15
            ax12.plot(minmax_array[i_id, it_d_, :kmax], zrange[:kmax], 'g', alpha=al, label='t='+str(it_d_))
            ax13.plot(minmax_array[i_id, it_t_, :kmax], zrange[:kmax], 'r', alpha=al, label='t='+str(it_t_))
        ax12.plot(minmax_array[i_id, it_d[i_id], :kmax], zrange[:kmax], '--k')
        ax13.plot(minmax_array[i_id, it_t[i_id], :kmax], zrange[:kmax], '--k')
        ax11.legend()
        ax12.legend()
        ax13.legend()

        it_s_ = it_s[i_id] - delta_s[i_id]
        it_d_ = it_d[i_id] - delta_d[i_id]
        ax20.set_title('maxima time window: it_ws=' + str(it_s) + ', it_wd=' + str(it_d))
        ax20.plot(np.amax(minmax_array[i_id, :it_s_, :kmax], axis=0), zrange[:kmax], color='b', label='single, ' + ID)
        ax20.plot(np.amax(minmax_array[i_id, it_s_:it_d_, :kmax], axis=0), zrange[:kmax], '-', color='g', label='double')
        ax20.plot(np.amax(minmax_array[i_id, it_d_:, :kmax], axis=0), zrange[:kmax], '-', color='r', label='triple')
        ax21.plot(np.amax(minmax_array[i_id, :it_s_, :kmax], axis=0), zrange[:kmax], 'b', label='single, ' + ID)
        ax21.plot(np.amax(minmax_array[i_id, :it_s_ - 1, :kmax], axis=0), zrange[:kmax], ':b', label='single, ' + ID)
        ax21.plot(np.amax(minmax_array[i_id, :it_s_ + 1, :kmax], axis=0), zrange[:kmax], ':b', label='single, ' + ID)
        ax22.plot(np.amax(minmax_array[i_id, it_s_:it_d_, :kmax], axis=0), zrange[:kmax], 'g', label='double, ' + ID)
        ax22.plot(np.amax(minmax_array[i_id, it_s_:it_d_ - 1, :kmax], axis=0), zrange[:kmax], ':g',
                  label='double, ' + ID)
        ax22.plot(np.amax(minmax_array[i_id, it_s_:it_d_ + 1, :kmax], axis=0), zrange[:kmax], ':g',
                  label='double, ' + ID)
        ax23.plot(np.amax(minmax_array[i_id, it_d_:, :kmax], axis=0), zrange[:kmax], color='r', label='triple, ' + ID)
        ax23.plot(np.amax(minmax_array[i_id, it_d_ - 1:, :kmax], axis=0), zrange[:kmax], ':r', label='triple, ' + ID)
        ax23.plot(np.amax(minmax_array[i_id, it_d_ + 1:, :kmax], axis=0), zrange[:kmax], ':r', label='triple, ' + ID)
        ax21.plot(np.amax(minmax_array[i_id, :it_s_, :kmax], axis=0), zrange[:kmax], '--k')
        ax22.plot(np.amax(minmax_array[i_id, it_s_:it_d_, :kmax], axis=0), zrange[:kmax], '--k')
        ax23.plot(np.amax(minmax_array[i_id, it_d_:, :kmax], axis=0), zrange[:kmax], '--k')

        ax00.legend(loc='best')
        fig.suptitle(var_name + ', ' + ID)
        plt.savefig(os.path.join(path_out_figs, fig_name_ + '_' + ID + '.png'))
        plt.close(fig)
    return



def plot_profiles_crosssections_test(var_name, ref_array, type_min_max, id_list, zrange, kmax,
                                     filename_in, fig_name_):
    delta_s = [7, 21, 36]  # 10: 7, 15: 21
    delta_d = [1, 4, 4]  # 10: 1, 15: 4
    for i_id, ID in enumerate(id_list[:]):
        fig, axis = plt.subplots(3, 4, figsize=(18, 15), sharex='col', sharey='all')
        ax00 = axis[0, 0]
        ax01 = axis[0, 1]
        ax02 = axis[0, 2]
        ax03 = axis[0, 3]
        ax10 = axis[1, 0]
        ax11 = axis[1, 1]
        ax12 = axis[1, 2]
        ax13 = axis[1, 3]
        ax20 = axis[2, 0]
        ax21 = axis[2, 1]
        ax22 = axis[2, 2]
        ax23 = axis[2, 3]

        fullpath_in = os.path.join(path_root, ID, filename_in)
        print(fullpath_in)
        file = nc.Dataset(fullpath_in, 'r', format='NETCDF4')
        max_arr_s = file.groups[var_name].variables[type_min_max + '_s'][:, :kmax]
        max_arr_d = file.groups[var_name].variables[type_min_max + '_d'][:, :kmax]
        max_arr_t = file.groups[var_name].variables[type_min_max + '_t'][:, :kmax]
        i_tw_s = np.int(file.groups['w'].variables['time_single'][:] / dt_fields)
        i_tw_d = np.int(file.groups['w'].variables['time_double'][:] / dt_fields)
        # i_tw_t = i_tw_t_[i_id]
        i_tw_t = i_tw_d + 7

        ax00.set_title('maxima over time windows')
        ax00.plot(ref_array[i_id, 0, :kmax], zrange[:kmax], color='b', label='single, ' + ID)
        ax00.plot(ref_array[i_id, 1, :kmax], zrange[:kmax], color='g', label='double')
        ax00.plot(ref_array[i_id, 2, :kmax], zrange[:kmax], color='r', label='triple')
        for i, delta in enumerate(np.arange(0, 50, 1)):
            it_s = i_tw_s - delta
            if it_s > 0:
                if type_min_max == 'max':
                    max_s = np.amax(max_arr_s[:it_s, :], axis=0)
                elif type_min_max == 'min':
                    max_s = np.amin(max_arr_s[:it_s, :], axis=0)
                al = 1. - np.double(delta) / (1 + i_tw_s)
                # print('-----------', it_s, al)
                ax01.plot(max_s, zrange[:kmax], color='b', alpha=al, label='t<'+str(it_s))
                delta_max = delta
        ax01.plot(ref_array[i_id, 0, :kmax], zrange[:kmax], '--k', label='single, ' + ID)
        ax01.set_title('it_s=' + str(i_tw_s) + ',...,' + str(i_tw_s - delta_max))
        ax02.set_title('it_d>=' + str(i_tw_s) + ', delta=' + str(delta))
        ax03.set_title('it_t>=' + str(i_tw_d) + ', delta=' + str(delta))
        for i, delta in enumerate(np.arange(0, 10, 1)):
            it_s = i_tw_s - delta
            it_d = i_tw_d - delta
            if type_min_max == 'max':
                max_d = np.amax(max_arr_d[it_s:it_d, :], axis=0)
                max_t = np.amax(max_arr_t[it_d:, :], axis=0)
            elif type_min_max == 'min':
                max_d = np.amin(max_arr_d[it_s:it_d, :], axis=0)
                max_t = np.amin(max_arr_t[it_d:, :], axis=0)
            al = 1. - np.double(delta) / 15
            ax02.plot(max_d, zrange[:kmax], color='g', alpha=al, label=str(it_s)+'<=t<'+str(it_d))
            ax03.plot(max_t, zrange[:kmax], color='r', alpha=al, label=str(it_d)+'<=t')
        ax02.plot(ref_array[i_id, 1, :kmax], zrange[:kmax], '--k', label='double')
        ax03.plot(ref_array[i_id, 2, :kmax], zrange[:kmax], '--k', label='triple')

        ax10.set_title('instantaneous maxima')
        ax10.plot(max_arr_s[i_tw_s, :kmax], zrange[:kmax], 'b', label='single, ' + ID)
        ax10.plot(max_arr_d[i_tw_d, :kmax], zrange[:kmax], 'g', label='double')
        ax10.plot(max_arr_t[i_tw_t, :kmax], zrange[:kmax], 'r', label='triple')
        for i, delta in enumerate(np.arange(0, 50, 1)):
            it_s_ = i_tw_s - delta
            if it_s_ > 0:
                al = 1. - np.double(delta) / (1 + i_tw_s)
                ax11.plot(max_arr_s[it_s_, :kmax], zrange[:kmax], 'b', alpha=al, label='t='+str(it_s_))
                delta_max = delta
        ax11.plot(max_arr_s[i_tw_s, :kmax], zrange[:kmax], 'k--')
        for i, delta in enumerate(np.arange(0, 10, 1)):
            it_d = i_tw_d - delta
            it_t = i_tw_t - delta
            al = 1. - np.double(delta) / 11
            ax12.plot(max_arr_d[it_d, :kmax], zrange[:kmax], 'g', alpha=al, label='t='+str(it_d))
            ax13.plot(max_arr_t[it_t, :kmax], zrange[:kmax], 'r', alpha=al, label='t='+str(it_t))
        ax11.plot(max_arr_s[i_tw_s, :kmax], zrange[:kmax], '--k')
        ax12.plot(max_arr_d[i_tw_d, :kmax], zrange[:kmax], '--k')
        ax13.plot(max_arr_t[i_tw_t, :kmax], zrange[:kmax], '--k')
        ax11.set_title('it_s=' + str(i_tw_s))
        ax12.set_title('it_d=' + str(i_tw_d))
        ax13.set_title('it_t=' + str(i_tw_t))
        ax11.legend()
        ax12.legend()
        ax13.legend()

        it_s = i_tw_s - delta_s[i_id]
        it_d = i_tw_d - delta_d[i_id]
        it_t = i_tw_t
        ax20.set_title('maxima time window')
        if type_min_max == 'max':
            ax20.plot(np.amax(max_arr_s[:it_s, :], axis=0), zrange[:kmax], color='b', label='single, ' + ID)
            ax20.plot(np.amax(max_arr_d[it_s:it_d, :], axis=0), zrange[:kmax], '-', color='g', label='double')
            ax20.plot(np.amax(max_arr_t[it_d:, :], axis=0), zrange[:kmax], '-', color='r', label='triple')
            ax21.plot(np.amax(max_arr_s[:it_s, :], axis=0), zrange[:kmax], 'b', label='single, ' + ID)
            ax21.plot(np.amax(max_arr_s[:it_s - 1, :], axis=0), zrange[:kmax], ':b', label='single, ' + ID)
            ax21.plot(np.amax(max_arr_s[:it_s + 1, :], axis=0), zrange[:kmax], ':b', label='single, ' + ID)
            ax22.plot(np.amax(max_arr_d[it_s:it_d, :], axis=0), zrange[:kmax], 'g', label='double, ' + ID)
            ax22.plot(np.amax(max_arr_d[it_s:it_d - 1, :], axis=0), zrange[:kmax], ':g', label='double, ' + ID)
            ax22.plot(np.amax(max_arr_d[it_s:it_d + 1, :], axis=0), zrange[:kmax], ':g', label='double, ' + ID)
            ax23.plot(np.amax(max_arr_t[it_d:, :], axis=0), zrange[:kmax], color='r', label='triple, ' + ID)
            ax23.plot(np.amax(max_arr_t[it_d - 1:, :], axis=0), zrange[:kmax], ':r', label='triple, ' + ID)
            ax23.plot(np.amax(max_arr_t[it_d + 1:, :], axis=0), zrange[:kmax], ':r', label='triple, ' + ID)
        elif type_min_max == 'min':
            ax20.plot(np.amin(max_arr_s[:it_s, :], axis=0), zrange[:kmax], color='b', label='single, ' + ID)
            ax20.plot(np.amin(max_arr_d[it_s:it_d, :], axis=0), zrange[:kmax], '-', color='g', label='double')
            ax20.plot(np.amin(max_arr_t[it_d:, :], axis=0), zrange[:kmax], '-', color='r', label='triple')
            ax21.plot(np.amin(max_arr_s[:it_s, :], axis=0), zrange[:kmax], 'b', label='single, ' + ID)
            ax21.plot(np.amin(max_arr_s[:it_s - 1, :], axis=0), zrange[:kmax], ':b', label='single, ' + ID)
            ax21.plot(np.amin(max_arr_s[:it_s + 1, :], axis=0), zrange[:kmax], ':b', label='single, ' + ID)
            ax22.plot(np.amin(max_arr_d[it_s:it_d, :], axis=0), zrange[:kmax], 'g', label='double, ' + ID)
            ax22.plot(np.amin(max_arr_d[it_s:it_d - 1, :], axis=0), zrange[:kmax], ':g', label='double, ' + ID)
            ax22.plot(np.amin(max_arr_d[it_s:it_d + 1, :], axis=0), zrange[:kmax], ':g', label='double, ' + ID)
            ax23.plot(np.amin(max_arr_t[it_d:, :], axis=0), zrange[:kmax], color='r', label='triple, ' + ID)
            ax23.plot(np.amin(max_arr_t[it_d - 1:, :], axis=0), zrange[:kmax], ':r', label='triple, ' + ID)
            ax23.plot(np.amin(max_arr_t[it_d + 1:, :], axis=0), zrange[:kmax], ':r', label='triple, ' + ID)
        ax21.plot(ref_array[i_id, 0, :kmax], zrange[:kmax], '--k', label='single, ' + ID)
        ax22.plot(ref_array[i_id, 1, :kmax], zrange[:kmax], '--k', label='double')
        ax23.plot(ref_array[i_id, 2, :kmax], zrange[:kmax], '--k', label='triple')
        ax21.set_title('t<'+ str(i_tw_s - delta_s[i_id]))
        ax22.set_title(str(i_tw_s - delta_s[i_id]) + '<=t<' + str(i_tw_d-delta_d[i_id]))
        ax23.set_title(str(i_tw_d - delta_d[i_id]) + '<=t')

        ax00.legend(loc='best')
        fig.suptitle(var_name + ', ' + ID)
        file.close()
        print('')
        plt.savefig(os.path.join(path_out_figs, fig_name_ + '_' + ID + '.png'))
    return



def plot_profiles_crosssections(var_name, max_array, max_array_tw, min_array, min_array_tw,
                  times_array, id_list, ncases, zrange, kmax, fig_name_):
    fig_name = fig_name_ + '_' + case_name + '.png'
    # fig_name = 'minmax_levels_all_' + var_name + '_' + case_name + '.png'
    fig, axis = plt.subplots(1, 3, figsize=(14, 10))
    ax0 = axis[0]
    ax1 = axis[1]
    ax2 = axis[2]
    for i_id, ID in enumerate(id_list):
        print('---id: '+ID, i_id, np.amax(max_array[i_id,0,:nz]))

        print(zrange)
        print(zrange.shape, max_array.shape)
        al = 1.-np.double(i_id+1)/(ncases+1)
        ax0.plot(max_array[i_id,0,:kmax], zrange[:kmax], color='b', alpha=al, label='single, '+ID)
        ax0.plot(max_array[i_id,1,:kmax], zrange[:kmax], '-', color='g', alpha=al, label='double')
        ax0.plot(max_array[i_id,2,:kmax], zrange[:kmax], '-', color='r', alpha=al, label='triple')
        ax1.plot(min_array[i_id,0,:kmax], zrange[:kmax], color='b', alpha=al, label='single, '+ID)
        ax1.plot(min_array[i_id,1,:kmax], zrange[:kmax], '-', color='g', alpha=al, label='double')
        ax1.plot(min_array[i_id,2,:kmax], zrange[:kmax], '-', color='r', alpha=al, label='triple')
        ax2.plot([1, 2], times_array[i_id, :], '-o', color=str(1.-al), label='times, ' + ID)
        ax2.fill_between(np.arange(1, 3), times_array[i_id, 0], times_array[i_id, 1],
                         facecolor=str(al), alpha=0.2, linewidth=0)
    for ax in axis[:2]:
        ax.set_xlabel(var_name)
        ax.set_ylabel('height  [m]')
        ax.grid()
    ax2.set_xlabel('collision (#CPs)')
    ax2.set_ylabel('time [s]')
    ax2.grid()
    ax2.set_xlim(0.9,2.1)
    ax0.set_title('maxima')
    ax1.set_title('minima')
    ax0.legend(loc='upper left', bbox_to_anchor=(-0.1, -0.1),
               fancybox=False, shadow=False, ncol=len(id_list), fontsize=9)
    ax2.legend(loc='upper left', bbox_to_anchor=(0.1, -0.1),
               fancybox=False, shadow=False, ncol=1, fontsize=9)
    plt.suptitle('min/max for ' + var_name, fontsize=21)
    plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.9, hspace=0.4)
    plt.savefig(os.path.join(path_out_figs, fig_name))
    plt.close()



    fig_name = fig_name_ + '_' + case_name + '_tw.png'
    zrange = np.arange(nz) * dx[2]
    fig, axis = plt.subplots(1, 3, figsize=(14, 10))
    ax0 = axis[0]
    ax1 = axis[1]
    ax2 = axis[2]
    for i_id, ID in enumerate(id_list):
        print('---id: ' + ID, i_id, np.amax(max_array[i_id, 0, :nz]))

        print(zrange)
        print(zrange.shape, max_array.shape)
        al = 1. - np.double(i_id + 1) / (ncases + 1)
        ax0.plot(max_array_tw[i_id, 0, :kmax], zrange[:kmax], color='b', alpha=al, label='single, ' + ID)
        ax0.plot(max_array_tw[i_id, 1, :kmax], zrange[:kmax], '-', color='g', alpha=al, label='double')
        ax0.plot(max_array_tw[i_id, 2, :kmax], zrange[:kmax], '-', color='r', alpha=al, label='triple')
        ax1.plot(min_array_tw[i_id, 0, :kmax], zrange[:kmax], color='b', alpha=al, label='single, ' + ID)
        ax1.plot(min_array_tw[i_id, 1, :kmax], zrange[:kmax], '-', color='g', alpha=al, label='double')
        ax1.plot(min_array_tw[i_id, 2, :kmax], zrange[:kmax], '-', color='r', alpha=al, label='triple')
        ax2.plot([1, 2], times_array[i_id, :], '-o', color=str(1. - al), label='times, ' + ID)
        ax2.fill_between(np.arange(1, 3), times_array[i_id, 0], times_array[i_id, 1],
                         facecolor=str(al), alpha=0.2, linewidth=0)
    for ax in axis[:2]:
        ax.set_xlabel(var_name)
        ax.set_ylabel('height  [m]')
        ax.grid()
    ax2.set_xlabel('collision (#CPs)')
    ax2.set_ylabel('time [s]')
    ax2.grid()
    ax2.set_xlim(0.9, 2.1)
    ax0.set_title('maxima')
    ax1.set_title('minima')
    ax0.legend(loc='upper left', bbox_to_anchor=(-0.1, -0.1),
               fancybox=False, shadow=False, ncol=len(id_list), fontsize=9)
    ax2.legend(loc='upper left', bbox_to_anchor=(0.1, -0.1),
               fancybox=False, shadow=False, ncol=1, fontsize=9)
    plt.suptitle('min/max for ' + var_name, fontsize=21)
    plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.9, hspace=0.4)
    plt.savefig(os.path.join(path_out_figs, fig_name))
    plt.close()
    return


# --------------------------------------------------------------------------------
def set_input_parameters(args):
    global path_root, path_out_figs
    path_root = args.path_root
    path_out_figs = os.path.join(path_root, 'figs_crosssections')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)

    global case_name
    case_name = args.casename
    print('')
    print('casename: ' + case_name)
    dTh = np.int(args.dTh)
    print('dTh: ' + str(dTh))

    if case_name[:21] == 'ColdPoolDry_triple_3D':
        if dTh == 3:
            zstar = 2000
            rstar = 2000
        elif dTh == 5:
            zstar = 1000
            rstar = 1100
        else:
            print('dTh not defined !!!!')
        d_range = [10, 15, 20]
        # d_range = [10,15]
        id_list = []
        for dstar in d_range:
            id_list.append('dTh' + str(dTh) +'_z' + str(zstar) + '_r' + str(rstar) + '_d' + str(dstar) + 'km')
    print(id_list)
    print('')

    ''' determine time range '''
    global tmin, tmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = np.int(0)
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = np.int(7200)
    # times = np.arange(tmin, tmax+dt_fields, dt_fields)
    # times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
    #          and np.int(name[:-3]) >= tmin and np.int(name[:-3]) <= tmax]
    # times.sort()
    # print('times', times)
    # files = [str(t) + '.nc' for t in times]
    print('tmin: '+str(tmin))
    print('tmax: '+str(tmax))
    print('')

    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = 1
    krange = np.arange(kmax + 1)
    print('krange: ', krange)
    print ''

    return kmax, krange, id_list




def define_geometry(id_list):
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
    times = np.arange(tmin, tmax+dt_fields, dt_fields)
    print('times: '+str(times))

    return dt_fields, times
# ----------------------------------


if __name__ == '__main__':
    main()



