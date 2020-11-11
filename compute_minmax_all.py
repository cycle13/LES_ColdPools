import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import netCDF4 as nc
import argparse
import json as simplejson
import os
import sys

label_size = 12
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 18
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['grid.linewidth'] = 20
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.minor.size'] = 2.5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.minor.size'] = 2.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['ytick.minor.width'] = 1
plt.rcParams['pdf.fonttype'] = 42         # Output Type 3 (Type3) or Type 42 (TrueType)




def main():
    print('COMPUTING MIN MAX')
    # (A) domain maximum
    #
    # (B) maximum in cross-sections
    #   determine max v in yz-crosssections
    #   determine max u in xz-crosssections


    ''' set paths & parameters '''
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path_root")
    parser.add_argument("dTh", type=int)
    parser.add_argument("--zparams", nargs='+', type=int)
    parser.add_argument('--rparams', nargs='+', type=int)
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    dTh, z_params, r_params, ic_arr, jc_arr, times = set_input_parameters(args)
    nt = len(times)

    if len(z_params) != len(r_params):
        print('wrong number of parameters! ')
        # sys.exit()

    ''' (A) plot and output domain min/max '''
    var_list = ['w', 's', 'temperature', 'theta']
    fig_name = 'dTh' + str(dTh) + '_minmax_all.png'
    minmax_file_name = 'minmax_domain.nc'
    data_minmax_domain, kmax = plot_minmax_domain(dTh, z_params, r_params, var_list, dx, times,
                                       path_root, fig_name, path_out_figs, minmax_file_name)



    ''' (B) plot min/max from  azimuthally averaged fields '''
    # read in azimuthally averaged fields
    ng = len(r_params)
    zstar = z_params[0]

    id_list = []
    for i, rstar in enumerate(r_params):
        ID = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        id_list = np.append(id_list, ID)
    n_id = len(id_list)

    path_data = 'data_analysis'
    stats_file_name = 'stats_radial_averaged.nc'
    var_list = ['w', 's', 'temperature', 'v_rad']
    # test-file
    root = nc.Dataset(os.path.join(path_root, id_list[0], path_data, stats_file_name), 'r')
    rad = root.groups['stats'].variables['r'][:]
    nr = len(rad)
    rmax = nr-30
    print('radius: '+str(rad[:3]) + '...'+str(rad[-3:]))
    dims = root.groups['dimensions'].dimensions
    nk = dims['nz'].size
    root.close()
    del dims

    for var_name in var_list:
        var = np.zeros((n_id, nt, nr, nk))
        for i, ID in enumerate(id_list):
            fullpath_in = os.path.join(path_root, ID, path_data, stats_file_name)
            print(fullpath_in)
            root = nc.Dataset(fullpath_in, 'r')
            var[i, :,:,:] = root.groups['stats'].variables[var_name][:nt,:nr,:nk]
            root.close()
        var_max_arr = np.amax(var, axis=3)
        var_min_arr = np.amin(var, axis=3)
        var_max = np.amax(var_max_arr, axis=2)
        var_min = np.amax(var_min_arr, axis=2)

        figname = var_name + '_minmax_all.png'
        fig, axes = plt.subplots(2,1, figsize=(7, 12))
        ax1 = axes[0]
        ax2 = axes[1]

        for i, ID in enumerate(id_list):
            pass
            maxx = ax1.plot(times, var_max[i, :], 'o-', label=ID)
            minn = ax2.plot(times, var_min[i, :], 'o-', label=ID)
        ax1.legend(loc='best', fontsize=10)
        ax2.legend(loc='best', fontsize=10)
        ax1.set_xlim(0, times[-1])
        ax2.set_xlim(0, times[-1])
        ax1.set_title('max(' + var_name + ')')
        ax1.set_ylabel('max(' + var_name + ')')
        ax2.set_title('min(' + var_name + ')')
        ax2.set_ylabel('min(' + var_name + ')', fontsize=18)
        ax2.set_xlabel('time [s]', fontsize=18)
        fig.suptitle(var_name)
        plt.subplots_adjust(bottom=0.1, right=.9, left=0.15, top=0.9)#, hspace=0.2, wspace=0.25)
        fig.savefig(os.path.join(path_out_figs, figname))
        plt.close(fig)

        # plot for each id
        plot_minmax_radial(var_name, rad, var_max_arr, var_min_arr, rmax, times, id_list)


    # data_minmax_xz = plot_xz_minmax(dTh, z_params, r_params,
    #                            ic_arr, jc_arr, dx, times,
    #                            case_name, path_root, path_out_figs)






    return



# --------------------------------------------------------------------

def plot_minmax_radial(var_name, rad, var_max_arr, var_min_arr, rmax, times, id_list):
    cm = plt.cm.get_cmap('coolwarm')
    for i, ID in enumerate(id_list):
        print(ID)
        fullpath_out = os.path.join(path_root, ID, 'figs_minmax')
        if not os.path.exists(fullpath_out):
            os.mkdir(fullpath_out)
        figname = var_name + '_minmax_' + ID + '_rad.png'
        fig, axes = plt.subplots(2, 1, figsize=(6, 12))
        ax1 = axes[0]
        ax2 = axes[1]
        for it, t0 in enumerate(times):
            count_color = np.double(it) / len(times)
            maxx = ax1.plot(rad[:rmax], var_max_arr[i, it, :rmax], '-', color=cm(count_color), label='t=' + str(t0))
            minn = ax2.plot(rad[:rmax], var_min_arr[i, it, :rmax], '-', color=cm(count_color), label='t=' + str(t0))
        ax1.legend(loc='best', fontsize=10)
        ax1.set_xlim(0, rmax * dx[0])
        ax2.set_xlim(0, rmax * dx[0])
        ax1.set_title('max(' + var_name + ')')
        ax1.set_ylabel('max(' + var_name + ')')
        ax2.set_title('min(' + var_name + ')')
        ax2.set_xlabel('radius [m]', fontsize=18)
        fig.suptitle(var_name + '  (' + ID + ')')
        fig.savefig(os.path.join(fullpath_out, figname))
        plt.close(fig)
    return

# ----------------------------------

def plot_xz_minmax(dTh, z_params, r_params, ic_arr, jc_arr, dx, times,
                   case_name, path_root, path_out_figs):
    print('')
    print('compting min/max xz')
    var_list = ['u', 'w', 's', 'temperature']
    minmax = {}
    minmax['time'] = times
    for var_name in var_list:
        minmax[var_name] = {}
        minmax[var_name]['max'] = np.zeros(len(times), dtype=np.double)
        minmax[var_name]['min'] = np.zeros(len(times), dtype=np.double)

    ng = len(z_params)
    kmax = np.amax(z_params) + 2000./dx[2]
    for var_name in var_list:
        print('')
        print('xz: variable: ' + var_name)
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex='all', figsize=(5, 12))
        for istar in range(ng):
            zstar = z_params[istar]
            rstar = r_params[istar]
            ID = 'dTh' + str(dTh)+ '_z' + str(zstar) + '_r' + str(rstar)
            print('ID', ID)
            path_fields = os.path.join(path_root, ID, 'fields')
            print(path_fields)

            for it, t0 in enumerate(times):
                if var_name == 'theta':
                    s_var = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))
                    var = theta_s(s_var)
                else:
                    var = read_in_netcdf_fields(var_name, os.path.join(path_fields, str(t0) + '.nc'))
                minmax[var_name]['max'][it] = np.amax(var[:,jc_arr[0], :kmax])
                minmax[var_name]['min'][it] = np.amin(var[:,jc_arr[0], :kmax])
                del var
            maxx = ax1.plot(times, minmax[var_name]['max'][:], 'o-', label=ID)
            minn = ax2.plot(times, minmax[var_name]['min'][:], 'o-', label=ID)
        ax1.legend(loc='best', fontsize=10)
        ax2.legend(loc='best', fontsize=10)
        ax1.set_title('max(' + var_name + ')')
        ax1.set_ylabel('max(' + var_name + ')')
        ax2.set_title('min(' + var_name + ')')
        ax2.set_ylabel('min(' + var_name + ')')
        ax2.set_xlabel('time [s]')
        fig.suptitle('dTh=' + str(dTh))
        fig.savefig(os.path.join(path_out_figs, var_name + '_dTh' + str(dTh) + '_minmax_xz.png'))
        plt.close(fig)

    return minmax


# compute domain minimum and maximum of variables (s, temperature, w) for each timestep
def plot_minmax_domain(dTh, z_params, r_params,
                       var_list, dx, times,
                       path_root, save_name, path_out_figs, filename):
    print('computing min/max domain')

    # nth = geom_params.shape[0]
    ng = len(z_params)
    kmax = (np.amax(z_params) + 2000.) / dx[2]

    minmax = {}
    nt_ranges = np.zeros(ng)
    time_range = {}
    for istar in range(ng):
        zstar = z_params[istar]
        rstar = r_params[istar]
        ID = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)

        print('FILE: ', os.path.join(path_root, ID, 'data_analysis', filename))
        minmax_file = nc.Dataset(os.path.join(path_root, ID, 'data_analysis', filename), 'r', format='NETCDF4')
        #minmax_file = nc.Dataset(os.path.join(path_root, ID+'_dz2000_run2', 'data_analysis', filename), 'r', format='NETCDF4')
        ts_grp = minmax_file.groups['timeseries']

        minmax[ID] = {}
        nt_ranges[istar] = ts_grp.dimensions['nt'].size
        time_range[ID] = ts_grp.variables['time'][:]
        for var_name in var_list:
            minmax[ID][var_name] = {}
            minmax[ID][var_name]['max'] = ts_grp[var_name+'_max'][:]
            minmax[ID][var_name]['min'] = ts_grp[var_name+'_min'][:]
        minmax_file.close()

    for var_name in var_list:
        print('')
        print('variable: ' + var_name)
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex='all', figsize=(5, 12))
        for istar in range(ng):
            zstar = z_params[istar]
            rstar = r_params[istar]
            ID = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
            maxx = ax1.plot(time_range[ID], minmax[ID][var_name]['max'][:], 'o-', label=ID)
            minn = ax2.plot(time_range[ID], minmax[ID][var_name]['min'][:], 'o-', label=ID)
        ax1.legend(loc='best', fontsize=10)
        ax2.legend(loc='best', fontsize=10)
        ax1.set_title('max(' + var_name + ')')
        ax1.set_ylabel('max(' + var_name + ')')
        ax2.set_title('min(' + var_name + ')')
        ax2.set_ylabel('min(' + var_name + ')')
        ax2.set_xlabel('time [s]')
        fig.suptitle('dTh='+str(dTh))
        # fig.savefig(os.path.join(path_out_figs, var_name+'_dTh'+str(dTh)+'_minmax_all.png'))
        print('saving: ' + os.path.join(path_out_figs, var_name + '_' + save_name))
        fig.savefig(os.path.join(path_out_figs, var_name + '_' + save_name))
        plt.close(fig)

    return minmax, kmax



# ----------------------------------
def theta_s(s):
    # parameters from pycles
    T_tilde = 298.15
    sd_tilde = 6864.8
    cpd = 1004.0
    th_s = T_tilde * np.exp( (s - sd_tilde)/cpd )
    return th_s

# --------------------------------------------------------------------
def set_geom_parameters(dTh):

    if dTh == 1:
        z_params = [3465, 1730, 1155]  # run1
    elif dTh == 2:
        # z_params = [2450, 1225, 815]  # run1
        z_params = [500, 900, 1600, 1900, 2500]  # run2
        r_params_ = [1900, 1300, 900, 800, 600]  # run2
    elif dTh == 3:
        # z_params = [4000, 2000, 1500, 1000, 670, 500, 250] # run1
        z_params = [500, 1000, 1600, 2000, 2500]  # run2
        r_params_ = [1500, 1000, 700, 600, 500]  # run2
    elif dTh == 4:
        # z_params = [1730, 870, 430]     # run1
        z_params = [500, 900, 1600, 2000, 2500]  # run2
        r_params_ = [1300, 900, 600, 500, 400]  # run2


    try:
        r_params = r_params_
        del r_params_
    except:
        r_params = z_params[::-1]
    print('z*: ', z_params)
    print('r*: ', r_params)
    print('')

    return z_params, r_params


def set_input_parameters(args):

    global path_root, path_out_data, path_out_figs
    path_root = args.path_root
    path_out_data = os.path.join(path_root, 'data_analysis')
    if not os.path.exists(path_out_data):
        os.mkdir(path_out_data)
    path_out_figs = os.path.join(path_root, 'figs_minmax')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print('')
    print('paths: ')
    print('path data in: ' + path_root)
    print('path data out: ' + path_out_data)
    print('path figures:  ' + path_out_figs)
    print('')

    dTh = args.dTh
    z_params = args.zparams
    r_params = args.rparams
    print('z*: ', z_params)
    print('r*: ', r_params)

    global case_name
    case_name = args.casename
    ID0 = 'dTh' + str(dTh) + '_z' + str(z_params[0]) + '_r' + str(r_params[0])
    nml = simplejson.loads(open(os.path.join(path_root, ID0, case_name + '.in')).read())
    global nx, ny, nz, dx
    dx = np.ndarray(3, dtype=np.int)
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    # gw = nml['grid']['gw']


    try:
        print('(ic,jc) from nml')
        ic = nml['init']['ic']
        jc = nml['init']['jc']
        ic_arr = np.zeros(1)
        jc_arr = np.zeros(1)
        ic_arr[0] = ic
        jc_arr[0] = jc
    except:
        print('(ic,jc) NOT from nml')
        if case_name == 'ColdPoolDry_single_3D':
            ic = np.int(nx/2)
            jc = np.int(ny/2)
            ic_arr = np.zeros(1)
            jc_arr = np.zeros(1)
            ic_arr[0] = ic
            jc_arr[0] = jc
        else:
            print('ic, jc not defined')
    # try:
    #     marg = nml['init']['marg']
    #     print('marg from nml')
    # except:
    #     marg = 500.
    #     print('marg NOT from nml')

    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100
    times = np.arange(tmin, tmax + 100, 100)

    print('times', times)
    print('')

    return dTh, z_params, r_params, ic_arr, jc_arr, times


# --------------------------------------------------------------------

def dump_minmax_file(data, var_list, times, kmax,
                     dTh, z_params, r_params,
                     file_name, path_root):
    print(' ')
    print('-------- dump minmax data -------- ')
    # print(os.path.join(path_out_data, file_name))

    nt = len(times)
    ng = len(z_params)
    for istar in range(ng):
        zstar = z_params[istar]
        rstar = r_params[istar]
        ID = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        path_out_data = os.path.join(path_root, ID, 'data_analysis')

        print(os.path.join(path_out_data, file_name))
        rootgrp = nc.Dataset(os.path.join(path_out_data, file_name), 'w', format='NETCDF4')

        dims_grp = rootgrp.createGroup('dimensions')
        # dims_grp.createDimension('dx', dx[0])
        dims_grp.createDimension('nz', kmax)

        ts_grp = rootgrp.createGroup('timeseries')
        ts_grp.createDimension('nt', nt)
        var = ts_grp.createVariable('time', 'f8', ('nt'))
        var.unit = "s"
        var[:] = times

        for var_name in var_list:
            var = ts_grp.createVariable(var_name+'_max', 'f8', ('nt'))
            var[:] = data[ID][var_name]['max'][:]
            var = ts_grp.createVariable(var_name + '_min', 'f8', ('nt'))
            var[:] = data[ID][var_name]['min'][:]

        rootgrp.close()
    print('')
    return

# ----------------------------------



def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:,:,:]
    rootgrp.close()
    return data

if __name__ == '__main__':
    main()