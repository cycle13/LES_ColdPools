import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
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
    # (B) cross-sections
    #   determine max v in yz-crosssections
    #   determine max u in xz-crosssections
    #

    ''' set paths & parameters '''
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    files, times, krange, nml = set_input_parameters(args)

    id = os.path.split(path_in)[1]
    print ('id: ', id)

    var_list = ['w', 's', 'temperature', 'theta']
    var_list = ['w']
    # minmax_domain = plot_minmax_domain(var_list, id, times)


    # 2CP collision at t=1200s
    # 3CP collision at t=1400s
    time_bins = [0, 1000, 1300, 2100, times[-1]]
    time_bins = [0, 1000, 1300, times[-1]]
    it_bins = [np.argwhere(np.asarray(times) == a)[0][0] for a in time_bins]
    print('bins: ', it_bins)
    minmax_levels = plot_minmax_levels(var_list, id, times, time_bins, it_bins)
    plot_minmax_levels_binned(var_list, id, times, time_bins, it_bins)

    # minmax_xz = plot_xz_minmax(id, jc_arr, times)


    return



# ----------------------------------
# ----------------------------------


def plot_xz_minmax(id, jc_arr, times):
    print('')
    print('computing min/max xz')
    var_list = ['u', 'w', 's', 'temperature']
    minmax = {}
    minmax['time'] = times
    for var_name in var_list:
        minmax[var_name] = {}
        minmax[var_name]['max'] = np.zeros(len(times), dtype=np.double)
        minmax[var_name]['min'] = np.zeros(len(times), dtype=np.double)

    for var_name in var_list:
        print('')
        print('xz: variable: ' + var_name)
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex='all', figsize=(5, 12))
        for it, t0 in enumerate(times):
            var = read_in_netcdf_fields(var_name, os.path.join(path_fields, str(t0) + '.nc'))
            minmax[var_name]['max'][it] = np.amax(var[:,jc_arr[0],:])
            minmax[var_name]['min'][it] = np.amin(var[:,jc_arr[0],:])
            del var
        maxx = ax1.plot(times, minmax[var_name]['max'][:], 'o-', label=id)
        minn = ax2.plot(times, minmax[var_name]['min'][:], 'o-', label=id)
        ax1.legend(loc='best', fontsize=10)
        ax2.legend(loc='best', fontsize=10)
        ax1.set_title('max(' + var_name + ')')
        ax1.set_ylabel('max(' + var_name + ')')
        ax2.set_title('min(' + var_name + ')')
        ax2.set_ylabel('min(' + var_name + ')')
        ax2.set_xlabel('time [s]')
        fig.suptitle(id)
        fig.savefig(os.path.join(path_out_figs, var_name + '_dTh' + str(dTh) + '_minmax_xz.png'))
        plt.close(fig)

    return minmax



# compute domain minimum and maximum of variables (s, temperature, w) for each timestep
def plot_minmax_domain(var_list, id, times):
    print('computing min/max domain')

    minmax = {}
    minmax['time'] = times
    for var_name in var_list:
        minmax[var_name] = {}
        minmax[var_name]['max'] = np.zeros(len(times), dtype=np.double)
        minmax[var_name]['min'] = np.zeros(len(times), dtype=np.double)
    # # # minmax['w'] = {}
    # # # minmax['s'] = {}
    # # # minmax['temp'] = {}

    for var_name in var_list:
        print('')
        print('variable: ' + var_name)
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex='all', figsize=(5,12))
        for it, t0 in enumerate(times):
            if var_name == 'theta':
                s_var = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0)+'.nc'))
                var = theta_s(s_var)
            else:
                var = read_in_netcdf_fields(var_name, os.path.join(path_fields, str(t0)+'.nc'))
            minmax[var_name]['max'][it] = np.amax(var[:,:,:])
            minmax[var_name]['min'][it] = np.amin(var[:,:,:])
            del var
        maxx = ax1.plot(times, minmax[var_name]['max'][:], 'o-', label=id)
        minn = ax2.plot(times, minmax[var_name]['min'][:], 'o-', label=id)
        ax1.legend(loc='best', fontsize=10)
        ax2.legend(loc='best', fontsize=10)
        ax1.set_title('max('+var_name+')')
        ax1.set_ylabel('max('+var_name+')')
        ax2.set_title('min('+var_name+')')
        ax2.set_ylabel('min('+var_name+')')
        ax2.set_xlabel('time [s]')
        fig.suptitle(id)
        fig.savefig(os.path.join(path_out_figs, var_name+'_'+str(id)+'_minmax_all.png'))
        plt.close(fig)
        print('')

    return minmax



# compute minimum and maximum of variables (s, temperature, w) at each level for each timestep
def plot_minmax_levels(var_list, id, times, time_bins, it_bins):
    print('computing min/max at each level')

    cm = plt.cm.get_cmap('coolwarm')
    zrange = dx[2] * krange

    minmax = {}
    minmax['time'] = times
    for var_name in var_list:
        minmax[var_name] = {}
        minmax[var_name]['max'] = np.zeros((len(times), kmax), dtype=np.double)
        minmax[var_name]['min'] = np.zeros((len(times), kmax), dtype=np.double)
    # # # minmax['w'] = {}
    # # # minmax['s'] = {}
    # # # minmax['temp'] = {}
    for var_name in var_list:
        print('')
        print('variable: ' + var_name)
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey='all', figsize=(12,6))
        for it, t0 in enumerate(times):
            count_color = np.double(it) / len(times)
            if var_name == 'theta':
                s_var = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0)+'.nc'))
                var = theta_s(s_var)
            else:
                var = read_in_netcdf_fields(var_name, os.path.join(path_fields, str(t0)+'.nc'))
            for k in range(kmax):
                minmax[var_name]['max'][it,k] = np.amax(var[:,:,k])
                minmax[var_name]['min'][it,k] = np.amin(var[:,:,k])
            del var
            minn = ax1.plot(minmax[var_name]['min'][it,:], zrange, '-',
                            color=cm(count_color), label='t='+str(it)+'s')
            maxx = ax2.plot(minmax[var_name]['max'][it,:], zrange, '-',
                            color=cm(count_color), label='t='+str(it)+'s')
        ax1.legend(loc='upper center', bbox_to_anchor=(2.75, 0.75),
                   fancybox=True, shadow=True, ncol=2, fontsize=10)
        fig.subplots_adjust(bottom=0.12, right=.75, left=0.1, top=0.95, wspace=0.25)
        ax1.set_xlim(np.amin(minmax[var_name]['min']), np.amax(minmax[var_name]['min']))
        ax2.set_xlim(np.amin(minmax[var_name]['max']), np.amax(minmax[var_name]['max']))
        ax1.set_xlabel('min('+var_name+')  [m/s]')
        ax2.set_xlabel('max('+var_name+')  [m/s]')
        ax1.set_ylabel('height z [m]')
        fig.suptitle(id)
        fig_name = var_name +'_'+str(id)+'_minmax_levels.png'
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)
        print('')

    return minmax


# compute minimum and maximum of variables (s, temperature, w) at each level for each timestep
def plot_minmax_levels_binned(var_list, id, times, time_bins, it_bins):
    print('computing binned min/max at each level')

    cm = plt.cm.get_cmap('coolwarm')
    cm2 = plt.cm.get_cmap('bwr')
    zrange = dx[2] * krange
    c1 = 0.
    c2 = 0.
    c3 = 0.

    minmax = {}
    minmax['time'] = times
    for var_name in var_list:
        minmax[var_name] = {}
        minmax[var_name]['max'] = np.zeros((len(times), kmax), dtype=np.double)
        minmax[var_name]['min'] = np.zeros((len(times), kmax), dtype=np.double)
    # # # minmax['w'] = {}
    # # # minmax['s'] = {}
    # # # minmax['temp'] = {}
    for var_name in var_list:
        print('')
        print('variable: ' + var_name)
        fig2, axes2 = plt.subplots(3, 2, sharey='all', figsize=(10, 15))
        for it, t0 in enumerate(times):
            count_color = np.double(it) / len(times)
            if var_name == 'theta':
                s_var = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))
                var = theta_s(s_var)
            else:
                var = read_in_netcdf_fields(var_name, os.path.join(path_fields, str(t0) + '.nc'))
            for k in range(kmax):
                minmax[var_name]['max'][it, k] = np.amax(var[:, :, k])
                minmax[var_name]['min'][it, k] = np.amin(var[:, :, k])
            del var
            if t0 < time_bins[1]:
                ax = axes2[0, :]
                c = c1 / (it_bins[1] - it_bins[0])
                c1 += 1.
            elif t0 < time_bins[2]:
                ax = axes2[1, :]
                c = c2 / (it_bins[2] - it_bins[1])
                c2 += 1.
            elif t0 < time_bins[3]:
                ax = axes2[2, :]
                c = c3 / (it_bins[3] - it_bins[2])
                print('c', c)
                c3 += 1.
            else:
                continue

            ax[0].plot(minmax[var_name]['min'][it, :], zrange, '-',
                       color=cm(c), label='t=' + str(it) + 's')
            ax[1].plot(minmax[var_name]['max'][it, :], zrange, '-',
                       color=cm(c), label='t=' + str(it) + 's')

        fig2.subplots_adjust(bottom=0.05, right=.85, left=0.1, top=.95, wspace=0.25)
        for i in range(3):
            axes2[i, 0].set_xlim(np.amin(minmax[var_name]['min']), np.amax(minmax[var_name]['min']))
            axes2[i, 1].set_xlim(np.amin(minmax[var_name]['max']), np.amax(minmax[var_name]['max']))
            axes2[i, 0].set_ylabel('height z [m]')
            axes2[i, 1].set_title('t=' + str(time_bins[i]) + ' - ' + str(time_bins[i + 1]) + 's', fontsize=18)
        axes2[0, 0].set_title('single CP', fontsize=18)
        axes2[1, 0].set_title('2-CP collision', fontsize=18)
        axes2[2, 0].set_title('3-CP collision', fontsize=18)
        for i in range(2):
            axes2[2, i].set_xlabel('min(' + var_name + ')  [m/s]')
            axes2[2, i].set_xlabel('max(' + var_name + ')  [m/s]')
        fig2.suptitle('min/max of ' + var_name, fontsize=24)
        fig_name2 = var_name + '_' + str(id) + '_minmax_binned.png'
        fig2.savefig(os.path.join(path_out_figs, fig_name2))
        plt.close(fig2)
        print('')

    return


# ----------------------------------
def theta_s(s):
    # parameters from pycles
    T_tilde = 298.15
    sd_tilde = 6864.8
    cpd = 1004.0
    th_s = T_tilde * np.exp( (s - sd_tilde)/cpd )
    return th_s

# ----------------------------------
def set_input_parameters(args):
    ''' setting parameters '''
    print''
    print'paths: '
    global path_in, path_out_data, path_out_figs, path_fields
    path_in = args.path
    if os.path.exists(os.path.join(path_in, 'fields')):
        path_fields = os.path.join(path_in, 'fields')
    elif os.path.exists(os.path.join(path_in, 'fields_k120')):
        path_fields = os.path.join(path_in, 'fields_k120')
    print path_in
    path_out_data = os.path.join(path_in, 'data_analysis')
    if not os.path.exists(path_out_data):
        os.mkdir(path_out_data)
    path_out_figs = os.path.join(path_in, 'figs_minmax')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print path_out_data
    print path_out_figs
    print path_fields

    global case_name
    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path_in, case_name + '.in')).read())
    global nx, ny, nz, dx
    dx = np.ndarray(3, dtype=np.int)
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']


    try:
        print('(ic,jc) from nml')
        ic = nml['init']['ic']
        jc = nml['init']['jc']
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
    try:
        marg = nml['init']['marg']
        print('marg from nml')
    except:
        marg = 500.
        print('marg NOT from nml')


    global krange, kmax
    if args.kmin:
        kmin = np.int(args.kmin)
    else:
        kmin = 1
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = 1
    krange = np.arange(kmin, kmax)

    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100

    ''' file range '''
    files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc']
    if len(files[0]) <= 7:  # 100.nc
        files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc'
                 and np.int(name[:-3]) >= tmin and np.int(name[:-3]) <= tmax]
        times = [np.int(name[:-3]) for name in files]
        times.sort()
        for it, t0 in enumerate(times):
            files[it] = str(t0) + '.nc'
    else:  # 100_k120.nc
        files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc'
                 and np.int(name[:-8]) >= tmin and np.int(name[:-8]) <= tmax]
        times = [np.int(name[:-8]) for name in files]
        times = times.sort()
        for it, t0 in enumerate(times):
            files[it] = str(t0) + files[0][3:]

    print('')
    print('files', files)
    print('len', len(files))
    print('')
    print('times', times)
    print('')
    print('krange', krange)
    print('')
    return files, times, krange, nml
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