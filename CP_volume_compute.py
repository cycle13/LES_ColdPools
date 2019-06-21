import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import netCDF4 as nc
import argparse
import json as simplejson
import os
import time


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35

def main():
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')
    cm_fall = plt.cm.get_cmap('winter')
    cm_summer = plt.cm.get_cmap('spring')

    times, nml, ID = set_input_parameters(args)



    ''' create output file '''
    filename_CP_vol = 'CP_volume_' + ID + '.nc'
    create_output_file(filename_CP_vol, times)



    ''' (1) from CP height (threshold on entropy)'''
    s_crit = 0.5
    filename_CP_height = 'CP_height_' + ID + '_sth' + str(s_crit) + '.nc'
    fullpath_in = os.path.join(path, 'data_analysis', filename_CP_height)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    time_CP_height = rootgrp.groups['timeseries'].variables['time'][:]
    CP_height_radav = rootgrp.groups['stats'].variables['CP_height_rad'][:, :]
    r_range = rootgrp.groups['stats'].variables['r'][:]
    rootgrp.close()

    nt = len(time_CP_height)
    CP_vol = np.zeros(nt)
    for ir, r in enumerate(r_range):
        CP_vol[:] += r_range[ir]*CP_height_radav[:,ir]
    dr = dx[0]
    CP_vol = 2*np.pi*dr*CP_vol

    # plotting
    plot_CP_volume(CP_vol, CP_height_radav, time_CP_height, r_range, ID)

    # output CP volume
    s_crit = 5e-1
    var_name = 'CP_vol_sth' + str(s_crit)
    dump_CP_height(CP_vol, var_name, filename_CP_vol)


    return

# ----------------------------------------------------------------------
def plot_CP_volume(CP_vol, CP_height_radav, time_CP_height, r_range, ID):
    rad = [np.where(CP_height_radav[it, :] == np.amax(CP_height_radav[it, :]))[0][0]
           for it, t0 in enumerate(time_CP_height)]

    rmax_plot = 8e3
    fig, axis = plt.subplots(1, 3, sharex='none', figsize=(18, 5))
    ax0 = axis[0]
    ax1 = axis[1]
    ax2 = axis[2]
    time_list = np.arange(0, 3700, 600)

    for it_, t0 in enumerate(time_list):
        if t0 <= tmin and t0 <= tmax:
            it = np.where(time_CP_height == t0)[0][0]
            ax0.plot(r_range, CP_height_radav[it, :], label='h_max=' + str(np.amax(CP_height_radav[it, :])))
            ax1.plot([t0, t0], [0, 1e4], '.25', linewidth=0.5)
            ax2.plot([t0, t0], [0, 1e10], '.25', linewidth=0.5)
            ax2.plot([0, tmax], [CP_vol[it], CP_vol[it]], '.25', linewidth=0.5, label='V=' + str(CP_vol[it]))

    ca = 'k'
    ax1.plot(time_CP_height, np.amax(CP_height_radav, axis=1), color=ca, label='max(CP height)')
    ax1.set_ylabel('max(CP height) [m]', color='black')
    ax1b = ax1.twinx()
    cb = '0.4'
    ax1b.plot(time_CP_height, rad, color=cb, label='radius of max(CP height)')
    ax1b.set_ylabel('radius of max(CP height) [m]', color=cb)
    ax1b.tick_params(axis='y', labelcolor=cb)
    ax2.plot(time_CP_height, CP_vol)

    ax0.legend(loc='best')
    ax2.legend(loc='best')

    ax0.set_title('CP height')
    ax1.set_title('max. CP height and corresp. radius')
    ax2.set_title('CP volume')
    ax0.set_xlabel('radius [m]')
    ax1.set_xlabel('time  [s]')
    ax2.set_xlabel('time  [s]')
    ax0.set_ylabel('CP height [m]')
    ax2.set_ylabel('CP volume [m3]')
    ax0.set_xlim(0, rmax_plot)
    ax1.set_xlim(0, tmax)
    ax2.set_xlim(0, tmax)
    ax1.set_ylim(0, 1e3)
    ax2.set_ylim(0, np.amax(CP_vol))
    plt.suptitle(ID)
    fig.tight_layout()
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.04, top=0.8, wspace=0.25)
    fig.savefig(os.path.join(path_out_figs, 'CP_vol_radav.png'))
    plt.close(fig)

    return


def dump_CP_height(CP_vol, var_name, filename_CP_vol):
    rootgrp = nc.Dataset(os.path.join(path_out_data, filename_CP_vol), 'r+')
    tsgrp = rootgrp.groups['timeseries']
    print var_name
    if not var_name in tsgrp.variables.keys():
        var = tsgrp.createVariable(var_name, 'f8', ('nt'))
        var.long_name = 'CP volume from CP_height_sth'
        var.units = "m"
    else:
        var = tsgrp.variables[var_name]
    var[:] = CP_vol[:]
    rootgrp.close()

    return

# ----------------------------------------------------------------------
def create_output_file(filename, times):
    if os.path.exists(os.path.join(path_out_data, filename)):
        print('CP_vol stats-file already existing! Not creating new file')
        print ''
        return
    else:
        nt = len(times)
        print('create output file in: ', path_out_data)
        print('size: ', nz, nt)
        print ''

        rootgrp = nc.Dataset(os.path.join(path_out_data, filename), 'w', format='NETCDF4')

        ts_grp = rootgrp.createGroup('timeseries')
        ts_grp.createDimension('nt', nt)
        var = ts_grp.createVariable('time', 'f8', ('nt'))
        var.units = "s"
        var[:] = times

        rootgrp.close()

        return
# ----------------------------------------------------------------------
def set_input_parameters(args):
    print ''' setting parameters '''
    global path, path_fields, path_out_figs, path_out_data
    path = args.path
    path_fields = os.path.join(path, 'fields')
    path_out_figs = os.path.join(path, 'figs_CP_volume')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    path_out_data = os.path.join(path, 'data_analysis')
    if not os.path.exists(path_out_data):
        os.mkdir(path_out_data)

    global case_name
    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    ID = os.path.basename(path[:])
    if ID == '':
        ID = os.path.basename(path[:-1])
    global nx, ny, nz, dx, dV, gw
    dx = np.ndarray(3, dtype=np.int)
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx[0]*dx[1]*dx[2]

    # global kmax
    # if args.kmax:
    #     kmax = np.int(args.kmax)+1
    # else:
    #     kmax = nz

    global tmin, tmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100

    ''' time range '''
    times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
             and tmin <= np.int(name[:-3]) <= tmax]
    times.sort()

    print('')
    print('times', times)
    print('')
    # print('kmax ', kmax, 'nx ', nx)
    # print('')

    return times, nml, ID

# ----------------------------------------

if __name__ == '__main__':
    main()
