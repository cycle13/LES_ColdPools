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
    parser.add_argument("path_root")
    parser.add_argument("dTh", type=int)
    parser.add_argument("--zparams", nargs='+', type=int)
    parser.add_argument('--rparams', nargs='+', type=int)
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    dTh, z_params, r_params, ic_arr, jc_arr, marg, times = set_input_parameters(args)

    if len(z_params) != len(r_params):
        print('wrong number of parameters! ')
        sys.exit()

    # minmax_domain = plot_domain_minmax(dTh, z_params, r_params, dx, times,
    #                                    case_name, path_root, path_out_figs, path_out_data)
    #

    minmax_xz = plot_xz_minmax(dTh, z_params, r_params,
                               ic_arr, jc_arr, dx, times,
                               case_name, path_root, path_out_figs, path_out_data)


    return



# ----------------------------------
# ----------------------------------


def plot_xz_minmax(dTh, z_params, r_params, ic_arr, jc_arr, dx, times,
                   case_name, path_root, path_out_figs, path_out_data):
    print('')
    print('COMPUTING MIN MAX (XZ)')
    # var_list = ['u', 'w', 's', 'temperature']
    var_list = ['u']
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
            id = 'dTh'+str(dTh)+'_z'+str(zstar)+'_r'+str(rstar)
            print('id', id)
            path_fields = os.path.join(path_root, id, 'fields')
            print(path_fields)
            for it, t0 in enumerate(times):
                var = read_in_netcdf_fields(var_name, os.path.join(path_fields, str(t0) + '.nc'))
                minmax[var_name]['max'][it] = np.amax(var[:,jc_arr[0],:kmax])
                minmax[var_name]['min'][it] = np.amin(var[:,jc_arr[0],:kmax])
                del var
            maxx = ax1.plot(times, minmax[var_name]['max'][:], 'o-', label=id)
            minn = ax2.plot(times, minmax[var_name]['min'][:], 'o-', label=id)
        ax1.legend(loc='best', fontsize=10)
        ax2.legend(loc='best', fontsize=10)
        ax1.set_title('max(' + var_name + ')')
        ax1.set_ylabel('max(' + var_name + ')')
        ax2.set_title('min(' + var_name + ')')
        ax2.set_xlabel('time [s]')
        ax2.set_ylabel('max(' + var_name + ')')
        fig.suptitle('dTh=' + str(dTh))
        fig.savefig(os.path.join(path_out_figs, var_name + '_dTh' + str(dTh) + '_minmax_xy.png'))
        plt.close(fig)

    return minmax



# compute domain minimum and maximum of variables (s, temperature, w) for each timestep
def plot_domain_minmax(dTh, z_params, r_params, dx, times,
                       case_name, path_root, path_out_figs, path_out_data):
    print('COMPUTING MIN MAX')

    var_list = ['w', 's', 'temperature']
    minmax = {}
    minmax['time'] = times
    for var_name in var_list:
        minmax[var_name] = {}
        minmax[var_name]['max'] = np.zeros(len(times), dtype=np.double)
        minmax[var_name]['min'] = np.zeros(len(times), dtype=np.double)
    # # # minmax['w'] = {}
    # # # minmax['s'] = {}
    # # # minmax['temp'] = {}

    # nth = geom_params.shape[0]
    ng = len(z_params)
    kmax = np.amax(z_params) + 2000. / dx[2]

    for var_name in var_list:
        print('variable: ' + var_name)
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex='all', figsize=(5,12))
        for istar in range(ng):
            zstar = z_params[istar]
            rstar = r_params[istar]
            id = 'dTh'+str(dTh)+'_z'+str(zstar)+'_r'+str(rstar)
            print('id', id)
            path_fields = os.path.join(path_root, id, 'fields')
            print(path_fields)
            for it, t0 in enumerate(times):
                var = read_in_netcdf_fields(var_name, os.path.join(path_fields, str(t0)+'.nc'))
                minmax[var_name]['max'][it] = np.amax(var[:,:,:kmax])
                minmax[var_name]['min'][it] = np.amin(var[:,:,:kmax])
                del var
            maxx = ax1.plot(times, minmax[var_name]['max'][:], 'o-', label=id)
            minn = ax2.plot(times, minmax[var_name]['min'][:], 'o-', label=id)
        ax1.legend(loc='best', fontsize=10)
        ax2.legend(loc='best', fontsize=10)
        ax1.set_title('max('+var_name+')')
        ax1.set_ylabel('max('+var_name+')')
        ax2.set_title('min('+var_name+')')
        ax2.set_xlabel('time [s]')
        ax2.set_ylabel('max('+var_name+')')
        fig.suptitle('dTh='+str(dTh))
        fig.savefig(os.path.join(path_out_figs, var_name+'_dTh'+str(dTh)+'_minmax_all.png'))
        plt.close(fig)
        print('')

    return minmax



# ----------------------------------
def set_input_parameters(args):

    global path_root, path_out_data, path_out_figs
    path_root = args.path_root
    path_out_data = os.path.join(path_root, 'data_analysis')
    if not os.path.exists(path_out_data):
        os.mkdir(path_out_data)
    path_out_figs = os.path.join(path_root, 'figs_minmax')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print path_root
    print path_out_data
    print path_out_figs


    dTh = args.dTh
    z_params = args.zparams
    r_params = args.rparams
    print('z*: ', z_params)
    print('r*: ', r_params)

    global case_name
    case_name = args.casename
    id0 = 'dTh' + str(dTh) + '_z' + str(z_params[0]) + '_r' + str(r_params[0])
    nml = simplejson.loads(open(os.path.join(path_root, id0, case_name + '.in')).read())
    global nx, ny, nz, dx
    dx = np.arange(3)
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

    return dTh, z_params, r_params, ic_arr, jc_arr, marg, times


# ----------------------------------
def theta_s(s):
    # parameters from pycles
    T_tilde = 298.15
    sd_tilde = 6864.8
    cpd = 1004.0
    th_s = np.exp( (s - sd_tilde)/cpd )
    return th_s


def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:,:,:]
    rootgrp.close()
    return data

if __name__ == '__main__':
    main()