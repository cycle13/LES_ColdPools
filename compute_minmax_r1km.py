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
    print('COMPUTING MIN MAX r=1km')
    # (A) domain maximum
    #
    # (B) maximum in cross-sections
    #   determine max v in yz-crosssections
    #   determine max u in xz-crosssections
    #

    ''' set paths & parameters '''
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path_root")
    parser.add_argument("dTh_params", nargs=3, type=int)
    parser.add_argument("zparams", nargs=3, type=int)
    parser.add_argument('rparams', nargs=3, type=int)
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    dTh_params, z_params, r_params, ic_arr, jc_arr, marg, times = set_input_parameters(args)

    if len(z_params) != len(r_params):
        print('wrong number of parameters! ')
        sys.exit()


    # plot for r=1km
    save_name = 'r1km_minmax_all.png'
    dTh_params = [2, 3, 4]
    plot_domain_minmax_r1km(dTh_params, z_params, r_params,
                            times, case_name, path_root, save_name, path_out_figs)



    return



# ----------------------------------
# ----------------------------------




# compute domain minimum and maximum of variables (s, temperature, w) for each timestep
def plot_domain_minmax_r1km(dTh_params, z_params, r_params, times,
                       case_name, path_root, save_name, path_out_figs):
    print('compting min/max domain')

    var_list = ['w', 's', 'temperature', 'theta']
    # var_list = ['theta']
    minmax = {}
    minmax['time'] = times
    for var_name in var_list:
        minmax[var_name] = {}
        minmax[var_name]['max'] = np.zeros(len(times), dtype=np.double)
        minmax[var_name]['min'] = np.zeros(len(times), dtype=np.double)
    # # # minmax['w'] = {}
    # # # minmax['s'] = {}
    # # # minmax['temp'] = {}

    kmax = 1e4

    for var_name in var_list:
        print('')
        print('variable: ' + var_name)
        # (1) r = 1km, varying dTh
        zstar_max = 0
        rstar_max = 0
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex='all', figsize=(5, 12))
        for istar, dTh in enumerate(dTh_params):
            zstar = z_params[istar]
            rstar = r_params[istar]
            zstar_max = np.maximum(zstar_max, zstar)
            rstar_max = np.maximum(rstar_max, rstar)
            id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
            path_fields = os.path.join(path_root, id, 'fields')
            for it, t0 in enumerate(times):
                if var_name == 'theta':
                    s_var = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0)+'.nc'))
                    var = theta_s(s_var)
                else:
                    var = read_in_netcdf_fields(var_name, os.path.join(path_fields, str(t0)+'.nc'))
                minmax[var_name]['max'][it] = np.amax(var[:,:,:kmax])
                minmax[var_name]['min'][it] = np.amin(var[:,:,:kmax])
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
        fig.suptitle('r = 1km')
        fig.savefig(os.path.join(path_out_figs, var_name + '_' + save_name))
        plt.close(fig)

    return minmax



# ----------------------------------
def theta_s(s):
    # parameters from pycles
    T_tilde = 298.15
    sd_tilde = 6864.8
    cpd = 1004.0
    th_s = T_tilde * np.exp( (s - sd_tilde)/cpd )
    return th_s

# ----------------------------------
# ----------------------------------
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
    print path_root
    print path_out_data
    print path_out_figs

    dTh_params = args.dTh_params
    z_params = args.zparams
    r_params = args.rparams
    print('dTh: ', dTh_params)
    print('z*: ', z_params)
    print('r*: ', r_params)

    global case_name
    case_name = args.casename
    id0 = 'dTh' + str(dTh_params[0]) + '_z' + str(z_params[0]) + '_r' + str(r_params[0])
    nml = simplejson.loads(open(os.path.join(path_root, id0, case_name + '.in')).read())
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

    return dTh_params, z_params, r_params, ic_arr, jc_arr, marg, times


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