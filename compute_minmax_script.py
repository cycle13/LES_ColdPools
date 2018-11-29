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

def compute_minmax(thermo_params, geom_params,
                   kmin, kmax, times, case_name, path_root, path_out_figs, path_out_data):
# def compute_minmax(dTh_range, zstar_range, rstar_range,
#                    kmin, kmax, times, path_root, case_name):
    print('COMPUTING MIN MAX')
    # (A) full fields
    # determine time range
    # determine domain (max height, evtl. lateral boundaries) for analysis
    # loop over time
    #   read in fields
    #   compute theta
    #   determine max of w, s, temperature, theta
    #
    # (B) cross-sections
    # determin CP center
    # extract crosssection through CP center
    # loop over time
    #   determine max v in yz-crosssections
    #   determine max u in xz-crosssections
    #

    # ''' set paths & parameters '''
    # # Parse information from the command line
    # parser = argparse.ArgumentParser(prog='LES_CP')
    # parser.add_argument("casename")
    # parser.add_argument("path")
    # # parser.add_argument("zparams")
    # parser.add_argument('--zparams', nargs='+', type=int)
    # parser.add_argument("--rparams", nargs='+')
    # parser.add_argument("--tmin")
    # parser.add_argument("--tmax")
    # parser.add_argument("--kmin")
    # parser.add_argument("--kmax")
    # args = parser.parse_args()

    # times = set_inpout_parameters()

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
    #
    nth = geom_params.shape[0]
    ng = geom_params.shape[2]
    for var_name in var_list:
        print('')
        print('variable: ' + var_name)
        dTh = thermo_params[0]
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex='all', figsize=(5,12))
        for istar in range(ng):
            id = 'dTh'+str(dTh)+'_z'+str(geom_params[0,0,istar])+'_r'+str(geom_params[0,1,istar])
            print('id', id)
            path_fields = os.path.join(path_root, id, 'fields')
            print(path_fields)
            for it, t0 in enumerate(times):
                var = read_in_netcdf_fields(var_name, os.path.join(path_fields, str(t0)+'.nc'))
                minmax[var_name]['max'][it] = np.amax(var)
                minmax[var_name]['min'][it] = np.amin(var)
            maxx = ax1.plot(times, minmax[var_name]['max'][:], 'o-', label=id)
            minn = ax2.plot(times, minmax[var_name]['min'][:], 'o-', label=id)
        ax1.legend(loc='best')
        ax2.legend(loc='best')
        ax1.set_title('max('+var_name+')')
        ax1.set_ylabel('max('+var_name+')')
        ax2.set_title('min('+var_name+')')
        ax2.set_xlabel('time [s]')
        ax2.set_ylabel('max('+var_name+')')
        fig.suptitle('dTh='+str(dTh))
        fig.savefig(os.path.join(path_out_figs, var_name+'_dTh'+str(dTh)+'_minmax_all.png'))
        plt.close(fig)

    return



# ----------------------------------
# ----------------------------------
def set_inpout_parameters():
    print '''setting parameters '''
    global path_root, path_out_data, path_out_figs
    path_root = args.path
    path_out_data = os.path.join(path_root, 'data_analysis')
    if not os.path.exists(path_out_data):
        os.mkdir(path_out_data)
    path_out_figs = os.path.join(path_root, 'figs_minmax')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print('paths:')
    print path_root
    print path_out_data
    print path_out_figs

    # global case_name
    # case_name = args.casename
    # nml = simplejson.loads(open(os.path.join(path_in, case_name + '.in')).read())
    # global nx, ny, nz, dx, dy, dz
    # nx = nml['grid']['nx']
    # ny = nml['grid']['ny']
    # nz = nml['grid']['nz']
    # dx = nml['grid']['dx']
    # dy = nml['grid']['dy']
    # dz = nml['grid']['dz']
    # gw = nml['grid']['gw']

    global zstars, rstars
    zstars = args.zparams
    rstars = args.rparams
    print('zstars:', zstars)
    print('rstars:', rstars)

    global kmin, kmax
    if args.kmin:
        kmin = np.int(args.kmin)
    else:
        kmin = 1
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = 1

    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100
    dt = 100
    times = np.arange(tmin, tmax+dt, dt)

    # ''' file range '''
    # files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc']
    # if len(files[0]) <= 7:  # 100.nc
    #     files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc'
    #              and np.int(name[:-3]) >= tmin and np.int(name[:-3]) <= tmax]
    #     times = [np.int(name[:-3]) for name in files]
    #     times.sort()
    #     for it, t0 in enumerate(times):
    #         files[it] = str(t0) + '.nc'
    # else:  # 100_k120.nc
    #     files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc'
    #              and np.int(name[:-8]) >= tmin and np.int(name[:-8]) <= tmax]
    #     times = [np.int(name[:-8]) for name in files]
    #     times = times.sort()
    #     for it, t0 in enumerate(times):
    #         files[it] = str(t0) + files[0][3:]
    #
    # print('')
    # print('files', files)
    # print('len', len(files))
    # print('')
    print('times', times)
    # print('')
    print('krange', kmin, kmax)
    print('')

    return times
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