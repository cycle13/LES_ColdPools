import argparse
import json as simplejson
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--timemin")
    parser.add_argument("--timemax")
    args = parser.parse_args()

    case_name = args.casename
    path_in = args.path
    if os.path.exists(os.path.join(path_in, 'fields')):
        path_fields = os.path.join(path_in, 'fields')
    elif os.path.exists(os.path.join(path_in, 'fields_k120')):
        path_fields = os.path.join(path_in, 'fields_k120')
    path_out = os.path.join(path_in, 'convergence')
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    global nx, ny, nz, dx, dy, dz
    nml = simplejson.loads(open(os.path.join(path_in, case_name + '.in')).read())
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    gw = nml['grid']['gw']


    if args.timemin:
        time_min = np.int(args.timemin)
    else:
        time_min = 100
    if args.timemax:
        time_max = np.int(args.timemax)
    else:
        time_max = 10000

    ''' determine file range '''
    files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc']
    if len(files[0])<7:     # 100.nc
        files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc'
             and np.int(name[:-3])>=time_min and np.int(name[:-3])<=time_max]
        times = [np.int(name[:-3]) for name in files]
        nz_ = nz
    else:       # 100_k120.nc
        files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc'
                 and np.int(name[:-8]) >= time_min and np.int(name[:-8]) <= time_max]
        times = [np.int(name[:-8]) for name in files]
        nz_ = read_in_netcdf_fields('u', os.path.join(path_fields, file)).shape[2]
    files.sort(key=len)
    times.sort()
    print('')
    print('files', files)
    print('')
    print('times', times)
    print('')


    create_conv_field(path_out)

    var_list = ['u', 'v', 'w']
    vel = np.ndarray((3, nx, ny, nz_))
    for it, file in enumerate(files):
        t0 = times[it]
        print('t', t0)
        vel[0, :, :, :] = read_in_netcdf_fields('u', os.path.join(path_fields, file))
        vel[1, :, :, :] = read_in_netcdf_fields('v', os.path.join(path_fields, file))
        vel[2, :, :, :] = read_in_netcdf_fields('w', os.path.join(path_fields, file))

        conv = compute_convergence(conv, vel)
        save_conv_field(path_out)
    return

# ----------------------------------
def compute_convergence(vel):
    conv = np.ndarray((nx, ny, nz_))
    i = 10
    j = 10
    k = 10
    
    return


def create_conv_field(path_out):
    # create file, such that at all time steps the conv field can be added
    return

def save_conv_field(path_out):
    return

#  ----------------------------------
def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:,:,:]
    rootgrp.close()
    return data

if __name__ == '__main__':
    main()