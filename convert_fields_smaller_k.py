import netCDF4 as nc
import argparse
import json as simplejson
import os
import numpy as np
import sys


def main():

    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    # parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--kmax")
    parser.add_argument("--kmin")
    parser.add_argument("--k0")
    args = parser.parse_args()


    # case_name = args.casename
    path_in = args.path
    # nml = simplejson.loads(open(os.path.join(path_in, case_name + '.in')).read())

    path_fields = os.path.join(path_in, 'fields')
    if args.kmax:
        k_max = np.int(args.kmax)
    else:
        k_max = 120     # leads to about 20% reduction of size for fields of original size 200x200x150
    if args.kmin:
        k_min = np.int(args.kmin)
        path_out = os.path.join(path_in, 'fields_k' + str(k_min) + '_' + str(k_max))
    else:
        k_min = k_max
        path_out = os.path.join(path_in, 'fields_k' + str(k_max))
    if args.k0:
        k0 = np.int(args.k0)
        path_out = os.path.join(path_in, 'fields_k' + str(k0))
    else:
        k0 = -1
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    print('')
    print('path in', path_in)
    print('path out', path_out)
    print('')
    print('k0', k0)
    print('kmin', k_min)
    print('kmax', k_max)
    print('')


    # (1) for each time/file in path_fields do...


    # (2) read in original fields file
    #       >> field_keys
    #       >> how to access dimensions of variables?
    #       >> data

    # (3) create new fields file, e.g. 100_k_kmax.nc

    # (4) save data[:,:,:k_max] in new variable of new fields file

    files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc']
    print(files)
    print('')

    ''' reduce number of vertical levels; leave all variables and other dimensions '''
    # convert_file_forall_variables(files, path_fields, path_out, k_min, k_max)

    # ''' reduce number of vertical levels and output file with only varialbes given in var_list '''
    # var_list = ['u', 'v', 'w']
    # convert_file_for_varlist(var_list, files, path_fields, path_out, k_min, k_max)

    ''' output file with one level of one variable for all times '''
    var_list = ['u', 'v']
    for var in var_list:
        if k0 >= 0:
            convert_file_for_singlevariable(var, files, path_fields, path_out, k0)
    # var_list = ['u', 'v']
    # for var in var_list:
    #     file_for_simple_array(var, files, path_fields, path_out, k0)
    return



def convert_file_for_singlevariable(var, files, path_fields, path_out, k0):
    print('var', var)
    times = [np.int(name[:-3]) for name in files]
    times.sort()
    print('times', times)
    files = [str(t) + '.nc' for t in times]
    nt = len(times)
    print('nt:', nt)
    t_ini = times[0]
    t_end = times[-1]
    # file_name = var + '_merged_t' +str(t_ini) + '_t' + str(t_end) + '_k' + str(k0) + '.nc'
    file_name = var + '_merged.nc'
    fullpath_out = os.path.join(path_out,file_name)
    print('filename', file_name)

    if os.path.exists(fullpath_out):
        print('')
        print('file ' + fullpath_out + ' already exists! ')
        print('')
    else:

        # read in test fields file
        fullpath_in = os.path.join(path_fields, files[0])
        rootgrp_in = nc.Dataset(fullpath_in, 'r')
        # field_keys = rootgrp_in.groups['fields'].variables.keys()
        # dims_keys = rootgrp_in.groups['fields'].dimensions.keys()
        dims = rootgrp_in.groups['fields'].dimensions
        nx = dims['nx'].size
        ny = dims['ny'].size
        nz = dims['nz'].size
        rootgrp_in.close()

        rootgrp_out = nc.Dataset(fullpath_out, 'w', format='NETCDF4')
        rootgrp_out.createDimension('time', None)
        # rootgrp_out.createDimension('time', nt)
        rootgrp_out.createDimension('nx', nx)
        rootgrp_out.createDimension('ny', ny)
        time_out = rootgrp_out.createVariable('time', 'f8', ('time',))
        time_out.long_name = 'Time'
        time_out.units = 's'
        time_out[:] = times
        var_out = rootgrp_out.createVariable(var, 'f8', ('time', 'nx', 'ny'))

        for it,file in enumerate(files):
            t0 = file[:-3]
            print('file: ', file)
            fullpath_in = os.path.join(path_fields, file)
            rootgrp_in = nc.Dataset(fullpath_in, 'r')
            data = rootgrp_in.groups['fields'].variables[var][:,:,:]
            var_out[it, :, :] = data[:, :, k0]

        rootgrp_out.close()




    return



def convert_file_for_varlist(var_list, files, path_fields, path_out, k_min, k_max):
    for file in files:
        t0 = file[:-3]
        print('')

        # (2) read in original fields file
        fullpath_in = os.path.join(path_fields, file)
        rootgrp_in = nc.Dataset(fullpath_in, 'r')
        dims = rootgrp_in.groups['fields'].dimensions
        nx = dims['nx'].size
        ny = dims['ny'].size
        nz = dims['nz'].size

        nz_new = k_max - k_min
        print('reducing 3D output fields:')
        print('>> from ' + str(nz) + ' levels to ' + str(k_max) + ' levels')

        fullpath_out = os.path.join(path_out, str(1000000 + np.int(t0)) + '_kmin' + str(k_min) + '_kmax' + str(k_max-1) + '.nc')
        if not os.path.exists(fullpath_out):
            print('t=' + str(t0) + ' not existing')
            # create file
            rootgrp_out = nc.Dataset(fullpath_out, 'w', format='NETCDF4')
            rootgrp_out.createDimension('time', None)
            rootgrp_out.createDimension('nx', nx)
            rootgrp_out.createDimension('ny', ny)
            rootgrp_out.createDimension('nz', nz_new)

            for var in var_list:
                print(var)
                data = rootgrp_in.groups['fields'].variables[var][:,:,:]
                print('var: ', data[:, :, k_min:k_max].shape, nx, ny, nz_new)
                # write_field(fullpath_out, var, data[:,:,k_min:k_max])
                var = rootgrp_out.createVariable(var, 'f8', ('nx', 'ny', 'nz'))
                var[:, :, :] = data[:, :, k_min:k_max]

            # rootgrp_out = nc.Dataset(fullpath_out, 'r+')
            time_out = rootgrp_out.createVariable('time','f8',('time',))
            time_out.long_name = 'Time'
            time_out.units = 's'
            time_out[:] = t0
            rootgrp_out.close()

        else:
            print('')
            print('file '+fullpath_out + ' already exists! ')
            print('')

        rootgrp_in.close()

    return


# _______________________________________________________
def convert_file_forall_variables(files, path_fields, path_out, k_min, k_max):
    for file in files:
        t0 = file[:-3]
        print('')

        # (2) read in original fields file
        fullpath_in = os.path.join(path_fields, file)
        rootgrp_in = nc.Dataset(fullpath_in, 'r')
        field_keys = rootgrp_in.groups['fields'].variables.keys()
        dims_keys = rootgrp_in.groups['fields'].dimensions.keys()
        dims = rootgrp_in.groups['fields'].dimensions
        nx = dims['nx'].size
        ny = dims['ny'].size
        nz = dims['nz'].size

        if k_min == k_max:
            print('reducing 3D output fields:')
            print('>> from '+str(nz) + ' levels to '+str(k_max) + ' levels')
            print('')

            # (3)
            nz_new = k_max
            fullpath_out = os.path.join(path_out, str(t0)+'_k'+str(k_max)+'.nc')
            if not os.path.exists(fullpath_out):
                print('t='+str(t0) + ' not existing')
                create_file(fullpath_out, nx, ny, nz_new)

                # (4)
                for var in field_keys:
                    print(var)
                    data = rootgrp_in.groups['fields'].variables[var][:,:,:]
                    write_field(fullpath_out, var, data[:,:,:k_max])
            else:
                print('')
                print('file '+fullpath_out + ' already exists! ')
                print('')
        else:

            nz_new = k_max - k_min
            print('nz_new: ', nz_new)
            fullpath_out = os.path.join(path_out, str(t0) + '_kmin' + str(k_min) + '_kmax' + str(k_max-1) + '.nc')
            if not os.path.exists(fullpath_out):
                print('t=' + str(t0) + ' not existing')
                create_file(fullpath_out, nx, ny, nz_new)

                for var in field_keys:
                    print(var)
                    data = rootgrp_in.groups['fields'].variables[var][:,:,:]
                    print('var: ', data[:,:,k_min:k_max].shape, nx, ny, nz_new)
                    write_field(fullpath_out, var, data[:,:,k_min:k_max])
            else:
                print('')
                print('file '+fullpath_out + ' already exists! ')
                print('')

        rootgrp_in.close()
    return
# _______________________________________________________

def create_file(fname, nx, ny, nz):
    rootgrp = nc.Dataset(fname, 'w', format='NETCDF4')
    rootgrp.createDimension('time', None)
    fieldgrp = rootgrp.createGroup('fields')
    fieldgrp.createDimension('nx', nx)
    fieldgrp.createDimension('ny', ny)
    fieldgrp.createDimension('nz', nz)

    rootgrp.close()
    return


def write_field(fname, f, data):
    rootgrp = nc.Dataset(fname, 'r+')
    fields = rootgrp.groups['fields']
    var = fields.createVariable(f, 'f8', ('nx', 'ny', 'nz'))
    var[:, :, :] = data
    rootgrp.close()


if __name__ == '__main__':
    main()
