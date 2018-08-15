import netCDF4 as nc
import argparse
import json as simplejson
import os
import numpy as np


def main():

    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path")
    args = parser.parse_args()


    case_name = args.casename
    path_in = args.path
    nml = simplejson.loads(open(os.path.join(path_in, case_name + '.in')).read())

    path_fields = os.path.join(path_in, 'fields')
    k_max = 120     # leads to about 20% reduction of size for fields of original size 200x200x150
    path_out = os.path.join(path_in, 'fields_k'+str(k_max))
    if not os.path.exists(path_out):
        os.mkdir(path_out)


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

        rootgrp_in.close()

    return


def create_file(fname, nx, ny, nz):
    rootgrp = nc.Dataset(fname, 'w', format='NETCDF4')
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
