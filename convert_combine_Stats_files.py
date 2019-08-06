import netCDF4 as nc
import argparse
import json as simplejson
import os
import numpy as np
import sys

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("t_end")
    parser.add_argument("--n_restart")
    args = parser.parse_args()
    case_name = args.casename
    path_root = args.path

    nml = simplejson.loads(open(os.path.join(path_root, case_name + '.in')).read())
    dt = nml['stats_io']['frequency']
    t_max = np.int(args.t_end)
    it_max = np.int(np.double(t_max)/dt) + 1

    if args.n_restart:
        num_restart = np.int(args.n_restart)
    else:
        num_restart = 1
    filename_out = 'Stats.' + case_name + '.all.nc'
    fullpath_out = os.path.join(path_root, 'stats', filename_out)
    root_out = nc.Dataset(fullpath_out, 'w', format='NETCDF4')
    print('')
    for n in range(num_restart+1):
        if n == 0:
            filename_in = 'Stats.'+case_name+'.nc'
        else:
            filename_in = 'Stats.' + case_name + '.Restart_' + str(n-1) + '.nc'
        print(filename_in)
        fullpath_in = os.path.join(path_root, 'stats', filename_in)
        root_in = nc.Dataset(fullpath_in, 'r')
        groups = root_in.groups.keys()
        t_in = root_in.groups['timeseries'].variables['t'][0]
        t_fi = root_in.groups['timeseries'].variables['t'][-1]
        it_in = np.int(t_in / dt)
        it_fi = np.int(t_fi / dt)
        print('    time: ' + str(t_in), t_fi, it_in, it_fi)
        print(root_in.groups.keys())


        for grp_name in groups:
            print('')
            print('------' + grp_name)
            if not grp_name in root_out.groups.keys():
                print('create group '+grp_name)
                grp_out = root_out.createGroup(grp_name)
            else:
                grp_out = root_out.groups[grp_name]
            grp_in = root_in.groups[grp_name]
            dims = grp_in.dimensions.keys()
            dim_list = []

            for dim in dims:
                print('dims: ', dim, grp_in.dimensions[dim].size)
                dim_list.append(dim)
                if not dim in grp_out.dimensions.keys():
                    print('create dimensions ' + dim)
                    if dim == 't':
                        grp_out.createDimension(dim, it_max)
                    else:
                        grp_out.createDimension(dim, grp_in.dimensions[dim].size)
            print('..', dim_list)

            var_list = grp_in.variables.keys()
            for var_name in var_list:
                print('---' + var_name)
                if not var_name in grp_out.variables.keys():
                    print('create variable ' + var_name)
                    if grp_name == 'profiles':
                        if var_name in ['z', 'z_half']:
                            var_out = grp_out.createVariable(var_name, 'f8', ('z'))
                        elif var_name in ['t']:
                            var_out = grp_out.createVariable(var_name, 'f8', ('t'))
                        else:
                            var_out = grp_out.createVariable(var_name, 'f8', ('t', 'z'))
                    elif grp_name == 'reference':
                        if var_name == 'z_half':
                            var_out = grp_out.createVariable(var_name, 'f8', ('z_full'))
                        else:
                            var_out = grp_out.createVariable(var_name, 'f8', ('z'))
                    else:
                        var_out = grp_out.createVariable(var_name, 'f8', dim_list)
                else:
                    var_out= grp_out.variables[var_name]
                print('created ' + str(var_out.shape))


                if grp_name == 'profiles':
                    print('output profiles-variable ' + var_name)
                    if var_name in ['z', 'z_half']:
                        var_in = grp_in.variables[var_name][:]
                        print('z', var_out.shape, var_in.shape)
                        var_out[:] = np.array(var_in)[:]
                    elif var_name == 't':
                        var_in = grp_in.variables[var_name][:]
                        var_out[it_in:it_fi+1] = np.array(var_in)[:]
                    else:
                        var_in = grp_in.variables[var_name][:,:]
                        var_out[it_in:it_fi+1,:] = np.array(var_in)[:,:]
                        print(var_in.shape, var_out.shape)
                elif grp_name == 'reference':
                    var_in = grp_in.variables[var_name][:]
                    print(var_in.shape, var_out.shape)
                    var_out[:] = var_in[:]
                elif grp_name == 'timeseries':
                    var_in = grp_in.variables[var_name][:]
                    print(var_in.shape, var_out.shape)
                    var_out[it_in:it_fi+1] = var_in[:]


            print('')
            print('---------------------------------------------------')
        root_in.close()




    root_out.close()

    return


# _______________________________________________________

if __name__ == '__main__':
    main()