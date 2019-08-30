import netCDF4 as nc
import argparse
import json as simplejson
import os
import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
import time



from convert_fields_smaller_k import convert_file_for_varlist_horsection

label_size = 10
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.labelsize'] = 15


def main():
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    path = set_input_output_parameters(args)
    path_out_figs = os.path.join(path, 'figs_massflux')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)

    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100
    times = [np.int(name[:-3]) for name in os.listdir(os.path.join(path, 'fields')) if name[-2:] == 'nc'
             and tmin <= np.int(name[:-3]) <= tmax]
    times.sort()
    nt = len(times)
    print('times: ' + str(times))
    print('nt:', nt)
    files = [str(t) + '.nc' for t in times]
    print(files)
    print('')

    zlevel = 1000
    k0 = np.int(zlevel/dx[0])
    print('considering mass flux at level: z='+str(zlevel) + ', k='+str(k0))
    convert_file_for_varlist_horsection(['w'], times, files, os.path.join(path, 'fields'), os.path.join(path, 'fields_merged'), k0)

    field_k0 = nc.Dataset(os.path.join(path, 'fields_merged', 'fields_allt_xy_k'+str(k0)+'.nc'), 'r')
    t_merged = field_k0.variables['time'][:]
    print('time merged field: ', t_merged, t_merged[0], t_merged[-1])
    if t_merged[0] > tmin or t_merged[-1] < tmax:
        print('removing file: ', os.path.join(path, 'fields_merged', 'fields_allt_xy_k'+str(k0)+'.nc'))
        os.remove(os.path.join(path, 'fields_merged', 'fields_allt_xy_k'+str(k0)+'.nc'))
        convert_file_for_varlist_horsection(['w'], times, files, os.path.join(path, 'fields'),
                                            os.path.join(path, 'fields_merged'), k0)
    mf = np.sum(field_k0.variables['w'], axis=0)
    field_k0.close()
    print('massflux: ', mf.shape)


    fig_name = 'massflux_k'+str(k0)+'.png'
    fig, axis = plt.subplots((1,2), figsize=(10,5))
    ax0 = axis[0]
    ax0.contourf(mf)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)



    # zlevel = 1000
    # k0 = np.int(zlevel/dz)
    # #[mf] = kg * m^-2 * s^-1
    # # [mf] = [rho] * [w] = kg * m^-3 * m * s^-1 = kg * m^-2 * s
    # try:
    #     stats_file = nc.Dataset(os.path.join(path, 'stat', 'Stats.'+case_name+'.nc'), 'r')
    # except:
    #     stats_file = nc.Dataset(os.path.join(path, 'Stats.'+case_name+'.nc'), 'r')
    # rho0 = stats_file.groups['reference'].variables['rho0_full'][k0]
    # stats_file.close()
    # massflux = np.zeros((nx,ny), dtype=np.double)
    # for it,t0 in times:
    #     fullpath_in = os.path.join(path, 'fields', str(t0)+'.nc')
    #     field = nc.Dataset(fullpath_in, 'r')
    #     w = field.groups['fields'].variables['w'][:,:,k0]
    #     field.close()
    #     massflux[:,:] += rho0*w[:,:]


    return

# ----------------------------------------------------------------------

def set_input_output_parameters(args):
    print('')
    print('--- set input parameters ---')
    path = args.path

    global case_name
    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dy, dz, dV, gw
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = np.zeros(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx[0] * dx[1] * dx[2]

    return path
# _______________________________________________________

if __name__ == '__main__':
    main()