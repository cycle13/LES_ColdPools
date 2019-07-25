import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 12


def main():
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path1")
    parser.add_argument("path2")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    kmax, krange = set_input_parameters(args)
    dt_fields, times = define_geometry()
    nt = len(times)
    # ncases = len(id_list)
    path_list = [path1, path2]
    id_list = ['neutral', 'stable']
    ncases = len(path_list)
    print('')

    var_name = 'w'
    filename_in = os.path.join('data_analysis', 'stats_crosssections.nc')
    times_array = np.ndarray((ncases, 2), dtype=np.double)
    min_array = np.ndarray((ncases, 3, kmax), dtype=np.double)
    max_array = np.ndarray((ncases, 3, kmax), dtype=np.double)
    #     print('--- '+ID)
    for ip, path in enumerate(path_list):
        print(ip, path)
        print(os.path.join(path, filename_in))
        fullpath_in = os.path.join(path, filename_in)
        file = nc.Dataset(fullpath_in, 'r', format='NETCDF4')
        min_array[ip, 0, :] = file.groups[var_name].variables['min_single'][:kmax]
        min_array[ip, 1, :] = file.groups[var_name].variables['min_double'][:kmax]
        min_array[ip, 2, :] = file.groups[var_name].variables['min_triple'][:kmax]
        max_array[ip, 0, :] = file.groups[var_name].variables['max_single'][:kmax]
        max_array[ip, 1, :] = file.groups[var_name].variables['max_double'][:kmax]
        max_array[ip, 2, :] = file.groups[var_name].variables['max_triple'][:kmax]
        times_array[ip, 0] = file.groups[var_name].variables['time_single'][:]
        times_array[ip, 1] = file.groups[var_name].variables['time_double'][:]
        file.close()


    fig_name = 'minmax_levels_' + var_name + '_' + case_name + '.png'
    zrange = np.arange(kmax) * dx[2]
    fig, axis = plt.subplots(1, 3, figsize=(14, 10))
    ax0 = axis[0]
    ax1 = axis[1]
    ax2 = axis[2]
    for ip in range(ncases):
        al = np.double(ip+1)/(ncases+1)
        ax0.plot(max_array[ip,0,:kmax], zrange, color='b', alpha=al, label='single, '+id_list[ip])
        ax0.plot(max_array[ip,1,:kmax], zrange, '-', color='g', alpha=al, label='double')
        ax0.plot(max_array[ip,2,:kmax], zrange, '-', color='r', alpha=al, label='triple')
        ax1.plot(min_array[ip,0,:kmax], zrange, color='b', alpha=al, label='single, '+id_list[ip])
        ax1.plot(min_array[ip,1,:kmax], zrange, '-', color='g', alpha=al, label='double')
        ax1.plot(min_array[ip,2,:kmax], zrange, '-', color='r', alpha=al, label='triple')
        ax2.plot([1, 2], times_array[ip, :], '-o', color=str(1-al), label='times, ' +id_list[ip])
        ax2.fill_between(np.arange(1, 3), times_array[ip, 0], times_array[ip, 1],
                          facecolor=str(al), alpha=0.2, linewidth=0)
    for ax in axis[:2]:
        ax.set_xlabel(var_name)
        ax.set_ylabel('height  [m]')
        ax.grid()
    ax2.set_xlabel('collision (#CPs)')
    ax2.set_ylabel('time [s]')
    ax2.grid()
    ax2.set_xlim(0.9,2.1)
    ax0.set_title('maxima')
    ax1.set_title('minima')
    ax0.legend(loc='upper left', bbox_to_anchor=(-0.1, -0.1),
               fancybox=False, shadow=False, ncol=2, fontsize=9)
    ax2.legend(loc='upper left', bbox_to_anchor=(0.1, -0.1),
               fancybox=False, shadow=False, ncol=1, fontsize=9)
    plt.suptitle('min/max for ' + var_name, fontsize=21)
    plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.9, hspace=0.4)
    plt.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)


    fig_name = 'times_' + var_name + '_' + case_name + '.png'
    zrange = np.arange(kmax) * dx[2]
    fig, axis = plt.subplots(1, 2, figsize=(9, 10))
    ax0 = axis[0]
    ax1 = axis[1]
    for ip in range(ncases):
        al = np.double(ip + 1) / (ncases + 1)
        ax0.plot([1,2],times_array[ip, :], '-o', color='b', alpha=al, label='single, ' + str(ip))
        ax1.plot(times_array[ip, 1], color='g', alpha=al, label='double')
    for ax in axis:
        ax.set_xlabel('collision (#CPs)')
        ax.set_ylabel('time [s]')
        ax.grid()
    # ax0.set_title('maxima')
    # ax1.set_title('minima')
    ax0.legend(loc='upper left', bbox_to_anchor=(-0.1, -0.1),
               fancybox=False, shadow=False, ncol=3, fontsize=9)
    plt.suptitle('min/max for ' + var_name, fontsize=21)
    plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.9, hspace=0.4)
    plt.savefig(os.path.join(path_out_figs, fig_name))
    plt.close()

    return



def set_input_parameters(args):
    global path1, path2, path_out_figs
    path1 = args.path1
    path2 = args.path2
    path_out_figs = os.path.join(path1, 'figs_crosssections_neutral_vs_stable')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print('')
    print('path out: ')
    print(path_out_figs)
    print('')

    global case_name
    case_name = args.casename
    print('')
    print('casename: ' + case_name)
    if case_name == 'ColdPoolDry_triple_3D':
        dTh = 3
        zstar = 2000
        rstar = 2000
        # d_range = [10, 15, 20]
        d_range = [10]
        # id_list = []
        # for dstar in d_range:
        #     id_list.append('dTh' + str(dTh) +'_z' + str(zstar) + '_r' + str(rstar) + '_d' + str(dstar) + 'km')
    # print(id_list)
    print('')

    ''' determine time range '''
    global tmin, tmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = np.int(0)
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = np.int(7200)
    print('tmin: '+str(tmin))
    print('tmax: '+str(tmax))
    print('')

    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = 5
    krange = np.arange(kmax + 1)
    print('krange: ', krange)
    print ''

    return kmax, krange
    # return kmax, krange, id_list




def define_geometry():
    print 'define geometry'
    global nx, ny, nz, dx, gw
    nml = simplejson.loads(open(os.path.join(path1, case_name + '.in')).read())
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = np.ndarray(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']
    global dt_fields
    dt_fields = np.int(nml['fields_io']['frequency'])
    times = np.arange(tmin, tmax+dt_fields, dt_fields)
    print('times: '+str(times))

    return dt_fields, times
# ----------------------------------


if __name__ == '__main__':
    main()



