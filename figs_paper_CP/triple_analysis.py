import netCDF4 as nc
import argparse
import json as simplejson
import os
import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
import matplotlib.patches as patch
import time


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

    path, nml, times = define_parameters(args)
    ic_arr, jc_arr, icoll, jcoll = define_geometry(nml)
    path_out_figs = os.path.join(path, 'figs_massflux')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    path_fields = os.path.join(path, 'fields')
    path_figs_out = os.path.join(path, 'CP_rim')
    if not os.path.exists(path_figs_out):
        os.mkdir()
    kmax = 1


    var_name = 'w'

    fig, axis = plt.subplots(3, 1, figsize=(12, 12))
    for it, t0 in enumerate(times):
        print('--- t0: '+str(t0)+ ' ---')
        rootgrp = nc.Dataset(os.path.join(path_fields, str(t0) + '.nc'), 'r')
        var = rootgrp.groups['fields'].variables[var_name][:, jcoll, :kmax]
        rootgrp.close()

        ax2.plot(x_arr, var[:, jc1], color=cm(count_color), label='t=' + str(file[:-3]))
        ax3.plot(x_arr, var[:, j0], color=cm(count_color), label='t=' + str(file[:-3]))
    plt.suptitle(var_name + ' (z=' + str(k0 * dx[2]) + 'm, k=' + str(k0) + ')', fontsize=21)
    plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.95, hspace=0.26)
    plt.savefig(os.path.join(path_out, 'xz_plane_' + var_name + '_z' + str(np.int(k0 * dx[2])) + 'm.png'))
    plt.close()
    return




# _______________________________________________________

# ----------------------------------------------------------------------

def define_parameters(args):
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

    global dt_fields, tmin, tmax
    dt_fields = nml['fields_io']['frequency']
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 0
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100
    # times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
    #          and np.int(name[:-3]) >= tmin and np.int(name[:-3]) <= tmax]
    # times.sort()
    times = np.arange(tmin, tmax+dt_fields, dt_fields)
    print('times', times)
    print('tmin, tmax', tmin, tmax)


    return path, nml, times
# _______________________
# _______________________________________________________

def define_geometry(nml):
    global rstar

    '''--- define geometry ---'''
    if case_name == 'ColdPoolDry_double_2D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        isep = 4 * irstar
        ic1 = np.int(nx / 3)
        ic2 = ic1 + isep
        jc1 = np.int(ny / 2)
        jc2 = jc1
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name[:21] == 'ColdPoolDry_single_3D':
        rstar = nml['init']['r']
        try:
            print('(ic,jc) from nml')
            ic = np.int(nml['init']['ic'])
            jc = np.int(nml['init']['jc'])
        except:
            print('(ic,jc) NOT from nml')
            ic = np.int(nx / 2)
            jc = np.int(ny / 2)
        ic_arr = [ic]
        jc_arr = [jc]
    elif case_name == 'ColdPoolDry_double_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        isep = 4 * irstar
        jsep = 0
        ic1 = np.int(np.round((nx + 2 * gw) / 3)) - gw
        jc1 = np.int(np.round((ny + 2 * gw) / 2)) - gw
        ic2 = ic1 + isep
        jc2 = jc1 + jsep
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_triple_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        try:
            d = nml['init']['d']
        except:
            d = 1e3
        i_d = np.int(np.round(d / dx[0]))
        idhalf = np.int(np.round(i_d / 2))
        a = np.int(np.round(i_d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        r_int = np.int(np.sqrt(3.) / 6 * i_d)  # radius of inscribed circle
        # point of 3-CP collision (ic, jc)
        ic = np.int(np.round(nx / 2))
        jc = np.int(np.round(ny / 2))
        ic1 = ic - r_int
        ic2 = ic1
        ic3 = ic + (a - r_int)
        jc1 = jc - idhalf
        jc2 = jc + idhalf
        jc3 = jc
        ic_arr = np.asarray([ic1, ic2, ic3])
        jc_arr = np.asarray([jc1, jc2, jc3])
        print(ic1, ic2, ic3)
        print(nx, ny, i_d, idhalf)
        print(rstar, r_int, ic)
        print('')

    return ic_arr, jc_arr, ic, jc




if __name__ == '__main__':
    main()