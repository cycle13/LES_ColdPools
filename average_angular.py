import numpy as np
#import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os
import sys

def main():

    ''' set paths & parameters '''
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    times, nml = set_input_parameters(args)
    ic_arr, jc_arr = define_geometry(nml)

    r_field = np.zeros((nx, ny), dtype=np.int)

    ic = np.int(nx/2)
    jc = np.int(ny/2)
    irange = np.minimum(nx-ic, ic)
    jrange = np.minimum(ny-jc, jc)
    rmax = np.int(np.ceil(np.sqrt(irange**2 + jrange**2)))

    # compute radius
    for i in range(irange):
        for j in range(jrange):
            r_field[ic+i,jc+j] = np.round(np.sqrt(i**2+j**2))
            r_field[ic-i,jc+j] = r_field[ic+i,jc+j]
            r_field[ic-i,jc-j] = r_field[ic+i,jc+j]
            r_field[ic+i,jc-j] = r_field[ic+i,jc+j]

    t0 = tmin
    k0 = 0
    var_list = ['w', 's', 'phi']
    fullpath_in = os.path.join(path, 'fields', str(t0) + '.nc')
    #rootgrp = nc.Dataset(fullpath_in, 'r')
    #for var_name in var_list:
    #    var = rootgrp.groups['fields'].variables[var_name][:,:,:]
    #    var_av = compute_average_var(var[:,:,k0], rmax, r_field)
    #rootgrp.close()
    data_dict = read_in_vars(fullpath_in, var_list)
    for var_name in var_list:
        var_av = compute_average_var(data_dict[var_name][:,:,k0], rmax, r_field)

    ## test field without reading in data
    #var1 = np.ones((nx, ny))  # should be from 3D LES fields, read in
    #var1_av = compute_average_var(var1, rmax, r_field)

    return


# _______________________________

def read_in_vars(fullpath_in, var_list):

    var_dict = {}

    rootgrp = nc.Dataset(fullpath_in, 'r')
    for var_name in var_list:
        var = rootgrp.groups['fields'].variables[var_name]
        data = var[:, :, :]
        var_dict[var_name] = data
    rootgrp.close()

    return var_dict


# _______________________________

def compute_average_var(var1, rmax, r_field):

    var1_av = np.zeros(rmax, dtype=np.double)
    count = np.zeros(rmax, dtype=np.int)
    for i in range(nx):
        for j in range(ny):
            r = r_field[i, j]
            count[r] += 1
            var1_av[r] += var1[i, j]

    for r in range(rmax):
        if count[r] > 0:
            var1_av[r] /= count[r]

    return var1_av
# _______________________________

# _______________________________
def set_input_parameters(args):
    print ''' setting parameters '''
    global path, path_fields #path_out_data, path_out_figs
    path = args.path
    path_fields = os.path.join(path, 'fields')

    # path_out_data = os.path.join(path, 'data_analysis')
    # if not os.path.exists(path_out_data):
    #     os.mkdir(path_out_data)
    # path_out_figs = os.path.join(path, 'figs_minmax')
    # if not os.path.exists(path_out_figs):
    #     os.mkdir(path_out_figs)
    print ''
    print 'paths:'
    #print path_out_data
    #print path_out_figs
    print path_fields
    print ''

    global case_name
    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx
    dx = np.ndarray(3, dtype=np.int)
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']


    global ic_arr, jc_arr
    try:
        print('(ic,jc) from nml')
        ic = nml['init']['ic']
        jc = nml['init']['jc']
    except:
        print('(ic,jc) NOT from nml')
        if case_name == 'ColdPoolDry_single_3D':
            ic = np.int(nx/2)
            jc = np.int(ny/2)
            ic_arr = [ic]
            jc_arr = [jc]
        else:
            print('ic, jc not defined')

    global kmax
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = nx

    global tmin, tmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100

    ''' file range '''
    # times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
    #          and np.int(name[:-3]) >= tmin and np.int(name[:-3]) <= tmax]
    times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
             and tmin <= np.int(name[:-3]) <= tmax]
    times.sort()
    # files = [str(t) + '.nc' for t in times]

    print('')
    # print('files', files, len(files))
    print('times', times)
    print('')
    print('kmax ', kmax, 'nx ', nx)
    print('')

    return times, nml
# _______________________________


# def define_geometry(nml, files):
def define_geometry(nml):
    a = nml['grid']['nx']
    '''--- define geometry ---'''
    global rstar
    if case_name == 'ColdPoolDry_double_2D':
        rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        # zstar = nml['init']['h']
        isep = 4 * irstar
        ic1 = np.int(nx / 3)
        ic2 = ic1 + isep
        jc1 = np.int(ny / 2)
        jc2 = jc1
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_single_3D':
        rstar = nml['init']['r']
        # irstar = np.int(np.round(rstar / dx))
        # zstar = nml['init']['h']
        dTh = nml['init']['dTh']
        ic = np.int(nx / 2)
        jc = np.int(ny / 2)
        # xc = Gr.x_half[ic + Gr.dims.gw]  # center of cold-pool
        # yc = Gr.y_half[jc + Gr.dims.gw]  # center of cold-pool
        ic_arr = [ic]
        jc_arr = [jc]
    elif case_name == 'ColdPoolDry_double_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        # zstar = nml['init']['h']
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
        irstar = np.int(np.round(rstar / dx))
        d = np.int(np.round(ny / 2))
        dhalf = np.int(np.round(ny / 4))
        a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        ic1 = np.int(np.round(a / 2))  # + gw
        ic2 = ic1
        ic3 = ic1 + np.int(np.round(a))
        jc1 = np.int(np.round(d / 2))  # + gw
        jc2 = jc1 + d
        jc3 = jc1 + np.int(np.round(d / 2))
        ic_arr = [ic1, ic2, ic3]
        jc_arr = [jc1, jc2, jc3]

        isep = dhalf

    # ''' --- auxiliary arrays (since no Grid.pyx) ---'''
    # global nx_, ny_, nz_
    # # test file:
    # var = read_in_netcdf_fields('u', os.path.join(path_fields, files[0]))
    # [nx_, ny_, nz_] = var.shape
    #
    # x_half = np.empty((nx_), dtype=np.double, order='c')
    # y_half = np.empty((ny_), dtype=np.double, order='c')
    # z_half = np.empty((nz_), dtype=np.double, order='c')
    # count = 0
    # for i in xrange(nx_):
    #     x_half[count] = (i + 0.5) * dx
    #     count += 1
    # count = 0
    # for j in xrange(ny_):
    #     y_half[count] = (j + 0.5) * dy
    #     count += 1
    # count = 0
    # for i in xrange(nz_):
    #     z_half[count] = (i + 0.5) * dz
    #     count += 1
    #
    # return x_half, y_half, z_half

    return ic_arr, jc_arr

# _______________________________


if __name__ == '__main__':
    main()
