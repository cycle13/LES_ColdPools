import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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
    # parser.add_argument("--t0")
    parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    times, nml = set_input_parameters(args)
    ic_arr, jc_arr = define_geometry(nml)
    ic1 = ic_arr[0]
    jc1 = jc_arr[0]


    t0 = 100
    for t0 in times:
        fullpath_in = os.path.join(path_fields, str(t0)+'.nc')
        print fullpath_in
        rootgrp = nc.Dataset(fullpath_in, 'r')
        phi = rootgrp.groups['fields'].variables['phi'][:,:,:]
        s = rootgrp.groups['fields'].variables['s'][:,:,:]
        rootgrp.close()

        kmax = 80

        fig_name = 'test_fig_t'+str(t0)+'.png'
        fig, axes = plt.subplots(2, 3, figsize=(12, 10))
        ax = axes[0,0]
        cf = ax.contourf(s[:,:,0].T)
        plt.colorbar(cf, ax=ax, shrink=0.8)
        plt.plot([ic1, ic1], [0, ny], 'k-')
        plt.plot([0, nx], [jc_arr[2], jc_arr[2]], 'k--')
        ax = axes[0,1]
        cf = ax.contourf(s[ic1, :, :].T)
        plt.colorbar(cf, ax=ax, shrink=0.8)
        ax = axes[0,2]
        cf = ax.contourf(phi[:, jc1, :kmax].T)
        plt.colorbar(cf, ax=ax, shrink=0.8)
        ax.set_title('j=jc1')
        ax = axes[1,0]
        cf = ax.contourf(phi[:, :, 0].T)
        plt.colorbar(cf, ax=ax, shrink=0.8)
        ax.set_title('z='+str(0*dx[2])+'m')
        ax = axes[1,1]
        k0 = 10
        cf = ax.contourf(phi[:, :, k0].T)
        ax.set_title('z='+str(k0*dx[2])+'m')
        plt.colorbar(cf, ax=ax, shrink=0.8)
        ax = axes[1,2]
        cf = ax.contourf(phi[ic1, :, :kmax].T)
        plt.colorbar(cf, ax=ax, shrink=0.8)
        ax.set_title('i=ic1')
        plt.savefig(os.path.join(path_out_figs, fig_name))
        plt.close()


        fig_name = 'test_fig2_t'+str(t0)+'.png'
        fig, axes = plt.subplots(1, 5, figsize=(25, 5))
        krange = np.arange(0, 5)
        for k0 in krange:
            ax = axes[k0]
            cf = ax.contourf(phi[:, :, k0*5].T)
            plt.colorbar(cf, ax=ax, shrink=0.8)
            ax.set_title('z='+str(k0*5*dx[2])+'m')
        plt.savefig(os.path.join(path_out_figs, fig_name))
        plt.close()

    return

# _______________________________
# _______________________________
def set_input_parameters(args):
    print ''' setting parameters '''

    global path, path_fields, path_out_figs
    path = args.path

    path_fields = os.path.join(path, 'fields')
    path_out_figs = os.path.join(path, 'figs_phi')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print ''
    print 'paths:'
    print path_fields
    print path_out_figs
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

    global kmax
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = nz

    global tmin, tmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100

    # ''' time range '''
    times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
             and tmin <= np.int(name[:-3]) <= tmax]
    times.sort()

    print('')
    print('times', times)
    print('')
    print('tmin, tmax', tmin, tmax)
    print('kmax ', kmax, 'nx ', nx)
    print('')

    return times, nml

# _______________________________

def define_geometry(nml):
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
        # dTh = nml['init']['dTh']
        try:
            print('(ic,jc) from nml')
            ic = nml['init']['ic']
            jc = nml['init']['jc']
        except:
            print('(ic,jc) NOT from nml')
            ic = np.int(nx / 2)
            jc = np.int(ny / 2)
        ic_arr = [ic]
        jc_arr = [jc]
        # xc = Gr.x_half[ic + Gr.dims.gw]  # center of cold-pool
        # yc = Gr.y_half[jc + Gr.dims.gw]  # center of cold-pool
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
        irstar = np.int(np.round(rstar / dx[0]))
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
    return ic_arr, jc_arr

# _______________________________
# _______________________________


if __name__ == '__main__':
    main()