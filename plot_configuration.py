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
import sys

plt.rcParams['lines.linewidth'] = 3

def main():
    print('PLOT CONFIGURATION')
    # (A) domain maximum
    #
    # (B) cross-sections
    #   determine max v in yz-crosssections
    #   determine max u in xz-crosssections
    #

    ''' set paths & parameters '''
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path_root")
    parser.add_argument("dTh", type=int)
    parser.add_argument("--zparams", nargs='+', type=int)
    parser.add_argument('--rparams', nargs='+', type=int)
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()
    #
    dTh, z_params, r_params, ic_arr, jc_arr, marg, times = set_input_parameters(args)

    # if len(z_params) != len(r_params):
    #     print('wrong number of parameters! ')
    #     sys.exit()
    #
    # # plot contourfigure of entropy and centers of CPs (ic_arr, jc_arr)
    # test_configuration(dTh, z_params, r_params, ic_arr, jc_arr, marg, path_root, path_out_figs)

    # plot Gaussian envelops of configurations
    dTh_params = [2, 3, 4]
    plot_envelope(dTh_params, ic_arr, path_out_figs)

    return



# ----------------------------------
def plot_envelope(dTh_params, ic_arr, path_out_figs):
    ''' geometry '''

    # global x_half, y_half, z_half
    x_half = np.empty((nx + 2 * gw), dtype=np.double, order='c')
    z_half = np.empty((nz + 2 * gw), dtype=np.double, order='c')
    count = 0
    for i in xrange(-gw, nx + gw):
        x_half[count] = (i + 0.5) * dx[0]
        count += 1
    count = 0
    for i in xrange(-gw, nz + gw):
        z_half[count] = (i + 0.5) * dx[2]
        count += 1
    # ''' geometry '''
    ic = np.int(ic_arr[0])
    xc = x_half[ic]

    fig, (axes) = plt.subplots(1, 4, figsize=(20,5))
    ax0 = axes[0]
    ax1 = axes[1]
    ax2 = axes[2]
    ax3 = axes[3]
    # dTh_params = [2, 3, 4]

    z_params = [1225, 1000, 870]
    r_params = z_params
    n_params = len(dTh_params)
    for istar in range(n_params):
        zstar = z_params[istar]
        rstar = r_params[istar]
        irstar = np.int(rstar/dx[0])
        x_half_aux = (xc+rstar) * np.ones(shape=x_half.shape)
        for i in range(ic-irstar, ic+irstar):
            x_half_aux[i] = x_half[i]
        z_max = zstar * (np.cos((x_half_aux-xc) / rstar * np.pi / 2) ** 2)
        imin = ic - irstar - 10
        imax = ic + irstar + 10
        ax0.plot(x_half[imin:imax], z_max[imin:imax], label='dTh'+str(dTh_params[istar]))
    ax0.set_title('r = 1km')


    for n,dTh in enumerate(dTh_params):
        print('n', n, dTh)
        if dTh == 4:
            z_params = [1730, 870, 430]
            # r_params=[430, 870, 1730]
        elif dTh == 3:
            z_params = [4000, 2000, 1500, 1000, 670, 500, 250]
        elif dTh == 2:
            z_params = [2450, 1225, 815]
        elif dTh == 1:
            z_params = [3465, 1730, 1155]
        r_params = z_params[::-1]
        print('dTh: ', dTh)
        print('z*: ', z_params)
        print('r*: ', r_params)

        n_params = len(z_params)
        for istar in range(n_params):
            zstar = z_params[istar]
            rstar = r_params[istar]
            irstar = np.int(rstar / dx[0])
            x_half_aux = (xc + rstar) * np.ones(shape=x_half.shape)
            for i in range(ic - irstar, ic + irstar):
                x_half_aux[i] = x_half[i]
            z_max = zstar * (np.cos((x_half_aux - xc) / rstar * np.pi / 2) ** 2)
            ax = axes[n + 1]
            imin = ic - irstar - 10
            imax = ic + irstar + 10
            ax.plot(x_half[imin:imax], z_max[imin:imax], label='z*'+str(zstar)+', r*'+str(rstar))
            ax.set_title('dTh='+str(dTh))


    irstar = 5000 / dx[0]
    imin = ic - irstar
    imax = ic + irstar
    xmin = x_half[imin]
    xmax = x_half[imax]
    for ax in axes:
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(0, 4000)
        ax.legend(fontsize=9)
    plt.tight_layout
    plt.savefig(os.path.join(path_out_figs, 'envelope.png'))
    plt.close(fig)


    return
# ----------------------------------


def test_configuration(dTh, z_params, r_params, ic_arr, jc_arr, marg, path_root, path_out_figs):
    print('test configuration: dTh='+str(dTh))
    t0 = 0
    ng = len(z_params)
    # path_out = os.path.join(path_out_figs, 'config')
    # print(path_out)
    # if not os.path.exists(path_out):
    #     os.mkdir(path_out)

    for istar in range(ng):
        zstar = z_params[istar]
        rstar = r_params[istar]
        kstar = zstar / dx[2]
        print('z*, r*: ', zstar, rstar)

        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        path_fields = os.path.join(path_root, id, 'fields')
        var = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))
        nx, ny, nz = var.shape
        imin = 120
        imax = nx - imin
        jmin = imin
        jmax = ny - jmin
        kmax = 50
        x_ = dx[0]*np.arange(imin, imax)
        y_ = dx[1]*np.arange(jmin, jmax)
        z_ = dx[2]*np.arange(0, kmax)
        # X, Y = np.meshgrid(x_, y_, indexing='ij')

        plt.figure()
        plt.contourf(x_, y_, var[imin:imax, jmin:jmax, 0].T, origin='lower')
        plt.colorbar()
        plt.plot([x_[ic_arr[0]-imin], x_[ic_arr[0]-imin]], [y_[0], y_[-1]], 'k')
        plt.plot([x_[ic_arr[0]-imin]-rstar, x_[ic_arr[0]-imin]-rstar], [y_[0], y_[-1]], 'k--')
        plt.plot([x_[ic_arr[0]-imin]-rstar-marg, x_[ic_arr[0]-imin]-rstar-marg], [y_[0], y_[-1]], 'k:')
        plt.plot([x_[ic_arr[0]-imin]+rstar, x_[ic_arr[0]-imin]+rstar], [y_[0], y_[-1]], 'k--')
        plt.plot([x_[ic_arr[0]-imin]+rstar+marg, x_[ic_arr[0]-imin]+rstar+marg], [y_[0], y_[-1]], 'k:')
        plt.plot([x_[0], x_[-1]], [y_[jc_arr[0]-jmin], y_[jc_arr[0]-jmin]], 'k')
        plt.plot([x_[0], x_[-1]], [y_[jc_arr[0]-jmin]-rstar, y_[jc_arr[0]-jmin]-rstar], 'k--')
        plt.plot([x_[0], x_[-1]], [y_[jc_arr[0]-jmin]+rstar, y_[jc_arr[0]-jmin]+rstar], 'k--')
        plt.tight_layout()
        plt.savefig(os.path.join(path_out_figs, 'config_' + id + '.png'))
        plt.close()


        plt.figure()
        plt.contourf(x_, z_, var[imin:imax, jc_arr[0], :kmax].T)
        plt.colorbar()
        plt.plot([x_[0],x_[-1]],[z_[kstar], z_[kstar]], 'w', linewidth=2)
        plt.plot([x_[0],x_[-1]],[z_[kstar]+marg, z_[kstar]+marg], 'w', linewidth=1)
        plt.plot([x_[ic_arr[0] - imin], x_[ic_arr[0] - imin]], [z_[0], z_[-1]], 'k', linewidth=2)
        plt.plot([x_[ic_arr[0] - imin] - rstar, x_[ic_arr[0] - imin] - rstar], [z_[0], z_[-1]], 'k--')
        plt.plot([x_[ic_arr[0] - imin] - rstar - marg, x_[ic_arr[0] - imin] - rstar - marg], [z_[0], z_[-1]], 'k:')
        plt.plot([x_[ic_arr[0] - imin] + rstar, x_[ic_arr[0] - imin] + rstar], [z_[0], z_[-1]], 'k--')
        plt.plot([x_[ic_arr[0] - imin] + rstar + marg, x_[ic_arr[0] - imin] + rstar + marg], [z_[0], z_[-1]], 'k:')
        # plt.plot([ic_arr[0] - imin, ic_arr[0] - imin], [0, 30], 'k')
        # plt.tight_layout()
        plt.savefig(os.path.join(path_out_figs, 'config_xz_' + id + '.png'))
        plt.close()
    return


# ----------------------------------
# ----------------------------------
def set_input_parameters(args):

    global path_root, path_out_data, path_out_figs
    path_root = args.path_root
    path_out_data = os.path.join(path_root, 'data_analysis')
    if not os.path.exists(path_out_data):
        os.mkdir(path_out_data)
    path_out_figs = os.path.join(path_root, 'figs_config')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print('path root: ', path_root)
    print('path out data: ', path_out_data)
    print('path out figs: ', path_out_figs)


    dTh = args.dTh
    if args.zparams:
        z_params = args.zparams
    else:
        if dTh == 4:
            z_params = [1730, 870, 430]
        elif dTh == 3:
            z_params = [4000, 2000, 1500, 1000, 670, 500, 250]
        elif dTh == 2:
            z_params = [2450, 1225, 815]
        elif dTh == 1:
            z_params = [3465, 1730, 1155]
    if args.rparams:
        r_params = args.rparams
    else:
        r_params = z_params[::-1]
    print('z*: ', z_params)
    print('r*: ', r_params)

    global case_name
    case_name = args.casename
    id0 = 'dTh' + str(dTh) + '_z' + str(z_params[0]) + '_r' + str(r_params[0])
    nml = simplejson.loads(open(os.path.join(path_root, id0, case_name + '.in')).read())
    global nx, ny, nz, dx, gw
    dx = np.arange(3)
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']


    try:
        print('(ic,jc) from nml')
        ic = nml['init']['ic']
        jc = nml['init']['jc']
    except:
        print('(ic,jc) NOT from nml')
        if case_name == 'ColdPoolDry_single_3D':
            ic = np.int(nx/2)
            jc = np.int(ny/2)
            ic_arr = np.zeros(1)
            jc_arr = np.zeros(1)
            ic_arr[0] = ic
            jc_arr[0] = jc
        else:
            print('ic, jc not defined')
    try:
        marg = nml['init']['marg']
        print('marg from nml')
    except:
        marg = 500.
        print('marg NOT from nml')

    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100
    times = np.arange(tmin, tmax + 100, 100)

    print('times', times)
    print('')

    return dTh, z_params, r_params, ic_arr, jc_arr, marg, times




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