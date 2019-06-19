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

label_size = 15
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['grid.linewidth'] = 20
# plt.rcParams['xtick.major.size'] = 8.5
# plt.rcParams['xtick.minor.size'] = 5
# plt.rcParams['ytick.major.size'] = 8.5
# plt.rcParams['ytick.minor.size'] = 5
# plt.rcParams['xtick.major.width'] = 2
# plt.rcParams['xtick.minor.width'] = 1.5
# plt.rcParams['ytick.major.width'] = 2
# plt.rcParams['ytick.minor.width'] = 1.5
plt.rcParams['pdf.fonttype'] = 42         # Output Type 3 (Type3) or Type 42 (TrueType)

def main():
    print('PLOT CONFIGURATION')
    print('')
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
    # parser.add_argument("--tmin")
    # parser.add_argument("--tmax")
    args = parser.parse_args()

    colorlist = ['seagreen', 'navy', 'firebrick', 'grey', 'k', 'r', 'r']

    dTh = args.dTh
    z_params, r_params = set_geom_parameters(dTh)
    if args.zparams:
        z_params = args.zparams
    if args.rparams:
        r_params = args.rparams
    ic_arr, jc_arr, marg = set_input_parameters(args, dTh, z_params, r_params)
    print('')
    print('dTh: ', dTh)
    print('z*: ', z_params)
    print('r*: ', r_params)
    print('')

    # # if len(z_params) != len(r_params):
    # #     print('wrong number of parameters! ')
    # #     sys.exit()
    # #
    # # # plot contourfigure of entropy and centers of CPs (ic_arr, jc_arr)
    # # test_configuration(dTh, z_params, r_params, ic_arr, jc_arr, marg, path_root, path_out_figs)
    # #
    # # plot Gaussian envelops of configurations
    # dTh_params = [2, 3, 4]
    # plot_envelope_run1(dTh_params, z_params, r_params, ic_arr, colorlist, path_out_figs)

    dTh = 5
    # id_ref = 'dTh3_z1000_r1000'
    # path_ref = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/dTh3_z1000_r1000'
    dTh_ref = 3
    r_star_ref = 1000
    z_star_ref = 1000
    plot_envelope_run5(dTh, z_params, r_params, dTh_ref, z_star_ref, r_star_ref, ic_arr,
                       colorlist, path_out_figs)

    return



# ----------------------------------
def plot_envelope(dTh_params, z_params, r_params, ic_arr, colorlist, path_out_figs):
    ''' geometry '''

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

    zstar_max = 0
    rstar_max = 0
    fig, (axes) = plt.subplots(1, 4, figsize=(20,5))

    # (1) r = 1km, varying dTh
    ax0 = axes[0]
    for ith, dTh in enumerate(dTh_params):
        # z_params, r_params = set_geom_parameters(dTh)
        istar = 0
        while (z_params[istar] < 1200):
            istar+=1
        zstar = z_params[istar-1]
        rstar = r_params[istar-1]
        irstar = np.int(rstar/dx[0])
        zstar_max = np.maximum(zstar_max, zstar)
        rstar_max = np.maximum(rstar_max, rstar)
        x_half_aux = (xc+rstar) * np.ones(shape=x_half.shape)
        for i in range(ic-irstar, ic+irstar):
            x_half_aux[i] = x_half[i]
        z_max = zstar * (np.cos((x_half_aux-xc) / rstar * np.pi / 2) ** 2)
        imin = ic - irstar - 10
        imax = ic + irstar + 10
        ax0.plot(x_half[imin:imax], z_max[imin:imax], color=colorlist[ith],
                 label='dTh'+str(dTh)+': z*'+str(zstar)+', r*'+str(rstar))
    ax0.set_title('r = 1km')
    ax0.set_ylabel('z   [m]')


    # (2) for all dTh, varying z*, r*
    for n,dTh in enumerate(dTh_params):
        z_params, r_params = set_geom_parameters(dTh)
        print('dTh: ', dTh)
        print('z*: ', z_params)
        print('r*: ', r_params)

        n_params = len(z_params)
        for istar in range(n_params):
            zstar = z_params[istar]
            rstar = r_params[istar]
            irstar = np.int(rstar / dx[0])
            zstar_max = np.maximum(zstar_max, zstar)
            rstar_max = np.maximum(rstar_max, rstar)
            x_half_aux = (xc + rstar) * np.ones(shape=x_half.shape)
            for i in range(ic - irstar, ic + irstar):
                x_half_aux[i] = x_half[i]
            z_max = zstar * (np.cos((x_half_aux - xc) / rstar * np.pi / 2) ** 2)
            ax = axes[n + 1]
            imin = ic - irstar - 10
            imax = ic + irstar + 10
            ax.plot(x_half[imin:imax], z_max[imin:imax], color=colorlist[istar],
                    label='z*'+str(zstar)+', r*'+str(rstar))
            ax.set_title('dTh='+str(dTh))


    irstar = 5000 / dx[0]
    imin = ic - irstar
    imax = ic + irstar
    xmin = x_half[imin]
    xmax = x_half[imax]
    xmin = x_half[ic] - rstar_max - 500
    xmax = x_half[ic] + rstar_max + 500
    for ax in axes:
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(0, zstar_max + 200)
        ax.set_xlabel('x  [m]')
        ax.legend()
    plt.tight_layout
    plt.savefig(os.path.join(path_out_figs, 'envelope.png'))
    plt.close(fig)

    return



def plot_envelope_run5(dTh, z_params, r_params, dTh_ref, zstar_ref, rstar_ref, ic_arr,
                       colorlist, path_out_figs):
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

    zstar_max = 0
    rstar_max = 0
    n_params = len(r_params)


    # (1) varying z*, r*
    # print('dTh: ', dTh)
    # print('z*: ', z_params)
    # print('r*: ', r_params)

    fig, (axes) = plt.subplots(1, 2, figsize=(16,5))
    for istar in range(n_params):
        zstar = z_params[0]
        rstar = r_params[istar]
        irstar = np.int(rstar / dx[0])
        zstar_max = np.maximum(zstar_max, zstar)
        rstar_max = np.maximum(rstar_max, rstar)
        x_half_aux = (xc + rstar) * np.ones(shape=x_half.shape)
        for i in range(ic - irstar, ic + irstar):
            x_half_aux[i] = x_half[i]
        z_max = zstar * (np.cos((x_half_aux - xc) / rstar * np.pi / 2) ** 2)
        imin = ic - irstar - 10
        imax = ic + irstar + 10
        axes[0].plot(x_half[imin:imax], z_max[imin:imax],
                label='dTh'+str(dTh)+', z*'+str(zstar)+', r*'+str(rstar))
        axes[1].plot(x_half[imin:imax], z_max[imin:imax], color=colorlist[istar],
                label='dTh'+str(dTh)+', z*'+str(zstar)+', r*'+str(rstar))

    zstar = zstar_ref
    rstar = rstar_ref
    irstar = np.int(rstar / dx[0])
    zstar_max = np.maximum(zstar_max, zstar)
    rstar_max = np.maximum(rstar_max, rstar)
    x_half_aux = (xc + rstar) * np.ones(shape=x_half.shape)
    for i in range(ic - irstar, ic + irstar):
        x_half_aux[i] = x_half[i]
    z_max = zstar * (np.cos((x_half_aux - xc) / rstar * np.pi / 2) ** 2)
    imin = ic - irstar - 10
    imax = ic + irstar + 10
    axes[0].plot(x_half[imin:imax], z_max[imin:imax], 'k',
            label='dTh'+str(dTh)+', z*'+str(zstar_ref)+', r*'+str(rstar_ref))

    irstar = 5000 / dx[0]
    imin = ic - irstar
    imax = ic + irstar
    xmin = x_half[imin]
    xmax = x_half[imax]
    xmin = x_half[ic] - rstar_max - 500
    xmax = x_half[ic] + rstar_max + 500
    for ax in axes:
        x_ticks = ax.get_xticks() * 1e-3
        y_ticks = ax.get_yticks() * 1e-3
        ax.set_xticklabels(x_ticks)
        ax.set_yticklabels(y_ticks)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(0, zstar_max + 200)
        ax.set_title('dTh=' + str(dTh))
        ax.set_xlabel('x  [km]', fontsize=15)
        ax.set_ylabel('z  [km]', fontsize=15)
        # ax.legend()
    axes[0].legend(loc='upper right', bbox_to_anchor=(-0.2, 1.0),
                   fancybox=True, ncol=1)
    axes[1].legend(loc='upper left', bbox_to_anchor=(1.05, 1.0),
                   fancybox=True, ncol=1)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12, right=.85, left=0.2, top=0.9, wspace=0.25)
    plt.savefig(os.path.join(path_out_figs, 'envelope.png'))
    plt.close(fig)







    Lx = 6e3
    dx_ = 10
    nx_ = Lx/dx_
    x_arr = np.arange(0,Lx,dx_)
    ic = np.int(nx_ /2)
    xc = x_arr[ic]
    print nx_, ic, xc
    n_params = len(r_params)
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(16, 5))
    for istar in range(n_params):
        zstar = z_params[0]
        rstar = r_params[istar]
        irstar = np.int(np.double(rstar)/dx_)
        x_arr_ = np.array(x_arr, copy=True)
        x_arr_[:ic-irstar] = rstar + xc
        x_arr_[ic+irstar:] = rstar + xc
        zmax = zstar * np.cos( (x_arr_ - xc)/rstar * np.pi/2 )**2
        ax0.plot(x_arr-xc, zmax, label='dTh'+str(dTh)+', z*'+str(zstar)+', r*'+str(rstar))
        ax1.plot(x_arr-xc, zmax, label='dTh'+str(dTh)+', z*'+str(zstar)+', r*'+str(rstar),
                 color = colorlist[istar])
    x_arr_ = np.array(x_arr, copy=True)
    irstar = np.int(np.double(rstar_ref) / dx_)
    x_arr_[:ic - irstar] = rstar_ref + xc
    x_arr_[ic + irstar:] = rstar_ref + xc
    zmax = zstar_ref * np.cos((x_arr_ - xc) / rstar_ref * np.pi / 2) ** 2
    ax0.plot(x_arr-xc, zmax, 'k',
                 label='dTh' + str(dTh_ref) + ', z*' + str(zstar_ref) + ', r*' + str(rstar_ref))
    ax1.plot(x_arr-xc, zmax, 'k',
             label='dTh' + str(dTh_ref) + ', z*' + str(zstar_ref) + ', r*' + str(rstar_ref))

    ax0.set_title('dTh=' + str(dTh))
    ax1.set_title('dTh=' + str(dTh))
    ax0.set_xlabel('Radius r', fontsize=21)
    ax1.set_xlabel('Radius r', fontsize=21)
    ax0.set_ylabel('Height z', fontsize=21)
    ax1.set_ylabel('Height z', fontsize=21)
    ax0.set_ylim(0, zstar_max + 200)
    ax1.set_ylim(0, zstar_max + 200)


    x_ticks = [np.int(i*1e-3) for i in ax0.get_xticks()]
    # print ax0.get_yticks()
    # y_ticks = np.round(ax0.get_yticks()*1e-3,1)
    # print y_ticks
    y_ticks = [np.round(i*1e-3,1) for i in ax0.get_yticks()]
    print y_ticks
    ax0.set_xticklabels(x_ticks)
    ax1.set_xticklabels(x_ticks)
    ax0.set_yticklabels(y_ticks)
    ax1.set_yticklabels(y_ticks)

    ax0.legend(loc='upper right', bbox_to_anchor=(-0.2, 1.0),
                  fancybox=True, ncol=1)
    ax1.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0),
                   fancybox=True, ncol=1)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12, right=.85, left=0.2, top=0.9, wspace=0.25)
    plt.savefig(os.path.join(path_out_figs, 'envelope_smooth.png'))
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
def set_geom_parameters(dTh):

    if dTh == 1:
        z_params = [3465, 1730, 1155]  # run1
    elif dTh == 2:
        # z_params = [2450, 1225, 815]  # run1
        z_params = [500, 900, 1600, 1900, 2500]  # run2
        r_params_ = [1900, 1300, 900, 800, 600]  # run2
    elif dTh == 3:
        # z_params = [4000, 2000, 1500, 1000, 670, 500, 250] # run1
        z_params = [500, 1000, 1600, 2000, 2500]  # run2
        r_params_ = [1500, 1000, 700, 600, 500]  # run2
    elif dTh == 4:
        # z_params = [1730, 870, 430]     # run1
        z_params = [500, 900, 1600, 2000, 2500]  # run2
        r_params_ = [1300, 900, 600, 500, 400]  # run2
    elif dTh == 5:
        z_params = [1000, 1000, 1000, 1000] # run5
        r_params = [500, 1100, 1600, 2300]  # run5



    # try:
    #     r_params = r_params_
    #     del r_params_
    # except:
    #     r_params = z_params[::-1]
    # print('z*: ', z_params)
    # print('r*: ', r_params)
    # print('')

    return z_params, r_params





def set_input_parameters(args, dTh, z_params, r_params):

    global path_root
    path_root = args.path_root
    # global path_out_data
    # path_out_data = os.path.join(path_root, 'data_analysis')
    # if not os.path.exists(path_out_data):
    #     os.mkdir(path_out_data)
    global path_out_figs
    path_out_figs = os.path.join(path_root, 'figs_config')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print('path root: ', path_root)
    print('path out figs: ', path_out_figs)
    # print('path out data: ', path_out_data)
    print('')


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
        ic_arr = np.zeros(1)
        jc_arr = np.zeros(1)
        ic_arr[0] = ic
        jc_arr[0] = jc
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
        print('marg from nml: ', marg)
    except:
        marg = 200.
        print('marg NOT from nml: ', marg)

    # if args.tmin:
    #     tmin = np.int(args.tmin)
    # else:
    #     tmin = 100
    # if args.tmax:
    #     tmax = np.int(args.tmax)
    # else:
    #     tmax = 100
    # times = np.arange(tmin, tmax + 100, 100)
    # print('times', times)
    # print('')

    return ic_arr, jc_arr, marg





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