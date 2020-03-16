import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import netCDF4 as nc
import argparse
import json as simplejson
import os
import time

execfile('settings.py')
label_size = 15
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 21
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35

def main():
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='LES_CP')
    # parser.add_argument("casename")
    # parser.add_argument("path_root")
    # parser.add_argument("dTh", type=int)
    # parser.add_argument("--zparams", nargs='+', type=int)
    # parser.add_argument('--rparams', nargs='+', type=int)
    # parser.add_argument("--tmin")
    # parser.add_argument("--tmax")
    # args = parser.parse_args()

    path_out_figs = os.path.join('./figs_paper_CP')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print('')
    print('path figs out: ' + path_out_figs)
    print('')

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')
    cm_fall = plt.cm.get_cmap('winter')
    cm_summer = plt.cm.get_cmap('spring')

    # nml, dTh, z_params, r_params = set_input_parameters(args)
    # n_params = len(r_params)
    dTh_ref = 3
    rstar_ref = 1000
    zstar_ref = 1000
    id_ref = 'dTh3_z1000_r1000'
    path_ref = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/dTh3_z1000_r1000'

    path_run5 = '/nbi/ac/conv3/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/single_3D_noise/run5_PE_scaling_dx100m'
    dTh = 5
    r_params = [500, 1100, 1600, 2300]
    r_params_ = [500, 1000, 1100, 1600, 2300]
    z_params = [1000]
    PE_array = [0.5, 2, 4, 8]
    PE_array_ = 2.**np.arange(-1, 4)
    # print('PE: ' + str(PE_array))
    n_params = len(r_params)

    # --------------------------------------

    ''' envelope '''
    Lx = 6e3
    dx_ = 10
    nx_ = Lx / dx_
    x_arr = np.arange(0, Lx, dx_)
    ic = np.int(nx_ / 2)
    xc = x_arr[ic]
    # print nx_, ic, xc
    #
    zmax = np.zeros((n_params+1, nx_))
    for istar in range(n_params):
        zstar = z_params[0]
        rstar = r_params[istar]
        irstar = np.int(np.double(rstar) / dx_)
        x_arr_ = np.array(x_arr, copy=True)
        x_arr_[:ic - irstar] = rstar + xc
        x_arr_[ic + irstar:] = rstar + xc
        zmax[istar,:] = zstar * np.cos((x_arr_ - xc) / rstar * np.pi / 2) ** 2
    x_arr_ = np.array(x_arr, copy=True)
    irstar = np.int(np.double(rstar_ref) / dx_)
    x_arr_[:ic - irstar] = rstar_ref + xc
    x_arr_[ic + irstar:] = rstar_ref + xc
    zmax[-1, :] = zstar_ref * np.cos((x_arr_ - xc) / rstar_ref * np.pi / 2) ** 2




    fig_name = 'PE_scaling_run5.png'

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(16, 5))
    for istar in range(n_params):
        zstar = z_params[0]
        rstar = r_params[istar]
        ax0.plot(x_arr - xc, zmax[istar,:], label='dTh' + str(dTh) + ', z*' + str(zstar) + ', r*' + str(rstar))
    ax0.plot(x_arr - xc, zmax[-1], 'k',
             label='dTh' + str(dTh_ref) + ', z*' + str(zstar_ref) + ', r*' + str(rstar_ref))
    ax0.set_xlabel('Radius r', fontsize=21)
    ax0.set_ylabel('Height z', fontsize=21)
    ax0.set_ylim(0, z_params[-1] + 100)
    x_ticks = [np.int(i * 1e-3) for i in ax0.get_xticks()]
    y_ticks = [np.round(i * 1e-3, 1) for i in ax0.get_yticks()]
    ax0.set_xticklabels(x_ticks)
    ax0.set_yticklabels(y_ticks)

    ax0.legend(loc='upper right', bbox_to_anchor=(-0.2, 1.0),
               fancybox=True, ncol=1)

    ax1.plot(r_params_, PE_array_, 'k', linewidth=0.5)
    ax1.plot(r_params, PE_array, '--k', linewidth=0.5)
    for istar in range(n_params):
    # #     axes[0].plot(np.log2(PE_array[istar]), r_params_[istar], 'o', markersize=10, markeredgecolor='w', )
        ax1.plot(r_params[istar], PE_array[istar], 'o', markersize=10, markeredgecolor='w', )
    ax1.plot(r_params_[1], PE_array_[1], 'ko', markersize=10, markeredgecolor='w', )

    # ax1.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0),
    #            fancybox=True, ncol=1)
    x_ticks = [np.round(i * 1e-3,1) for i in ax1.get_xticks()]
    print x_ticks
    # ax1.set_xticklabels(x_ticks)
    ax1.set_xlabel('Radius r')
    ax1.set_ylabel(r'PE / PE$_0$')
    ax1.set_xlim(400,2500)
    ax1.set_ylim(0,8.5)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, right=.85, left=0.2, top=0.9, wspace=0.25)
    plt.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return


# --------------------------------------



# --------------------------------------
# --------------------------------------

def set_input_parameters(args):
    print('--- set input parameters ---')
    global case_name
    global path_root, path_out_figs
    global times

    path_root = args.path_root


    dTh = args.dTh
    z_params = args.zparams
    r_params = args.rparams
    print('z*: ', z_params)
    print('r*: ', r_params)

    case_name = args.casename
    id0 = 'dTh' + str(dTh) + '_z' + str(z_params[0]) + '_r' + str(r_params[0])
    nml = simplejson.loads(open(os.path.join(path_root, id0, case_name + '.in')).read())
    global nx, ny, nz, dx, dV, gw
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = np.zeros(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx[0] * dx[1] * dx[2]

    ''' determine file range '''
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = np.int(100)
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = np.int(100)
    times = np.arange(tmin, tmax + 100, 100)
    # times = [np.int(name[:-3]) for name in files]
    times.sort()
    print('times', times)

    return nml, dTh, z_params, r_params



# _____________________________________________________________________
if __name__ == '__main__':
    main()
