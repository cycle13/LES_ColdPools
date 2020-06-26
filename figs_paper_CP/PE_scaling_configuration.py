import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import netCDF4 as nc
import argparse
import json as simplejson
import os
import time

label_size = 12
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 15
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35

execfile('settings.py')

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    # parser.add_argument("path_root")
    parser.add_argument("case")
    # parser.add_argument("dTh", type=int)
    # parser.add_argument("--zparams", nargs='+', type=int)
    # parser.add_argument('--rparams', nargs='+', type=int)
    # parser.add_argument("--tmin")
    # parser.add_argument("--tmax")
    args = parser.parse_args()

    dTh, z_params, r_params = set_input_parameters(args)
    n_params = len(r_params)
    dTh_ref = 3
    zstar_ref = 1000
    rstar_ref = 1000
    id_ref = 'dTh'+str(dTh_ref)+'_z'+str(zstar_ref)+'_r'+str(rstar_ref)
    # id_ref = 'dTh3_z1000_r1000'
    path_ref = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/dTh3_z1000_r1000'

    path_out_figs = os.path.join(path_root, 'figs_configuration')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print('')
    print('path figs out: ')
    print('   ' + path_out_figs)
    print('')

    PE_array = 2. ** np.arange(-1, 4)
    print('PE: ' + str(PE_array))
    fig_name = 'PE_array_envelope.png'
    plot_PE_array_enevelop(dTh, r_params, z_params, n_params,
                           dTh_ref, zstar_ref, PE_array,
                           path_out_figs, fig_name)

    return

# --------------------------------------

def plot_PE_array_enevelop(dTh, r_params, z_params, n_params,
                           dTh_ref, zstar_ref, PE_array,
                           path_out_figs, fig_name):
    # dTh_ref = 3
    # rstar_ref = 1000
    # zstar_ref = 1000
    zstar = z_params[0]
    zstar_max = np.maximum(np.amax(z_params), zstar_ref)

    Lx = 5.e3
    dx_ = 10
    nx_ = np.int(Lx / dx_)
    ic = nx_ / 2
    xc = (ic + 0.5) * dx_
    x_arr = np.arange(0, Lx, dx_)
    fig, (ax0, ax1) = plt.subplots(1, 2, sharex='none', figsize=(12, 5))
    for istar, rstar in enumerate(r_params):
        irstar = np.int(np.double(rstar) / dx_)
        x_arr_ = np.array(x_arr, copy=True)
        x_arr_[:ic - irstar] = rstar + xc
        x_arr_[ic + irstar:] = rstar + xc
        zmax = zstar * np.cos((x_arr_ - xc) / rstar * np.pi / 2) ** 2
        ax0.plot(x_arr - xc, zmax, color=colorlist5[istar], label='dTh' + str(dTh) + ', z*' + str(zstar) + ', r*' + str(rstar))
    # x_arr_ = np.array(x_arr, copy=True)
    # irstar = np.int(np.double(rstar_ref) / dx_)
    # x_arr_[:ic - irstar] = rstar_ref + xc
    # x_arr_[ic + irstar:] = rstar_ref + xc
    # zmax = zstar_ref * np.cos((x_arr_ - xc) / rstar_ref * np.pi / 2) ** 2
    # ax0.plot(x_arr - xc, zmax, 'k',
    #          label='dTh' + str(dTh_ref) + ', z*' + str(zstar_ref) + ', r*' + str(rstar_ref))

    r_params_ = np.zeros(len(r_params) + 1)
    print('.........')
    print(r_params)
    print(PE_array)
    print('....', len(r_params), np.shape(r_params_), np.shape(PE_array))
    # PE_array_ = np.zeros(len(r_params))
    # for istar, rstar in enumerate(r_params):
    #     if istar == 0:
    #         r_params_[istar] = r_params[istar]
    #         PE_array_[istar] = PE_array[istar]
    #     elif istar == 1:
    #         r_params_[1] = 1000.
    #     if istar >= 1:
    #         r_params_[istar + 1] = r_params[istar]
    #         PE_array_[istar] = PE_array[istar + 1]
    ax1.plot(r_params_, PE_array, 'k', linewidth=0.5)
    # ax1.plot(r_params, PE_array_, '--k', linewidth=0.5)
    # for istar in range(n_params + 1):
    #     ax1.plot(r_params_[istar], PE_array[istar], 'o', markersize=10, markeredgecolor='w',
    #              label='dTh' + str(dTh) + ', z*' + str(zstar) + ', r*' + str(rstar))
    # ax1.plot(r_params_[1], PE_array[1], 'ko', markersize=10, markeredgecolor='w',
    #          label='dTh' + str(dTh_ref) + ', z*' + str(zstar_ref) + ', r*' + str(rstar_ref))
    #
    # ax0.set_xlim(-Lx / 2, Lx / 2)
    # ax0.set_ylim(0, zstar_max + 100)
    # ax1.set_xlim(400, 2500)
    # ax1.set_ylim(0, 8.5)
    #
    # x_ticks = [np.int(i * 1e-3) for i in ax0.get_xticks()]
    # y_ticks = [np.round(i * 1e-3, 1) for i in ax0.get_yticks()]
    # ax0.set_xticklabels(x_ticks)
    # ax0.set_yticklabels(y_ticks)
    # x_ticks = [np.int(i * 1e-3) for i in ax1.get_xticks()]
    # ax1.set_xticklabels(x_ticks)
    #
    # ax1.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0),
    #            fancybox=True, ncol=1)
    #
    ax0.set_xlabel('Radius R [km]')
    ax0.set_ylabel('Height z  [km]')
    ax1.set_xlabel('Radius R [km]')
    ax1.set_ylabel('PE / PE_ref')
    # fig.tight_layout()
    plt.subplots_adjust(bottom=0.12, right=.85, left=0.07, top=0.9, wspace=0.2)
    fig.savefig(os.path.join(path_root, path_out_figs, fig_name))
    plt.close(fig)
    return
# --------------------------------------

def set_input_parameters(args):
    print('--- set input parameters ---')
    global case_name
    global path_root
    # global times

    if args.case == 'run5':
        path_root = '/nbi/ac/conv4/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/single_3D_noise/run5_PE_scaling_dx100m'
    elif args.case == 'case6':
        path_root = '/nbi/ac/conv4/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/single_3D_noise/run6_PE_scaling_dx50m'
    print('path root: ' + path_root)
    root_name = os.path.basename(path_root)
    print(root_name[:4])


    if root_name[:4] == 'run5':
        dTh = 5
        z_params = [1000]
        r_params = [500, 1100, 1600, 2000, 2300]

    # dTh = args.dTh
    # z_params = args.zparams
    # r_params = args.rparams
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
    # if args.tmin:
    #     tmin = np.int(args.tmin)
    # else:
    #     tmin = np.int(100)
    # if args.tmax:
    #     tmax = np.int(args.tmax)
    # else:
    #     tmax = np.int(100)
    # times = np.arange(tmin, tmax + 100, 100)
    # # times = [np.int(name[:-3]) for name in files]
    # times.sort()
    # print('times', times)

    global ic_arr, jc_arr
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
            ic = np.int(nx / 2)
            jc = np.int(ny / 2)
            ic_arr = np.zeros(1)
            jc_arr = np.zeros(1)
            ic_arr[0] = ic
            jc_arr[0] = jc
        else:
            print('ic, jc not defined')

    return dTh, z_params, r_params


# _____________________________________________________________________
if __name__ == '__main__':
    main()