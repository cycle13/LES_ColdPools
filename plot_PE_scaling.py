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

def main():
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

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')
    cm_fall = plt.cm.get_cmap('winter')
    cm_summer = plt.cm.get_cmap('spring')

    nml, dTh, z_params, r_params = set_input_parameters(args)
    n_params = len(r_params)
    id_ref = 'dTh3_z1000_r1000'
    path_ref = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/dTh3_z1000_r1000'

    aux = np.arange(-1, 4)
    PE_array = 2.**aux
    print('PE: ' + str(PE_array))

    # --------------------------------------

    ''' CP height '''
    fig_name = 'CP_height.png'
    # path_out = 'figs_CP_height'
    path_out = path_out_figs
    plot_CP_height(dTh, z_params, r_params, n_params, PE_array, id_ref, path_ref, path_out_figs, fig_name)


    ''' PE array '''
    fig_name = 'PE_array.png'
    plot_PE_array(r_params, PE_array, n_params, fig_name)

    fig_name = 'PE_array_envelope.png'
    plot_PE_array_enevelop(dTh, r_params, z_params, n_params, PE_array, fig_name)


    return


# --------------------------------------

def plot_CP_height(dTh, z_params, r_params, n_params, PE_array, id_ref, path_ref, path_out, fig_name):
    min_CP_height = np.zeros(n_params + 1)
    for istar in range(n_params):
        zstar = z_params[0]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        sth = 0.5
        file_name = 'CP_height_' + id +'_sth' + str(sth) + '.nc'
        fullpath_in = os.path.join(path_root, id, 'data_analysis', file_name)
        print('')
        print(fullpath_in)
        print('')

        file = nc.Dataset(fullpath_in, 'r')
        # CP_height_2d = file.groups['fields_2D'].variables['CP_height_2d'][:,:,:]
        CP_height = file.groups['timeseries'].variables['CP_height'][:]
        if istar == 0:
            min_CP_height[istar] = np.amin(CP_height)
        elif istar >= 1:
            min_CP_height[istar+1] = np.amin(CP_height)
        file.close()
    file_name = 'CP_height_' + id_ref + '_sth' + str(sth) + '.nc'
    fullpath_in = os.path.join(path_ref, 'data_analysis', file_name)
    file = nc.Dataset(fullpath_in, 'r')
    CP_height = file.groups['timeseries'].variables['CP_height'][:]
    min_CP_height[1] = np.amin(CP_height)
    file.close()

    fig, axes = plt.subplots(1, 2, sharex='none', figsize=(12, 5))
    axes[0].plot(PE_array, min_CP_height)
    # axes[0].plot(min_CP_height)
    fig.tight_layout()
    fig.savefig(os.path.join(path_root, path_out, fig_name))
    plt.close(fig)

    return

# --------------------------------------

def plot_PE_array(r_params, PE_array, n_params, fig_name):
    r_params_ = np.zeros(len(r_params) + 1)
    PE_array_ = np.zeros(len(r_params))
    fig, axes = plt.subplots(1, 2, sharex='none', figsize=(12, 5))
    for istar, r in enumerate(r_params):
        if istar == 0:
            r_params_[istar] = r_params[istar]
            PE_array_[istar] = PE_array[istar]
        elif istar == 1:
            r_params_[1] = 1000.
        if istar >= 1:
            r_params_[istar + 1] = r_params[istar]
            PE_array_[istar] = PE_array[istar + 1]
    axes[0].plot(np.log2(PE_array), r_params_, 'k', linewidth=0.5)
    axes[0].plot(np.log2(PE_array_), r_params, '--k', linewidth=0.5)
    axes[1].plot(r_params_, PE_array, 'k', linewidth=0.5)
    axes[1].plot(r_params, PE_array_, '--k', linewidth=0.5)
    for istar in range(n_params + 1):
        axes[0].plot(np.log2(PE_array[istar]), r_params_[istar], 'o', markersize=10, markeredgecolor='w', )
        axes[1].plot(r_params_[istar], PE_array[istar], 'o', markersize=10, markeredgecolor='w', )
    axes[0].plot(np.log2(PE_array[1]), r_params_[1], 'ko', markersize=10, markeredgecolor='w', )
    axes[1].plot(r_params_[1], PE_array[1], 'ko', markersize=10, markeredgecolor='w', )

    axes[0].set_xlim(-1.2, 3.2)
    axes[0].set_ylim(400, 2500)
    axes[1].set_xlim(400, 2500)
    axes[1].set_ylim(0, 8.5)
    axes[0].set_xlabel('log2(PE / PE_ref)')
    axes[0].set_ylabel('radius')
    axes[1].set_xlabel('radius')
    axes[1].set_ylabel('PE / PE_ref')
    fig.tight_layout()
    fig.savefig(os.path.join(path_root, path_out_figs, fig_name))
    plt.close(fig)
    return



def plot_PE_array_enevelop(dTh, r_params, z_params, n_params, PE_array, fig_name):
    dTh_ref = 3
    rstar_ref = 1000
    zstar_ref = 1000
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
        ax0.plot(x_arr - xc, zmax, label='dTh' + str(dTh) + ', z*' + str(zstar) + ', r*' + str(rstar))
    x_arr_ = np.array(x_arr, copy=True)
    irstar = np.int(np.double(rstar_ref) / dx_)
    x_arr_[:ic - irstar] = rstar_ref + xc
    x_arr_[ic + irstar:] = rstar_ref + xc
    zmax = zstar_ref * np.cos((x_arr_ - xc) / rstar_ref * np.pi / 2) ** 2
    ax0.plot(x_arr - xc, zmax, 'k',
             label='dTh' + str(dTh_ref) + ', z*' + str(zstar_ref) + ', r*' + str(rstar_ref))

    r_params_ = np.zeros(len(r_params) + 1)
    PE_array_ = np.zeros(len(r_params))
    for istar, rstar in enumerate(r_params):
        if istar == 0:
            r_params_[istar] = r_params[istar]
            PE_array_[istar] = PE_array[istar]
        elif istar == 1:
            r_params_[1] = 1000.
        if istar >= 1:
            r_params_[istar + 1] = r_params[istar]
            PE_array_[istar] = PE_array[istar + 1]
    ax1.plot(r_params_, PE_array, 'k', linewidth=0.5)
    ax1.plot(r_params, PE_array_, '--k', linewidth=0.5)
    for istar in range(n_params + 1):
        ax1.plot(r_params_[istar], PE_array[istar], 'o', markersize=10, markeredgecolor='w',
                 label='dTh' + str(dTh) + ', z*' + str(zstar) + ', r*' + str(rstar))
    ax1.plot(r_params_[1], PE_array[1], 'ko', markersize=10, markeredgecolor='w',
             label='dTh' + str(dTh_ref) + ', z*' + str(zstar_ref) + ', r*' + str(rstar_ref))

    ax0.set_xlim(-Lx / 2, Lx / 2)
    ax0.set_ylim(0, zstar_max + 100)
    ax1.set_xlim(400, 2500)
    ax1.set_ylim(0, 8.5)

    x_ticks = [np.int(i * 1e-3) for i in ax0.get_xticks()]
    y_ticks = [np.round(i * 1e-3, 1) for i in ax0.get_yticks()]
    ax0.set_xticklabels(x_ticks)
    ax0.set_yticklabels(y_ticks)
    x_ticks = [np.int(i * 1e-3) for i in ax1.get_xticks()]
    ax1.set_xticklabels(x_ticks)

    ax1.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0),
               fancybox=True, ncol=1)

    ax0.set_xlabel('Radius  [km]')
    ax0.set_ylabel('Height z  [km]')
    ax1.set_xlabel('Radius  [km]')
    ax1.set_ylabel('PE / PE_ref')
    fig.tight_layout()
    plt.subplots_adjust(bottom=0.12, right=.85, left=0.07, top=0.9, wspace=0.2)
    fig.savefig(os.path.join(path_root, path_out_figs, fig_name))
    plt.close(fig)
    return

# --------------------------------------

def set_input_parameters(args):
    print('--- set input parameters ---')
    global case_name
    global path_root, path_out_figs
    global times

    path_root = args.path_root
    path_out_figs = os.path.join(path_root, 'figs_comparison')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print('')
    print('path figs out: ')
    print('   ' + path_out_figs)
    print('')

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
            ic = np.int(nx/2)
            jc = np.int(ny/2)
            ic_arr = np.zeros(1)
            jc_arr = np.zeros(1)
            ic_arr[0] = ic
            jc_arr[0] = jc
        else:
            print('ic, jc not defined')

    return nml, dTh, z_params, r_params



# _____________________________________________________________________
if __name__ == '__main__':
    main()