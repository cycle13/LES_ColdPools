import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--casename")
    parser.add_argument("--path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    # parser.add_argument("--k0", nargs = '+', type = int)
    parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    global path_fields, path_out
    if args.path:
        path = args.path
    else:
        # path = '/Users/bettinameyer/polybox/ClimatePhysics/Copenhagen/Projects/LES_ColdPool/' \
        #        'triple_3D_noise/Out_CPDry_triple_dTh2K/'
        path = '/nbi/ac/cond1/meyerbe/ColdPools/triple_3D_noise/Out_CPDry_triple_Th3K/run1/'
    path_mask = os.path.join(path, 'fields_cp_rim')
    path_out = os.path.join(path, 'figs_cp_rim')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    global case_name
    if args.casename:
        case_name = args.casename
    else:
        case_name = 'ColdPoolDry_triple_3D'

    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = tmin
    timerange = np.arange(tmin, tmax + 100, 100)
    nt = len(timerange)

    if args.kmin:
        kmin = np.int(args.kmin)
    else:
        kmin = 5
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = 5
    krange = np.arange(kmin, kmax+1, 1)

    perc = 95
    global nx_, ny_, nk
    global ic, jc, ishift, jshift, id, jd
    mask_file_name = 'rimmask_perc' + str(perc) + 'th' + '_t' + str(np.int(timerange[0])) + '.nc'
    print(path_mask, mask_file_name)
    rootgrp = nc.Dataset(os.path.join(path_mask, mask_file_name), 'r')
    mask = rootgrp.groups['fields'].variables['mask'][:, :, :]
    # make sure that only non-zero values in krange if nc-files contain this level
    contr = rootgrp.groups['profiles'].variables['k_dumped'][:]
    zrange_ = rootgrp.groups['profiles'].variables['zrange'][:]
    krange_ = rootgrp.groups['profiles'].variables['krange'][:]
    ic = np.int(rootgrp.groups['description'].variables['ic'][:])
    jc = np.int(rootgrp.groups['description'].variables['jc'][:])
    ishift = np.int(rootgrp.groups['description'].variables['ishift'][:])
    jshift = np.int(rootgrp.groups['description'].variables['jshift'][:])
    [nx_, ny_, nk] = mask.shape
    id = np.int(nx_ / 2)
    jd = np.int(ny_ / 2)
    rootgrp.close()
    nk_ = np.int(np.sum(contr))
    nk = len(krange_)
    krange = np.zeros(nk_)
    zrange = np.zeros(nk_)
    krange[:] = krange_[0:nk_]
    zrange[:] = zrange_[0:nk_]
    del contr, krange_, zrange_
    del mask
    print('ic,jc', ic, jc, 'shifts', ishift, jshift)

    t0 = tmin
    k0 = kmin


    plot_mask_vel(path_fields, k0, t0)

    return
# ----------------------------------

def plot_mask_vel(path_in, k0, t0):
    u = read_in_netcdf_fields('u', os.path.join(path_in, str(t0) + '.nc'))
    v = read_in_netcdf_fields('v', os.path.join(path_in, str(t0) + '.nc'))
    w = read_in_netcdf_fields('w', os.path.join(path_in, str(t0) + '.nc'))

    # define mask
    u_ = np.roll(np.roll(u[:, :, k0], ishift, axis=0), jshift, axis=1)[ic + ishift - id:ic + ishift + id,
         jc + jshift - jd:jc + jshift + jd]
    v_ = np.roll(np.roll(v[:, :, k0], ishift, axis=0), jshift, axis=1)[ic + ishift - id:ic + ishift + id,
         jc + jshift - jd:jc + jshift + jd]
    w_roll = np.roll(np.roll(w[:, :, k0], ishift, axis=0), jshift, axis=1)  # nbi
    # w_roll = np.roll(w[:, :, k0], [ishift, jshift], [0, 1])     # on Laptop gives the same as above
    w_ = w_roll[ic + ishift - id:ic + ishift + id, jc + jshift - jd:jc + jshift + jd]

    plt.figure(figsize=(16, 6))
    # plt.subplot(2, 4, 2)
    # max = np.maximum(np.abs(np.amin(u_)), np.amax(u_))
    # plt.contourf(u_[:, :].T, cmap=cm_bwr, vmin=-max, vmax=max)
    # plt.colorbar()
    # plt.subplot(2, 4, 3)
    # max = np.maximum(np.abs(np.amin(v_)), np.amax(v_))
    # plt.contourf(v_[:, :].T, cmap=cm_bwr, vmin=-max, vmax=max)
    # plt.colorbar()
    # plt.subplot(2, 4, 4)
    # max = np.maximum(np.abs(np.amin(w_)), np.amax(w_))
    # plt.contourf(w_[:, :].T, cmap=cm_bwr, vmin=-max, vmax=max)
    # plt.colorbar()
    #
    # plt.subplot(2, 4, 5)
    # plt.contourf(mask_[:, :].T)
    # plt.colorbar()
    # plt.subplot(2, 4, 6)
    # max = np.maximum(np.abs(np.amin(u_)), np.amax(u_))
    # plt.contourf(u_mask[:, :].T, cmap=cm_bwr, vmin=-max, vmax=max)
    # plt.colorbar()
    # plt.subplot(2, 4, 7)
    # max = np.maximum(np.abs(np.amin(v_)), np.amax(v_))
    # plt.contourf(v_mask[:, :].T, cmap=cm_bwr, vmin=-max, vmax=max)
    # plt.colorbar()
    # plt.subplot(2, 4, 8)
    # max = np.maximum(np.abs(np.amin(w_)), np.amax(w_))
    # plt.contourf(w_mask[:, :].T, cmap=cm_bwr, vmin=-max, vmax=max)
    # plt.colorbar()
    plt.savefig(os.path.join(path_out, 'comparison_rim_k' + str(k0) + '.png'))
    plt.close()

    del u, v, w
    return


# ----------------------------------

def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:, :, :]
    rootgrp.close()
    return data

if __name__ == '__main__':
    main()