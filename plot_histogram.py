import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import netCDF4 as nc
import os
import argparse
import json as simplejson


def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--casename")
    parser.add_argument("--path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    # parser.add_argument("--kmin")
    # parser.add_argument("--kmax")
    args = parser.parse_args()

    global cm_bwr, cm_grey, cm_vir
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_vir = plt.cm.get_cmap('viridis')
    cm_grey = plt.cm.get_cmap('gist_gray_r')

    global path_out, path_fields, path_rim
    if args.path:
        path = args.path
    else:
        path = '/Users/bettinameyer/polybox/ClimatePhysics/Copenhagen/Projects/LES_ColdPool/' \
               'triple_3D_noise/Out_CPDry_triple_dTh2K/'
        # path = '/nbi/ac/cond1/meyerbe/ColdPools/triple_3D_noise/Out_CPDry_triple_Th3K/run1/'
    # path_in = os.path.join(path, 'fields_cp_rim')
    path_fields = os.path.join(path, 'fields')
    path_out = os.path.join(path, 'figs_hist')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    path_in = os.path.join(path, 'fields_cp_rim')

    global case_name
    if args.casename:
        case_name = args.casename
    else:
        case_name = 'ColdPoolDry_triple_3D'
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dy, dz, dV
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx * dy * dz
    global kmax
    kmax = 100

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
    print('times:', timerange)


    for t0 in timerange:
        # # :::
        # mask_file_name = 'rimmask_perc' + str(perc) + 'th' + '_t100.nc'
        # rootgrp = nc.Dataset(os.path.join(path_rim, mask_file_name), 'r')
        # mask = rootgrp.groups['fields'].variables['mask'][:, :, :]
        # # make sure that only non-zero values in krange if nc-files contain this level
        # contr = rootgrp.groups['profiles'].variables['k_dumped'][:]
        # zrange_ = rootgrp.groups['profiles'].variables['zrange'][:]
        # krange_ = rootgrp.groups['profiles'].variables['krange'][:]
        # ic = np.int(rootgrp.groups['description'].variables['ic'][:])
        # jc = np.int(rootgrp.groups['description'].variables['jc'][:])
        # ishift = np.int(rootgrp.groups['description'].variables['ishift'][:])
        # jshift = np.int(rootgrp.groups['description'].variables['jshift'][:])
        # [nx_, ny_, nk] = mask.shape
        # id = np.int(nx_ / 2)
        # jd = np.int(ny_ / 2)
        # rootgrp.close()
        # w_ = np.roll(np.roll(w[:, :, :], ishift, axis=0), jshift, axis=1)[i0, jc + jshift - jd:jc + jshift + jd, :nz_]
        # # :::

        w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0) + '.nc'))
        percentiles = np.arange(90,101,1)
        plot_hist(w, percentiles, t0)


    return



#----------------------------------------------------------------------
def plot_hist(data, perc_range, t0):
    n_bins = 1e2
    min = np.amin(data)
    max = np.amax(data)
    bins = np.linspace(min, max, n_bins)
    print('min', min, 'max', max)

    perc = np.zeros(len(perc_range))
    for i,p in enumerate(perc_range):
        perc[i] = np.percentile(data[:,:,:kmax], p)

    fig, axes = plt.subplots(2,1, sharex=True, figsize=(12,10))
    plt.subplot(211)
    n, bins, patches = plt.hist(data[:,:,:kmax].flatten(), bins)
    max_hist = np.amax(n)
    print('n', np.amax(n), n.shape)
    # n, bins, patches = plt.hist(data, bins, normed=False, facecolor='red', alpha=0.4, label='w')
    for i,p in enumerate(perc_range):
        plt.plot([perc[i], perc[i]], [0, max_hist], 'k', linewidth=1, label=str(p)+'th percentile')

    # histogram on log scale.
    # Use non-equal bin sizes, such that they look equal on log scale.
    n_bins = 50
    bins_pos = np.linspace(0,max,n_bins)
    logbins = np.logspace(np.log10(bins_pos[1]), np.log10(bins_pos[-1]), len(bins_pos))
    plt.subplot(212)
    plt.hist(data[:,:,:kmax].flatten(), bins=logbins)
    plt.xscale('log')
    max_hist = 4e4
    for i,p in enumerate(perc_range):
        plt.plot([perc[i], perc[i]], [0, max_hist], 'k', linewidth=3, label=str(p)+'th percentile')
    plt.legend(loc='right')#, bbox_to_anchor=(0, 0))

    plt.xlabel('w  [m/s]')
    plt.suptitle('Histogram of w' + '   (t=' + str(t0) + ')')
    plt.savefig(os.path.join(path_out, 'hist_w_t' + str(t0) + '.png'))
    plt.close()

    return


# ----------------------------------

def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:,:,:]
    rootgrp.close()
    return data


if __name__ == '__main__':
    main()
