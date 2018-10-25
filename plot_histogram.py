import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import netCDF4 as nc
import os
import argparse
import json as simplejson
from scipy.optimize import curve_fit


def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--casename")
    parser.add_argument("--path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    # parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    global cm_bwr, cm_grey, cm_vir
    cm_bwr = plt.cm.get_cmap('bwr')
    try:
        cm_vir = plt.cm.get_cmap('viridis')
    except:
        cm_vir = plt.cm.get_cmap('jet')
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
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = nz

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
        percentiles_l = np.append(np.arange(95,99,1), np.arange(99,100,0.5))
        percentiles_s = np.append(np.arange(0,1,0.2), np.arange(1,4,1))
        # plot_hist_all_levels(w, percentiles_l, t0)

        plot_hist_twolevels(w, percentiles_l, percentiles_s, 50, 20, t0)
        # plot_hist_twolevels_log(w, percentiles_l, nz, 50, t0)

    # plot w along crosssection
    perc = 95
    global ic, jc
    mask_file_name = 'rimmask_perc' + str(perc) + 'th' + '_t100.nc'
    print(mask_file_name)
    rootgrp = nc.Dataset(os.path.join(path_in, mask_file_name), 'r')
    mask = rootgrp.groups['fields'].variables['mask'][:, :, :]
    # # make sure that only non-zero values in krange if nc-files contain this level
    contr = rootgrp.groups['profiles'].variables['k_dumped'][:]
    # zrange_ = rootgrp.groups['profiles'].variables['zrange'][:]
    krange_ = rootgrp.groups['profiles'].variables['krange'][:]
    ic = np.int(rootgrp.groups['description'].variables['ic'][:])
    jc = np.int(rootgrp.groups['description'].variables['jc'][:])
    ishift = np.int(rootgrp.groups['description'].variables['ishift'][:])
    jshift = np.int(rootgrp.groups['description'].variables['jshift'][:])
    [nx_, ny_, nk] = mask.shape
    id = np.int(nx_ / 2)
    jd = np.int(ny_ / 2)
    rootgrp.close()

    d = np.int(np.round(ny / 2))
    dhalf = np.int(np.round(ny / 4))
    a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
    ic1 = np.int(np.round(a / 2))# + gw
    jc1 = np.int(dhalf)# + gw  # np.int(np.round(d/2) + Gr.dims.gw)

    nk_ = np.int(np.sum(contr))
    # nk = len(krange_)
    krange = np.zeros(nk_)
    krange[:] = krange_[0:nk_]
    # zrange = np.zeros(nk_)
    # zrange[:] = zrange_[0:nk_]

    # w cross-section and fit curve
    # s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
    for t0 in timerange:
        i0 = id
        w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0) + '.nc'))
    #     w_ = np.roll(np.roll(w[:, :, :], ishift, axis=0), jshift, axis=1)[i0, jc + jshift - jd:jc + jshift + jd, :]
    #     yrange = np.arange(0,ny_)
    #     # plt.figure(figsize=(10,8))
    #     # for ik,k in enumerate(krange):
    #     #     k = np.int(k)
    #     #     plt.plot(yrange, w_[:,k], label='k='+str(k))
    #     # plt.legend()
    #     # plt.savefig(os.path.join(path_out, 'w_crosssection_yz_t' + str(t0) + '.png'))
    #     # plt.close()
    #
    #     # # plt.figure(figsize=(10, 8))
    #     # s = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))
    #     # fig, ax1 = plt.subplots(figsize=(10, 8))
    #     # ax1.contourf(s[ic,:,:50].T, alpha=0.4)
    #     # ax1.plot([jc, jc], [0, 49], 'k')
    #     # ax1.set_ylabel('k', color='b')
    #     # ax2 = ax1.twinx()
    #     # yrange = np.arange(0, ny)
    #     # for ik,k in enumerate(krange):
    #     #     k = np.int(k)
    #     #     ax2.plot(yrange, w[ic, :,k], label='k='+str(k))
    #     # ax2.set_ylabel('w')
    #     # ax2.set_ylim(-5,5)
    #     # plt.legend(loc='best', bbox_to_anchor=(0, 1))
    #     # fig.tight_layout()
    #     # plt.savefig(os.path.join(path_out, 'w_crosssection_yz_t' + str(t0) + '.png'))
    #     # plt.close()
    #
        # fit curve
        k0 = 5
        ymin = jc
        ymax = jc + 31 + (t0-500)/100 * np.int(300./dy)
        ymax_ = ymax + 20
        yrange_ = np.arange(ymin, ymax_)
        s = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))
        fig, ax1 = plt.subplots(figsize=(10, 8))
        ax1.contourf(s[ic, ymin:ymax_, :50].T, alpha=0.4)
        ax1.plot([jc - ymin, jc - ymin], [0, 49], 'k')
        ax1.set_xlabel('y  [-]')
        ax1.set_ylabel('k  [-]')
        ax2 = ax1.twinx()
        ax2.set_ylabel('w')
        for ik,k in enumerate(krange):
            k = np.int(k)
            if k == k0:
                ax2.plot(yrange_-ymin, w[ic, ymin:ymax_,k], label='k='+str(k), linewidth=3)
            else:
                ax2.plot(yrange_-ymin, w[ic, ymin:ymax_,k], label='k='+str(k))
        print('-- fit curve --')
        xdata = np.arange(0, ymax - ymin)
        ydata = w[ic, ymin:ymax, k0]
        try:
            popt2, pcov = curve_fit(fit_func, xdata, ydata)
            ax2.plot(yrange_ - ymin, fit_func(yrange_ - ymin, *popt2), 'kx-', label='fit', linewidth=1)
            ax2.plot(yrange - ymin, fit_func(yrange - ymin, *popt2), 'k-', label='fitting range', linewidth=3)
        except:
            pass
        yrange = np.arange(ymin, ymax)
        ax2.set_ylim(-5,5)
        plt.legend(loc='best', bbox_to_anchor=(0, 1))
        fig.tight_layout()
        plt.savefig(os.path.join(path_out, 'fit_interior_w_crosssection_yz_t' + str(t0) + '.png'))
        plt.close()

    # ymin = jc
    # ymax = jc+30
    # k0 = 5
    # xdata = np.arange(0,30)
    # ydata = w[ic, ymin:ymax, k0]
    # popt, pcov = curve_fit(fit_func, xdata, ydata)
    # fig, ax1 = plt.subplots(figsize=(10, 8))
    # ax1.plot(xdata, ydata, 'o', label='data')
    # ax2 = ax1.twinx()
    # ax2.plot(xdata, fit_func(xdata, *popt))
    # plt.legend()
    # plt.savefig(os.path.join(path_out, 'fitting.png'))
    # plt.close()

    return



#----------------------------------------------------------------------
def plot_hist_all_levels(data, perc_range, t0):
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
    plt.ylabel('n')

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
    plt.ylabel('log_10(n)')
    plt.suptitle('Histogram of w' + '   (t=' + str(t0) + r',  k$\in$[0,' + str(kmax)+'])')
    plt.savefig(os.path.join(path_out, 'hist_w_kmax'+str(kmax)+'_t' + str(t0) + '.png'))
    plt.close()

    return



def plot_hist_twolevels(data, perc_range_l, perc_range_s, kmax1, kmax2, t0):
    n_bins = 1e2
    min = np.amin(data)
    max = np.amax(data)
    bins = np.linspace(min, max, n_bins)
    print('min', min, 'max', max)

    perc1 = np.zeros(len(perc_range_l))
    perc2 = np.zeros(len(perc_range_l))
    perc3 = np.zeros(len(perc_range_l))
    for i,p in enumerate(perc_range_l):
        perc1[i] = np.percentile(data[:,:,:nz], p)
        perc2[i] = np.percentile(data[:,:,:kmax1], p)
        perc3[i] = np.percentile(data[:,:,:kmax2], p)
    perc1s = np.zeros(len(perc_range_s))
    perc2s = np.zeros(len(perc_range_s))
    perc3s = np.zeros(len(perc_range_s))
    for i,p in enumerate(perc_range_s):
        perc1s[i] = np.percentile(data[:,:,:nz], p)
        perc2s[i] = np.percentile(data[:,:,:kmax1], p)
        perc3s[i] = np.percentile(data[:,:,:kmax2], p)

    fig, axes = plt.subplots(3,1, sharex=True, figsize=(15,25))
    # fig, axes = plt.subplots(3,1, figsize=(15,25))
    plt.subplot(311)
    n, bins_, patches = plt.hist(data[:,:,:nz].flatten(), bins, label='kmax=nz ('+str(nz)+')')
    max_hist = np.amax(n)
    n, bins_, patches = plt.hist(data[:,:,:kmax1].flatten(), bins=bins,
                                 facecolor='y', alpha=0.4, label='kmax='+str(kmax1))
    n, bins_, patches = plt.hist(data[:,:,:kmax2].flatten(), bins=bins,
                                 facecolor='r', alpha=0.4, label='kmax='+str(kmax2))
    max_hist = np.maximum(np.amax(n), max_hist)
    print('n', np.amax(n), n.shape)
    # n, bins, patches = plt.hist(data, bins, normed=False, facecolor='red', alpha=0.4, label='w')
    for i,p in enumerate(perc_range_l):
        plt.plot([perc1[i], perc1[i]], [0, max_hist], 'k', linewidth=2, label=str(p)+'th percentile')
        plt.plot([perc2[i], perc2[i]], [0, max_hist], 'k--', linewidth=2)
        plt.plot([perc3[i], perc3[i]], [0, max_hist], 'k:', linewidth=2)
    plt.ylim([0,4e5])
    plt.title('all w: max='+str(np.round(max,2))+', min='+str(np.round(min,2)))
    plt.ylabel('n')
    plt.legend(loc='upper left')  # , bbox_to_anchor=(0, 0))

    plt.subplot(312)
    plt.hist(data[:,:,:nz].flatten(), bins, label='kmax=nz ('+str(nz)+')')
    plt.hist(data[:,:,:kmax1].flatten(), bins, facecolor='y', alpha=0.4, label='kmax='+str(kmax1))
    plt.hist(data[:,:,:kmax2].flatten(), bins, facecolor='r', alpha=0.4, label='kmax='+str(kmax2))
    for i,p in enumerate(perc_range_l):
        plt.plot([perc1[i], perc1[i]], [0, max_hist], 'k', linewidth=2, label=str(p)+'th percentile')
        plt.plot([perc2[i], perc2[i]], [0, max_hist], 'k--', linewidth=2)
        plt.plot([perc3[i], perc3[i]], [0, max_hist], 'k:', linewidth=2)
    for i,p in enumerate(perc_range_s):
        plt.plot([perc1s[i], perc1s[i]], [0, max_hist], 'g', linewidth=2, label=str(p)+'th percentile')
        plt.plot([perc2s[i], perc2s[i]], [0, max_hist], 'g--', linewidth=2)
        plt.plot([perc3s[i], perc3s[i]], [0, max_hist], 'g:', linewidth=2)
    plt.yscale('log', nonposy='clip')
    plt.title('all w')
    plt.ylabel('log(n)')
    plt.legend(loc='upper left', fontsize=10)#, bbox_to_anchor=(0, 0))

    plt.subplot(313)
    # perc_range = np.arange(0,10,1)
    perc1 = np.zeros(len(perc_range_l))
    perc2 = np.zeros(len(perc_range_l))
    for i, p in enumerate(perc_range_l):
        perc1[i] = np.percentile(np.abs(data[:, :, :kmax1]), p)
        perc2[i] = np.percentile(np.abs(data[:, :, :kmax2]), p)
    bins_pos = np.arange(0, np.maximum(np.abs(min), max), 1e-1)
    plt.hist(np.abs(data[:,:,:nz]).flatten(), bins=bins_pos)
    plt.hist(np.abs(data[:,:,:kmax1]).flatten(), bins=bins_pos, facecolor='y', alpha=0.4)
    plt.hist(np.abs(data[:,:,:kmax2]).flatten(), bins=bins_pos, facecolor='r', alpha=0.4)
    for i,p in enumerate(perc_range_l):
        plt.plot([perc1[i], perc1[i]], [0, max_hist], 'k', linewidth=2, label=str(p)+'th percentile')
        plt.plot([perc2[i], perc2[i]], [0, max_hist], 'k--', linewidth=2)
    plt.yscale('log', nonposy='clip')
    plt.title('abs(w)')
    plt.ylabel('log(n)')

    plt.xlabel('w  [m/s]')
    plt.suptitle('Histogram of w' + '   (t=' + str(t0) + '), max=' +str(np.amax(data)))
    plt.savefig(os.path.join(path_out, 'hist_w_comp_t' + str(t0) + '.png'))
    plt.close()

    return




# def plot_hist_twolevels_log(data, perc_range, kmax1, kmax2, t0):
#     min = np.floor(np.amin(data))
#     max = np.amax(data)
#     print('min', min, 'max', max)
#     kmax1 = nz
#
#     perc1 = np.zeros(len(perc_range))
#     perc2 = np.zeros(len(perc_range))
#     for i, p in enumerate(perc_range):
#         perc1[i] = np.percentile(data[:, :, :kmax1], p)
#         perc2[i] = np.percentile(data[:, :, :kmax2], p)
#
#     # fig, axes = plt.subplots(2,1, sharex=True, figsize=(12,10))
#     fig, axes = plt.subplots(2,1, figsize=(12,10))
#
#     plt.subplot(211)
#     n_bins = 50
#     bins = np.linspace(min, max, n_bins)
#     bins_pos = np.arange(0, max, 1e-1)
#     bins_neg = np.arange(-min, 0, -1e-1)
#     print(bins_pos)
#     print(bins_neg)
#     logbins_pos = np.logspace(np.log10(bins_pos[1]), np.log10(bins_pos[-1]), len(bins_pos))
#     logbins_neg = -np.logspace(np.log10(bins_neg[1]), np.log10(bins_neg[-1]), len(bins_neg))
#     print(logbins_pos)
#     print(logbins_neg)
#
#
#     n1, bins_, patches = plt.hist(data[:,:,:kmax1].flatten(), bins=logbins_pos, label='kmax=nz ('+str(kmax1)+')')
#     # n2, bins_, patches = plt.hist(data[:,:,:kmax1].flatten(), bins=logbins_neg)
#     # max_hist = np.maximum(np.amax(n1), np.amax(n2))
#     plt.hist(data[:,:,:kmax2].flatten(), bins=logbins_pos, facecolor='red', alpha=0.4, label='kmax='+str(kmax2))
#     plt.hist(data[:,:,:kmax2].flatten(), bins=logbins_neg, facecolor='red', alpha=0.4)
#     # plt.xscale('log')
#
#     # max_hist = np.maximum(np.amax(n), max_hist)
#     # print('n', np.amax(n), n.shape)
#     # # n, bins, patches = plt.hist(data, bins, normed=False, facecolor='red', alpha=0.4, label='w')
#     # for i,p in enumerate(perc_range):
#     #     plt.plot([perc1[i], perc1[i]], [0, max_hist], 'k', linewidth=1, label=str(p)+'th percentile')
#     #     plt.plot([perc2[i], perc2[i]], [0, max_hist], 'k--', linewidth=1)
#     # plt.title('all w')
#     # plt.ylabel('n')
#     # plt.legend(loc='upper left')  # , bbox_to_anchor=(0, 0))
#     #
#     # histogram on log scale.
#     # Use non-equal bin sizes, such that they look equal on log scale.
#
#     logbins = np.logspace(np.log10(bins_pos[1]), np.log10(bins_pos[-1]), len(bins_pos))
#     plt.subplot(212)
#     n, bins, patches = plt.hist(data[:,:,:kmax1].flatten(), bins=logbins)
#     max_hist = np.amax(n)
#     plt.hist(data[:,:,:kmax2].flatten(), bins=logbins, facecolor='red', alpha=0.4)
#     plt.xscale('log')
#     max_hist = 4e4
#     for i,p in enumerate(perc_range):
#         if i == 0:
#             plt.plot([perc1[i], perc1[i]], [0, max_hist], 'k', linewidth=1, label = 'level 1 (kmax='+str(kmax1)+')')
#             plt.plot([perc2[i], perc2[i]], [0, max_hist], 'k--', linewidth=1, label = 'level 2 (kmax='+str(kmax2)+')')
#         else:
#             plt.plot([perc1[i], perc1[i]], [0, max_hist], 'k', linewidth=1)
#             plt.plot([perc2[i], perc2[i]], [0, max_hist], 'k--', linewidth=2)
#     plt.legend(loc='upper right')#, bbox_to_anchor=(0, 0))
#     plt.title('w>0')
#     plt.xlabel('w  [m/s]')
#     plt.ylabel('log_10(n)')
#     plt.suptitle('Histogram of w' + '   (t=' + str(t0) + ')')
#     plt.savefig(os.path.join(path_out, 'log_hist_w_comp_t' + str(t0) + '.png'))
#     plt.close()
#
#     return

# ----------------------------------

def fit_func(x, a, b, c):
    return a*np.exp(b * x) + c



def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:,:,:]
    rootgrp.close()
    return data


if __name__ == '__main__':
    main()
