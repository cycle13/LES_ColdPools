
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import netCDF4 as nc
import argparse
import json as simplejson
import os


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    # parser.add_argument("casename")
    # parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    # parser.add_argument("--kmin")
    # parser.add_argument("--kmax")
    args = parser.parse_args()

    global path, path_fields, path_figs
    path = '/nbi/ac/cond2/meyerbe/ColdPools_dry/3D_sfc_fluxes_off/' \
           'triple_3D_noise/old_config/dTh3K_z2000_r2000_triple/'
    path_single = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run4/dTh3_z1000_r1000'
    path_fields = os.path.join(path, 'fields')
    path_figs = '/nbi/ac/cond2/meyerbe/figs_EGU_2019/'
    if not os.path.exists(path_figs):
        os.mkdir(path_figs)

    case_name = 'ColdPoolDry_triple_3D'
    print('casename: ' + case_name)

    ''' determine file range '''
    if args.tmin:
        time_min = np.int(args.tmin)
    else:
        time_min = np.int(100)
    if args.tmax:
        time_max = np.int(args.tmax)
    else:
        time_max = np.int(10000)
    times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
             and np.int(name[:-3]) >= time_min and np.int(name[:-3]) <= time_max]
    times.sort()
    print('times', times)
    files = [str(t) + '.nc' for t in times]
    print('files', files)

    cm_cw = plt.cm.get_cmap('coolwarm')
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_gray = plt.cm.get_cmap('gray')
    cm_inf = plt.cm.get_cmap('gnuplot')
    cm_blues = plt.cm.get_cmap('seismic')
    cm_blues_r = plt.cm.get_cmap('seismic_r')
    cm_spec = plt.cm.get_cmap('spectral')
    # cm_spec = plt.cm.get_cmap('plasma')
    # cm_spec = plt.cm.get_cmap('plasma')

    # # plots contourfigures
    # plt.rcParams['xtick.labelsize'] = 1
    # plt.rcParams['ytick.labelsize'] = 1
    # plt.rcParams['axes.linewidth'] = 2
    # k0 = 0
    # time_range = [0, 900, 1000, 1100, 1200, 1400, 1500, 1600, 2400]
    # for it, t0 in enumerate(time_range):
    #
    #     temp = read_in_netcdf_fields('temperature', os.path.join(path_fields, str(t0) + '.nc'))
    #     s0 = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))
    #     th0 = theta_s(s0)
    #     if it == 0:
    #         # min = np.floor(np.amin(th0[:,:,k0]))
    #         min = np.amin(th0[:,:,k0])
    #         # max = np.ceil(np.amax(th0[:,:,k0]))
    #         max = (np.amax(th0[:,:,k0]))
    #     w0 = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0) + '.nc'))
    #
    #     plt.figure(figsize=(12, 12))
    #     plt.imshow(th0[:, :, k0].T, cmap = cm_cw, vmin = min, vmax = max)#, origin="lower")
    #     # plt.colorbar()
    #     cbar = plt.colorbar(shrink=0.75, ticks=np.arange(np.floor(min), np.ceil(max), 1), aspect=15)
    #     cbar.ax.tick_params(labelsize=18)
    #     # plt.xlabel('Lx', fontsize=24)
    #     # plt.ylabel('Ly', fontsize=24)
    #     plt.tight_layout()
    #     plt.savefig(os.path.join(path_figs, 'xz_plane_th_k' + str(np.int(k0)) + 'm_t'+str(t0)+'.png'))
    #     plt.close()
    #
    #     plt.figure(figsize=(12, 12))
    #     plt.imshow(w0[:, :, k0].T, cmap = cm_bwr, origin="lower")
    #     # plt.xlabel('y')
    #     # plt.ylabel('x')
    #     plt.savefig(os.path.join(path_figs, 'xz_plane_w_k' + str(np.int(k0)) + 'm_t'+str(t0)+'.png'))
    #     plt.close()



    # '''plot profile w_max'''
    # plot_profile_wmax()


    '''plot max(x) crosssections'''
    x_half, y_half, z_half = define_geometry(case_name, files, path)
    krange_ = [1, 4, 10, 15]
    cm = cm_gray
    # cm = cm_blues_r
    # time_bins = [700, 1100, 1500, 1900]
    # plot_xz_crosssections_multilevel(time_bins, 'w', jc_arr[2], krange_, ic_arr, jc_arr, cm,
    #                                  times, path_figs, path_fields)
    '''plot CP height'''
    # time_bins = [time_min, 1100, 1500, time_max]
    # path = os.path.join(path, 'figs_CP_height', 'CP_height_dTh3K_triple_sth0.5.nc')
    # CP_top_file = nc.Dataset(path, 'r')
    # # rootgrp = nc.Dataset(fullpath_in, 'r')
    # CP_top_time= CP_top_file.groups['timeseries'].variables['time'][:]
    # CP_top_2d = CP_top_file.groups['fields_2D'].variables['CP_height_2d'][:,:,:]
    # plot_xz_crosssections_CP_height(CP_top_2d, CP_top_time, time_bins, 'w', jc_arr[2], krange_, ic_arr, jc_arr, cm,
    #                                  times, path_figs, path_fields)
    # CP_top_file.close()

    '''plot xy-field'''
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    plt.rcParams['axes.labelsize'] = 24
    path_out = os.path.join(path_figs, 'xy_snapshots')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    timerange = np.append(np.arange(200, 1000, 200), 1400)
    timerange = [200, 400, 600, 800, 1000, 1400, 1800, 2400]
    k0 = 0
    dx = 25
    var_name = 'temperature'
    for it, t0 in enumerate(timerange):
        var = nc.Dataset(os.path.join(path_single, 'fields', str(t0)+'.nc'), 'r').groups['fields'].variables[var_name][:,:,k0]
        if it == 0:
            min = np.amin(var)
            max = np.amax(var)
        print min, max
        min = 297
        plt.figure(figsize=(12,10))
        # im = plt.imshow(var.T, origin='lower', cmap =cm_inf, vmin=min, vmax=max, extend='both')
        # im = plt.contourf(var.T, levels=np.linspace(min, max, 1e2),
        #                   cmap =cm_gray, vmin=min, vmax=max, extend='both')
        # cbar = plt.colorbar(im, shrink=0.75, ticks=np.arange(297, 300, 1), aspect=15)
        # im = plt.pcolor(var.T, cmap =cm_gray, norm=colors.LogNorm(vmin=min, vmax=max))
        im = plt.pcolor(var.T, norm=colors.colors.PowerNorm(gamma=1./2.), cmap=cm_gray)
        cbar = plt.colorbar(im, shrink=0.75, aspect=15)
        plt.axis('equal')
        cbar.ax.tick_params(labelsize=21)
        plt.xlim(0, 800)
        plt.ylim(0, 800)
        ax = plt.gca()
        plt.locator_params(nbins=8)
        x_ticks = [np.int((ti) * dx * 1e-3) for ti in ax.get_xticks()]
        y_ticks = [np.int((ti) * dx * 1e-3) for ti in ax.get_yticks()]
        ax.set_xticklabels(x_ticks)
        ax.set_yticklabels(y_ticks)
        # plt.tight_layout()
        plt.xlabel('x  [km]')
        plt.xlabel('y  [km]')

        plt.savefig(os.path.join(path_out, str(t0) + 's.png'))
        plt.close()

    return



def plot_xz_crosssections_multilevel(time_bins, var_name, j0, krange, ic_arr, jc_arr, cm,
                                     times, path_out, path_fields):
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    # plt.rcParams['title.fontsize'] = 18
    plt.rcParams['axes.labelsize'] = 24
    plt.rcParams['xtick.direction']='out'
    plt.rcParams['ytick.direction']='out'

    cm_cw = plt.cm.get_cmap('coolwarm')
    cm_gray = plt.cm.get_cmap('gray')
    cm_blue_r = plt.cm.get_cmap('seismic_r')
    cm_blue = plt.cm.get_cmap('seismic')

    ic1 = ic_arr[0]
    ic2 = ic_arr[1]
    jc1 = jc_arr[0]
    jc2 = jc_arr[1]

    nax = 5
    fig, axis = plt.subplots(nax, 1, figsize=(8, 18), sharex='all')
    ax1 = axis[0]
    s0 = read_in_netcdf_fields('s', os.path.join(path_fields, str(times[-1])+'.nc'))
    ax1.imshow(s0[:, :, 1].T, origin='lower', cmap=cm_cw)
    ax1.plot([0, nx], [j0, j0], 'k', linewidth=2)
    ax1.set_xlim(0, nx)
    ax1.set_ylim(0, ny)

    print time_bins
    for it, t0 in enumerate(times):
        print 'it', it
        if t0 > time_bins[0] and t0 <= time_bins[-1]:
            var = read_in_netcdf_fields(var_name, os.path.join(path_fields, str(t0)+'.nc'))

            if t0 <= time_bins[1]:
                cm = cm_gray
                count_color = (np.double(t0)-time_bins[0]) / (time_bins[1]- time_bins[0]) * 0.8
            elif t0 <= time_bins[2]:
                cm = cm_blue_r
                count_color = (np.double(t0)-time_bins[1]) / (time_bins[2]- time_bins[1]) * 0.5 + 0.5
            elif t0 <= time_bins[3]:
                cm = cm_blue_r
                count_color = (np.double(t0)-time_bins[2]) / (time_bins[3]- time_bins[2]) * 0.5

            for i in range(nax-1):
                ax = axis[i+1]
                ax.plot([0, nx], [0, 0], '0.5', linewidth=1)
                ax.plot(var[:, j0, krange[i]], color=cm(count_color), label='t=' + str(t0))
                textstr = 'z=' + str(np.int(dz * krange[i])) + 'm'
                props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='white')
                ax.text(0.73, 0.9, textstr, transform=ax.transAxes, fontsize=24,
                        verticalalignment='top', bbox=props)
                ax.set_ylim([-3,9])
                ax.set_ylabel('w  [m/s]')
                ax.set_xlabel(r'Ly')
            # axis[-1].set_xlabel(r'y  [km]')
            # x_ticks = [np.int(ti * dx * 1e-3) for ti in ax.get_xticks()]
            # ax.set_xticklabels(x_ticks)


    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
                   fancybox=True, shadow=True, ncol=3, fontsize=12)
    plt.subplots_adjust(bottom=0.1, right=.95, left=0.15, top=.99, wspace=0.5, hspace=0.3)
    # plt.subplots_adjust(bottom=0.1, right=.95, left=0.15, top=.99, wspace=0.5, hspace=0.1)
    plt.savefig(os.path.join(path_out, 'xz_multilevel_'+var_name+'.png'))
    plt.close()
    return




def plot_xz_crosssections_CP_height(CP_top, CP_top_time, time_bins, var_name, j0, krange, ic_arr, jc_arr, cm,
                                     times, path_out, path_fields):
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    # plt.rcParams['title.fontsize'] = 18
    plt.rcParams['axes.labelsize'] = 24
    plt.rcParams['xtick.direction']='out'
    plt.rcParams['ytick.direction']='out'

    cm_cw = plt.cm.get_cmap('coolwarm')
    cm_gray = plt.cm.get_cmap('gray')
    cm_blue_r = plt.cm.get_cmap('seismic_r')
    cm_blue = plt.cm.get_cmap('seismic')

    ic1 = ic_arr[0]
    ic2 = ic_arr[1]
    jc1 = jc_arr[0]
    jc2 = jc_arr[1]

    nax = 3
    fig, axis = plt.subplots(nax, 1, figsize=(8, 18), sharex='none')
    ax1 = axis[0]
    ax2 = axis[1]
    s0 = read_in_netcdf_fields('temperature', os.path.join(path_fields, str(times[-1])+'.nc'))
    im = ax1.imshow(s0[:, :, 1].T, origin='lower', cmap=cm_cw, vmin=296.8, vmax=298.6)
    x_ticks = [np.int(ti * dx * 1e-3 + 10) for ti in ax1.get_xticks()]
    ax1.set_xticklabels(x_ticks)
    ax1.set_yticklabels(x_ticks)
    cbar = plt.colorbar(im, ax=ax1, shrink=0.75, ticks=np.arange(297, 299, 1), aspect=15)
    cbar.ax.tick_params(labelsize=21)
    ax1.set_title('temperature', fontsize=24)
    it = np.where(CP_top_time == times[-1])[0][0]
    min = 0
    max = 2200
    im = ax2.imshow(CP_top[it, :, :].T, origin='lower', cmap=cm_cw, vmin=min, vmax=max)
    cbar = plt.colorbar(im, ax=ax2, shrink=0.75, ticks=np.arange(min, max, 500), aspect=15)
    cbar.ax.tick_params(labelsize=21)
    ax2.set_title('CP top', fontsize=24)
    x_ticks = [np.int(ti * dx * 1e-3 + 10) for ti in ax2.get_xticks()]
    ax2.set_xticklabels(x_ticks)
    ax2.set_yticklabels(x_ticks)
    ax1.plot([0, nx], [j0, j0], 'k', linewidth=2)
    ax2.plot([0, nx], [j0, j0], 'k', linewidth=1)
    ax1.set_xlim(0, nx)
    ax1.set_ylim(0, ny)
    ax2.set_xlim(0, nx)
    ax2.set_ylim(0, ny)

    print time_bins
    for it, t0 in enumerate(times):
        it_CP = np.where(CP_top_time == t0)[0][0]
        print 'it', t0, it, it_CP, CP_top_time[it_CP]
        if t0 > time_bins[0] and t0 <= time_bins[-1]:
            # var = read_in_netcdf_fields(var_name, os.path.join(path_fields, str(t0)+'.nc'))

            if t0 <= time_bins[1]:
                cm = cm_gray
                count_color = (np.double(t0)-time_bins[0]) / (time_bins[1]- time_bins[0]) * 0.8
            elif t0 <= time_bins[2]:
                cm = cm_blue_r
                count_color = (np.double(t0)-time_bins[1]) / (time_bins[2]- time_bins[1]) * 0.5 + 0.5
            elif t0 <= time_bins[3]:
                cm = cm_blue
                count_color = (np.double(t0)-time_bins[2]) / (time_bins[3]- time_bins[2]) * 0.45 + 0.5

            for i in range(nax-2):
                ax = axis[i+2]
                ax.plot([0, nx], [0, 0], '0.5', linewidth=1)
                ax.plot(CP_top[it_CP, :, j0], color=cm(count_color), label='t=' + str(t0))
                ax.set_ylim([-0.,2300])
                ax.set_ylabel('CP top  [m]')
                ax.set_xlabel(r'Ly')
            axis[-1].set_xlabel(r'y  [km]')
            x_ticks = [np.int(ti * dx * 1e-3) for ti in ax.get_xticks()]
            ax.set_xticklabels(x_ticks)


    plt.legend(loc='upper center', bbox_to_anchor=(0.4, -0.2),
                   fancybox=True, shadow=True, ncol=5, fontsize=12)
    plt.subplots_adjust(bottom=0.1, right=.95, left=0.15, top=.99, wspace=0.5, hspace=0.3)
    # plt.subplots_adjust(bottom=0.1, right=.95, left=0.15, top=.99, wspace=0.5, hspace=0.1)
    plt.savefig(os.path.join(path_out, 'CP_top.png'))
    plt.close()
    return



def plot_profile_wmax():
    # plot profile w_max
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    plt.rcParams['axes.labelsize'] = 28
    plt.rcParams['lines.linewidth'] = 3
    # plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    single_3K = [1.4, 2.3, 2.8, 3.1, 2.9, 2.4, 2.1, 1.8, 1.6, 1.5, 1.4,
                 1.3, 1.1, 1.2, 1., 0.8, 0.8, 0.8, 0.9, 0.9, 0.8]
    # double_3K = [3,5.3,7,8.1,8.6,8.5,8,7.5,6.7,6,5,
    #              3.8,2.7,2,2,1.8,
    #              1.2,1,1.1,1,0.6]
    # triple_3K = [2.3,4.1,5.9,6.88,7.5,8,8.1,7.7,6.9,7.1,7.4,
    #              7.6,7.6,7.1,5.9,5.7,
    #              5.7,5.3,4.5,3.4,2.3]
    double_3K = [3, 5.3, 7, 8.1, 8.6, 8.5, 7.9, 7.4, 6.8, 5.9, 4.9,
                 3.8, 2.6, 2.2, 1.9, 1.7,
                 1.3, 1.2, 1.2, 1.1, 0.9]
    triple_3K = [2.5, 4.5, 5.8, 6.8, 7.5, 8, 8.1, 7.7, 6.9, 7.1, 7.4, 7.6,
                 7.6, 7.1, 5.9, 5.7, 5.7,
                 5.3, 4.5, 3.4, 2.3]
    height_range = np.arange(0, 2100, 100)
    c1 = 'navy'
    c2 = 'gray'
    c3 = 'coral'
    # path = '/nbi/ac/cond2/meyerbe/ColdPools_dry/3D_sfc_fluxes_off/triple_3D_noise/dTh3K_triple'
    plt.figure(figsize=(6, 9))
    plt.plot(single_3K, height_range, '-o', label='undisturbed front', color=c2)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22),
               fancybox=True, shadow=False, ncol=1, fontsize=18)
    plt.ylabel('height z  [m]')
    plt.xlabel('max(w)  [m/s]')
    plt.xlim([0, 9])
    plt.subplots_adjust(bottom=0.1, right=.95, left=0.3, top=.84, wspace=0.25)
    plt.savefig(os.path.join(path_figs, 'w_max_profile_single.png'))
    # plt.show()
    #
    plt.figure(figsize=(6, 9))
    plt.plot(single_3K, height_range, '-o', label='undisturbed front', color=c2)
    plt.plot(double_3K, height_range, '-o', label='2 CP collison', color=c1)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22),
               fancybox=True, shadow=False, ncol=1, fontsize=18)
    plt.ylabel('height z  [m]')
    plt.xlabel('max(w)  [m/s]')
    plt.xlim([0, 9])
    # plt.subplots_adjust(bottom=0.2, right=.95, left=0.2, top=.95, wspace=0.25)
    plt.subplots_adjust(bottom=0.1, right=.95, left=0.3, top=.84, wspace=0.25)
    plt.savefig(os.path.join(path_figs, 'w_max_profile_single_double.png'))
    # plt.show()

    plt.figure(figsize=(6, 10))
    plt.plot(single_3K, height_range, '-o', label='undisturbed front', color=c2)
    plt.plot(double_3K, height_range, '-o', label='2 CP collison', color=c1)
    plt.plot(triple_3K, height_range, '-o', label='3 CP collison', color=c3)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13),
               fancybox=True, shadow=False, ncol=1, fontsize=18)
    plt.ylabel('height z  [m]')
    plt.xlabel('max(w)  [m/s]')
    plt.xlim([0, 9])
    plt.subplots_adjust(bottom=0.23, right=.95, left=0.3, top=.95, wspace=0.25)
    # plt.subplots_adjust(bottom=0.2, right=.95, left=0.2, top=.95, wspace=0.25)
    plt.savefig(os.path.join(path_figs, 'w_max_profile_single_double_triple.png'))
    # plt.show()
    #
    plt.figure(figsize=(6, 9))
    plt.plot(single_3K, height_range, '-o', label='undisturbed front', color=c2)
    plt.plot(double_3K, height_range, '-o', label='2 CP collison', color=c1)
    plt.plot(triple_3K, height_range, '-o', label='3 CP collison', color=c3)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22),
               fancybox=True, shadow=False, ncol=1, fontsize=18)
    plt.ylabel('height z  [m]')
    plt.xlabel('max(w)  [m/s]')
    plt.xlim([0, 9])
    plt.subplots_adjust(bottom=0.1, right=.95, left=0.3, top=.84, wspace=0.25)
    # plt.subplots_adjust(bottom=0.2, right=.95, left=0.2, top=.95, wspace=0.25)
    plt.savefig(os.path.join(path_figs, 'w_max_profile_single_double_triple.png'))
    # plt.show()

    return





# ----------------------------------
# ----------------------------------
def define_geometry(case_name, files, path):
    print 'define geometry'
    global nx, ny, nz, dx, dy, dz, gw
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    gw = nml['grid']['gw']

    # set coordinates for plots
    # (a) double 3D
    global ic_arr, jc_arr
    global isep

    if case_name == 'ColdPoolDry_single_3D':
        rstar = nml['init']['r']
        irstar = np.int(np.round(rstar / dx))
        zstar = nml['init']['h']
        try:
            ic = nml['init']['ic']
            jc = nml['init']['jc']
        except:
            ic = np.int(nx/2)
            jc = np.int(ny/2)
        ic_arr = np.zeros(1)
        jc_arr = np.zeros(1)
        ic_arr[0] = ic
        jc_arr[0] = jc
    # (b) double 2D
    elif case_name == 'ColdPoolDry_double_2D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        isep = 4 * irstar
        ic1 = np.int(nx / 3)  # np.int(Gr.dims.ng[0] / 3)
        ic2 = ic1 + isep
        jc1 = np.int(ny / 2)
        jc2 = jc1
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_double_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        # zstar = nml['init']['h']
        isep = 4 * irstar
        jsep = 0
        # ic1 = np.int(np.round((nx + 2 * gw) / 3)) - gw
        ic1 = np.int(np.round((nx - gw) / 3)) + 1
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

    print('ic, jc: ', ic_arr, jc_arr)

    ''' --- auxiliary arrays (since no Grid.pyx) ---'''
    global nx_, ny_, nz_
    # test file:
    var = read_in_netcdf_fields('u', os.path.join(path_fields, files[0]))
    [nx_, ny_, nz_] = var.shape

    x_half = np.empty((nx_), dtype=np.double, order='c')
    y_half = np.empty((ny_), dtype=np.double, order='c')
    z_half = np.empty((nz_), dtype=np.double, order='c')
    count = 0
    for i in xrange(nx_):
        x_half[count] = (i + 0.5) * dx
        count += 1
    count = 0
    for j in xrange(ny_):
        y_half[count] = (j + 0.5) * dy
        count += 1
    count = 0
    for i in xrange(nz_):
        z_half[count] = (i + 0.5) * dz
        count += 1

    return x_half, y_half, z_half
# ----------------------------------
def theta_s(s):
    # parameters from pycles
    T_tilde = 298.15
    sd_tilde = 6864.8
    cpd = 1004.0
    th_s = T_tilde * np.exp( (s - sd_tilde)/cpd )
    return th_s

# ----------------------------------
def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    # shape = var.shape
    # data = np.ndarray(shape = var.shape)
    data = var[:,:,:]
    rootgrp.close()
    return data

if __name__ == '__main__':
    main()

