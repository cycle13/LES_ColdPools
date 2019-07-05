import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os
import sys


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35
# plt.rcParams['savefig.edgecolor'] = 'white'


def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')

    nml = set_input_parameters(args)
    # nt = len(times)

    # threshold on vorticity: vort_c
    vort_c = 1.5e-2


    # read in vort-field
    file_name = 'field_vort_yz.nc'
    rootgrp = nc.Dataset(os.path.join(path_in, file_name), 'r')
    grp_descr =rootgrp.groups['description'].variables
    ic_arr = grp_descr['ic_arr'][:]
    jc_arr = grp_descr['jc_arr'][:]

    grp_fields = rootgrp.groups['fields']
    print grp_fields.dimensions.keys()
    nt = grp_fields.dimensions['time'].size
    nx = grp_fields.dimensions['nx'].size
    ny = grp_fields.dimensions['ny'].size
    nz = grp_fields.dimensions['nz'].size
    print('shape fields: ', nx, ny, nz)
    vort_yz = grp_fields.variables['vort_yz'][:,:,:]
    time_range = grp_fields.variables['time'][:]
    rootgrp.close()


    # plot crosssections to see threshold
    print''
    print''
    fig_name = 'vort_threshold.png'
    timesample = np.arange(1, 37, 5)
    print timesample
    krange = [0,1,2,3]
    cm = plt.cm.get_cmap('coolwarm')
    fig, axes = plt.subplots(len(krange), 1, figsize=(10, 10), sharex='all')
    for ik, k0 in enumerate(krange):
        ax = axes[ik]
        ax.plot([0, nx],[0,0], 'k', linewidth=0.5)
        ax.plot([0, nx],[vort_c, vort_c], 'k', linewidth=0.5)
        ax.plot([0, nx],[-vort_c, -vort_c], 'k', linewidth=0.5)
        ax.fill_between([0, nx], -vort_c, vort_c, color='gray', alpha=0.4)
        for count, it in enumerate(timesample):
            count_color = np.double(count)/len(timesample)
            ax.plot(vort_yz[it, :, k0], color=cm(count_color), label='t='+str(time_range[it]))
        ax.set_title('k='+str(k0))
    axes[-1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=8)
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, wspace=0.25)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    fig_name = 'vort_threshold_var.png'
    fig, axes = plt.subplots(len(krange), 1, figsize=(10, 10), sharex='all')
    for ik, k0 in enumerate(krange):
        ax = axes[ik]
        ax.plot([0, nx], [0, 0], 'k', linewidth=0.5)
        ax.plot([0, nx], [vort_c**2, vort_c**2], 'k', linewidth=0.5)
        # ax.plot([0, nx], [-vort_c**2, -vort_c**2], 'k', linewidth=0.5)
        ax.fill_between([0, nx], 0, vort_c**2, color='gray', alpha=0.4)
        for count, it in enumerate(timesample):
            count_color = np.double(count) / len(timesample)
            ax.plot(vort_yz[it, :, k0]**2, color=cm(count_color), label='t=' + str(time_range[it]))
        ax.set_title('k=' + str(k0))
        ax.set_ylim(0,10*vort_c**2)
    axes[-1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=8)
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, wspace=0.25)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    print''
    print''



    ''' FILTER FIELD '''
    import scipy as sp
    import scipy.ndimage
    sigma_x = 1.
    sigma_y = 1.
    n= 200
    m = 1000
    x = np.zeros((n, n))
    i = np.random.choice(range(0, n*n), size=m)
    x[i/n, i%n] = 1.
    x_conv = sp.ndimage.filters.gaussian_filter(x, [sigma_x, sigma_y], mode='constant')

    for it in timesample:
        vort_mask = sp.ndimage.filters.gaussian_filter(vort_yz[it,:,:], [sigma_x, sigma_y])
        fig_name = 'convolution_test_t'+str(it)+'.png'
        fig, axis = plt.subplots(3, 1, figsize=(10, 10), sharex='all')
        cf = axis[0].imshow(vort_yz[it,:,:].T, origin='lower')
        plt.colorbar(cf, ax=axis[0])
        cf = axis[1].imshow(vort_mask.T, origin='lower')
        plt.colorbar(cf, ax=axis[1])
        cf = axis[2].imshow((vort_mask-vort_yz[it,:,:]).T, origin='lower')
        plt.colorbar(cf, ax=axis[2])
        plt.suptitle('sigma_x='+str(sigma_x)+', sigma_y='+str(sigma_y))
        plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, wspace=0.1)
        plt.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)

    for sigma in np.linspace(0.2,1.6,8):
        fig_name = 'convoultion_vort_threshold_sigma'+str(sigma)+'_var.png'
        fig, axes = plt.subplots(len(krange), 1, figsize=(10, 10), sharex='all')
        for ik, k0 in enumerate(krange):
            ax = axes[ik]
            ax.plot([0, nx], [0, 0], 'k', linewidth=0.5)
            ax.plot([0, nx], [vort_c ** 2, vort_c ** 2], 'k', linewidth=0.5)
            ax.fill_between([0, nx], 0, vort_c ** 2, color='gray', alpha=0.4)
            for count, it in enumerate(timesample):
                vort_mask = sp.ndimage.filters.gaussian_filter(vort_yz[it, :, :], [sigma, sigma])
                count_color = np.double(count) / len(timesample)
                ax.plot(vort_mask[:, k0] ** 2, color=cm(count_color), label='t=' + str(time_range[it]))
            ax.set_title('k=' + str(k0))
            ax.set_ylim(0, 10 * vort_c ** 2)
        axes[-1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=8)
        plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, wspace=0.25)
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)



    ''' COMPUTE CIRCULATION '''
    circ = np.zeros(nt)
    area = np.zeros(nt)
    area_ = np.zeros(nt)
    vort_mean = np.zeros(nt)    # mean
    vort_var = np.zeros(nt)     # variance
    vort_max = np.zeros(nt)
    vort_min = np.zeros(nt)
    vort_max_cond = np.zeros(nt)    # max conditioned on vort>vort_c
    vort_mean2 = 0.0            # to compute variance
    for it, t0 in enumerate(time_range):
        print('--- time: t0='+str(t0)+' ---')
        # location of maximum / minimum
        vort_max[it] = np.amax(vort_yz[it,:,:])
        vort_min[it] = np.amin(vort_yz[it,:,:])
        vort_max_i = np.unravel_index(np.argmax(vort_yz[it,:,:]), (ny,nz))
        vort_min_i = np.unravel_index(np.argmin(vort_yz[it,:,:]), (ny,nz))

        jmin = 50
        jmax = 350
        kmax = 20

        # computing circulation
        # !!!!!!!! need to define limits of search area?! (dynamically as it moves with time, maybe read in radius)

        # vort_mask = np.masked_where(...)
        vort_mask = np.zeros((jmax-jmin, kmax), dtype = np.int)
        for j in range(jmin, jmax):
            for k in range(kmax):
                if vort_yz[it, j,k]**2 >= vort_c**2:
                    circ[it] += vort_yz[it, j,k]
                    vort_mask[j-jmin, k] = 1
                    vort_mean[it] += vort_yz[it,j,k]
                    vort_mean2 += vort_yz[it,j,k]*vort_yz[it,j,k]
                    if vort_yz[it, j, k] > vort_max_cond[it]:
                        vort_max_cond[it] = vort_yz[it, j, k]
        area[it] = np.sum(vort_mask)
        if area[it] > 0:
            vort_mean[it] /= area[it]
            vort_mean2 /= area[it]
            vort_var[it] = vort_mean2 - vort_mean[it]*vort_mean[it]
        else:
            vort_mean[it] =  0.0
            vort_var[it] = 0.0
        area[it] = np.sum(vort_mask)*dx[0]*dx[1]
        area_[it] = circ[it] / vort_mean[it]*dx[0]*dx[1]

        # plotting circulation
        lvls = np.linspace(vort_min[it], vort_max[it], 1e3)
        fig, axes = plt.subplots(1, 2, figsize=(15, 4))
        ax = axes[0]
        ax.contourf(vort_yz[it, :, :kmax].T, levels=lvls)
        ax.plot(vort_max_i[0], vort_max_i[1], 'kx', markersize=20)
        ax = axes[1]
        ax.imshow(vort_mask[vort_max_i[0]-50-jmin:vort_max_i[0]+25-jmin,:].T, origin='lower')
        ax.set_title('circulation: C='+str(circ[it]))
        plt.suptitle('t='+str(t0)+'s')
        # plt.tight_layout
        plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, wspace=0.1)
        fig_name = 'vort_max_t' + str(np.int(t0)) + '.png'
        plt.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)


    # plotting time series
    fig_name = 'vort_threshold_vortc' + str(vort_c)+'.png'
    ncol = 4
    fig, axes = plt.subplots(1, ncol, sharey='none', figsize=(5 * ncol, 5))
    count_color = 2 * np.double(it) / len(time_range)
    ax = axes[0]
    ax.plot(time_range, circ, '-o', label='circulation')
    # ax.legend(fontsize=10)
    ax.set_xlabel('time')
    ax.set_title('circulation')
    ax = axes[1]
    ax.plot(time_range, area, '-o', label='area vort>crit')
    ax.plot(time_range, area_, '-o', label='area=circ/<vort>')
    ax.legend(loc='best')
    ax.set_title('area of vort>crit')
    ax.set_xlabel('time')
    ax.set_ylabel('area  [m2]')
    ax = axes[2]
    ax.plot(time_range, np.sqrt(area/(np.pi)), '-o', label='effective radius')
    ax.set_title('effective radius')
    ax.set_ylabel('radius [m]')
    ax.set_xlabel('time')
    ax.set_ylabel('r  [m]')
    ax = axes[3]
    ax.plot(time_range, vort_mean, '-o', label='mean vorticity')
    ax.plot(time_range, vort_max, '-o', label='max vorticity')
    ax.plot(time_range, -vort_min, '-o', label='-min vorticity')
    ax.plot(time_range, vort_max_cond, '-o', label='max vorticity cond')
    ax.plot(time_range, vort_var*1e2, '-o', label='variance vorticity (*1e2)')
    ax.plot([time_range[0],time_range[-1]],[vort_c, vort_c], 'k-')
    plt.legend(loc='best')
    ax.set_title('mean vorticity')
    ax.set_xlabel('time')
    ax.set_ylabel('vort  [1/s]')
    plt.tight_layout()
    plt.suptitle('vort_c = '+str(vort_c)+'1/s')
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)


    file_name = 'Stats_vorticity_rim.nc'
    create_statsfile(time_range, file_name)
    dump_timeseries('area', area, time_range, file_name)
    dump_timeseries('circulation', circ, time_range, file_name)
    dump_timeseries('vort_mean', vort_mean, time_range, file_name)
    dump_timeseries('vort_max', vort_max, time_range, file_name)
    dump_timeseries('vort_min', vort_min, time_range, file_name)

    return

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def create_statsfile(time_range, file_name):
    if os.path.exists(os.path.join(path_out_data, file_name)):
        print''
        print('Stats-file already existing, not dumping')
        # rootgrp = nc.Dataset(os.path.join(path_out_data, file_name), 'r+')    # append data
        print''
        return
    else:
        rootgrp = nc.Dataset(os.path.join(path_out_data, file_name), 'w')
        ts_grp = rootgrp.createGroup('timeseries')
        ts_grp.createDimension('nt', len(time_range))
        var = ts_grp.createVariable('time', 'f8', ('nt'))
        var[:] = time_range[:]
    return


def dump_timeseries(var_name, var, time_range, file_name):
    rootgrp = nc.Dataset(os.path.join(path_out_data, file_name), 'r+')
    ts_grp = rootgrp.groups['timeseries']
    try:
        var_ = ts_grp.createVariable(var_name, 'f8', ('nt'))
    except:
        var_ = ts_grp.variables[var_name][:]
    var_[:] = var[:]
    rootgrp.close()
    return

# def dump_timeseries(circ, area, time_range, file_name):
#     rootgrp = nc.Dataset(os.path.join(path_out_data, file_name), 'r+')
#     ts_grp = rootgrp.groups['timeseries']
#     var = ts_grp.createVariable('circ', 'f8', ('nt'))
#     var[:] = circ[:]
#     var = ts_grp.createVariable('area', 'f8', ('nt'))
#     var[:] = area[:]
#
#     rootgrp.close()
#
#     return

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

def set_input_parameters(args):
    print('--- set input parameters ---')
    global case_name, times
    global path, path_in, path_fields
    global path_out_figs, path_out_data

    path = args.path
    path_in = os.path.join(path, 'fields_vorticity')
    path_fields = os.path.join(path, 'fields')
    path_out_figs = os.path.join(path, 'figs_vorticity')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    path_out_data = os.path.join(path, 'fields_vorticity')
    if not os.path.exists(path_out_data):
        os.mkdir(path_out_data)

    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
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

    ''' determine time range '''
    global tmin, tmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = tmin
    # times = np.arange(tmin, tmax + 100, 100)
    # nt = len(times)
    # print('timerange', times)
    print('tmin ', tmin)
    print('tmax ', tmax)
    print ''

    # return nml, times
    return nml

# ----------------------------------
# ----------------------------------

if __name__ == '__main__':
    main()
