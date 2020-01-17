import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os
import matplotlib.patches as mpatches



label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35
plt.rcParams['text.usetex'] = 'true'

execfile('settings.py')


def main():
    path_single = '/nbi/ac/conv4/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/single_3D_noise/run5_PE_scaling_dx100m'
    path_single_dx50 = '/nbi/ac/conv4/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/single_3D_noise/run6_PE_scaling_dx50m'
    path_double = '/nbi/ac/coag/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/double_3D'
    path_triple = '/nbi/ac/cond2/meyerbe/ColdPools/3D_sfc_fluxes_off/triple_3D/'
    case_name_single = 'ColdPoolDry_single_3D'
    case_name_double = 'ColdPoolDry_double_3D'
    case_name_triple = 'ColdPoolDry_triple_3D'

    # path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single/'
    path_out_figs = '/nbi/home/meyerbe/paper_CP_single/'
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)

    print('path figs: ' + path_out_figs)
    print('')

    dTh = 5
    zstar = 1000
    # rstar_range = [1100, 2000]
    rstar = 1100
    # rstar = 2000
    print('zstar: '+str(zstar))
    print('rstar: '+str(rstar))
    print('')
    if rstar == 1100:
        d_range = [10, 12, 15]
        t_2CP = [1100, 1500, 2300]
        t_3CP = [1500, 2200, 3300]
    elif rstar == 2000:
        d_range = [10, 15, 20]
    # d_range = [15, 20]
    id_list_s = ['dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)]
    id_list_d = []
    id_list_t = []
    for dstar in d_range:
        id_list_d.append('dTh'+str(dTh)+'_z'+str(zstar)+'_r'+str(rstar)+'_d'+str(dstar) + 'km')
        id_list_t.append('dTh'+str(dTh)+'_z'+str(zstar)+'_r'+str(rstar)+'_d'+str(dstar) + 'km')
    ncases = len(id_list_t)
    print('ID-list: ', id_list_d)
    print('ID-list: ', id_list_t)
    print('')

    define_geometry(case_name_single, case_name_double, case_name_triple,
                    path_single, path_double, path_triple,
                    id_list_s, id_list_d, id_list_t)

    ''' determine max. height '''
    zmax = 5000.
    kmax = np.int(zmax/dx[2])
    tmax_s = 1200.
    itmax_s = np.int(tmax_s/dt_fields)
    print('zmax: '+str(zmax))
    print('kmax: '+str(kmax))
    print('tmax s: '+str(tmax_s))
    print('')

    ''' determine time range '''
    global tmin, tmax
    tmin = np.int(0)
    tmax = np.int(7200)
    tmax = np.int(3600)
    tmax = np.int(3500)  # bcs. of tracer statistics
    times = np.arange(tmin, tmax + dt_fields, dt_fields)
    nt = len(times)
    print('tmin: ' + str(tmin))
    print('tmax: ' + str(tmax))
    print('nt  : ' + str(nt))
    # print('times: ' + str(times))
    print('')


    ''' determine sampling subdomain '''
    # for single CP: search maximum within circle of radius CP has at 2-CP / 3-CP collision time
    # tracer statistics
    k0 = 0
    # times = np.arange(0, 3600, 100)
    # nt = len(times)
    fullpath_in = os.path.join(path_single, id_list_s[0], 'tracer_k' + str(k0), 'output')
    n_tracers = get_number_tracers(fullpath_in)
    n_cps = get_number_cps(fullpath_in)
    dist_av = np.zeros((nt))
    for it, t0 in enumerate(times):
        print('it, t0', it, t0)
        cp_id = 2
        dist_av[it] = get_radius(fullpath_in, it, cp_id, n_tracers, n_cps)
    r_av = dist_av * dx[0]
    rad_2CP = np.empty(3)
    rad_3CP = np.empty(3)
    delta_s = 6.e2/dx[0]
    for d in range(len(d_range)):
        rad_2CP[d] = r_av[np.int(t_2CP[d]/dt_fields)]+delta_s
        rad_3CP[d] = r_av[np.int(t_3CP[d]/dt_fields)]+delta_s
    [xs,ys] = nx_s[:2]*.5

    delta_d = np.asarray([2.e3/dx[0],6.e2/dx[1]])
    [xd,yd] = nx_d[:2]*.5-delta_d*.5
    rect_double = mpatches.Rectangle((xd, yd), delta_d[0], delta_d[1], linewidth=1, edgecolor='grey',
                                     facecolor='none')
    delta_t = 6.e2/dx[0]
    [xt, yt] = nx_t[:2]*.5-delta_t*.5
    rect_triple = mpatches.Rectangle((xt, yt), delta_t, delta_t, linewidth=1, edgecolor='grey',
                                     facecolor='none')



    ''' testing '''
    k0 = 0
    if k0 == 0:
        w_max = 3.
    elif k0 == 10:
        w_max = 1.
    # lvls = np.arange(-w_max, w_max, .25)
    lvls = np.linspace(-4, 4, 10)
    for d, dstar in enumerate(d_range):
        fig_name = 'collisions_subdomains_k'+str(k0)+'_d'+str(dstar)+ 'km_test.png'
        fig, axis = plt.subplots(3, 2, figsize=(11, 15))
        fullpath_in = os.path.join(path_single, id_list_s[0], 'fields', str(t_2CP[d])+'.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:,:,k0]
        root.close()
        cf = axis[0,0].contourf(w.T, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[0,0], shrink=0.8)
        fullpath_in = os.path.join(path_single, id_list_s[0], 'fields', str(t_3CP[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        cf = axis[0,1].contourf(w.T, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[0, 1], shrink=0.8)
        fullpath_in = os.path.join(path_double, id_list_d[d], 'fields', str(t_2CP[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        # [nx_d,ny_d] = w.shape
        root.close()
        cf = axis[1,0].contourf(w, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[1, 0], shrink=0.8)
        fullpath_in = os.path.join(path_double, id_list_d[d], 'fields', str(t_3CP[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        cf = axis[1,1].contourf(w, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[1, 1], shrink=0.8)
        fullpath_in = os.path.join(path_triple, id_list_t[d], 'fields', str(t_2CP[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        # [nx_t, ny_t] = w.shape
        cf = axis[2,0].contourf(w.T, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[2, 0], shrink=0.8)
        fullpath_in = os.path.join(path_triple, id_list_t[d], 'fields', str(t_3CP[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        cf = axis[2,1].contourf(w.T, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[2, 1], shrink=0.8)

        circle1 = plt.Circle((xs, ys), rad_2CP[d], fill=False, color='lime', linewidth=1)
        circle2 = plt.Circle((xs, ys), rad_3CP[d], fill=False, color='lime', linewidth=1)
        axis[0,0].add_artist(circle1)
        axis[0,1].add_artist(circle2)
        # ic = np.int(nx_d[0]*.5)
        # jc = np.int(nx_d[1]*.5)
        # di = 20
        # dj = 5
        # rect_double = mpatches.Rectangle((ic - di, jc - dj), 2 * di, 2 * dj, linewidth=1, edgecolor='grtey', facecolor='none')
        # axis[1,0].plot(ic, jc, 'o', markersize=20)
        axis[1,0].plot(xd,yd, 'o', markersize=20)
        print('HAAAAAAAAAAAAAAAAAA', xd, yd)
        # ic = np.int(nx_t*0.5)
        # jc = np.int(ny_t*0.5)
        # di = 20
        # rect_triple = mpatches.Rectangle((ic - di, jc - di), 2 * di, 2* di, linewidth=1, edgecolor='grey', facecolor='none')
        axis[1,0].add_patch(rect_double)
        axis[2,0].add_patch(rect_triple)

        axis[0,0].set_title('t='+str(t_2CP[d])+'s')
        axis[0,1].set_title('t='+str(t_3CP[d])+'s')
        for ax in axis[1,:].flat:
            ax.set_xlim(100,nx_d[1]-100)
            ax.set_ylim(200,nx_d[0]-200)
        for ax in axis[2, :].flat:
            if d == 0:
                delta = 80
            else:
                delta = 200
            ax.set_xlim(delta, nx_t[0]-delta)
            ax.set_ylim(delta, nx_t[1]-delta)
        for ax in axis.flat:
            ax.set_aspect('equal')
        plt.subplots_adjust(bottom=0.05, right=.95, left=0.05, top=0.95, hspace=0.05)
        plt.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)






    # ''' (A) plot from local (unaveraged) min / max'''
    # filename = 'minmax_domain_noaverage.nc'
    #
    # ''' compute domain min/max values '''
    # # w_min_s, w_max_s, th_min_s, th_max_s, z, z_half \
    # #         = compute_domain_max(path_single, id_list_s[0], case_name_single, kmax, times, nt)
    # # path_out = os.path.join(path_single, id_list_s[0], 'data_analysis')
    # # dump_minmax_file(w_min_s, w_max_s, th_min_s, th_max_s, z, z_half, kmax, times, filename, path_out)
    # # w_min_s, w_max_s, th_min_s, th_max_s, z, z_half \
    # #         = compute_domain_max(path_single_dx50, id_list_s[0], case_name_single, kmax, times, nt)
    # # path_out = os.path.join(path_single_dx50, id_list_s[0], 'data_analysis')
    # # dump_minmax_file(w_min_s, w_max_s, th_min_s, th_max_s, z, z_half, kmax, times, filename, path_out)
    # for d,dstar in enumerate(d_range):
    # # for d,dstar in enumerate([12]):
    #     w_min_d, w_max_d, th_min_d, th_max_d, z, z_half \
    #         = compute_domain_max(path_double, id_list_d[d], case_name_double, kmax, times, nt)
    #     path_out = os.path.join(path_double, id_list_d[d], 'data_analysis')
    #     if not os.path.exists(path_out):
    #         os.mkdir(path_out)
    #     dump_minmax_file(w_min_d, w_max_d, th_min_d, th_max_d, z, z_half, kmax, times, filename, path_out)
    #     path_out = os.path.join(path_triple, id_list_t[d], 'data_analysis')
    #     if not os.path.exists(path_out):
    #         os.mkdir(path_out)
    #     w_min_t, w_max_t, th_min_t, th_max_t, z, z_half \
    #         = compute_domain_max(path_triple, id_list_t[d], case_name_triple, kmax, times, nt)
    #     dump_minmax_file(w_min_t, w_max_t, th_min_t, th_max_t, z, z_half, kmax, times, filename, path_out)
    #
    # ''' read in min/max values '''
    # fig_name = 'collisions_minmax_profiles_domain_unaveraged_rstar'+str(rstar)+'.png'
    # zmax_plot = 3000.
    # kmax_plot = np.int(zmax_plot/dx[2])
    # fig, axis = plt.subplots(1, 3, figsize=(14, 10))
    # ax0 = axis[0]
    # ax1 = axis[1]
    # ax2 = axis[2]
    # for d,dstar in enumerate(d_range):
    # # for d,dstar in enumerate([10, 15]):
    #     print('.... d: '+str(dstar))
    #     al = 1.-d*1./len(d_range)
    #     path = os.path.join(path_single, id_list_s[0], 'data_analysis')
    #     print(path)
    #     w_max_s, th_min_s, z, z_half = read_in_minmax(kmax_plot, path, filename)
    #     path = os.path.join(path_double, id_list_d[d], 'data_analysis')
    #     print(path)
    #     w_max_d, th_min_d, z, z_half = read_in_minmax(kmax_plot, path, filename)
    #     path = os.path.join(path_triple, id_list_t[d], 'data_analysis')
    #     print(path)
    #     w_max_t, th_min_t, z, z_half = read_in_minmax(kmax_plot, path, filename)
    #     if d > 0:
    #         lbl_s = ''
    #         lbl_d = ''
    #         lbl_t = ''
    #     else:
    #         lbl_s = 'single CP gust front'
    #         lbl_d = 'double CP collision'
    #         lbl_t = 'triple CP collision'
    #     ax0.plot(np.amax(w_max_s[:itmax_s,:], axis=0), z, color=colorlist3[0], alpha=al, label=lbl_s)
    #     ax0.plot(np.amax(w_max_d, axis=0), z, color=colorlist3[1], alpha=al, label=lbl_d)
    #     ax0.plot(np.amax(w_max_t, axis=0), z, color=colorlist3[2], alpha=al, label=lbl_t)
    #     ax1.plot(np.amin(th_min_s, axis=0), z_half, color=colorlist3[0], alpha=al, label=lbl_s)
    #     ax1.plot(np.amin(th_min_d, axis=0), z_half, color=colorlist3[1], alpha=al, label=lbl_d)
    #     ax1.plot(np.amin(th_min_t, axis=0), z_half, color=colorlist3[2], alpha=al, label=lbl_t)
    #
    # ax0.set_xlabel('max. w  [m/s]')
    # ax0.set_ylabel('height z  [m]')
    # ax1.set_xlabel(r'min. $\theta$ [K]')
    # # # ax2.grid()
    # # # ax2.set_xlim(0.9, 2.1)
    # # # ax0.set_title('maxima')
    # # # ax1.set_title('minima')
    # ax0.legend(loc='upper left', bbox_to_anchor=(0.35, .95),
    #            fancybox=False, shadow=False, ncol=1, fontsize=12)
    # # # ax2.legend(loc='upper left', bbox_to_anchor=(0.1, -0.1),
    # # #            fancybox=False, shadow=False, ncol=1, fontsize=9)
    # # # plt.suptitle('min/max for ' + var_name, fontsize=21)
    # plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.9, hspace=0.4)
    # plt.savefig(os.path.join(path_out_figs, fig_name))
    # plt.close(fig)


    #
    #
    # ''' (B) plot from averaged min / max'''
    # # - single - azimuthal average
    # # - double - average along collision line
    # # - triple - average about a few grid points at center
    #
    #
    # # # fig_name = 'collisions_minmax_profiles_domain.png'
    # # # zrange = np.arange(nz) * dx[2]
    # # # fig, axis = plt.subplots(1, 3, figsize=(14, 5))
    # # # ax0 = axis[0]
    # # # ax1 = axis[1]
    # # # for i_id, ID in enumerate(id_list):
    # # #     al = np.double(i_id + 1) / (ncases + 1)
    # # #     ax0.plot(ts_max_w[i_id, :nz], zrange, color='b', alpha=al, label='single, ' + ID)
    # # #     ax1.plot(ts_min_th[i_id, :nz], zrange, color='b', alpha=al, label='single, ' + ID)
    # # # ax0.legend(loc='upper left', bbox_to_anchor=(-0.1, -0.1),
    # # #            fancybox=False, shadow=False, ncol=3, fontsize=9)
    # # # # ax1.legend(loc='upper left', bbox_to_anchor=(0.1, -0.1),
    # # # #            fancybox=False, shadow=False, ncol=1, fontsize=9)
    # # # ax0.set_ylabel('height [km]')
    # # # ax0.set_xlabel('max(w) [m/s]')
    # # # ax1.set_xlabel(r'min(potential temperature  ) [K]')
    # # # for ax in axis.flatten():
    # # #     y_ticks = [np.int(ti * 1e-3) for ti in ax.get_yticks()]
    # # #     ax.set_yticklabels(y_ticks)
    # # # plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.9, hspace=0.4)
    # # # plt.savefig(os.path.join(path_out_figs, fig_name))
    # # # plt.close()
    # # #
    # # # time_bins = np.ndarray((ncases,3), dtype=np.double)
    # # # time_bins[0,:] = [1000, 1400, times[-1]]
    # # # time_bins[1,:] = [4500, 5900, times[-1]]
    # # # # time_bins[1,:] = [2300, 3200, times[-1]]
    # # # # time_bins[2,:] = [4500, 5900, times[-1]]
    # # # time_bins_ = {}
    # # # time_bins_['d10km'] = [1000, 1400, times[-1]]
    # # # time_bins_['d15km'] = [2300, 3200, times[-1]]
    # # # time_bins_['d20km'] = [4500, 5900, times[-1]]
    # # #
    # # # fig_name = 'collisions_minmax_profiles.png'
    # # # zrange = np.arange(nz) * dx[2]
    # # # fig, axis = plt.subplots(1, 3, figsize=(14, 5))
    # # # ax0 = axis[0]
    # # # ax1 = axis[1]
    # # # for i_id, ID in enumerate(id_list):
    # # #     print(ID[-5:])
    # # #     # it_s = np.int(times_array_w[i_id, 0]/dt_fields)
    # # #     # it_d = np.int(times_array_w[i_id, 1]/dt_fields)
    # # #     # it_s = np.int(time_bins[i_id, 0]/dt_fields)
    # # #     # it_d = np.int(time_bins[i_id, 1]/dt_fields)
    # # #     it_s = np.int(time_bins_[ID[-5:]][0]/dt_fields)
    # # #     it_d = np.int(time_bins_[ID[-5:]][1]/dt_fields)
    # # #     # it_t = times_array_w[i_id, 2]
    # # #     print('it single, double: ', it_s, it_d)
    # # #     al = np.double(i_id + 1) / (ncases + 1)
    # # #     aux = np.amax(prof_max_w[i_id, :it_s, :nz], axis=0)
    # # #     print('', prof_max_w.shape, aux.shape, nz, zrange.shape)
    # # #     ax0.plot(np.amax(prof_max_w[i_id, :it_s, :nz], axis=0), zrange, color='b', alpha=al, label='single, ' + ID)
    # # #     # ax0.plot(prof_max_w[i_id, :nz], zrange, '-', color='g', alpha=al, label='double')
    # # #     # ax0.plot(prof_max_w[i_id, :nz], zrange, '-', color='r', alpha=al, label='triple')
    # # # #     ax1.plot(prof_min_th[i_id, :nz], zrange, color='b', alpha=al, label='single, ' + ID)
    # # # #     # ax1.plot(prof_min_th[i_id, :nz], zrange, '-', color='g', alpha=al, label='double')
    # # # #     # ax1.plot(prof_min_th[i_id, :nz], zrange, '-', color='r', alpha=al, label='triple')
    # # # # ax0.legend(loc='upper left', bbox_to_anchor=(-0.1, -0.1),
    # # # #            fancybox=False, shadow=False, ncol=3, fontsize=9)
    # # # # # ax1.legend(loc='upper left', bbox_to_anchor=(0.1, -0.1),
    # # # # #            fancybox=False, shadow=False, ncol=1, fontsize=9)
    # # # ax0.set_ylabel('height [km]')
    # # # ax0.set_xlabel('max(w) [m/s]')
    # # # ax1.set_xlabel(r'min(potential temperature  ) [K]')
    # # # # ax1.set_xlim(250,350)
    # # # for ax in axis.flatten():
    # # #     y_ticks = [np.int(ti * 1e-3) for ti in ax.get_yticks()]
    # # #     ax.set_yticklabels(y_ticks)
    # # # plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.9, hspace=0.4)
    # # # plt.savefig(os.path.join(path_out_figs, fig_name))
    # # # plt.close()

    return


# --------------------------------------------------------------------
def compute_domain_max(path, ID, casename, kmax, times, nt):
    w_min = 999.9 * np.ones((nt, kmax), dtype=np.double)
    w_max = 999.9 * np.ones((nt, kmax), dtype=np.double)
    th_min = 999.9 * np.ones((nt, kmax), dtype=np.double)
    th_max = 999.9 * np.ones((nt, kmax), dtype=np.double)

    root = nc.Dataset(os.path.join(path, ID, 'stats', 'Stats.'+casename+'.nc'))
    zrange = root.groups['profiles'].variables['z'][:kmax]
    zrange_half = root.groups['profiles'].variables['z_half'][:kmax]
    root.close()


    for it, t0 in enumerate(times):
        fullpath_in = os.path.join(path, ID, 'fields', str(t0) + '.nc')
        print(fullpath_in)
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, :kmax]
        w_min[it] = np.amin(np.amin(w, axis=0), axis=0)
        w_max[it] = np.amax(np.amax(w, axis=0), axis=0)
        del w
        s = root.groups['fields'].variables['s'][:, :, :kmax]
        s_min = np.amin(np.amin(s, axis=0), axis=0)
        s_max = np.amax(np.amax(s, axis=0), axis=0)
        del s
        th_min[it] = thetas_c(s_min, 0)
        th_max[it] = thetas_c(s_max, 0)
        root.close()
        del s_min, s_max
    return w_min, w_max, th_min, th_max, zrange, zrange_half


def dump_minmax_file(w_min, w_max, th_min, th_max, z, z_half, kmax, times, filename, path_out):
    print('dumping: ', filename, path_out)
    # output for each CP:
    # - min, max (timeseries)
    # - CP height (field; max=timeseries)
    # - (ok) CP rim (field)
    nt = len(times)
    print('create output file: ', os.path.join(path_out, filename))
    print('size: ', nz, nt)

    rootgrp = nc.Dataset(os.path.join(path_out, filename), 'w', format='NETCDF4')

    ts_grp = rootgrp.createGroup('timeseries')
    ts_grp.createDimension('nt', nt)
    var = ts_grp.createVariable('time', 'f8', ('nt'))
    var.units = "s"
    var[:] = times
    ts_grp.createDimension('nz', kmax)
    var = ts_grp.createVariable('z', 'f8', ('nz'))
    var.units = "m"
    var[:] = z
    var = ts_grp.createVariable('z_half', 'f8', ('nz'))
    var.units = "m"
    var[:] = z_half
    #
    var = ts_grp.createVariable('th_max', 'f8', ('nt', 'nz'))
    var.long_name = 'domain max potential temperature'
    var.units = "K"
    var[:,:] = th_max[:,:]
    var = ts_grp.createVariable('th_min', 'f8', ('nt', 'nz'))
    var.long_name = 'domain min potential temperature'
    var.units = "K"
    var[:,:] = th_min[:,:]

    var = ts_grp.createVariable('w_max', 'f8', ('nt', 'nz'))
    var.long_name = 'domain max vertical velocity'
    var.units = "m/s"
    var[:,:] = w_max[:,:]
    var = ts_grp.createVariable('w_min', 'f8', ('nt', 'nz'))
    var.long_name = 'domain min vertical velocity'
    var.units = "m/s"
    var[:,:] = w_min[:,:]

    # field_grp = rootgrp.createGroup('fields_2D')
    # field_grp.createDimension('nt', nt)
    # field_grp.createDimension('nx', nx)
    # field_grp.createDimension('ny', ny)
    # var = field_grp.createVariable('CP_height_2d', 'f8', ('nt', 'nx', 'ny'))
    # var.units = "m"
    # var = field_grp.createVariable('CP_height_gradient_2d', 'f8', ('nt', 'nx', 'ny'))
    # var.units = "m"
    # var = field_grp.createVariable('w_max_2d', 'f8', ('nt', 'nx', 'ny'))
    # var.units = "m/s"
    # var = field_grp.createVariable('w_max_height_2d', 'f8', ('nt', 'nx', 'ny'))
    # var.units = "m"

    rootgrp.close()
    print ''
    return

def read_in_minmax(kmax, path_out, filename):
    root = nc.Dataset(os.path.join(path_out, filename), 'r')
    z = root.groups['timeseries'].variables['z'][:kmax]
    z_half = root.groups['timeseries'].variables['z_half'][:kmax]
    w_max = root.groups['timeseries'].variables['w_max'][:,:kmax]
    th_min = root.groups['timeseries'].variables['th_min'][:,:kmax]
    root.close()
    return w_max, th_min, z, z_half
# --------------------------------------------------------------------
# --------------------------------------------------------------------

# ---------------------------- TRACER STATISTICS -----------------------

def get_number_tracers(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in + '/coldpool_tracer_out.txt', 'r')
    lines = f.readlines()
    count = 0
    # while CP age is 0 and CP ID is 1
    cp_age = int(lines[count].split()[0])
    cp_ID = int(lines[count].split()[3])
    print('cp_age', cp_age)
    while (cp_age == 1 and cp_ID == 1):
        count += 1
        cp_age = int(lines[count].split()[0])
        cp_ID = int(lines[count].split()[3])
    n_tracers = count
    f.close()

    return n_tracers



def get_number_cps(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in + '/coldpool_tracer_out.txt', 'r')
    lines = f.readlines()
    cp_number = int(lines[-1].split()[3])
    f.close()

    return cp_number



def get_radius(fullpath_in, t0, cp_id, n_tracers, n_cps):
    f = open(fullpath_in + '/coldpool_tracer_out.txt', 'r')
    # f = open(DIR+EXPID+'/'+child+'/output/irt_tracks_output_pure_sort.txt', 'r')
    lines = f.readlines()
    dist = []

    count = t0 * n_cps * n_tracers + (cp_id - 1)*n_tracers
    # while CP age is 0 and CP ID is cp_id
    timestep = int(lines[count].split()[0])
    cp_ID = int(lines[count].split()[3])
    # print(timestep, cp_ID)
    while (timestep-1 == t0 and int(lines[count].split()[3])==cp_id):
        columns = lines[count].split()
        dist.append(float(columns[8]))
        count += 1
        timestep = int(lines[count].split()[0])
    f.close()
    r_av = np.average(dist)

    return r_av

# ----------------------------------------------------------------------
def define_geometry(case_name_single, case_name_double, case_name_triple,
                    path_single, path_double, path_triple, id_list_s, id_list_d, id_list_t):
    print 'define geometry'


    global nx_s, nx_d, nx_t
    nx_s = np.empty(3, dtype=np.int)
    nx_d = np.empty(3, dtype=np.int)
    nx_t = np.empty(3, dtype=np.int)
    dx_s = np.empty(3, dtype=np.int)
    dx_d = np.empty(3, dtype=np.int)
    dx_t = np.empty(3, dtype=np.int)

    nml = simplejson.loads(open(os.path.join(path_single, id_list_s[0], case_name_single+ '.in')).read())
    nx_s[0] = nml['grid']['nx']
    nx_s[1] = nml['grid']['ny']
    nx_s[2] = nml['grid']['nz']
    # dx_s[0] = nml['grid']['dx']
    # dx_s[1] = nml['grid']['dy']
    # dx_s[2] = nml['grid']['dz']
    # dt_fields_s = np.int(nml['fields_io']['frequency'])
    nml = simplejson.loads(open(os.path.join(path_double, id_list_d[0], case_name_double+ '.in')).read())
    nx_d[0] = nml['grid']['nx']
    nx_d[1] = nml['grid']['ny']
    nx_d[2] = nml['grid']['nz']
    dx_d[0] = nml['grid']['dx']
    dx_d[1] = nml['grid']['dy']
    dx_d[2] = nml['grid']['dz']
    dt_fields_d = np.int(nml['fields_io']['frequency'])
    nml = simplejson.loads(open(os.path.join(path_triple, id_list_t[0], case_name_triple+ '.in')).read())
    nx_t[0] = nml['grid']['nx']
    nx_t[1] = nml['grid']['ny']
    nx_t[2] = nml['grid']['nz']
    dx_t[0] = nml['grid']['dx']
    dx_t[1] = nml['grid']['dy']
    dx_t[2] = nml['grid']['dz']
    dt_fields_t = np.int(nml['fields_io']['frequency'])
    print('')
    print('nx single: '+str(nx_s))
    print('nx double: '+str(nx_d))
    print('nx triple: '+str(nx_t))
    print('dx single: ' + str(dx_s))
    print('dx double: ' + str(dx_d))
    print('dx triple: ' + str(dx_t))
    # print('dt fields single: ' + str(dt_fields_s))
    print('dt fields double: ' + str(dt_fields_d))
    print('dt fields triple: ' +str(dt_fields_t))
    print('')

    global nx, ny, nz, dx
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = np.ndarray(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    global dt_fields
    dt_fields = np.maximum(dt_fields_d, dt_fields_t)

    return
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def cpm_c(qt):
    cpd = 1004.0
    cpv = 1859.0
    return (1.0-qt) * cpd + qt * cpv

def thetas_c(s, qt):
    T_tilde = 298.15
    sd_tilde = 6864.8
    sv_tilde = 10513.6
    return T_tilde*np.exp((s-(1.0-qt)*sd_tilde - qt*sv_tilde)/cpm_c(qt))

# ----------------------------------




if __name__ == '__main__':
    print('')
    main()

