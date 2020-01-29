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
    path_out_figs = '/nbi/home/meyerbe/paper_CP/'
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print('path figs: ' + path_out_figs)
    print('')

    dTh = 5
    zstar = 1000
    rstar = 1100
    print('zstar: '+str(zstar))
    print('rstar: '+str(rstar))
    print('')
    if rstar == 1100:
        d_range = [10, 12, 15]
        t_ini = [600, 600, 600]
        t_2CP = [1100, 1500, 2400]
        t_3CP = [1500, 2200, 3300]
        t_final = [2100, 2800, 3900]
    elif rstar == 2000:
        d_range = [10, 15, 20]
        t_ini = [400, 400, 400]
        t_2CP = [600, 800, 2200]
        t_3CP = [800, 1100, 3200]
        t_final = [1800, 2500, 3500]
    id_list_s = ['dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)]
    id_list_d = []
    id_list_t = []
    for dstar in d_range:
        id_list_d.append('dTh'+str(dTh)+'_z'+str(zstar)+'_r'+str(rstar)+'_d'+str(dstar) + 'km')
        id_list_t.append('dTh'+str(dTh)+'_z'+str(zstar)+'_r'+str(rstar)+'_d'+str(dstar) + 'km')
    print('ID-list: ', id_list_d)
    print('ID-list: ', id_list_t)
    print('')

    define_geometry(case_name_single, case_name_double, case_name_triple,
                    path_single, path_double, path_triple,
                    id_list_s, id_list_d, id_list_t)

    ''' determine max. height '''
    zmax = 5000.
    kmax = np.int(zmax/dx[2])
    print('zmax: '+str(zmax))
    print('kmax: '+str(kmax))
    print('')

    ''' determine time range '''
    global tmin, tmax
    tmin = np.int(0)
    tmax = np.int(np.amax(t_final)+2*dt_fields)  # bcs. of tracer statistics
    times = np.arange(tmin, tmax + dt_fields, dt_fields)
    nt = len(times)
    print('tmin: ' + str(tmin))
    print('tmax: ' + str(tmax))
    print('nt  : ' + str(nt))
    print('')


    ''' determine sampling subdomain '''
    print('define sampling domain')
    # for single CP: search maximum within circle of radius CP has at 2-CP / 3-CP collision time
    # single: from tracer statistics
    k0 = 0
    fullpath_in = os.path.join(path_single, id_list_s[0], 'tracer_k' + str(k0), 'output')
    n_tracers = get_number_tracers(fullpath_in)
    n_cps = get_number_cps(fullpath_in)
    dist_av = np.zeros((nt))
    for it, t0 in enumerate(times):
        cp_id = 2
        dist_av[it] = get_radius(fullpath_in, it, cp_id, n_tracers, n_cps)
    r_av = dist_av * dx[0]
    rad_1CP_ini = np.empty(3)
    rad_2CP_ini = np.empty(3)
    rad_3CP_ini = np.empty(3)
    rad_3CP_end = np.empty(3)
    delta_s = 6.e2
    for d in range(len(d_range)):
        rad_1CP_ini[d] = r_av[np.int(t_ini[d]/dt_fields)]
        rad_2CP_ini[d] = r_av[np.int(t_2CP[d]/dt_fields)]
        rad_3CP_ini[d] = r_av[np.int(t_3CP[d]/dt_fields)]
        rad_3CP_end[d] = r_av[np.int(t_final[d]/dt_fields)]
    [xs,ys] = nx_s[:2]*.5
    # double
    delta_d = np.asarray([1.e3, 8.e3]/dx[0])
    # triple
    delta_t = np.int(2.e3/dx[0])

    # plotting limits
    if rstar == 1100:
        lim_single = [0,0,0]
        lim_double = [[100, 200], [100, 250], [100, 220]]
        lim_triple = [120, 220, 300]
    elif rstar == 2000:
        lim_single = [0, 0, 0]
        lim_double = [[100, 280], [80, 280], [80, 220]]
        lim_triple = [200, 250, 250]

    # plot_CPs_at_times(rstar, xs, ys, delta_s, delta_d, delta_t, lim_single, lim_double, lim_triple,
    #                   d_range, t_ini, t_2CP, t_3CP, t_final,
    #                   rad_1CP_ini, rad_2CP_ini, rad_3CP_ini, rad_3CP_end,
    #                   id_list_s, id_list_d, id_list_t,
    #                   path_single, path_double, path_triple, path_out_figs)
    print('')










    ''' (A) plot from local (unaveraged) min / max in subdomains'''
    ''' compute min/max values in subdomains (circle, rectangles) '''
    path_out = os.path.join(path_single, id_list_s[0], 'data_analysis')
    filename = 'minmax_subdomains_noaverage.nc'

    #w_min, w_max, th_min, th_max, s_min, s_max, z, z_half = \
    #    compute_subdomains_max_single(path_single, id_list_s[0], case_name_single,
    #                                      delta_s,
    #                                      kmax, times, nt,
    #                                      t_ini[0], t_2CP, t_3CP, t_final, r_av, d_range, path_out_figs)
    #dump_minmax_file(w_min, w_max, th_min, th_max, s_min, s_max,
    #                 z, z_half, kmax, times, filename, path_out)
    #
    #
    # for d,dstar in enumerate(d_range):
    # #     w_min_d, w_max_d, th_min_d, th_max_d, s_min_d, s_max_d, z, z_half \
    # #         = compute_subdomains_max_double(path_double, id_list_d[d], case_name_double,
    # #                                         d, dstar, rstar,
    # #                                         r_av, delta_d,
    # #                                         times, nt, t_2CP[d], t_3CP[d],
    # #                                         kmax, path_out_figs)
    # #     path_out = os.path.join(path_double, id_list_d[d], 'data_analysis')
    # #     if not os.path.exists(path_out):
    # #         os.mkdir(path_out)
    # #     dump_minmax_file(w_min_d, w_max_d, th_min_d, th_max_d, s_min_d, s_max_d,
    # #                      z, z_half, kmax, times, filename, path_out)
    # #
    #     w_min_t, w_max_t, th_min_t, th_max_t, s_min_t, s_max_t, z, z_half \
    #         = compute_subdomains_max_triple(path_triple, id_list_t[d], case_name_triple,
    #                                         d, dstar,
    #                                         times, nt, t_3CP[d], t_final[d], delta_t,
    #                                         kmax, path_out_figs)
    #     path_out = os.path.join(path_triple, id_list_t[d], 'data_analysis')
    #     if not os.path.exists(path_out):
    #         os.mkdir(path_out)
    #     dump_minmax_file(w_min_t, w_max_t, th_min_t, th_max_t, s_min_t, s_max_t, z, z_half, kmax, times, filename, path_out)


    # plot min/max in each subdomain for all times
    plot_minmax_timeseries_subdomains(rstar, d_range, id_list_s, id_list_d, id_list_t,
                                     t_final,
                                     path_single, path_double, path_triple,
                                     filename, path_out_figs)


    print(path_double)
    # plot min/max in each subdomain for time windows
    plot_minmax_local_subdomain(rstar, d_range, id_list_s, id_list_d, id_list_t,
                               t_ini, t_2CP, t_3CP, t_final,
                               path_single, path_double, path_triple,
                               filename, path_out_figs)



    ''' (B) plot from local (unaveraged) min / max in total domain'''
    filename = 'minmax_domain_noaverage.nc'
    ''' compute domain min/max values '''
    w_min_s, w_max_s, th_min_s, th_max_s, s_min_s, s_max_s, z, z_half \
            = compute_domain_max(path_single, id_list_s[0], case_name_single, kmax, times, nt)
    path_out = os.path.join(path_single, id_list_s[0], 'data_analysis')
    dump_minmax_file(w_min_s, w_max_s, th_min_s, th_max_s, s_min_s, s_max_s, z, z_half, kmax, times, filename, path_out)
    # w_min_s, w_max_s, th_min_s, th_max_s, z, z_half \
    #         = compute_domain_max(path_single_dx50, id_list_s[0], case_name_single, kmax, times, nt)
    # path_out = os.path.join(path_single_dx50, id_list_s[0], 'data_analysis')
    # dump_minmax_file(w_min_s, w_max_s, th_min_s, th_max_s, z, z_half, kmax, times, filename, path_out)
    for d,dstar in enumerate(d_range):
    #     w_min_d, w_max_d, th_min_d, th_max_d, s_min_d, s_max_d, z, z_half \
    #         = compute_domain_max(path_double, id_list_d[d], case_name_double, kmax, times, nt)
    #     path_out = os.path.join(path_double, id_list_d[d], 'data_analysis')
    #     if not os.path.exists(path_out):
    #         os.mkdir(path_out)
    #     dump_minmax_file(w_min_d, w_max_d, th_min_d, th_max_d, s_min_d, s_max_d, z, z_half, kmax, times, filename, path_out)
    #     path_out = os.path.join(path_triple, id_list_t[d], 'data_analysis')
    #     if not os.path.exists(path_out):
    #         os.mkdir(path_out)
        w_min_t, w_max_t, th_min_t, th_max_t, s_min_t, s_max_t, z, z_half \
            = compute_domain_max(path_triple, id_list_t[d], case_name_triple, kmax, times, nt)
        dump_minmax_file(w_min_t, w_max_t, th_min_t, th_max_t, s_min_t, s_max_t, z, z_half, kmax, times, filename, path_out)

    plot_minmax_timeseries_domain(rstar, d_range, id_list_s, id_list_d, id_list_t,
                                  t_final,
                                  path_single, path_double, path_triple,
                                  filename, path_out_figs)




    ''' (C) plot from averaged min / max'''
    # - single - azimuthal average
    # - double - average along collision line
    # - triple - average about a few grid points at center


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
def compute_subdomains_max_triple(path, ID, casename,
                                  d, dstar,
                                  times, nt, t_ini, t_fi, delta_t,
                                  kmax, path_out_figs):
    print('triple: compute max in subdomain')
    # times = [0, ..., t_final+2*dt_fields], len(times)=nt
    w_min = np.zeros((nt, kmax), dtype=np.double)
    w_max = np.zeros((nt, kmax), dtype=np.double)
    th_min = np.ones((nt, kmax), dtype=np.double)
    th_max = np.ones((nt, kmax), dtype=np.double)
    s_min = np.ones((nt, kmax), dtype=np.double)
    s_max = np.ones((nt, kmax), dtype=np.double)

    root = nc.Dataset(os.path.join(path, ID, 'stats', 'Stats.' + casename + '.nc'))
    zrange = root.groups['profiles'].variables['z'][:kmax]
    zrange_half = root.groups['profiles'].variables['z_half'][:kmax]
    root.close()

    if not os.path.exists(os.path.join(path_out_figs, 'figs_collision_test')):
        os.mkdir(os.path.join(path_out_figs, 'figs_collision_test'))

    # it_ini = np.int(t_ini / dt_fields)
    # it_fi = np.int(t_fi / dt_fields)
    # ic = np.int(nx_t[d][0] * .5)
    # jc = np.int(nx_t[d][1] * .5)
    for it, t0 in enumerate(times):
        print('--- t: ', it, t0)
        fullpath_in = os.path.join(path, ID, 'fields', str(t0) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, :kmax]
        s = root.groups['fields'].variables['s'][:, :, :kmax]
        root.close()

        xt = np.int(nx_t[d][0]*.5 - delta_t*.5)
        yt = np.int(nx_t[d][1]*.5 - delta_t*.5)
        #print('triple: ', xt, yt, delta_t)
        w_min[it, :] = np.amin(np.amin(w[xt:xt+delta_t,yt:yt+delta_t,:], axis=0), axis=0)
        w_max[it, :] = np.amax(np.amax(w[xt:xt+delta_t,yt:yt+delta_t,:], axis=0), axis=0)
        s_min[it,:] = np.amin(np.amin(s[xt:xt+delta_t,yt:yt+delta_t,:], axis=0), axis=0)
        s_max[it,:] = np.amax(np.amax(s[xt:xt+delta_t,yt:yt+delta_t,:], axis=0), axis=0)
        th_min[it, :] = thetas_c(s_min[it,:], 0)
        th_max[it, :] = thetas_c(s_max[it,:], 0)
    return w_min, w_max, th_min, th_max, s_min, s_max, zrange, zrange_half



def compute_subdomains_max_double(path, ID, casename, d, dstar, rstar,
                                  r_av, delta_d,
                                  times, nt, t_ini, t_fi,
                                  kmax, path_out_figs):
    print('double: compute max in subdomain')
    w_min = np.zeros((nt, kmax), dtype=np.double)
    w_max = np.zeros((nt, kmax), dtype=np.double)
    th_min = np.zeros((nt, kmax), dtype=np.double)
    th_max = np.zeros((nt, kmax), dtype=np.double)
    s_min = np.zeros((nt, kmax), dtype=np.double)
    s_max = np.zeros((nt, kmax), dtype=np.double)

    root = nc.Dataset(os.path.join(path, ID, 'stats', 'Stats.' + casename + '.nc'))
    zrange = root.groups['profiles'].variables['z'][:kmax]
    zrange_half = root.groups['profiles'].variables['z_half'][:kmax]
    root.close()

    if not os.path.exists(os.path.join(path_out_figs, 'figs_collision_test')):
        os.mkdir(os.path.join(path_out_figs, 'figs_collision_test'))

    it_ini = np.int(t_ini / dt_fields)
    it_fi = np.int(t_fi[d] / dt_fields)
    ic = np.int(nx_d[d][0] * .5)
    jc = np.int(nx_d[d][1] * .5)
    for it, t0 in enumerate(times):
        # it = it_+it_ini
        print('--- t: ', it, t0)
        fullpath_in = os.path.join(path, ID, 'fields', str(t0) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, :kmax]
        s = root.groups['fields'].variables['s'][:, :, :kmax]
        root.close()
        # rad = (delta_s + r_av[np.int(t0 / dt_fields)]) / dx[0]
        # mask = np.ma.masked_less(r_field, rad)
        # w_masked = r_mask.mask[:, :, np.newaxis] * w
        # s_masked = r_mask.mask[:, :, np.newaxis] * s
        rad2 = r_av[it]**2
        d2 = (dstar*1.e3/2)**2
        #print('double: ', rad2, d2)
        if rad2 >= d2:
            delta_d[1] = np.floor(2.e3/dx[0]) + 2./dx[0]*np.sqrt(rad2 - d2)
            [xd, yd] = [ic, jc] - delta_d * .5
            w_min[it, :] = np.amin(np.amin(w[xd:xd + delta_d[0], yd:yd + delta_d[1], :], axis=0), axis=0)
            w_max[it, :] = np.amax(np.amax(w[xd:xd + delta_d[0], yd:yd + delta_d[1], :], axis=0), axis=0)
            s_min[it, :] = np.amin(np.amin(s[xd:xd + delta_d[0], yd:yd + delta_d[1], :], axis=0), axis=0)
            s_max[it, :] = np.amax(np.amax(s[xd:xd + delta_d[0], yd:yd + delta_d[1], :], axis=0), axis=0)
            th_min[it, :] = thetas_c(s_min[it, :], 0)
            th_max[it, :] = thetas_c(s_max[it, :], 0)
        else:
            [xd,yd] = [ic,jc]
            delta_d[1] = 1
            w_min[it, :] = 0
            w_max[it, :] = 0
            s_min[it, :] = np.amin(np.amin(s, axis=0), axis=0)
            s_max[it, :] = np.amax(np.amax(s, axis=0), axis=0)
            th_min[it, :] = thetas_c(s_min[it, :], 0)
            th_max[it, :] = thetas_c(s_max[it, :], 0)
        #print('double: ', ic,jc, xd, yd, delta_d)
        #print('...', np.amax(w_max[it,:]), np.amin(s_min[it,:]), np.amin(th_min[it,:]))

    return w_min, w_max, th_min, th_max, s_min, s_max, zrange, zrange_half



def compute_subdomains_max_single(path, ID, casename,
                           delta_s, kmax, times, nt,
                           t_ini, t_2CP, t_3CP, t_final,
                           r_av, d_range, path_out_figs):
    print('single: compute max in subdomain')
    w_min = np.zeros((nt, kmax), dtype=np.double)
    w_max = np.zeros((nt, kmax), dtype=np.double)
    th_min = np.zeros((nt, kmax), dtype=np.double)
    th_max = np.zeros((nt, kmax), dtype=np.double)
    s_min = np.zeros((nt, kmax), dtype=np.double)
    s_max = np.zeros((nt, kmax), dtype=np.double)

    root = nc.Dataset(os.path.join(path, ID, 'stats', 'Stats.'+casename+'.nc'))
    zrange = root.groups['profiles'].variables['z'][:kmax]
    zrange_half = root.groups['profiles'].variables['z_half'][:kmax]
    root.close()

    if not os.path.exists(os.path.join(path_out_figs, 'figs_collision_test')):
        os.mkdir(os.path.join(path_out_figs, 'figs_collision_test'))

    root = nc.Dataset(os.path.join(path, ID, 'fields_v_rad', 'v_rad.nc'), 'r')
    r_field = root.variables['r_field'][:,:]
    root.close()

    it_ini = np.int(t_ini/dt_fields)
    for it, t0 in enumerate(times):
        print('--- t: ', it, t0)
        fullpath_in = os.path.join(path, ID, 'fields', str(t0) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, :kmax]
        s = root.groups['fields'].variables['s'][:, :, :kmax]
        root.close()
        rad = (delta_s + r_av[np.int(t0 / dt_fields)]) / dx[0]
        r_mask_l = np.ma.masked_less(r_field, rad)
        r_mask_g = np.ma.masked_greater_equal(r_field, rad)
        w_masked = r_mask_l.mask[:,:,np.newaxis]*w
        s_masked = (1+9*r_mask_g.mask[:,:,np.newaxis])*s
        w_min[it,:] = np.amin(np.amin(w_masked, axis=0), axis=0)
        w_max[it,:] = np.amax(np.amax(w_masked, axis=0), axis=0)
        s_min[it,:] = np.amin(np.amin(s_masked, axis=0), axis=0)
        s_max[it,:] = np.amax(np.amax(s_masked, axis=0), axis=0)
        th_min[it, :] = thetas_c(s_min[it,:], 0)
        #print('single: smin', s_min.shape, s_min[it,:], th_min[it,5])
        th_max[it, :] = thetas_c(s_max[it,:], 0)


        if it < 100:
            for k0 in [0,5,10,15,20]:
                w_max_loc_all = np.unravel_index(np.argmax(w[:, :, k0], axis=None), nx_s[:2])
                w_max_loc = np.unravel_index(np.argmax(w_masked[:, :, k0], axis=None), nx_s[:2])
                s_min_loc_all = np.unravel_index(np.argmin(s[:, :, k0], axis=None), nx_s[:2])
                s_min_loc = np.unravel_index(np.argmin(s_masked[:, :, k0], axis=None), nx_s[:2])
                # print('loc: ', w_max_loc, w_max_loc_all, s_min_loc, s_min_loc_all)
                fig_name = 'w_masked_single_t'+str(t0)+'_k'+str(k0)+'.png'
                fig, axis = plt.subplots(2,5, figsize=(16,7))
                ax0 = axis[0,0]
                ax1 = axis[0,1]
                ax2 = axis[0,2]
                ax3 = axis[0,3]
                ax4 = axis[0,4]
                ax0.set_title('radius')
                cf = ax0.imshow(r_field)
                plt.colorbar(cf, ax=ax0, shrink=0.7)
                ax1.set_title('radius masked (r>=r_av[t])')
                cf = ax1.imshow(r_mask_l)
                plt.colorbar(cf, ax=ax1, shrink=0.7)
                ax2.set_title('mask')
                cf = ax2.imshow(r_mask_l.mask)
                plt.colorbar(cf, ax=ax2, shrink=0.7)
                ax3.set_title('w')
                cf = ax3.imshow(w[:,:,k0])
                plt.colorbar(cf, ax=ax3, shrink=0.7)
                ax3.plot(w_max_loc_all[0], w_max_loc_all[1],'x', color='k', markersize=10)
                ax3.plot(w_max_loc[0], w_max_loc[1],'ok', markersize=5)
                ax4.set_title('w*mask')
                cf = ax4.imshow(w_masked[:,:,k0])
                plt.colorbar(cf, ax=ax4, shrink=0.7)

                ax2 = axis[1, 2]
                ax3 = axis[1, 3]
                ax4 = axis[1, 4]
                ax2.set_title('mask g')
                cf = ax2.imshow(r_mask_g.mask)
                plt.colorbar(cf, ax=ax2, shrink=0.7)
                ax3.set_title('s')
                cf = ax3.imshow(s[:, :, k0], vmin=np.amin(s[:,:,k0]), vmax=np.amax(s[:,:,k0]))
                plt.colorbar(cf, ax=ax3, shrink=0.7)
                ax3.plot(s_min_loc_all[0], s_min_loc_all[1], 'kx', markersize=10)
                ax3.plot(s_min_loc[0], s_min_loc[1], 'ok', markersize=5)
                ax4.set_title('s*mask')
                cf = ax4.imshow(s_masked[:, :, k0], vmin=6864, vmax=np.amax(s[:,:,k0]))
                plt.colorbar(cf, ax=ax4, shrink=0.7, extend='max')

                for ax in axis.flat:
                    ax.set_xlim(nx_s[0]*.5-rad-10,nx_s[0]*.5+rad+10)
                    ax.set_ylim(nx_s[0]*.5-rad-10,nx_s[0]*.5+rad+10)
                plt.subplots_adjust(bottom=0.1, right=.95, left=0.05, top=0.9, hspace=0.4, wspace=0.2)
                plt.savefig(os.path.join(path_out_figs, 'figs_collision_test', fig_name))
                plt.close()

        del w_masked, s_masked
        del w, s
        # del s_min, s_max


    fig,axis = plt.subplots(1,3)
    for it, t0 in enumerate(times):
        if t0>t_ini:
            if t0<=t_2CP[0]:
                axis[0].plot(w_max[it, :], zrange)
                axis[0].set_title('t='+str(t_ini) + '-'+str(t_2CP[0]))
            if t0<=t_2CP[1]:
                axis[1].plot(w_max[it,:],zrange)
                axis[1].set_title('t='+str(t_ini) + '-'+str(t_2CP[1]))
            if t0<=t_2CP[2]:
                axis[2].plot(w_max[it,:],zrange)
                axis[2].set_title('t='+str(t_ini) + '-'+str(t_2CP[2]))
    for ax in axis:
        ax.set_xlim(0,4)
    plt.savefig(os.path.join(path_out_figs, 'figs_collision_test', 'w_max_test.png'))

    return w_min, w_max, th_min, th_max, s_min, s_max, zrange, zrange_half
# --------------------------------------------------------------------
def compute_domain_max(path, ID, casename, kmax, times, nt):
    w_min = np.ones((nt, kmax), dtype=np.double)
    w_max = np.ones((nt, kmax), dtype=np.double)
    th_min = np.ones((nt, kmax), dtype=np.double)
    th_max = np.ones((nt, kmax), dtype=np.double)
    s_min = np.ones((nt, kmax), dtype=np.double)
    s_max = np.ones((nt, kmax), dtype=np.double)

    root = nc.Dataset(os.path.join(path, ID, 'stats', 'Stats.'+casename+'.nc'))
    zrange = root.groups['profiles'].variables['z'][:kmax]
    zrange_half = root.groups['profiles'].variables['z_half'][:kmax]
    root.close()

    for it, t0 in enumerate(times):
        fullpath_in = os.path.join(path, ID, 'fields', str(t0) + '.nc')
        print(fullpath_in)
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, :kmax]
        w_min[it,:] = np.amin(np.amin(w, axis=0), axis=0)
        w_max[it,:] = np.amax(np.amax(w, axis=0), axis=0)
        del w
        s = root.groups['fields'].variables['s'][:, :, :kmax]
        s_min[it,:] = np.amin(np.amin(s, axis=0), axis=0)
        s_max[it,:] = np.amax(np.amax(s, axis=0), axis=0)
        del s
        th_min[it,:] = thetas_c(s_min[it,:], 0)
        th_max[it,:] = thetas_c(s_max[it,:], 0)
        root.close()
        # del s_min, s_max
    return w_min, w_max, th_min, th_max, s_min, s_max, zrange, zrange_half


def dump_minmax_file(w_min, w_max, th_min, th_max, s_min, s_max,
                     z, z_half, kmax, times, filename, path_out):
    print('dumping: ', filename, path_out)
    # output for each CP:
    # - min, max (timeseries)
    # - CP height (field; max=timeseries)
    # - (ok) CP rim (field)
    nt = len(times)
    print('create output file: ', os.path.join(path_out, filename))
    #print('size: ', nz, nt)

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
    var.long_name = 'max potential temperature'
    var.units = "K"
    var[:,:] = th_max[:,:]
    var = ts_grp.createVariable('th_min', 'f8', ('nt', 'nz'))
    var.long_name = 'min potential temperature'
    var.units = "K"
    var[:,:] = th_min[:,:]

    var = ts_grp.createVariable('s_max', 'f8', ('nt', 'nz'))
    var.long_name = 'max entropy'
    var.units = "K"
    var[:, :] = s_max[:, :]
    var = ts_grp.createVariable('s_min', 'f8', ('nt', 'nz'))
    var.long_name = 'min entropy'
    var.units = "K"
    var[:, :] = s_min[:, :]

    var = ts_grp.createVariable('w_max', 'f8', ('nt', 'nz'))
    var.long_name = 'domain max vertical velocity'
    var.units = "m/s"
    var[:,:] = w_max[:,:]
    var = ts_grp.createVariable('w_min', 'f8', ('nt', 'nz'))
    var.long_name = 'domain min vertical velocity'
    var.units = "m/s"
    var[:,:] = w_min[:,:]

    rootgrp.close()
    print ''
    return


def read_in_minmax(kmax, path_out, filename):
    root = nc.Dataset(os.path.join(path_out, filename), 'r')
    time = root.groups['timeseries'].variables['time'][:]
    z = root.groups['timeseries'].variables['z'][:kmax]
    z_half = root.groups['timeseries'].variables['z_half'][:kmax]
    w_max = root.groups['timeseries'].variables['w_max'][:,:kmax]
    th_min = root.groups['timeseries'].variables['th_min'][:,:kmax]
    s_min = root.groups['timeseries'].variables['s_min'][:,:kmax]
    root.close()
    return w_max, th_min, s_min, z, z_half, time

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
    nx_d = []
    nx_t = []
    dx_s = np.empty(3, dtype=np.int)
    dx_d = []
    dx_t = []

    nml = simplejson.loads(open(os.path.join(path_single, id_list_s[0], case_name_single+ '.in')).read())
    nx_s[0] = nml['grid']['nx']
    nx_s[1] = nml['grid']['ny']
    nx_s[2] = nml['grid']['nz']
    dx_s[0] = nml['grid']['dx']
    dx_s[1] = nml['grid']['dy']
    dx_s[2] = nml['grid']['dz']
    dt_fields_s = np.int(nml['fields_io']['frequency'])
    for d, ID in enumerate(id_list_d):
        nml = simplejson.loads(open(os.path.join(path_double, ID, case_name_double+ '.in')).read())
        nx_d.append(np.empty(3, dtype=np.int))
        nx_d[d][0] = nml['grid']['nx']
        nx_d[d][1] = nml['grid']['ny']
        nx_d[d][2] = nml['grid']['nz']
        dx_d.append(np.empty(3, dtype=np.int))
        dx_d[d][0] = nml['grid']['dx']
        dx_d[d][1] = nml['grid']['dy']
        dx_d[d][2] = nml['grid']['dz']
    dt_fields_d = np.int(nml['fields_io']['frequency'])
    for d, ID in enumerate(id_list_t):
        nml = simplejson.loads(open(os.path.join(path_triple, ID, case_name_triple+ '.in')).read())
        nx_t.append(np.empty(3, dtype=np.int))
        nx_t[d][0] = nml['grid']['nx']
        nx_t[d][1] = nml['grid']['ny']
        nx_t[d][2] = nml['grid']['nz']
        dx_t.append(np.empty(3, dtype=np.int))
        dx_t[d][0] = nml['grid']['dx']
        dx_t[d][1] = nml['grid']['dy']
        dx_t[d][2] = nml['grid']['dz']
    dt_fields_t = np.int(nml['fields_io']['frequency'])
    print('')
    print('nx single: '+str(nx_s))
    print('nx double: '+str(nx_d))
    print('nx triple: '+str(nx_t))
    print('dx single: ' + str(dx_s))
    print('dx double: ' + str(dx_d))
    print('dx triple: ' + str(dx_t))
    print('dt fields single: ' + str(dt_fields_s))
    print('dt fields double: ' + str(dt_fields_d))
    print('dt fields triple: ' + str(dt_fields_t))
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
def plot_CPs_at_times(rstar, xs, ys, delta_s, delta_d, delta_t, lim_single, lim_double, lim_triple,
                      d_range, t_ini, t_2CP, t_3CP, t_final,
                      rad_1CP_ini, rad_2CP_ini, rad_3CP_ini, rad_3CP_end,
                      id_list_s, id_list_d, id_list_t,
                      path_single, path_double, path_triple, path_out_figs):
    ''' testing '''
    k0 = 0
    if k0 == 0:
        w_max = 3.
    elif k0 == 10:
        w_max = 1.
    lvls = np.arange(-w_max, w_max, .25)
    lvls = np.arange(-w_max, w_max, .5)
    # lvls = np.linspace(-4, 4, 10)
    for d, dstar in enumerate(d_range):
        print('plotting CPs: d='+str(dstar), d)
        print('times: ', t_ini[d], t_2CP[d], t_3CP[d], t_final[d])
        fig_name = 'collisions_subdomains_k'+str(k0)+'_'+id_list_s[0] + '_d'+str(dstar)+ 'km.png'
        fig, axis = plt.subplots(3, 4, figsize=(20, 15))
        fullpath_in = os.path.join(path_single, id_list_s[0], 'fields', str(t_ini[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        cf = axis[0, 0].contourf(w.T, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[0, 0], shrink=0.8)
        fullpath_in = os.path.join(path_single, id_list_s[0], 'fields', str(t_2CP[d])+'.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:,:,k0]
        root.close()
        cf = axis[0,1].contourf(w.T, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[0,1], shrink=0.8)
        fullpath_in = os.path.join(path_single, id_list_s[0], 'fields', str(t_3CP[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        cf = axis[0,2].contourf(w.T, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[0, 2], shrink=0.8)
        fullpath_in = os.path.join(path_single, id_list_s[0], 'fields', str(t_final[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        cf = axis[0, 3].contourf(w.T, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[0, 3], shrink=0.8)

        fullpath_in = os.path.join(path_double, id_list_d[d], 'fields', str(t_ini[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        cf = axis[1, 0].contourf(w, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[1, 0], shrink=0.8)
        fullpath_in = os.path.join(path_double, id_list_d[d], 'fields', str(t_2CP[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        cf = axis[1, 1].contourf(w, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[1, 1], shrink=0.8)
        fullpath_in = os.path.join(path_double, id_list_d[d], 'fields', str(t_3CP[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        cf = axis[1, 2].contourf(w, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[1, 2], shrink=0.8)
        fullpath_in = os.path.join(path_double, id_list_d[d], 'fields', str(t_final[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        cf = axis[1, 3].contourf(w, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[1, 3], shrink=0.8)

        fullpath_in = os.path.join(path_triple, id_list_t[d], 'fields', str(t_ini[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        cf = axis[2, 0].contourf(w.T, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[2, 0], shrink=0.8)
        fullpath_in = os.path.join(path_triple, id_list_t[d], 'fields', str(t_2CP[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        cf = axis[2, 1].contourf(w.T, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[2, 1], shrink=0.8)
        fullpath_in = os.path.join(path_triple, id_list_t[d], 'fields', str(t_3CP[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        cf = axis[2, 2].contourf(w.T, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[2, 2], shrink=0.8)
        fullpath_in = os.path.join(path_triple, id_list_t[d], 'fields', str(t_final[d]) + '.nc')
        root = nc.Dataset(fullpath_in, 'r')
        w = root.groups['fields'].variables['w'][:, :, k0]
        root.close()
        cf = axis[2, 3].contourf(w.T, levels=lvls, cmap=cm_bwr, extend='both')
        plt.colorbar(cf, ax=axis[2, 3], shrink=0.8)

        circle0 = plt.Circle((xs, ys), (rad_1CP_ini[d]+delta_s)/dx[0], fill=False, color='k', linewidth=1)
        circle1 = plt.Circle((xs, ys), (rad_2CP_ini[d]+delta_s)/dx[0], fill=False, color='k', linewidth=1)
        circle2 = plt.Circle((xs, ys), (rad_3CP_ini[d]+delta_s)/dx[0], fill=False, color='k', linewidth=1)
        circle3 = plt.Circle((xs, ys), (rad_3CP_end[d]+delta_s)/dx[0], fill=False, color='k', linewidth=1)
        axis[0,0].add_artist(circle0)
        axis[0,1].add_artist(circle1)
        axis[0,2].add_artist(circle2)
        axis[0,3].add_artist(circle3)
        ic = np.int(nx_d[d][0]*.5)
        jc = np.int(nx_d[d][1]*.5)
        # delta_d[1] = np.floor(1.e3/dx[0]) + 2./dx[0]*np.sqrt(rad_1CP_ini[d]**2 - (dstar*1.e3)**2/4)
        # print('delta_d', rad_1CP_ini[d], (dstar*1.e3)/2, rad_1CP_ini[d]**2, (dstar*1.e3)**2/4, rad_1CP_ini[d]**2 - (dstar*1.e3)**2/4)
        # [xd, yd] = [ic,jc] - delta_d*.5
        # rect_double0 = mpatches.Rectangle((yd, xd), delta_d[1], delta_d[0], linewidth=2, edgecolor='k', facecolor='none')
        rad2 = rad_2CP_ini[d]**2
        d2 = (dstar*1.e3/2)**2
        if rad2 >= d2:
            delta_d[1] = np.floor(2.e3/dx[0]) + 2./dx[0]*np.sqrt(rad_2CP_ini[d]**2 - (dstar*1.e3)**2/4)
            [xd, yd] = [ic,jc] - delta_d*.5
            rect_double1 = mpatches.Rectangle((yd, xd), delta_d[1], delta_d[0], linewidth=2, edgecolor='k', facecolor='none')
            axis[1,1].add_patch(rect_double1)
        delta_d[1] = np.floor(2.e3/dx[0]) + 2./dx[0]*np.sqrt(rad_3CP_ini[d]**2 - (dstar*1.e3/2)**2)
        [xd, yd] = [ic,jc] - delta_d*.5
        rect_double2 = mpatches.Rectangle((yd, xd), delta_d[1], delta_d[0], linewidth=2, edgecolor='k', facecolor='none')
        delta_d[1] = np.floor(2.e3/dx[0]) + 2./dx[0] * np.sqrt(rad_3CP_end[d]**2 - (dstar*1.e3/2)**2)
        [xd, yd] = [ic, jc] - delta_d * .5
        rect_double3 = mpatches.Rectangle((yd, xd), delta_d[1], delta_d[0], linewidth=2, edgecolor='k', facecolor='none')
        axis[1,2].add_patch(rect_double2)
        axis[1,3].add_patch(rect_double3)

        [xt, yt] = nx_t[d][:2] * .5 - delta_t * .5
        rect_triple0 = mpatches.Rectangle((xt, yt), delta_t, delta_t, linewidth=1, edgecolor='k', facecolor='none')
        rect_triple1 = mpatches.Rectangle((xt, yt), delta_t, delta_t, linewidth=1, edgecolor='k', facecolor='none')
        rect_triple2 = mpatches.Rectangle((xt, yt), delta_t, delta_t, linewidth=1, edgecolor='k', facecolor='none')
        rect_triple3 = mpatches.Rectangle((xt, yt), delta_t, delta_t, linewidth=1, edgecolor='k', facecolor='none')
        axis[2, 0].add_patch(rect_triple0)
        axis[2, 1].add_patch(rect_triple1)
        axis[2, 2].add_patch(rect_triple2)
        axis[2, 3].add_patch(rect_triple3)

        axis[0,0].set_title('t='+str(t_ini[d])+'s')
        axis[0,1].set_title('t='+str(t_2CP[d])+'s')
        axis[0,2].set_title('t='+str(t_3CP[d])+'s')
        axis[0,3].set_title('t='+str(t_final[d])+'s')
        for ax in axis[0,:].flat:
            ax.set_xlim(lim_single[d], nx_s[0]-lim_single[d])
            ax.set_ylim(lim_single[0], nx_s[1]-lim_single[d])
        for ax in axis[1,:].flat:
            ax.set_xlim(lim_double[d][0],nx_d[d][1]-lim_double[d][0])
            ax.set_ylim(lim_double[d][1],nx_d[d][0]-lim_double[d][1])
        for ax in axis[2, :].flat:
            ax.set_xlim(lim_triple[d], nx_t[d][0]-lim_triple[d])
            ax.set_ylim(lim_triple[d], nx_t[d][1]-lim_triple[d])
        for ax in axis.flat:
            ax.set_aspect('equal')
        plt.subplots_adjust(bottom=0.05, right=.95, left=0.05, top=0.95, hspace=0.1, wspace=0.1)
        print('saving: ', fig_name)
        plt.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)

    return

# ----------------------------------------------------------------------
def plot_minmax_timeseries_domain(rstar, d_range, id_list_s, id_list_d, id_list_t,
                                      t_final,
                                      path_single, path_double, path_triple,
                                      filename, path_out_figs):
    for d, dstar in enumerate(d_range):
        fig_name = 'collisions_minmax_allltimes_domain_unaveraged_rstar' + str(rstar) + '_d'+str(dstar)+'km.png'
        zmax_plot = 3000.
        kmax_plot = np.int(zmax_plot / dx[2])
        path = os.path.join(path_single, id_list_s[0], 'data_analysis')
        w_max_s, th_min_s, s_min_s, z, z_half, t_s = read_in_minmax(kmax_plot, path, filename)
        #path = os.path.join(path_double, id_list_d[d], 'data_analysis')
        #w_max_d, th_min_d, s_min_d, z, z_half, t_d = read_in_minmax(kmax_plot, path, filename)
        path = os.path.join(path_triple, id_list_t[d], 'data_analysis')
        w_max_t, th_min_t, s_min_t, z, z_half, t_t = read_in_minmax(kmax_plot, path, filename)

        fig, axis = plt.subplots(2, 4, figsize=(14, 12), sharey='none')
        # maxw = np.amax(w_max_s)+.1
        # #maxw = np.maximum(np.amax(w_max_s), np.amax(w_max_d))+.1
        # print('time single: ', t_s, w_max_s.shape)
        axis[0, 0].plot(t_s, np.amax(w_max_s[:, :], axis=1), 'o-', color=colorlist3[0], label='single CP gust front')
        axis[0, 0].plot(t_t, np.amax(w_max_t[:, :], axis=1), 'o-', color=colorlist3[2], label='triple CP collision')
        axis[1, 0].plot(t_s, np.amin(s_min_s[:, :], axis=1), 'o-', color=colorlist3[0], label='single CP gust front')
        axis[1, 0].plot(t_t, np.amin(s_min_t[:, :], axis=1), 'o-', color=colorlist3[2], label='triple CP collision')
        for ax in axis[0, 1:].flat:
            ax.plot([0., maxw], [1000, 1000], 'k-', linewidth=0.5)
        for ax in axis[1, 1:].flat:
            ax.plot([298, 300.1], [1000, 1000], 'k-', linewidth=0.5)
        for it,t0 in enumerate(range(0, t_final[d], dt_fields)):
            lbl = 't='+str(t0)+'s'
            cl = t0*1./t_final[d]

            axis[0, 1].plot(w_max_s[it, :kmax_plot], z[:kmax_plot], color=cm(cl), label=lbl)
            #axis[0, 2].plot(w_max_d[it, :kmax_plot], z[:kmax_plot], color=cm(cl), label=lbl)
            axis[0, 3].plot(w_max_t[it, :kmax_plot], z[:kmax_plot], color=cm(cl), label=lbl)

            axis[1, 1].plot(th_min_s[it, :kmax_plot], z[:kmax_plot], color=cm(cl), label=lbl)
            #axis[1, 2].plot(th_min_d[it, :kmax_plot], z[:kmax_plot], color=cm(cl), label=lbl)
            axis[1, 3].plot(th_min_t[it, :kmax_plot], z[:kmax_plot], color=cm(cl), label=lbl)

        axis[0, 1].set_title('single CP')
        axis[0, 2].set_title('double CP, collision line')
        axis[0, 3].set_title('triple CP, collision point')
        for ax in axis[:,1:].flat:
            ax.set_ylabel('height z  [m]')
        for ax in axis[0,1:].flat:
            ax.set_xlabel('max(w)')
            ax.set_xlim(-0.1, maxw)
        for ax in axis[1,1:].flat:
            ax.set_xlim(298,300.1)
            ax.set_xlabel('min(theta)')
        axis[0,0].set_xlabel('time [s]')
        axis[1,0].set_xlabel('time [s]')
        axis[0,0].set_ylabel('max(w)')
        axis[1,0].set_ylabel('min(s)')

        for ax in axis[:, 2].flat:
            ax.axis('off')

        axis[0, 0].legend(loc=1, fontsize=12)
        axis[0, 3].legend(loc='upper left', bbox_to_anchor=(1, 1.),
                   fancybox=False, shadow=False, ncol=1, fontsize=12)
        plt.subplots_adjust(bottom=0.1, right=.9, left=0.05, top=0.95, hspace=0.2, wspace=0.1)
        plt.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)
    return




def plot_minmax_timeseries_subdomains(rstar, d_range, id_list_s, id_list_d, id_list_t,
                                      t_final,
                                      path_single, path_double, path_triple,
                                      filename, path_out_figs):
    for d, dstar in enumerate(d_range):
        fig_name = 'collisions_minmax_alltimes_subdomain_unaveraged_rstar' + str(rstar) + '_d'+str(dstar)+'km.png'
        zmax_plot = 3000.
        kmax_plot = np.int(zmax_plot / dx[2])
        path = os.path.join(path_single, id_list_s[0], 'data_analysis')
        w_max_s, th_min_s, s_min_s, z, z_half, t_s = read_in_minmax(kmax_plot, path, filename)
        #path = os.path.join(path_double, id_list_d[d], 'data_analysis')
        #w_max_d, th_min_d, s_min_d, z, z_half, t_d = read_in_minmax(kmax_plot, path, filename)
        path = os.path.join(path_triple, id_list_t[d], 'data_analysis')
        w_max_t, th_min_t, s_min_t, z, z_half, t_t = read_in_minmax(kmax_plot, path, filename)

        fig, axis = plt.subplots(2, 4, figsize=(14, 12), sharey='none')
        maxw = np.amax(w_max_s)+.1
        #maxw = np.maximum(np.amax(w_max_s), np.amax(w_max_d))+.1
        print('time single: ', t_s, w_max_s.shape)
        axis[0, 0].plot(t_s, np.amax(w_max_s[:, :], axis=1), 'o-', color=colorlist3[0], label='single CP')
        axis[0, 0].plot(t_t, np.amax(w_max_t[:, :], axis=1), 'o-', color=colorlist3[2], label='triple CP')
        axis[1, 0].plot(t_s, np.amin(s_min_s[:, :], axis=1), 'o-', color=colorlist3[0], label='single CP')
        axis[1, 0].plot(t_t, np.amin(s_min_t[:, :], axis=1), 'o-', color=colorlist3[2], label='triple CP')
        for ax in axis[0, 1:].flat:
            ax.plot([0., maxw], [1000, 1000], 'k-', linewidth=0.5)
        for ax in axis[1, 1:].flat:
            ax.plot([298, 300.1], [1000, 1000], 'k-', linewidth=0.5)
        for it,t0 in enumerate(range(0, t_final[d], dt_fields)):
            lbl = 't='+str(t0)+'s'
            cl = t0*1./t_final[d]

            axis[0, 1].plot(w_max_s[it, :kmax_plot], z[:kmax_plot], color=cm(cl), label=lbl)
            #axis[0, 2].plot(w_max_d[it, :kmax_plot], z[:kmax_plot], color=cm(cl), label=lbl)
            axis[0, 3].plot(w_max_t[it, :kmax_plot], z[:kmax_plot], color=cm(cl), label=lbl)

            axis[1, 1].plot(th_min_s[it, :kmax_plot], z[:kmax_plot], color=cm(cl), label=lbl)
            #axis[1, 2].plot(th_min_d[it, :kmax_plot], z[:kmax_plot], color=cm(cl), label=lbl)
            axis[1, 3].plot(th_min_t[it, :kmax_plot], z[:kmax_plot], color=cm(cl), label=lbl)

        axis[0, 1].set_title('single CP')
        axis[0, 2].set_title('double CP, collision line')
        axis[0, 3].set_title('triple CP, collision point')
        for ax in axis[:,1:].flat:
            ax.set_ylabel('height z  [m]')
        for ax in axis[0,1:].flat:
            ax.set_xlabel('max(w)')
            ax.set_xlim(-0.1, maxw)
        for ax in axis[1,1:].flat:
            ax.set_xlim(298,300.1)
            ax.set_xlabel('min(theta)')
        axis[0,0].set_xlabel('time [s]')
        axis[1,0].set_xlabel('time [s]')
        axis[0,0].set_ylabel('max(w)')
        axis[1,0].set_ylabel('min(s)')

        for ax in axis[:, 2].flat:
            ax.axis('off')

        axis[0, 0].legend(loc=1, fontsize=12)
        axis[0, 3].legend(loc='upper left', bbox_to_anchor=(1, 1.),
                   fancybox=False, shadow=False, ncol=1, fontsize=12)
        plt.subplots_adjust(bottom=0.1, right=.9, left=0.05, top=0.95, hspace=0.2, wspace=0.1)
        plt.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)
    return



def plot_minmax_local_subdomain(rstar, d_range, id_list_s, id_list_d, id_list_t,
                                t_ini, t_2CP, t_3CP, t_final,
                                path_single, path_double, path_triple,
                                filename, path_out_figs):
    ''' read in min/max values '''
    print('plotting min / max in subdomain')
    fig_name = 'collisions_minmax_profiles_subdomain_unaveraged_rstar'+str(rstar)+'.png'
    zmax_plot = 3000.
    kmax_plot = np.int(zmax_plot/dx[2])
    fig, axis = plt.subplots(2, 4, figsize=(14, 12), sharey='all')

    path = os.path.join(path_single, id_list_s[0], 'data_analysis')
    w_max_s, th_min_s, s_min_s, z, z_half, t_s = read_in_minmax(kmax_plot, path, filename)

    for d,dstar in enumerate(d_range):
        print('.... d: '+str(dstar))
        path = os.path.join(path_double, id_list_d[d], 'data_analysis')
        print(os.path.join(path, filename))
        w_max_d, th_min_d, s_min_d, z, z_half, t_d = read_in_minmax(kmax_plot, path, filename)
        path = os.path.join(path_triple, id_list_t[d], 'data_analysis')
        w_max_t, th_min_t, s_min_t, z, z_half, t_t = read_in_minmax(kmax_plot, path, filename)

        al = 1.-d*1./(len(d_range)+1)
        # w_max_ss, th_min_ss, z, z_half, t_ = read_in_minmax(kmax_plot, path, filename_ss)
        # w_max_sd, th_min_sd, z, z_half, t_ = read_in_minmax(kmax_plot, path, filename_sd)
        # w_max_st, th_min_st, z, z_half, t_ = read_in_minmax(kmax_plot, path, filename_st)
        #     path = os.path.join(path_double, id_list_d[d], 'data_analysis')
        #     w_max_d, th_min_d, z, z_half, t_ = read_in_minmax(kmax_plot, path, filename)
        #     path = os.path.join(path_triple, id_list_t[d], 'data_analysis')
        #     w_max_t, th_min_t, z, z_half, t_ = read_in_minmax(kmax_plot, path, filename)
        if d > 0:
            lbl_s = ''
            lbl_d = ''
            lbl_t = ''
        else:
            lbl_s = 'single CP gust front'
            lbl_d = 'double CP collision'
            lbl_t = 'triple CP collision'

        it_ini = np.int(t_ini[d]/dt_fields)
        it_2CP = np.int(t_2CP[d]/dt_fields)
        it_3CP = np.int(t_3CP[d]/dt_fields)
        it_final = np.int(t_final[d]/dt_fields)

        ax = axis[0,0]
        ax.plot(np.amax(w_max_s[it_ini:,:], axis=0), z, color=colorlist3[0], alpha=al, label=lbl_s)
        ax.plot(np.amax(w_max_d[it_ini:,:], axis=0), z, color=colorlist3[1], alpha=al, label=lbl_d)
        ax.plot(np.amax(w_max_t[it_ini:,:], axis=0), z, color=colorlist3[2], alpha=al, label=lbl_t)
        ax = axis[0,1]
        ax.plot(np.amax(w_max_s[it_ini:it_2CP,:], axis=0), z, color=colorlist3[0], alpha=al, label=lbl_s)
        ax.plot(np.amax(w_max_d[it_ini:it_2CP,:], axis=0), z, color=colorlist3[1], alpha=al, label=lbl_d)
        ax.plot(np.amax(w_max_t[it_ini:it_2CP,:], axis=0), z, color=colorlist3[2], alpha=al, label=lbl_t)
        ax = axis[0,2]
        ax.plot(np.amax(w_max_s[it_2CP:it_3CP,:], axis=0), z, color=colorlist3[0], alpha=al, label=lbl_s)
        ax.plot(np.amax(w_max_d[it_2CP:it_3CP,:], axis=0), z, color=colorlist3[1], alpha=al, label=lbl_d)
        ax.plot(np.amax(w_max_t[it_2CP:it_3CP,:], axis=0), z, color=colorlist3[2], alpha=al, label=lbl_t)
        ax = axis[0,3]
        ax.plot(np.amax(w_max_s[it_3CP:it_final,:], axis=0), z, color=colorlist3[0], alpha=al, label=lbl_s)
        ax.plot(np.amax(w_max_d[it_3CP:it_final,:], axis=0), z, color=colorlist3[1], alpha=al, label=lbl_d)
        ax.plot(np.amax(w_max_t[it_3CP:it_final,:], axis=0), z, color=colorlist3[2], alpha=al, label=lbl_t)

        ax = axis[1, 0]
        ax.plot(np.amax(th_min_s[it_ini:, :], axis=0), z_half, color=colorlist3[0], alpha=al, label=lbl_s)
        ax.plot(np.amax(th_min_d[it_ini:, :], axis=0), z_half, color=colorlist3[1], alpha=al, label=lbl_d)
        ax.plot(np.amax(th_min_t[it_ini:, :], axis=0), z_half, color=colorlist3[2], alpha=al, label=lbl_t)
        ax = axis[1, 1]
        ax.plot(np.amax(th_min_s[it_ini:it_2CP, :], axis=0), z_half, color=colorlist3[0], alpha=al, label=lbl_s)
        ax.plot(np.amax(th_min_d[it_ini:it_2CP, :], axis=0), z_half, color=colorlist3[1], alpha=al, label=lbl_d)
        ax.plot(np.amax(th_min_t[it_ini:it_2CP, :], axis=0), z_half, color=colorlist3[2], alpha=al, label=lbl_t)
        ax = axis[1, 2]
        ax.plot(np.amax(th_min_s[it_2CP:it_3CP, :], axis=0), z_half, color=colorlist3[0], alpha=al, label=lbl_s)
        ax.plot(np.amax(th_min_d[it_2CP:it_3CP, :], axis=0), z_half, color=colorlist3[1], alpha=al, label=lbl_d)
        ax.plot(np.amax(th_min_t[it_2CP:it_3CP, :], axis=0), z_half, color=colorlist3[2], alpha=al, label=lbl_t)
        ax = axis[1, 3]
        ax.plot(np.amax(th_min_s[it_3CP:it_final, :], axis=0), z_half, color=colorlist3[0], alpha=al, label=lbl_s)
        ax.plot(np.amax(th_min_d[it_3CP:it_final, :], axis=0), z_half, color=colorlist3[1], alpha=al, label=lbl_d)
        ax.plot(np.amax(th_min_t[it_3CP:it_final, :], axis=0), z_half, color=colorlist3[2], alpha=al, label=lbl_t)



    # #     ax00.plot(np.amax(w_max_d, axis=0), z, color=colorlist3[1], alpha=al, label=lbl_d)
    # #     ax00.plot(np.amax(w_max_t, axis=0), z, color=colorlist3[2], alpha=al, label=lbl_t)
    # #     ax01.plot(np.amin(th_min_ss, axis=0), z_half, color=colorlist3[0], alpha=al, label=lbl_s)
    # #     ax01.plot(np.amin(th_min_d, axis=0), z_half, color=colorlist3[1], alpha=al, label=lbl_d)
    # #     ax01.plot(np.amin(th_min_t, axis=0), z_half, color=colorlist3[2], alpha=al, label=lbl_t)
    # #     ax10.plot(np.amax(w_max_sd[:itmax_s, :], axis=0), z, color=colorlist3[0], alpha=al, label=lbl_s)
    # #     ax10.plot(np.amax(w_max_d, axis=0), z, color=colorlist3[1], alpha=al, label=lbl_d)
    # #     ax10.plot(np.amax(w_max_t, axis=0), z, color=colorlist3[2], alpha=al, label=lbl_t)
    # #     ax11.plot(np.amin(th_min_s, axis=0), z_half, color=colorlist3[0], alpha=al, label=lbl_s)
    # #     ax11.plot(np.amin(th_min_d, axis=0), z_half, color=colorlist3[1], alpha=al, label=lbl_d)
    # #     ax11.plot(np.amin(th_min_t, axis=0), z_half, color=colorlist3[2], alpha=al, label=lbl_t)

    for ax in axis[0,:].flat:
        ax.set_xlim(0,4)
        ax.set_xlabel('max. w  [m/s]')
        # ax.set_ylabel('height z  [m]')
    for ax in axis[1,:].flat:
        ax.set_xlim(299.2,300.1)
        ax.set_xlabel(r'min. $\theta$ [K]')
        # ax.set_ylabel('height z  [m]')
    for ax in axis[:,0].flat:
        ax.set_ylabel('height z  [m]')
    axis[0,0].set_title(r't $>$'+str(t_ini[d]))
    axis[0,1].set_title(str(t_ini[d])+r'$\leq$ t $<$'+str(t_2CP[d]))
    axis[0,2].set_title(str(t_2CP[d])+r'$\leq$ t $<$'+str(t_3CP[d]))
    axis[0,3].set_title(str(t_3CP[d])+r'$\leq$ t $<$'+str(t_final[d]))


    # # # # ax2.grid()
    axis[0,0].legend(loc='upper left', bbox_to_anchor=(0.12, .97),
               fancybox=False, shadow=False, ncol=1, fontsize=12)
    # # # ax2.legend(loc='upper left', bbox_to_anchor=(0.1, -0.1),
    # # #            fancybox=False, shadow=False, ncol=1, fontsize=9)
    # # # plt.suptitle('min/max for ' + var_name, fontsize=21)
    plt.subplots_adjust(bottom=0.1, right=.95, left=0.1, top=0.95, hspace=0.2)
    plt.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    return
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

