import netCDF4 as nc
import argparse
import json as simplejson
import os
import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import matplotlib.cbook as cbook
import time



# from convert_fields_smaller_k import convert_file_for_varlist_horsection





execfile('settings.py')


def main():
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    # parser.add_argument("casename")
    # parser.add_argument("path_3CP")
    parser.add_argument("--level")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    global case_name_3CP
    case_name_2CP = 'ColdPoolDry_double_3D'
    case_name_3CP = 'ColdPoolDry_triple_3D'

    dTh = 5
    zstar = 1000
    rstar = 2000
    sep = 12
    case = 'dTh'+str(dTh)+'_z'+str(zstar)+'_r'+str(rstar)+'_d'+str(sep)+'km'

    path_2CP = '/nbi/ac/coag/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/double_3D/'+case
    path_3CP = '/nbi/ac/cond2/meyerbe/ColdPools/3D_sfc_fluxes_off/triple_3D/'+case
    path_data = os.path.join(path_3CP, 'data_analysis')
    path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single'
    nml_2CP, nml_3CP, times, files = set_input_output_parameters(args, case_name_2CP, case_name_3CP, path_2CP, path_3CP)
    ic_arr_2CP, jc_arr_2CP, ic_arr_3CP, jc_arr_3CP, \
    ic_2CP, jc_2CP, ic_3CP, jc_3CP = define_geometry(nml_2CP, nml_3CP)
    print('CP geometry:')
    # print(ic_arr)
    # print(jc_arr)

    if args.level:
        zlevel = np.int(args.level)
    else:
        zlevel = 1000

    d = path_3CP[-4:-2]
    print('d='+d)
    print('')

    time_windows = {}
    time_windows['10'] = [600, 750, 3600]
    time_windows['12'] = [700, 1100, 3600]
    time_windows['15'] = [1200, 1600, 3600]
    time_windows['20'] = [1200, 1600, 3600]


    # read in flux (2D-field)
    filename_data_2CP = 'mass_flux_z'+str(zlevel)+'.nc'
    filename_data_3CP = 'mass_flux_z'+str(zlevel)+'.nc'
    root_in = nc.Dataset(os.path.join(path_3CP, 'data_analysis', filename_data_3CP), 'r')
    mass_flux_3CP = root_in.groups['fields_2D'].variables['mass_flux_2D'][:,:,:]
    mass_flux_pos_3CP = root_in.groups['fields_2D'].variables['mass_flux_2D_positive'][:,:,:]
    time_data_3CP = root_in.groups['timeseries'].variables['time'][:]
    root_in.close()
    root_in = nc.Dataset(os.path.join(path_2CP, 'data_analysis', filename_data_2CP), 'r')
    mass_flux_2CP = root_in.groups['fields_2D'].variables['mass_flux_2D'][:, :, :]
    mass_flux_pos_2CP = root_in.groups['fields_2D'].variables['mass_flux_2D_positive'][:, :, :]
    time_data_2CP = root_in.groups['timeseries'].variables['time'][:]
    root_in.close()

    # averaged mass_flux_3CP
    delta = 6
    itmax = 15
    MF_2CP = np.sum(mass_flux_2CP[:itmax, :, :], axis=0)
    MF_mean_2CP = np.mean(MF_2CP[ic_2CP - delta:ic_2CP + delta + 1, :], axis=0)
    MF_3CP = np.sum(mass_flux_3CP[:itmax, :, :], axis=0)
    MF_mean_3CP = np.mean(MF_3CP[:, jc_3CP - delta:jc_3CP + delta + 1], axis=1)



    fig_name = 'collisions_massflux_' + case + '.png'
    ncol = 3
    vmin = 1e-2
    vmax = 6.5
    fig, axis = plt.subplots(1, ncol, figsize=(ncol * 5, 5))
    ax = axis[0]
    pcm = ax.pcolormesh(np.arange(nx_2CP[1])-jc_2CP, np.arange(nx_2CP[0]), MF_2CP,
                        # norm = colors.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-max,vmax=max),
                        # norm = colors.PowerNorm(),
                        norm=colors.LogNorm(vmin=vmin, vmax=vmax), cmap=cm_bw)  # cmap='RdBu_r')
    plt.colorbar(pcm, ax=ax, extend='both')
    imin = 100
    jmin = 250
    jmax = nx_2CP[0] - jmin
    rect = mpatches.Rectangle((-imin, ic_2CP-delta), 200, 2*delta, fill=True,
                              linewidth=0, edgecolor='r', facecolor='lightyellow', alpha=0.2)
    ax.add_patch(rect)
    ax.set_xlim(-imin, imin)
    ax.set_ylim(jmin, jmax)
    ax.set_aspect('equal')

    ax = axis[1]
    ax.pcolormesh(np.arange(nx_3CP[0])-ic_arr_3CP[0], np.arange(nx_3CP[1]), MF_3CP.T,
                        # norm = colors.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-max,vmax=max),
                        # norm = colors.PowerNorm(),
                        norm=colors.LogNorm(vmin=vmin, vmax=vmax), cmap=cm_bw)  # cmap='RdBu_r')
    imin = 170
    imax = nx_3CP[0] - imin
    rect2 = mpatches.Rectangle((-100, jc_arr_3CP[2]-delta), 300, 2*delta, fill=True,
                               linewidth=0, edgecolor='r', facecolor='lightyellow', alpha=0.2)
    ax.add_patch(rect2)
    ax.set_xlim(-100, 200)
    ax.set_ylim(imin, imax)
    ax.set_aspect('equal')

    ax = axis[2]
    ax.plot(np.arange(nx_2CP[1]) - jc_2CP, MF_mean_2CP, label='2-CP', color=colorlist2[0])
    ax.plot(np.arange(nx_3CP[0]) - ic_arr_3CP[0], MF_mean_3CP, label='3-CP', color=colorlist2[1])
    ax.legend()
    ax.set_xlim(-150, 250)
    ax.set_ylim(-2,6)
    # imin = 0.5 * (np.maximum(nx_2CP[1], nx_3CP[0]) - np.minimum(nx_2CP[1], nx_3CP[0]))
    # if nx_2CP[1] >= nx_3CP[0]:
    #     # ax.plot(MF_mean_3CP)
    #     ax.plot(MF_mean_3CP)
    #     ax.plot(np.arange(nx_3CP[0])-ic_arr_3CP[0], MF_mean_3CP)
    #     # ax.plot(MF_mean_2CP[imin:nx_2CP[1]-imin])
    # else:
    #     ax.plot(MF_mean_3CP[imin:nx_3CP[0]-imin])
    #     # ax.plot(MF_mean_2CP)
    for ax in axis.flatten():
        ax.set_xlabel('x')
    axis[0].set_ylabel('y')
    axis[2].set_ylabel(r'Integrated Mass Flux  [kg/m$^2$]')
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, hspace=0.2, wspace=0.2)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)





    fig_name = 'collisions_massflux_2CP_3CP.png'
    itmax = 15

    times = time_windows[d]
    ncol = 3
    fig, axis = plt.subplots(2, ncol, figsize=(ncol * 5, 2*5))

    # cf = axis[0].contourf(MF_3CP.T, levels=lvls, cmap=cm_contourfigs, extend='both')
    lvls = np.log(np.linspace(0,6,12))
    ax = axis[0,0]
    cf = ax.contourf(MF_3CP.T, norm=colors.LogNorm(), levels=lvls, cmap=cm_bw)
    ax.set_title('tmax='+str(itmax*100))
    ax.plot([0,nx_3CP[0]],[jc_arr_3CP[2],jc_arr_3CP[2]],'r')
    rect = mpatches.Rectangle((200,jc_arr_3CP[2]-delta), 100, 2*delta, fill=False, linewidth=2, edgecolor='r', facecolor='white')
    ax.add_patch(rect)
    plt.colorbar(cf, ax=ax, shrink=0.8, extend='both')
    for i,ic_ in enumerate(ic_arr_3CP):
        ax.plot(ic_arr_3CP[i],jc_arr_3CP[i], 'o', markersize=10, label='i='+str(i))
    ax.legend()

    ax = axis[0,1]
    pcm = ax.pcolormesh(np.arange(nx_3CP[0]), np.arange(nx_3CP[1]), MF_3CP.T,
                             # norm = colors.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-max,vmax=max),
                             # norm = colors.PowerNorm(),
                             norm = colors.LogNorm(vmin=1e-2, vmax=1e1),
                             cmap=cm_bw)#cmap='RdBu_r')
    rect2 = mpatches.Rectangle((200,jc_arr_3CP[2]-delta), 100, 2*delta, fill=False, linewidth=2, edgecolor='r', facecolor='white')
    ax.add_patch(rect2)
    plt.colorbar(pcm, ax=ax, extend='both')

    ax = axis[0,2]
    imin = 0.5*(np.maximum(nx_2CP[1], nx_3CP[0]) - np.minimum(nx_2CP[1], nx_3CP[0]))
    # ax.plot(MF_mean_3CP)
    # ax.plot[]
    ax.plot(np.arange(nx_2CP[1])-jc_2CP, MF_mean_2CP)
    ax.plot(np.arange(nx_3CP[0])-ic_arr_3CP[0], MF_mean_3CP)
    # if nx_2CP[1] >= nx_3CP[0]:
    #     # ax.plot(MF_mean_3CP)
    #     ax.plot(MF_mean_3CP)
    #     ax.plot(np.arange(nx_3CP[0])-ic_arr_3CP[0], MF_mean_3CP)
    #     # ax.plot(MF_mean_2CP[imin:nx_2CP[1]-imin])
    # else:
    #     ax.plot(MF_mean_3CP[imin:nx_3CP[0]-imin])
    #     # ax.plot(MF_mean_2CP)
    imin = 180
    imax = nx_3CP[0] - imin
    for ax in axis[0,:2].flatten():
        ax.set_xlim(imin, imax)
        ax.set_ylim(imin, imax)
        ax.set_aspect('equal')

    ''' 2CP'''
    ax = axis[1, 0]
    ax.set_title('tmax=' + str(itmax * 100))
    cf = ax.contourf(MF_2CP.T, norm=colors.LogNorm(), levels=lvls, cmap=cm_bw)
    plt.colorbar(cf, ax=ax, shrink=0.8, extend='both')
    ax.plot([ic_2CP,ic_2CP],[0, nx_2CP[1]], 'r')
    rect = mpatches.Rectangle((ic_2CP-delta, jc_2CP-50), 2*delta, 100, fill=False, linewidth=2, edgecolor='r',
                              facecolor='white')
    ax.add_patch(rect)
    for i, ic_ in enumerate(ic_arr_2CP):
        ax.plot(ic_arr_2CP[i], jc_arr_2CP[i], 'o', markersize=10, label='i=' + str(i))
    ax.legend()

    ax = axis[1, 1]
    pcm = ax.pcolormesh(np.arange(nx_2CP[0]), np.arange(nx_2CP[1]), MF_2CP.T,
                        # norm = colors.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-max,vmax=max),
                        # norm = colors.PowerNorm(),
                        norm=colors.LogNorm(vmin=1e-2, vmax=1e1),
                        cmap=cm_bw)  # cmap='RdBu_r')
    plt.colorbar(pcm, ax=ax, extend='both')
    # rect2 = mpatches.Rectangle((200, jc_2CP - delta), 100, 2 * delta, fill=False, linewidth=2, edgecolor='r',
    #                            facecolor='white')
    # ax.add_patch(rect2)

    ax = axis[1, 2]
    ax.plot(MF_mean_2CP)

    imin = 200
    imax = nx_2CP[0] - imin
    jmin = 100
    jmax = nx_2CP[1] - jmin
    for ax in axis[1,:2].flatten():
        ax.set_xlim(imin, imax)
        ax.set_ylim(jmin, jmax)
        ax.set_aspect('equal')

    # axis[2].plot()
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, hspace=0.2, wspace=0.1)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    #
    #
    #
    #times = time_windows[d]
    # # fig_name = 'collisions_massflux_accumulated_test.png'
    # # imin = 150
    # # imax = nx - imin
    # # times = [6, 7, 8, 10, 11, 12, 14, 16, 20]
    # # ncol = len(times)
    # # print('mass flux 2d: ', mass_flux_3CP.shape)
    # # fig, axis = plt.subplots(3, ncol, figsize=(ncol * 5, 15), sharey='all')
    # # max = 5.
    # # lvls_cum = np.linspace(-max, max, 21)
    # # lvls = np.linspace(-2, 2, 41)
    # # for i,itmax in enumerate(times):
    # #     axis[0,i].set_title('t='+str(itmax*100))
    # #     cf = axis[0, i].contourf(np.sum(mass_flux_3CP[:itmax, :, :], axis=0).T, levels=lvls_cum, cmap=cm_contourfigs, extend='both')
    # #     # cf = axis[0, i].contourf(np.sum(mass_flux_3CP[:itmax, :, :], axis=0).T, norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
    # #     #                                       vmin=-1.0, vmax=1.0))
    # #     cb = plt.colorbar(cf, ax=axis[0, i], shrink=0.8)
    # #     cf = axis[1, i].contourf(mass_flux_3CP[itmax, :, :].T, levels=lvls, cmap=cm_contourfigs, extend='both')
    # #     cb = plt.colorbar(cf, ax=axis[1, i], shrink=0.8)
    # #     cf = axis[2, i].contourf(np.sum(mass_flux_pos_3CP[:itmax, :, :], axis=0).T, levels=lvls_cum, cmap=cm_contourfigs, extend='both')
    # #     # cf = axis[2, i].contourf(np.sum(mass_flux_pos_3CP[:itmax, :, :], axis=0).T, norm=colors.LogNorm(vmin=0.1, vmax=max))
    # #     cb = plt.colorbar(cf, ax=axis[2, i], shrink=0.8)
    # # for ax in axis.flatten():
    # #     ax.set_xlim(imin, imax)
    # #     ax.set_ylim(imin, imax)
    # #     ax.set_aspect('equal')
    # # plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, hspace=0.2, wspace=0.1)
    # # fig.savefig(os.path_3CP.join(path_out_figs, fig_name))
    # # plt.close(fig)
    # #
    # # fig_name = 'collisions_massflux_accumulated_testlog.png'
    # # imin = 150
    # # imax = nx - imin
    # # times = [6, 7, 8, 10, 11, 12, 14, 16, 20]
    # # ncol = len(times)
    # # print('mass flux 2d: ', mass_flux_3CP.shape)
    # # fig, axis = plt.subplots(3, ncol, figsize=(ncol * 5, 15), sharey='all')
    # # max = 5.
    # # lvls_cum = np.linspace(-max, max, 21)
    # # lvls = np.linspace(-2, 2, 41)
    # # for i, itmax in enumerate(times):
    # #     axis[0, i].set_title('t=' + str(itmax * 100))
    # #     # cf = axis[0, i].contourf(np.sum(mass_flux_3CP[:itmax, :, :], axis=0).T, levels=lvls_cum, cmap=cm_contourfigs,
    # #     #                          extend='both')
    # #     cf = axis[0, i].contourf(np.sum(mass_flux_3CP[:itmax, :, :], axis=0).T,
    # #                              norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
    # #                                                     vmin=-1.0, vmax=1.0))
    # #     cb = plt.colorbar(cf, ax=axis[0, i], shrink=0.8)
    # #     cf = axis[1, i].contourf(mass_flux_3CP[itmax, :, :].T, levels=lvls, cmap=cm_contourfigs, extend='both')
    # #     cb = plt.colorbar(cf, ax=axis[1, i], shrink=0.8)
    # #     # cf = axis[2, i].contourf(np.sum(mass_flux_pos_3CP[:itmax, :, :], axis=0).T, levels=lvls_cum, cmap=cm_contourfigs, extend='both')
    # #     cf = axis[2, i].contourf(np.sum(mass_flux_pos_3CP[:itmax, :, :], axis=0).T, norm=colors.LogNorm(vmin=0.1, vmax=max))
    # #     cb = plt.colorbar(cf, ax=axis[2, i], shrink=0.8)
    # # for ax in axis.flatten():
    # #     ax.set_xlim(imin, imax)
    # #     ax.set_ylim(imin, imax)
    # #     ax.set_aspect('equal')
    # # plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, hspace=0.2, wspace=0.1)
    # # fig.savefig(os.path_3CP.join(path_out_figs, fig_name))
    # # plt.close(fig)
    # #
    # #
    # # fig_name = 'collisions_massflux_accumulated.png'
    # # imin = 150
    # # imax = nx - imin
    # # times = [7, 11, 15]
    # # ncol = len(times)
    # # print('mass flux 2d: ', mass_flux_3CP.shape)
    # # fig, axis = plt.subplots(2, ncol, figsize=(ncol * 5, 10), sharey='all')
    # # max = 5.
    # # lvls_cum = np.linspace(-max, max, 21)
    # # lvls = np.linspace(-2, 2, 21)
    # # for i, itmax in enumerate(times):
    # #     axis[0, i].set_title('t=' + str(itmax * 100))
    # #     cf = axis[0, i].contourf(np.sum(mass_flux_3CP[:itmax, :, :], axis=0).T, levels=lvls_cum, cmap=cm_contourfigs,
    # #                              extend='both')
    # #     plt.colorbar(cf, ax=axis[0, i], shrink=0.8)
    # #     cf = axis[1, i].contourf(mass_flux_3CP[itmax, :, :].T, levels=lvls, cmap=cm_contourfigs, extend='both')
    # #     plt.colorbar(cf, ax=axis[1, i], shrink=0.8)
    # # for ax in axis.flatten():
    # #     ax.set_xlim(imin, imax)
    # #     ax.set_ylim(imin, imax)
    # #     ax.set_aspect('equal')
    # # plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, hspace=0.2, wspace=0.1)
    # # fig.savefig(os.path_3CP.join(path_out_figs, fig_name))
    # # plt.close(fig)
    #
    #

    return


# ----------------------------------------------------------------------

def set_input_output_parameters(args, case_name_2CP, case_name_3CP, path_2CP, path_3CP):
    print('')
    print('--- set input parameters ---')

    nml_2CP = simplejson.loads(open(os.path.join(path_2CP, case_name_2CP + '.in')).read())
    nml_3CP = simplejson.loads(open(os.path.join(path_3CP, case_name_3CP + '.in')).read())
    global nx_2CP, nx_3CP, dx, dV, gw
    nx_2CP = np.zeros(3, dtype=np.int)
    nx_3CP = np.zeros(3, dtype=np.int)
    nx_2CP[0] = nml_2CP['grid']['nx']
    nx_2CP[1] = nml_2CP['grid']['ny']
    nx_2CP[2] = nml_2CP['grid']['nz']
    nx_3CP[0] = nml_3CP['grid']['nx']
    nx_3CP[1] = nml_3CP['grid']['ny']
    nx_3CP[2] = nml_3CP['grid']['nz']
    dx = np.zeros(3, dtype=np.int)
    dx[0] = nml_3CP['grid']['dx']
    dx[1] = nml_3CP['grid']['dy']
    dx[2] = nml_3CP['grid']['dz']
    gw = nml_3CP['grid']['gw']
    dV = dx[0] * dx[1] * dx[2]
    print('grid 2CP: ', nx_2CP)
    print('grid 3CP: ', nx_3CP)


    global dt_fields, dt_stats, dt
    dt_fields = np.int(nml_3CP['fields_io']['frequency'])
    dt_stats = np.int(nml_3CP['stats_io']['frequency'])
    dt = np.int(np.maximum(dt_stats, dt_fields))

    global tmin, tmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100

    # times = [np.int(name[:-3]) for name in os.listdir(os.path.join(path_3CP, 'fields')) if name[-2:] == 'nc'
    #          and tmin <= np.int(name[:-3]) <= tmax]
    times = [np.int(name[:-3]) for name in os.listdir(os.path.join(path_3CP, 'fields')) if name[-2:] == 'nc'
             and tmin <= np.int(name[:-3]) <= tmax and np.mod(np.int(name[:-3]), dt) == 0]
    times.sort()
    nt = len(times)
    print('times: ' + str(times))
    print('tmin, tmax: ' + str(tmin) + ', ' + str(tmax))
    print('dt: ' + str(dt))
    print('nt:', nt)
    files = [str(t) + '.nc' for t in times]
    print('')

    return nml_2CP, nml_3CP, times, files
# _______________________________________________________

def define_geometry(nml_2CP, nml_3CP):
    global rstar

    '''--- define geometry ---'''
    rstar = nml_2CP['init']['r']
    sep = nml_2CP['init']['sep']
    isep = np.int(np.round(sep / dx[0]))
    jsep = 0
    ic_2CP = np.int(np.round(nx_2CP[0] / 2))
    jc_2CP = np.int(np.round(nx_2CP[1] / 2))
    ic1 = ic_2CP - np.int(np.round(isep / 2))
    jc1 = jc_2CP
    ic2 = ic1 + isep
    jc2 = jc1 + jsep
    ic_arr_2CP = [ic1, ic2]
    jc_arr_2CP = [jc1, jc2]

    try:
        rstar = nml_3CP['init']['r']
    except:
        rstar = 5000.0  # half of the width of initial cold-pools [m]
    try:
        d = nml_3CP['init']['d']
    except:
        d = 1e3
    i_d = np.int(np.round(d / dx[0]))
    idhalf = np.int(np.round(i_d / 2))
    a = np.int(np.round(i_d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
    r_int = np.int(np.sqrt(3.) / 6 * i_d)  # radius of inscribed circle
    # point of 3-CP collision (ic, jc)
    ic_3CP = np.int(np.round(nx_3CP[0] / 2))
    jc_3CP = np.int(np.round(nx_3CP[1] / 2))
    ic1 = ic_3CP - r_int
    ic2 = ic1
    ic3 = ic_3CP + (a - r_int)
    jc1 = jc_3CP - idhalf
    jc2 = jc_3CP + idhalf
    jc3 = jc_3CP
    ic_arr_3CP = [ic1, ic2, ic3]
    jc_arr_3CP = [jc1, jc2, jc3]
    print('')
    print('CP constellations: ')
    print('2CP: ', ic_2CP, jc_2CP, ic_arr_2CP, jc_arr_2CP)
    print('3CP: ', ic_3CP, jc_3CP, ic_arr_3CP, jc_arr_3CP)


    return ic_arr_2CP, jc_arr_2CP, ic_arr_3CP, jc_arr_3CP, ic_2CP, jc_2CP, ic_3CP, jc_3CP

# _______________________________________________________

if __name__ == '__main__':
    main()