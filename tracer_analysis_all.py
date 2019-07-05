import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os

label_size = 10
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['lines.markersize'] = 5

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path_root")
    parser.add_argument("dTh", type=int)
    parser.add_argument("--zparams_r1km", nargs='+', type=int)
    parser.add_argument("--rparams_r1km", nargs='+', type=int)
    parser.add_argument("--zparams", nargs='+', type=int)
    parser.add_argument('--rparams', nargs='+', type=int)
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--k0")
    args = parser.parse_args()

    # level of tracers
    if args.k0:
        k0 = np.int(args.k0)
    else:
        k0 = 0

    global cm_bwr, cm_grey, cm_vir, cm_hsv, cm_grey2
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')
    cm_grey2 = plt.cm.get_cmap('bone_r')
    cm_fall = plt.cm.get_cmap('winter')
    cm_summer = plt.cm.get_cmap('spring')

    dt_fields = 100

    dTh, z_params, r_params = set_input_parameters(args)
    id0 = 'dTh' + str(dTh) + '_z' + str(z_params[0]) + '_r' + str(r_params[0])
    fullpath_in = os.path.join(path_root, id0, 'tracer_k' + str(k0), 'output')
    n_tracers = get_number_tracers(fullpath_in)
    n_cps = get_number_cps(fullpath_in)
    print('number of CPs: ', n_cps)
    print('number of tracers per CP: ', n_tracers)
    print ''
    nt = len(times)

    krange = [0]
    nk = len(krange)
    k0 = 0

    # --------------------------------------
    ''' ---------------- for each dTh ---------------- '''
    ''' (a) read in data from tracer output (text-file)'''
    n_params = len(z_params)
    dist_av = np.zeros((n_params, nt, nk))
    r_av = np.zeros((n_params, nt, nk))
    drdt_av = np.zeros((n_params, nt, nk))
    U_rad_av = np.zeros((n_params, nt, nk))
    dU_rad_av = np.zeros((n_params, nt, nk))

    for istar in range(n_params):
        zstar = z_params[istar]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        print('id', id)
        fullpath_in = os.path.join(path_root, id, 'tracer_k'+str(k0), 'output')

        # read_in_txtfile(fullpath_in)
        for it, t0 in enumerate(times):
            print('---t0: '+str(t0)+'---', it)
            cp_id = 2
            # get_radius(fullpath_in, it, cp_id)
            dist_av[istar, it, k0], U_rad_av[istar, it, k0] = get_radius_vel(fullpath_in, it, cp_id, n_tracers, n_cps)
        r_av = dist_av * dx[0]
    for it, t0 in enumerate(times[:]):
        if it > 0:
            drdt_av[:,it,:] = 1./dt_fields * (r_av[:,it,:] - r_av[:,it-1,:])
            print('........ t0: '+str(t0)+'---', it, drdt_av[:,it,:])
    print ''
    print ''


    ''' (b) plot r_av, drdt_av, U_rad_av '''
    figname = 'CP_rim_dTh' + str(dTh) + '.png'
    title = 'CP rim (dTh=' + str(dTh) + ')'
    plot_dist_vel(r_av, drdt_av, U_rad_av,
                  [dTh], z_params, r_params, n_params, k0, title, figname)
    figname = 'comparison_drdt_Urad_dTh' + str(dTh) + '.png'
    plot_comparison_drdt_Urad(r_av, drdt_av, U_rad_av,
                              [dTh], z_params, r_params, n_params, k0, title, figname)





    ''' (c) plot normalized radius / velocity'''
    # (i) for vertical velocity from crosssection in 3D field
    # (ii) for azimuthally averaged vertical vleocity
    print('plotting normalized')
    for istar in range(n_params):
        zstar = z_params[istar]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        print('id', id)
        figname_norm = 'CP_rim_normalized_' + id + '.png'
        plot_vel_normalized(r_av, times, istar, k0, id, figname_norm)
        # figname_norm = 'CP_rim_normalized_' + id + '_av.png'
        # plot_vel_normalized_w_av(r_av, times, istar, k0, id, figname_norm)





    # # -----------------------------------------------
    # ''' ---------------- r = 1km ---------------- '''
    # # dTh_params = [2, 3, 4]
    # # z_params = args.zparams_r1km
    # # r_params = args.rparams_r1km
    # # print('r=1km')
    # # print('dTh: ', dTh_params)
    # # print('z*: ', z_params)
    # # print('r*: ', r_params)
    # # print ''
    # # n_params = len(dTh_params)
    # # '''(a) read in data from textfiles'''
    # # dist_av = np.zeros((n_params, nt, nk))
    # # r_av = np.zeros((n_params, nt, nk))
    # # drdt_av = np.zeros((n_params, nt, nk))
    # # U_rad_av = np.zeros((n_params, nt, nk))
    # # dU_rad_av = np.zeros((n_params, nt, nk))
    # # for istar in range(n_params):
    # #     dTh = dTh_params[istar]
    # #     zstar = z_params[istar]
    # #     rstar = r_params[istar]
    # #     id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
    # #     print('id', id)
    # #     fullpath_in = os.path.join(path_root, id, 'tracer_k'+str(k0), 'output')
    # #     for it, t0 in enumerate(times):
    # #         print('---t0: '+str(t0)+'---', it)
    # #         cp_id = 1
    # #         dist_av[istar, it, k0], U_rad_av[istar, it, k0] = get_radius_vel(fullpath_in, it, cp_id, n_tracers, n_cps)
    # #     r_av = dist_av * dx[0]
    # #     for it, t0 in enumerate(times[1:]):
    # #         drdt_av[:,it,:] = 1./dt_fields * (r_av[:,it,:] - r_av[:,it-1,:])
    # #     print''
    # #
    # # ''' (b) plot r_av, dtdt_av, U_rad_av'''
    # # # print('path', path_out_figs)
    # # # figname = 'CP_rim_r1km.png'
    # # # title = 'CP rim (r=1km)'
    # # # plot_dist_vel(r_av, drdt_av, U_rad_av, dTh_params, z_params, r_params, n_params, k0, title, figname)
    # # # figname = 'CP_rim_vel_r1km.png'
    # # # plot_vel(r_av, U_rad_av, dTh_params, z_params, r_params, n_params, k0, figname)
    #


    return


# ----------------------------------------------------------------------
def plot_vel_normalized_w_av(r_av, times, istar, k0, id, fig_name):
    # read in azimuthally averaged vertical velocity
    rootgrp = nc.Dataset(os.path.join(path_root, id, 'data_analysis', 'stats_radial_averaged.nc'))
    ts_grp = rootgrp.groups['timeseries']
    t_out = ts_grp.variables['time'][:]
    stats_grp = rootgrp.groups['stats']
    w = stats_grp.variables['w'][:,:,:] # dimensions (nt, nr, nz)
    nr = stats_grp.dimensions['nr'].size
    n_times = ts_grp.dimensions['nt'].size
    rootgrp.close()

    # defining geometry for single CP
    x_array = np.arange(0, nx) * dx[0]
    ic = nx / 2
    # r_array = x_array - dx[0] * ic
    r_array = np.arange(0,nr)*dx[0]

    # plotting parameters
    imin = 0
    imax = 80
    k0 = 0

    fig, axes = plt.subplots(2, 2, sharex='none', figsize=(20, 12))
    ax0 = axes[0, 0]
    ax1 = axes[1, 0]
    ax2 = axes[0, 1]
    ax3 = axes[1, 1]
    ax0.plot(r_av[istar, :, k0], times, '-o', label=id)
    ax0.set_xlim(r_array[imin], r_array[imax])
    ax0.grid()
    it = 0
    for t0 in t_out[1::2]:
        if t0 >= np.int(times[0]) and t0 <= times[-1]:
            count = np.double(t0) / times[-1]
            print '>> it', it, t0, count
            iR = r_av[istar, it, k0] / dx[0]
            ax1.plot(r_array[imin:imax], w[it, imin:imax, k0], '-', color=cm_grey2(count),
                     label='t=' + str(t0) + 's')
            ax1.plot(r_array[iR], w[it, iR, k0], 'd', color=cm_grey2(count), markersize=12)
            R_array = r_array / r_av[istar, it, k0]
            ax2.plot(R_array[imin:imax], w[it, imin:imax, k0], '-', color=cm_grey2(count), linewidth=3,
                     label='t=' + str(t0) + 's')
            ax3.plot(R_array[imin:imax], w[it, imin:imax, k0], '-', color=cm_grey2(count), linewidth=3,
                     label='t=' + str(t0) + 's')
            ax3.plot(R_array[iR], w[it, iR, k0], 'd', color=cm_grey(count), markersize=12)
            it += 1
    fig.suptitle(id)
    ax0.legend(loc=3)
    ax1.legend(loc='lower center', bbox_to_anchor=(0.5, 0.1),
               fancybox=True, shadow=True, ncol=6, fontsize=11)
    ax1.set_xlim(r_array[imin], r_array[imax])
    ax2.set_xlim([0., 1.5])
    ax3.set_xlim([0.75, 1.25])
    ax3.set_ylim([-1.1, 2.])
    ax0.set_xlabel('r_av  [m]')
    ax1.set_xlabel('r [m]')
    ax2.set_xlabel('r/r_av [m]')
    ax3.set_xlabel('r/r_av [m]')
    ax0.set_ylabel('times [s]')
    ax1.set_ylabel('w  [m/s]')
    ax2.set_ylabel('w  [m/s]')
    ax3.set_ylabel('w  [m/s]')
    ax2.set_title('w normalized by tracer radius')
    ax3.set_title('w normalized by tracer radius')
    ax1.grid()
    ax2.grid()
    ax3.grid()
    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    return


def plot_vel_normalized(r_av, times, istar, k0, id, fig_name):
    print('')
    print os.path.join(path_root, id, 'fields_merged', 'fields_allt_yz_i200.nc')
    # read in vertical velocity for all times
    rootgrp = nc.Dataset(os.path.join(path_root, id, 'fields_merged', 'fields_allt_yz_i200.nc'))
    w = rootgrp.variables['w'][:, :, :]
    t_out = rootgrp.variables['time'][:]
    rootgrp.close()

    # defining geometry for single CP
    x_array = np.arange(0, nx)*dx[0]
    ic = nx/2
    r_array = x_array - dx[0]*ic


    #plotting parameters
    imin = ic - 5
    imax = ic + 80

    fig, axes = plt.subplots(2, 2, sharex='none', figsize=(20, 12))
    ax0 = axes[0,0]
    ax1 = axes[1,0]
    ax2 = axes[0,1]
    ax3 = axes[1,1]
    ax0.plot(r_av[istar, :, k0], times, '-o', label=id)
    ax0.set_xlim(r_array[imin], r_array[imax])
    ax0.grid()
    it = 0
    for t0 in t_out[0::2]:
        if t0 >= np.int(times[0]) and t0 <= times[-1]:
            count = np.double(t0) / times[-1]
            print '>> it', it, t0, count
            iR = r_av[istar, it, k0] / dx[0] + ic
            ax1.plot(r_array[imin:imax], w[it, imin:imax, k0], '-', color=cm_grey2(count),
                     label='t=' + str(t0) + 's')
            ax1.plot(r_array[iR], w[it, iR, k0], 'd', color=cm_grey2(count), markersize=12)
            R_array = r_array / r_av[istar, it, k0]
            ax2.plot(R_array[imin:imax], w[it, imin:imax, k0], '-', color=cm_grey2(count), linewidth=3,
                     label='t=' + str(t0) + 's')
            ax3.plot(R_array[imin:imax], w[it, imin:imax, k0], '-', color=cm_grey2(count), linewidth=3,
                     label='t=' + str(t0) + 's')
            ax3.plot(R_array[iR], w[it, iR, k0], 'd', color=cm_grey(count), markersize=12)
            it += 1
    fig.suptitle(id)
    ax0.legend(loc=3)
    ax1.legend(loc='lower center', bbox_to_anchor=(0.5, 0.1),
               fancybox=True, shadow=True, ncol=6, fontsize=11)
    # ax2.legend()
    ax1.set_xlim(r_array[imin], r_array[imax])
    ax2.set_xlim([-0.2, 1.2])
    ax3.set_xlim([0.75, 1.2])
    ax3.set_ylim([-1.1, 2.])
    ax0.set_xlabel('r_av  [m]')
    ax1.set_xlabel('r [m]')
    ax2.set_xlabel('r/r_av [m]')
    ax3.set_xlabel('r/r_av [m]')
    ax0.set_ylabel('times [s]')
    ax1.set_ylabel('w  [m/s]')
    ax2.set_ylabel('w  [m/s]')
    ax3.set_ylabel('w  [m/s]')
    ax2.set_title('w normalized by tracer radius')
    ax3.set_title('w normalized by tracer radius')
    ax1.grid()
    ax2.grid()
    ax3.grid()
    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    return


# ----------------------------------------------------------------------
def plot_dist_vel(r_av, drdt_av, U_rad_av, dTh_params, z_params, r_params, n_params, k0,
                  title, fig_name):

    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, sharex='all', figsize=(18, 5))
    for istar in range(n_params):
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[istar]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        ax0.plot(times, r_av[istar, :, k0], 'o-', label=id)
        ax1.plot(times[1:], drdt_av[istar, 1:, k0], 'o-', label=id)
        ax2.plot(times, U_rad_av[istar, :, k0], 'o-', label=id)
    # ax0.set_title('r_av')
    # ax1.set_title('drdt_av')
    # ax2.set_title('U_av')
    ax0.set_xlabel('times [s]')
    ax1.set_xlabel('times [s]')
    ax2.set_xlabel('times [s]')
    ax0.set_ylabel('r_av  [m]')
    ax1.set_ylabel('drdt_av')
    ax2.set_ylabel('U_rad_av  [m/s]')
    ax1.legend()
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return

# ----------------------------------------------------------------------
def plot_comparison_drdt_Urad(r_av, drdt_av, U_rad_av,
                              dTh_params, z_params, r_params, n_params, k0,
                            title, fig_name):
    fig, axis = plt.subplots(2, n_params, sharex='all', sharey='row', figsize=(6*n_params, 10))
    print n_params
    colorlist = ['navy', 'blue', 'forestgreen']
    for istar in range(n_params):
        print istar
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[istar]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)

        ax2 = axis[0,istar].twinx()
        ax2.plot(times, r_av[istar, :, k0], '-o', color='gray',
                 linewidth=1, markersize=3, markeredgecolor='gray',
                 label=str(r_av[istar,:2,k0]))
        axis[0,istar].plot(times[:-1], drdt_av[istar, 1:, k0], 'o-', color=colorlist[0], label='dR/dt')
        axis[0,istar].plot(times[:], drdt_av[istar, :, k0], 'o-', color=colorlist[1], label='dR/dt')
        axis[0,istar].plot(times, U_rad_av[istar, :, k0], 'o-', color=colorlist[2], label='U_rad')
        axis[0,istar].set_title(id)
        axis[0,istar].legend(loc='best')
        axis[0,istar].set_ylabel('velocity [m/s]')

        axis[1, istar].plot(times[:-1], drdt_av[istar, 1:, k0]/U_rad_av[istar, :-1, k0],
                            color=colorlist[0], label='frac')
        axis[1, istar].plot(times[:], drdt_av[istar, :, k0]/U_rad_av[istar, :, k0], color=colorlist[1], label='frac')
        axis[1, istar].legend(loc='best')
        ax2.legend(loc='best')

        axis[1,istar].set_xlabel('time  [s]')
        axis[1,istar].set_ylabel('fraction [-]')

    fig.suptitle(title)
    fig.tight_layout()
    plt.subplots_adjust(bottom=0.08, right=.95, left=0.07, top=0.9, wspace=0.25)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return

# ----------------------------------------------------------------------
def plot_vel(r_av, U_rad_av, dTh_params, z_params, r_params, n_params, k0, fig_name):

    fig, (ax0, ax1) = plt.subplots(1, 2, sharex='all', figsize=(11, 5))
    for istar in range(n_params):
        if len(dTh_params) == 1:
            dTh = dTh_params[0]
        else:
            dTh = dTh_params[istar]
        zstar = z_params[istar]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        ax0.plot(times, r_av[istar, :, k0], 'o-', label=id)
        ax1.plot(times, U_rad_av[istar, :, k0], 'o-', label=id)
    ax0.set_xlabel('times [s]')
    ax1.set_xlabel('times [s]')
    ax0.set_ylabel('r_av')
    ax1.set_ylabel('U_rad_av')
    ax1.legend()
    fig.suptitle('CPs (r=1km)')
    fig.tight_layout()
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, wspace=0.25)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return

# ----------------------------------------------------------------------
def read_in_txtfile(fullpath_in):
    f = open(fullpath_in+'/coldpool_tracer_out.txt', 'r')
    # f = open(DIR+EXPID+'/'+child+'/output/irt_tracks_output_pure_sort.txt', 'r')
    lines = f.readlines()
    column =lines[0].split()
    # print column
    start_c = int(column[1])   #get the first timestep when tracking starts
    print 'precip starts at timestep', start_c
    column = lines[-1].split()
    end_c = int(column[1])
    print 'last timestep in precip tracking', end_c
    for line in lines:
        columns = line.split()
        tist = (int(columns[0]))    # timestep
        # age = (int(columns[1]))    # age = timestep - tstart
        tracerID = (int(columns[2]))
        cpID = (int(columns[3]))
        # position tracer: columns 4,5
        # rounded position tracer: columns 6, 7
        dist = (float(columns[8])) # distance of tracer from COG
        # angle = (float(columns[9])) # angle of vector (position tracer)-COG
        # u =(float(columns[10]))  # u- wind component
        # v =(float(columns[11]))  # v- wind component
        vel_rad =(float(columns[12]))  # radial wind component
        # vel_tang =(float(columns[13]))  # tangential wind component
        # dist_x = (float(columns[14])) # x-component of distance from COG
        # dist_y = (float(columns[15])) # x-component of distance from COG
        COGx =(float(columns[16]))  # center of CP (and initialized circle)
        COGy = (float(columns[17]))

        if tracerID == 1:
            print('t', tist)
            print('vel_rad', vel_rad, dist)

    f.close()

    return



def get_radius_vel(fullpath_in, t0, cp_id, n_tracers, n_cps):
    # print('in', fullpath_in)
    f = open(fullpath_in + '/coldpool_tracer_out.txt', 'r')
    # f = open(DIR+EXPID+'/'+child+'/output/irt_tracks_output_pure_sort.txt', 'r')
    lines = f.readlines()
    count = 0
    dist = []
    vel = []

    count = t0 * n_cps * n_tracers + (cp_id - 1)*n_tracers
    # while CP age is 0 and CP ID is cp_id
    timestep = int(lines[count].split()[0])
    cp_ID = int(lines[count].split()[3])
    # print(timestep, cp_ID)
    while (timestep-1 == t0 and int(lines[count].split()[3])==cp_id):
        columns = lines[count].split()
        dist.append(float(columns[8]))
        # vel.append(np.sqrt(float(columns[10])**2 + float(columns[11])**2))
        vel.append(float(columns[12]))
        count += 1
        timestep = int(lines[count].split()[0])
    f.close()
    r_av = np.average(dist)
    vel_av = np.average(vel)

    return r_av, vel_av


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
    # count = 0
    # # while CP age is 0 and CP ID is 1
    # while (int(lines[count].split()[0]) == 1):
    #     count += 1
    # cp_number = int(lines[count-1].split()[3])
    cp_number = int(lines[-1].split()[3])

    f.close()

    return cp_number


# ----------------------------------------------------------------------

def set_input_parameters(args):
    print('--- set input parameters ---')
    global case_name
    global path_root, path_out_figs
    global times

    path_root = args.path_root
    path_out_figs = os.path.join(path_root, 'figs_tracers')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print ''
    print 'path_root: ', path_root
    print 'path figs: ', path_out_figs
    print ''

    dTh = args.dTh
    if args.zparams:
        z_params = args.zparams
    else:
        if dTh == 1:
            z_params = [3465, 1730, 1155]  # run1
        elif dTh == 2:
            # z_params = [2450, 1225, 815]  # run1
            z_params = [500, 900, 1600, 1900, 2500]  # run2, run3
            r_params_ = [1900, 1300, 900, 800, 600]  # run2, run3
        elif dTh == 3:
            # z_params = [4000, 2000, 1500, 1000, 670, 500, 250] # run1
            z_params = [500, 1000, 1600, 2000, 2500]  # run2, run3
            r_params_ = [1500, 1000, 700, 600, 500]  # run2, run3
            z_params = [500, 1000, 2000]  # run2
            r_params_ = [1500, 1000, 600]  # run2
            z_params = [1000]  # run2
            r_params_ = [1000]  # run2
        elif dTh == 4:
            # z_params = [1730, 870, 430]     # run1
            z_params = [500, 900, 1600, 2000, 2500]  # run2, run3
            r_params_ = [1300, 900, 600, 500, 400]  # run2, run3
    if args.rparams:
        r_params = args.rparams
    else:
        try:
            r_params = r_params_
            del r_params_
        except:
            r_params = z_params[::-1]


    print('dTh: ', dTh)
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
    print ''


    return dTh, z_params, r_params


# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()
