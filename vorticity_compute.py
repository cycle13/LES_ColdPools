import numpy as np
import scipy.integrate as integrate  # for simpsons integration
import scipy.integrate as integrate  # for simpsons integration
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os
import sys
import time


'''
COMPUTE VORTICITY
compute_vorticity_yz: computes vorticity on yz-crosssection on Euclidian grid 
    > location: on grid faces, i.e. y- and z-faces of grid boxes
compute_vorticity_yz_staggered: computes vorticity on yz-crosssection on Euclidian grid
    (!! without interpolating v, w onto centered grid, i.e. vorticity location not well defined)
'''
# Comparison of computation from velocity fields on staggered grid vs. interpolated:
#   - big difference at lowest model level (change of sign)
#   - smaller difference at first model level
#   - minor differences at higher levels (k>=2)

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 12


def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')

    timerange, kmax, nml = set_input_output_parameters(args)
    n_CPs = define_geometry(case_name, nml)
    nx_half = np.int(nx_ / 2)
    ny_half = np.int(ny_ / 2)
    icshift = nx_half - 1
    jcshift = ny_half -1
    ic = ic_arr[0]
    jc = jc_arr[0]
    # x_arr = dx*np.arange(0,nx_)
    # y_arr = dy*np.arange(0,ny_)
    # z_arr = dz*np.arange(0,nz)
    # X_, Y_ = np.meshgrid(x_arr, y_arr)
    print('')

    # READ IN density profile
    try:
        rootgrp = nc.Dataset(os.path.join(path, 'stats', 'Stats.' + case_name + '.nc'))
    except:
        rootgrp = nc.Dataset(os.path.join(path, 'Stats.' + case_name + '.nc'))
    rho0 = rootgrp.groups['reference'].variables['rho0'][:]
    rho_unit = rootgrp.groups['reference'].variables['rho0'].units
    z_half = rootgrp.groups['reference'].variables['z'][:]
    rootgrp.close()



    # Note: interpolation for k < 0
    # - centered grid: (ir, jr, kr)
    # - velocities are on staggered grid (is,js,ks) = (ir+1/2, jr+1/2, kr+1/2)
    #           >> u[is,jr,kr] = vel[ir+1/2, jr, kr]
    #           >> v[ir,js,kr] = vel[ir, jr+1/2, kr]
    #           >> w[ir,jr,ks] = vel[ir, jr, kr+1/2]
    #           >> vel[is,js,ks] = vel[ir+1/2, jr+1/2, kr+1/2]
    # - therefore, the vertical velocity at the lowest output level ks=0
    #           is actually at the level kr=1/2 and so w[ks=0]!=0
    # - to guarantee zero vertical velocity at actual surfae, i.e. w[kr=0]=0,
    #           must have reflected field on vertical ghost points:
    #           >> w[kr=-1/2] = -w[ks=1/2]
    #           >> w[ks=-1] = -w[ks=0]
    #           >> w[kr=0] = 0.5*(w[kr=-1/2]+w[kr=1/2]) = 0.5*(w[ks=0]+w[ks=-1]) = 0

    ''' COMPUTE VORTICITY '''
    # create fields file
    fields_file_name = 'field_vort_yz.nc'
    create_vort_field_file(timerange, fields_file_name, kmax, n_CPs)
    for it, t0 in enumerate(timerange):
        print('--- time: t='+str(t0)+'s ---')
        # read in fields
        # s = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))
        # u = read_in_netcdf_fields('u', os.path.join(path_fields, str(t0) + '.nc'))
        # v = read_in_netcdf_fields('v', os.path.join(path_fields, str(t0) + '.nc'))
        # w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0) + '.nc'))
        # u_roll = np.roll(np.roll(u[:, :, :], ishift, axis=0), jshift, axis=1)
        # u_ = u_roll[ic - nx_half + ishift:ic + nx_half + ishift, jc - ny_half + jshift:jc + ny_half + jshift, :]
        # v_roll = np.roll(np.roll(v[:, :, :], ishift, axis=0), jshift, axis=1)
        # v_ = v_roll[ic - nx_half + ishift:ic + nx_half + ishift, jc - ny_half + jshift:jc + ny_half + jshift, :]
        # w_roll = np.roll(np.roll(w[:, :, :], ishift, axis=0), jshift, axis=1)
        # w_ = w_roll[ic - nx_half + ishift:ic + nx_half + ishift, jc - ny_half + jshift:jc + ny_half + jshift, :]
        # del u, v, w, u_roll, v_roll, w_roll

        fullpath_in = os.path.join(path_fields, str(t0) + '.nc')
        rootgrp = nc.Dataset(fullpath_in, 'r')
        grp_fields = rootgrp.groups['fields']
        nx = grp_fields.dimensions['nx'].size
        ny = grp_fields.dimensions['ny'].size
        nz = grp_fields.dimensions['nz'].size

        global gwz
        gwz = 1
        u = np.zeros((nx,ny,kmax+2*gwz), dtype=np.double)
        v = np.zeros((nx,ny,kmax+2*gwz), dtype=np.double)
        w = np.zeros((nx,ny,kmax+2*gwz), dtype=np.double)
        u[:,:,gwz:-gwz] = grp_fields.variables['u'][:,:,:kmax]
        v[:,:,gwz:-gwz] = grp_fields.variables['v'][:,:,:kmax]
        w[:,:,gwz:-gwz] = grp_fields.variables['w'][:,:,:kmax]
        rootgrp.close()
        # u[:,:,0] = -u[:,:,1]
        # v[:,:,0] = -v[:,:,1]
        w[:,:,0] = -w[:,:,1]

        u_ = np.roll(np.roll(u[:, :, :], ishift, axis=0), jshift,
                     axis=1)[ic - nx_half + ishift:ic + nx_half + ishift, jc - ny_half + jshift:jc + ny_half + jshift, :]
        v_ = np.roll(np.roll(v[:, :, :], ishift, axis=0), jshift,
                         axis=1)[ic - nx_half + ishift:ic + nx_half + ishift, jc - ny_half + jshift:jc + ny_half + jshift, :]
        w_ = np.roll(np.roll(w[:, :, :], ishift, axis=0), jshift,
                     axis=1)[ic - nx_half + ishift:ic + nx_half + ishift, jc - ny_half + jshift:jc + ny_half + jshift, :]



        ''' compute vorticity fields '''
        # compute and plot vorticity in yz-cross section
        vort_yz = compute_vorticity_yz(v[ic,:,:], w[ic,:,:], kmax)                      # interpolating velocity fields
        # compute and plot vorticity in xz-crosssection
        vort_xz = compute_vorticity_xz(u[:,jc,:], w[:,jc,:], kmax)


        ''' dump vorticity yz-field '''
        dump_vort_field(vort_yz, it, t0, fields_file_name)


        ''' PLOTTING '''
        # plot_configuration(u_, v_, icshift, jcshift, t0)
        # compare vort_yz and vort_xz
        comparison_vort_yz_vort_xz(vort_xz, vort_yz, kmax, t0)
        # plot fields
        plot_vorticity_field(vort_xz, vort_yz, t0, kmax)


        ''' Compare interpolation to staggered grid '''
        # vorticites without interpolation of velocity fields
        vort_yz_stag = compute_vorticity_yz_staggered(v[ic, :, :], w[ic, :, :],
                                                      kmax)  # NOT interpolating velocity fields
        vort_xz_stag = compute_vorticity_xz_staggered(u[:, jc, :], w[:, jc, :], kmax)
        # compare vorticities
        plot_comparison_vort_vort_stag(vort_yz, vort_yz_stag, jcshift, jc, kmax, t0)



    ''' compute vorticity statistics '''
    # create file for output time arrays
    stats_file_name = 'Stats_vorticity.nc'
    create_statistics_file(stats_file_name, path_out_fields, nt, timerange)
    add_statistics_variable('vort_yz_max', 's^-1', 'timeseries', stats_file_name, path_out_fields)
    add_statistics_variable('vort_yz_min', 's^-1', 'timeseries', stats_file_name, path_out_fields)
    add_statistics_variable('vort_yz_sum', 's^-1', 'timeseries', stats_file_name, path_out_fields)
    add_statistics_variable('vort_yz_env', 's^-1', 'timeseries', stats_file_name, path_out_fields)

    # read in vorticity field
    fields_file_name = 'field_vort_yz.nc'
    file = nc.Dataset(os.path.join(path_out_fields, fields_file_name), 'r')
    vort_yz = file.groups['fields'].variables['vort_yz'][:,:,:]
    file.close()
    # compute vorticity statistics
    vort_yz_max = np.amax(np.amax(vort_yz, axis=1), axis=1)
    vort_yz_min = np.amin(np.amin(vort_yz, axis=1), axis=1)
    vort_yz_sum = np.sum(np.sum(vort_yz[:, jc_arr[0]:, :], axis=1), axis=1)
    vort_yz_env = np.sum(vort_yz[:, nx-1, :], axis=1)  # profile outside of coldpool

    # output vorticity timeseries
    dump_statistics_variable(vort_yz_max, 'vort_yz_max', 'timeseries', stats_file_name, path_out_fields)
    dump_statistics_variable(vort_yz_min, 'vort_yz_min', 'timeseries', stats_file_name, path_out_fields)
    dump_statistics_variable(vort_yz_sum, 'vort_yz_sum', 'timeseries', stats_file_name, path_out_fields)
    dump_statistics_variable(vort_yz_env, 'vort_yz_env', 'timeseries', stats_file_name, path_out_fields)



    ''' PLOTTING '''
    figname = 'vort_yz_max_sum.png'
    plot_vorticity_timeseries(vort_yz_max, vort_yz_min, vort_yz_sum, vort_yz_env, timerange, figname)

    return





# ---------------------------------- PLOTTING ------------------------------------
def plot_vorticity_timeseries(vort_yz_max, vort_yz_min, vort_yz_sum, vort_yz_env, timerange, figname):
    fig, axes = plt.subplots(1,3, figsize=(15,5))
    ax1 = axes[0]
    ax1.plot(timerange, vort_yz_max, 'o-', label = 'max')
    ax1.plot(timerange, -vort_yz_min, 'o-', label='-min')
    ax1.legend()
    ax1.set_title('max vort_yz')
    ax1.set_xlabel('time t  [s]')
    ax1.set_ylabel('max vort_yz  [1/s]')
    ax2 = axes[1]
    ax2.plot(timerange, vort_yz_sum, '-o')
    ax2.set_title('sum vort_yz')
    ax2.set_xlabel('time t  [s]')
    ax2.set_ylabel('sum vort_yz  [1/s]')
    fig.savefig(os.path.join(path_out_figs, 'vort_yz_max_sum.png'))
    ax3 = axes[2]
    ax3.plot(timerange, vort_yz_env, '-o')
    ax3.set_title('sum vort_yz env')
    ax3.set_xlabel('time t  [s]')
    ax3.set_ylabel('sum vort_yz  [1/s]')
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, wspace=0.25)
    # fig.savefig(os.path.join(path_out_figs, 'vort_yz_max_sum.png'))
    fig.savefig(os.path.join(path_out_figs, figname))

    return


def plot_vorticity_field(vort_xz, vort_yz, t0, kmax):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,5))
    cf = ax1.imshow(vort_xz.T, origin='lower')
    plt.colorbar(cf, ax=ax1, shrink=0.65, aspect=12)
    ax1.set_title('vort_xz')
    cf = ax2.imshow(vort_yz.T, origin='lower')
    plt.colorbar(cf, ax=ax2, shrink=0.65, aspect=12)
    ax2.set_title('vort_yz')
    ax1.plot([jc_arr[0],jc_arr[0]],[1,ny_-2],'k')
    ax2.plot([ic_arr[0],ic_arr[0]],[1,ny_-2],'k')
    pltlim = 150
    ax1.set_xlim([ic_arr[0]-pltlim, ic_arr[0]+pltlim])
    ax1.set_ylim([0, kmax])
    ax2.set_xlim([jc_arr[0]-pltlim, jc_arr[0]+pltlim])
    ax2.set_ylim([0, kmax])
    ax1.set_ylabel('z  (dx='+str(dx[2]) + ')')
    ax2.set_ylabel('z  (dx='+str(dx[2]) + ')')
    ax1.set_xlabel('x  (dx='+str(dx[0]) + ')')
    ax2.set_xlabel('y  (dx='+str(dx[1]) + ')')
    fig.savefig(os.path.join(path_out_figs, 'vort_t' + str(t0) + 's.png'))
    plt.close(fig)

    return


def plot_configuration(u_, v_, icshift, jcshift, t0):
    ''' FIELDS / GEOMETRY '''
    plt.figure()
    plt.subplot(131)
    plt.imshow(u_[:,:,1].T, origin='lower')
    plt.plot([icshift,icshift],[1,ny_-2],'k')
    plt.plot([1,nx_-2],[jcshift,jcshift],'k')
    plt.xlim([0,nx_])
    plt.ylim([0,ny_])
    plt.subplot(132)
    plt.imshow(v_[:,:,1].T, origin='lower')
    plt.plot([icshift, icshift], [1, ny_ - 2], 'k')
    plt.plot([1, nx_ - 2], [jcshift, jcshift], 'k')
    plt.xlim([0, nx_])
    plt.ylim([0, ny_])
    plt.subplot(133)
    plt.imshow(w_[:,:,1].T, origin='lower')
    plt.plot([icshift, icshift], [1, ny_ - 2], 'k')
    plt.plot([1, nx_ - 2], [jcshift, jcshift], 'k')
    plt.xlim([0, nx_])
    plt.ylim([0, ny_])
    plt.savefig(os.path.join(path_out_figs, 'test_fields_t' + str(t0) + 's.png'))
    return



# ----------------------------------------------------------------------
def compute_vorticity_yz(v_, w_, kmax):
    # compute vorticity on staggered grid, i.e. on y- and k-faces of boxes
    # >> return vort_yz_stag[i, j_half, k_half], (where j_half, k_half are the faces of the boxes, i.e. the location of v and w)
    [ly, lz] = v_.shape
    vort_yz = np.zeros((ly, kmax+2*gwz), dtype=np.double)
    dyi = 1. / dx[1]
    dzi = 1. / dx[2]
    for j in range(ly-1):
        for k in range(kmax+2*gwz-1):
            vort_yz[j, k] = dyi * (w_[j+1, k] - w_[j, k]) - dzi * (v_[j, k+1] - v_[j, k])
    return vort_yz[:, gwz:kmax+gwz]


def compute_vorticity_xz(u_, w_, kmax):
    # compute vorticity on staggered grid, i.e. on x- and k-faces of boxes
    [lx, lz] = u_.shape
    vort_xz = np.zeros((lx, kmax+2*gwz), dtype=np.double)
    dxi = 1. / dx[0]
    dzi = 1. / dx[2]
    if lz < kmax+gwz:
        kmax = lz-gwz
    for i in range(1, lx - 1):
        for k in range(gwz, kmax+gwz):
            vort_xz[i, k] = dzi * (u_[i, k + 1] - u_[i, k]) - dxi * (w_[i+1, k] - w_[i, k])

    return vort_xz[:,gwz:kmax+gwz]



def compute_vorticity_yz_staggered(v_, w_, kmax):
    [ly, lz] = v_.shape
    if lz < kmax+gwz:
        kmax = lz-gwz
    vort_yz = np.zeros((ly, kmax+2*gwz), dtype=np.double)
    dyi2 = 1./(2*dx[1])
    dzi2 = 1./(2*dx[2])
    for j in range(1, ly - 1):
        # for k in range(1, lz-1):
        for k in range(gwz, kmax+gwz):
            vort_yz[j, k] = (w_[j + 1, k] - w_[j - 1, k]) * dyi2 \
                                - (v_[j, k + 1] - v_[j, k - 1]) * dzi2
    return vort_yz[:,gwz:kmax+gwz]


def compute_vorticity_xz_staggered(u_, w_, kmax):
    [lx, lz] = u_.shape
    if lz < kmax+gwz:
        kmax = lz-gwz
    # vort_xz = np.zeros(shape=u_.shape, dtype=np.double)
    vort_xz = np.zeros((lx, kmax+2*gwz), dtype=np.double)
    for j in range(1, lx - 1):
        # for k in range(1, lz-1):
        for k in range(gwz, kmax+gwz):
            vort_xz[j, k] = (u_[j, k + 1] - u_[j, k - 1]) / (2 * dx[2]) \
                            - (w_[j + 1, k] - w_[j - 1, k]) / (2 * dx[1])
    return vort_xz[:,gwz:kmax+gwz]



# ----------------------------------------------------------------------
# PLOTTING

def plot_comparison_vort_vort_stag(vort_yz, vort_yz_stag, jcshift, jc, kmax, t0):
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))

    ax = axes[0, 0]
    cf = ax.imshow(vort_yz.T, origin='lower')
    ax.set_title('vort_yz')
    ax.plot([jc, jc], [1, ny_ - 2], 'k')
    plt.colorbar(cf, ax=ax, shrink=0.6)

    ax = axes[0, 1]
    cf = ax.imshow(vort_yz_stag.T, origin='lower')
    ax.set_title('vort_yz_ staggered')
    ax.plot([jc, jc], [1, ny_ - 2], 'k')
    plt.colorbar(cf, ax=ax, shrink=0.6)

    for i in range(2):
        axes[0, i].set_xlim([jc - np.int(ny_ / 2), jc + np.int(ny_ / 2)])
        axes[0, i].set_ylim([0, kmax])

    for i in range(3):
        ax = axes[1, i]
        ax.plot(np.arange(0, ny_), vort_yz[jc - np.int(ny_ / 2):jc + np.int(ny_ / 2), i], label='vort_yz')
        ax.plot(np.arange(0, ny_), vort_yz_stag[jc - np.int(ny_ / 2):jc + np.int(ny_ / 2), i], label='vort_yz stag')
        ax.legend(loc='best', fontsize=8)
        ax.set_title('k='+str(i))
        ax.set_xlabel('y')
        ax.set_xlim(0,ny_)

    for i in range(3):
        ax = axes[2, i]
        ax.plot(np.arange(0, ny_), vort_yz[jc - np.int(ny_ / 2):jc + np.int(ny_ / 2), i] - vort_yz_stag[:ny_, i],
                label='diff')
        ax.set_title('vort_yz - vort_yz_ (k='+str(i)+')')
        ax.set_xlabel('y')
        ax.set_xlim(0,ny_)

    plt.suptitle('comparison vort_yz unstaggered vs. staggered')
    fig.savefig(os.path.join(path_out_figs, 'test_vort_t' + str(t0) + 's.png'))
    plt.close(fig)
    return


def comparison_vort_yz_vort_xz(vort_xz_, vort_yz_, kmax, t0):
    fig, axes = plt.subplots(3,1, figsize=(14,4), sharey='all')
    y_arr = np.arange(0,ny_)
    z_arr = np.arange(0,kmax)
    Y_, Z_ = np.meshgrid(y_arr, z_arr)
    plt.suptitle('t='+str(t0)+'s')
    ax1 = axes[0]
    a = ax1.imshow(vort_xz_[:,:].T, origin='lower')
    ax1.autoscale(False)
    ax1.set_title('vort_xz')
    ax2 = axes[1]
    b = ax2.imshow(vort_yz_[:,:].T, origin='lower')
    ax2.set_title('vort_yz')
    ax3 = axes[2]
    c = ax3.imshow(vort_xz_[:,:].T+vort_yz_[:,:].T, origin='lower')
    ax3.set_title('vort_xz - vort_yz')
    plt.colorbar(a, ax=ax1)
    plt.colorbar(b, ax=ax2)
    plt.colorbar(c, ax=ax3)
    plt.savefig(os.path.join(path_out_figs, 'vort_xz_yz_t'+str(t0)+'s.png'))
    plt.close()
    return

# ----------------------------------------------------------------------

def set_input_output_parameters(args):
    print('--- set input parameters ---')
    global path, path_in, path_out_figs, path_out_fields, path_stats, path_fields
    path = args.path
    path_in = os.path.join(path, 'fields_CP_rim')
    path_fields = os.path.join(path, 'fields')
    path_out_figs = os.path.join(path, 'figs_vorticity')
    path_out_fields = os.path.join(path, 'fields_vorticity')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    if not os.path.exists(path_out_fields):
        os.mkdir(path_out_fields)

    global case_name
    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dy, dz, dV, gw
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
    global nt
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

    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = 60
    print('nx, ny, nz', nx, ny, nz)
    print('times', timerange)


    return timerange, kmax, nml



def define_geometry(case_name, nml):
    print('--- define geometry ---')
    global x_half, y_half, z_half
    global nx_, ny_
    global shift, ishift, jshift
    global ic_arr, jc_arr
    global rstar, irstar, zstar, kstar

    x_half = np.empty((nx), dtype=np.double, order='c')
    y_half = np.empty((ny), dtype=np.double, order='c')
    z_half = np.empty((nz), dtype=np.double, order='c')
    count = 0
    for i in xrange(nx):
        x_half[count] = (i + 0.5) * dx[0]
        count += 1
    count = 0
    for j in xrange(ny):
        y_half[count] = (j + 0.5) * dx[1]
        count += 1
    count = 0
    for i in xrange(nz):
        z_half[count] = (i + 0.5) * dx[2]
        count += 1

    # set coordinates for plots
    if case_name == 'ColdPoolDry_single_3D':
        n_CPs = 1
        rstar = nml['init']['r']
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        kstar = np.int(np.round(zstar / dx[2]))
        try:
            ic = nml['init']['ic']
            jc = nml['init']['jc']
            # print('(ic,jc) from nml')
        except:
            ic = np.int(nx/2)
            jc = np.int(ny/2)
            # print('(ic,jc) NOT from nml')
        ic_arr = np.zeros(1)
        jc_arr = np.zeros(1)
        ic_arr[0] = ic
        jc_arr[0] = jc
    elif case_name == 'ColdPoolDry_double_2D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        kstar = np.int(np.round(zstar / dx[2]))
        isep = 4 * irstar
        ic1 = np.int(nx / 3)  # np.int(Gr.dims.ng[0] / 3)
        ic2 = ic1 + isep
        jc1 = np.int(ny / 2)
        jc2 = jc1
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_double_3D':
        n_CPs = 2
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        kstar = np.int(np.round(zstar / dx[2]))
        isep = 4 * irstar
        jsep = 0
        ic1 = np.int(nx / 3)
        jc1 = np.int(ny / 2)
        # ic2 = ic1 + isep
        # jc2 = jc1 + jsep
        ic = ic1
        jc = jc1

        # ic_arr = [ic1, ic2]
        # jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_triple_3D':
        n_CPs = 2
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        zstar = nml['init']['h']
        kstar = np.int(np.round(zstar / dx[2]))
        marg_i = 10  # width of margin
        # d = np.int(np.round(ny / 2))
        d = np.int(np.round((ny + gw) / 2))
        # d = np.int(np.round(10 * irstar)) # for r=1km, dTh=2K
        a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        ic1 = np.int(np.round(a / 2))
        # ic1 = 10 + np.int(np.round(a / 2)) + Gr.dims.gw # for r=1km
        ic2 = ic1
        ic3 = ic1 + np.int(np.round(a))
        jc1 = np.int(np.round(d / 2))
        # jc1 = np.int(np.round(d / 2) + gw)  # np.int(np.round(d/2) + Gr.dims.gw)
        jc2 = jc1 + d
        jc2 = jc1 + d
        jc3 = jc1 + np.int(np.round(d / 2))
        ic = ic1
        jc = jc1
        ic_arr = [ic1, ic2, ic3]
        jc_arr = [jc1, jc2, jc3]

    shift = 120
    nx_half = irstar + shift
    ny_half = irstar + shift
    ishift = np.max(nx_half - ic, 0)
    jshift = np.max(ny_half - jc, 0)
    nx_ = 2 * nx_half
    ny_ = 2 * ny_half


    print('rstar:     ' + str(rstar), irstar)
    print('ic,jc,id,jd', ic, jc, nx_half, ny_half)
    print('nx_,ny_    ', nx_, ny_)
    print('shift, ishift, jshift', shift, ishift, jshift)

    return n_CPs

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

def create_statistics_file(file_name, path, nt, timerange):#, nk, krange):
    print('-------- create statistics file -------- ')
    print(path + ', ' + file_name)
    print('')
    rootgrp = nc.Dataset(os.path.join(path, file_name), 'w', format='NETCDF4')
    ts_grp = rootgrp.createGroup('timeseries')
    ts_grp.createDimension('nt', nt)
    var = ts_grp.createVariable('time', 'f8', ('nt'))
    var.units = "s"
    var[:] = timerange[:]

    rootgrp.close()

    return

def add_statistics_variable(var_name, units, grp_name, file_name, path):
    rootgrp = nc.Dataset(os.path.join(path, file_name), 'r+', format='NETCDF4')
    try:
        grp = rootgrp.groups[grp_name]
    except:
        print 'except', grp_name
        grp = rootgrp.createGroup(grp_name)
    if grp_name == 'timeseries':
        var = grp.createVariable(var_name, 'f8', ('nt'))
        var.units = units
    rootgrp.close()
    return


def dump_statistics_variable(var_in, var_name, grp_name, file_name, path):
    rootgrp = nc.Dataset(os.path.join(path, file_name), 'r+', format='NETCDF4')
    grp = rootgrp.groups[grp_name]
    var = rootgrp.groups['timeseries'].variables['vort_yz_max']
    if grp_name == 'timeseries':
        var = grp.variables[var_name]
        var[:] = var_in[:]
    rootgrp.close()
    return

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

def create_vort_field_file(times, file_name, kmax, ncp):
    # create fields file
    # test field
    files = [name for name in os.listdir(path_fields) if name[-2:] == 'nc']
    fullpath_in = os.path.join(path_fields, files[0])
    rootgrp = nc.Dataset(fullpath_in, 'r')
    grp_fields = rootgrp.groups['fields']
    nx = grp_fields.dimensions['nx'].size
    ny = grp_fields.dimensions['ny'].size
    nz = grp_fields.dimensions['nz'].size
    rootgrp.close()

    rootgrp_fields = nc.Dataset(os.path.join(path_out_fields, file_name), 'w', format='NETCDF4')

    descr_grp = rootgrp_fields.createGroup('description')
    descr_grp.createDimension('ncp', ncp)
    var = descr_grp.createVariable('ic_arr', 'f8', 'ncp')
    var[:] = ic_arr
    var = descr_grp.createVariable('jc_arr', 'f8', 'ncp')
    var[:] = jc_arr
    var = descr_grp.createVariable('ishift', 'f8', )
    var[:] = ishift
    var = descr_grp.createVariable('jshift', 'f8', )
    var[:] = jshift

    fields_grp = rootgrp_fields.createGroup('fields')
    fields_grp.createDimension('time', None)
    fields_grp.createDimension('nx', nx)
    fields_grp.createDimension('ny', ny)
    # fields_grp.createDimension('nz', nz)
    fields_grp.createDimension('nz', kmax)
    fields_grp.createDimension('nz_ori', nz)

    time_out = fields_grp.createVariable('time', 'f8', ('time',))
    time_out.units = 's'
    # time_out[:] = times
    var = fields_grp.createVariable('vort_yz', 'f8', ('time', 'ny', 'nz'))
    var[:,:] = np.ones((len(times), ny, kmax))
    rootgrp_fields.close()
    return



def dump_vort_field(vort_yz, it, t0, file_name):
    rootgrp = nc.Dataset(os.path.join(path_out_fields, file_name), 'r+', format='NETCDF4')
    var = rootgrp.groups['fields'].variables['vort_yz']
    var[it, :,:] = vort_yz[:,:]
    time = rootgrp.groups['fields'].variables['time']
    time[it] =  t0
    rootgrp.close()
    return




# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:, :, :]
    rootgrp.close()
    return data

if __name__ == '__main__':
    main()