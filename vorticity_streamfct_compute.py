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
COMPUTE VORTICITY & STREAMFUNCTION
compute_vorticity_yz: computes vorticity on yz-crosssection on Euclidian grid 
    > location: on grid faces, i.e. y- and z-faces of grid boxes
compute_vorticity_yz_staggered: computes vorticity on yz-crosssection on Euclidian grid
    (!! without interpolating v, w onto centered grid, i.e. vorticity location not well defined)
 
'''

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


    ''' TEST FUNCTIONS '''
    # streamfunction_circular_flow()

    # create file for output time arrays
    stats_file_name = 'Stats_rotational.nc'
    create_statistics_file(stats_file_name, path_out_fields, nt, timerange)
    add_statistics_variable('vort_yz_max', 's^-1', 'timeseries', stats_file_name, path_out_fields)
    add_statistics_variable('vort_yz_min', 's^-1', 'timeseries', stats_file_name, path_out_fields)
    add_statistics_variable('vort_yz_sum', 's^-1', 'timeseries', stats_file_name, path_out_fields)
    add_statistics_variable('vort_yz_env', 's^-1', 'timeseries', stats_file_name, path_out_fields)

    # # create filds file
    fields_file_name = 'field_vort_yz.nc'
    create_vort_field_file(timerange, fields_file_name, kmax, n_CPs)

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
    vort_yz_max = np.zeros((len(timerange)))
    vort_yz_min = np.zeros((len(timerange)))
    vort_yz_sum = np.zeros((len(timerange)))
    vort_yz_env = np.zeros((len(timerange)))
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

        # del u, v, w


        # plot_configuration(u_, v_, icshift, jcshift)

        ''' VORTICITY '''
        # compute and plot vorticity in yz-cross section
        vort_yz = compute_vorticity_yz(v[ic,:,:], w[ic,:,:], kmax)
        vort_yz_ = compute_vorticity_yz(v_[icshift,:,:], w_[icshift,:,:], kmax)
        vort_yz_stag = compute_vorticity_yz_staggered(v[ic,:,:], w[ic,:,:], kmax)
        vort_yz_stag_ = compute_vorticity_yz_staggered(v_[icshift,:,:], w_[icshift,:,:], kmax)

        # compute and plot vorticity in xz-crosssection
        vort_xz_stag = compute_vorticity_xz(u[:,jc,:], w[:,jc,:], kmax)
        vort_xz_stag_ = compute_vorticity_xz(u_[:,jcshift,:], w_[:,jcshift,:], kmax)
        # print('vorticity', vort_yz.shape, vort_yz_.shape, vort_xz.shape, nx, ny, nz)


        # compare vorticities
        plot_comparison_vort_vort_stag(vort_yz, vort_yz_stag, vort_yz_stag_, jcshift, jc, kmax, t0)


        # compute vorticity statistics
        vort_yz_max[it] = np.amax(vort_yz)
        vort_yz_min[it] = np.amin(vort_yz)
        vort_yz_sum[it] = np.sum(vort_yz[jcshift:,:])
        vort_yz_env[it] = np.sum(vort_yz[jcshift+50:,:])       # profile outside of coldpool

        # dump vorticity_yz field
        dump_vort_field(vort_yz_stag, it, t0, fields_file_name)

        # compare vort_yz and vort_xz
        comparison_vort_yz_vort_xz(vort_xz_stag_, vort_yz_stag_, kmax, t0)

        # PLOTTING
        plot_vorticity_field(vort_xz_stag_, vort_yz, icshift, jcshift, t0)

    plot_vorticity_timeseries(vort_yz_max, vort_yz_min, vort_yz_sum, vort_yz_env, timerange)

    # output vorticity timeseries
    dump_statistics_variable(vort_yz_max, 'vort_yz_max', 'timeseries', stats_file_name, path_out_fields)
    dump_statistics_variable(vort_yz_min, 'vort_yz_min', 'timeseries', stats_file_name, path_out_fields)
    dump_statistics_variable(vort_yz_sum, 'vort_yz_sum', 'timeseries', stats_file_name, path_out_fields)
    dump_statistics_variable(vort_yz_env, 'vort_yz_env', 'timeseries', stats_file_name, path_out_fields)


    ''' STREAM FUNCTION '''
    # for it, t0 in enumerate(timerange):
    #     print('--- time: t=' + str(t0) + 's ---')
    #
    #     # compute 2D streamfunction in yz-crosssection
    #     psi_1 = compute_streamfunction_simps(v[ic,:,:], w[ic,:,:], t0)
    #     psi_2 = compute_streamfunction_simps2(v[ic,:,:], w[ic,:,:], t0)
    #     psi_3 = compute_streamfunction_simps3(v[ic,:,:], w[ic,:,:], t0)
    #     psi_4 = compute_streamfunction_simps4(v[ic,:,:], w[ic,:,:], rho0, t0)
    #     psi_5 = compute_streamfunction_simps5(v[ic,:,:], w[ic,:,:], rho0, 0, 0)
    #
    #     # compare two functions
    #     ymax = 120
    #     y = np.arange(0,ymax)
    #     z = np.arange(kmax)
    #     Y, Z = np.meshgrid(y, z,indexing='ij')
    #     print 'psi', psi_1.shape, psi_2.shape, y.shape, z.shape, Y.shape, Z.shape
    #     fig, axes = plt.subplots(2,3, figsize=(16,8))
    #     ax1 = axes[0,0]
    #     ax1.set_title('psi_1')
    #     max = np.maximum(np.amax(psi_1[:ymax, :kmax]),-np.amin(psi_1[:ymax, :kmax]))
    #     a = ax1.contourf(Y,Z,psi_1[:ymax, :kmax], levels=np.linspace(-max,max,1e2))
    #     plt.colorbar(a, ax=ax1)
    #     ax2 = axes[0,1]
    #     ax2.set_title('psi_2')
    #     max = np.maximum(np.amax(psi_2[:ymax, :kmax]), -np.amin(psi_2[:ymax, :kmax]))
    #     max = np.amax(psi_2[:ymax, :kmax])
    #     a = ax2.contourf(Y,Z,psi_2[:ymax, :kmax], levels=np.linspace(-max,max,1e2), extend = 'both')
    #     a.cmap.set_under('yellow')
    #     a.cmap.set_over('cyan')
    #     plt.colorbar(a, ax=ax2)
    #     ax3 = axes[0,2]
    #     ax3.set_title('psi_4')
    #     max = np.maximum(np.amax(psi_4[:ymax, :kmax]), -np.amin(psi_4[:ymax, :kmax]))
    #     a = ax3.contourf(Y,Z,psi_4[:ymax, :kmax], levels=np.linspace(-max,max,1e2))
    #     plt.colorbar(a, ax=ax3)
    #
    #     ax = axes[1,1]
    #     ax.set_title('w')
    #     var = w[ic,:ymax,:kmax]
    #     max = np.maximum(np.amax(var), -np.amin(var))
    #     a = ax.contourf(Y,Z,var[:ymax, :kmax], levels=np.linspace(-max,max,1e2), cmap=cm_bwr)
    #     speed_yz = np.sqrt(v[ic,:ymax,:kmax]*v[ic,:ymax,:kmax] + var*var)
    #     lw = 5 * speed_yz[:, :] / speed_yz[:, :].max()
    #     ax.streamplot(y,z,v[ic,:ymax,:kmax].T,var.T, color='k', linewidth=lw[:,:].T)
    #     ax3.streamplot(y,z,v[ic,:ymax,:kmax].T,var.T, color='k', linewidth=lw[:,:].T)
    #     plt.colorbar(a, ax=ax)
    #     ax = axes[1,2]
    #     ax.set_title('w')
    #     max = np.maximum(np.amax(psi_4[:ymax, :kmax]), -np.amin(psi_4[:ymax, :kmax]))
    #     a = ax.contourf(Y, Z, psi_4[:ymax, :kmax], levels=np.linspace(-max, max, 1e2))
    #     plt.colorbar(a, ax=ax)
    #     var = w[ic,:,:]
    #     max = np.maximum(np.amax(var[:ymax, :kmax]), -np.amin(var[:ymax, :kmax]))
    #     a = ax.contour(Y,Z,var[:ymax, :kmax], levels=np.linspace(-max,max,1e1+1), cmap=cm_bwr)
    #     # plt.colorbar(a, ax=ax)
    #
    #     ax = axes[1,0]
    #     ax.set_title('s')
    #     var = s[ic,:,:]
    #     max = np.amax(var[:ymax, :kmax])
    #     min = np.amin(var[:ymax, :kmax])
    #     a = ax.contourf(Y,Z,var[:ymax, :kmax], levels=np.linspace(min, max, 1e2))
    #     plt.colorbar(a, ax=ax)
    #
    #     fig.savefig(os.path.join(path_out_figs, 'psi_t'+str(t0)+'s.png'))
    #     plt.close(fig)
    #
    #     fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    #     plt.suptitle('t='+str(t0)+'s')
    #     ax1 = axes[0]
    #     ax1.set_title('v')
    #     max = np.maximum(np.amax(v[ic,20:20 + ymax, :kmax]), -np.amin(v[ic,20:20 + ymax, :kmax]))
    #     a = ax1.contourf(Y, Z, v[ic, 20:20 + ymax, :kmax], levels=np.linspace(-max, max, 1e2))
    #     plt.colorbar(a, ax=ax1)
    #     ax1 = axes[1]
    #     ax1.set_title('w')
    #     max = np.maximum(np.amax(w[ic,20:20 + ymax, :kmax]), -np.amin(w[ic,20:20 + ymax, :kmax]))
    #     a = ax1.contourf(Y, Z, w[ic, 20:20 + ymax, :kmax], levels=np.linspace(-max, max, 1e2))
    #     plt.colorbar(a, ax=ax1)
    #     ax2 = axes[2]
    #     ax2.set_title('vort_yz')
    #     print vort_yz.shape
    #     max = np.maximum(np.amax(vort_yz[20:20 + ymax, :kmax]), -np.amin(vort_yz[20:20 + ymax, :kmax]))
    #     a = ax2.contourf(Y, Z, vort_yz[20:20 + ymax, :kmax], levels=np.linspace(-max, max, 1e2))
    #     plt.colorbar(a, ax=ax2)
    #     ax2 = axes[3]
    #     ax2.set_title('psi_4')
    #     max = np.maximum(np.amax(psi_4[20:20 + ymax, :kmax]), -np.amin(psi_4[20:20 + ymax, :kmax]))
    #     a = ax2.contourf(Y, Z, psi_4[20:20 + ymax, :kmax], levels=np.linspace(-max, max, 1e2))
    #     plt.colorbar(a, ax=ax2)
    #     fig.savefig(os.path.join(path_out_figs, 'psi_2_t' + str(t0) + 's.png'))
    #     plt.close(fig)
    #
    #
    #     ym = np.int(3*ny_/3)
    #     ym = ny
    #     YM_, ZM_ = np.meshgrid(y_arr[:ym],z_arr[:kmax],indexing='ij')
    #     # compare streamfunction and streamlines
    #     fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,4))
    #     # # plt.imshow(psi[:,:40].T, origin='lower')
    #     a = ax1.contourf(YM_,ZM_,psi_1[:ym,:kmax])
    #     plt.colorbar(a, ax=ax1)
    #     # y_arr = np.arange(0,ny)
    #     # z_arr = np.arange(0,ny)
    #     ax1.streamplot(y_arr[:ym], z_arr[:kmax], v[ic,:ym,:kmax].T, w[ic,:ym,:kmax].T, density=1.5, color='k')
    #     # print '.....', y_arr.shape, z_arr.shape, ym, kmax, YM_.shape, ZM_.shape, v[ic,:ym,:kmax].shape, w[ic,:ym,:kmax].shape
    #     # ax1.streamplot(YM_,ZM_,v[ic,:ym,:kmax], w[ic,:ym,:kmax], density=1.5, color='k')
    #     b = ax2.contourf(y_arr[:ym], z_arr[:kmax], vort_yz[:ym,:kmax].T)
    #     ax2.streamplot(y_arr[:ym], z_arr[:kmax], v[ic,:ym,:kmax].T, w[ic,:ym,:kmax].T)
    #     plt.colorbar(b, ax=ax2)
    #     ax1.set_title('stream function (x(i-1)..x(i+1)')
    #     ax2.set_title('vorticity')
    #     plt.suptitle('t='+str(t0)+'s')
    #     plt.savefig(os.path.join(path_out_figs, 'psi_1_t'+str(t0)+'s.png'))
    #     plt.close()
    #
    #     print('')



    return



# ---------------------------------- PLOTTING ------------------------------------
def plot_vorticity_timeseries(vort_yz_max, vort_yz_min, vort_yz_sum, vort_yz_env, timerange):
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
    fig.savefig(os.path.join(path_out_figs, 'vort_yz_max_sum.png'))

    return


def plot_vorticity_field(vort_xz_, vort_yz_, icshift, jcshift, t0):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))
    ax1.imshow(vort_xz_.T, origin='lower')
    ax1.set_title('vort_xz')
    ax2.imshow(vort_yz_.T, origin='lower')
    ax2.set_title('vort_yz')
    ax1.plot([icshift,icshift],[1,ny_-2],'k')
    ax2.plot([icshift,icshift],[1,ny_-2],'k')
    # ax1.plot([1,nx_-2],[jcshift,jcshift],'k')
    # ax2.plot([1,nx_-2],[jcshift,jcshift],'k')
    ax1.set_xlim([0, nx_])
    ax1.set_ylim([0, ny_])
    ax2.set_xlim([0, nx_])
    ax2.set_ylim([0, ny_])
    fig.savefig(os.path.join(path_out_figs, 'vort_t' + str(t0) + 's.png'))
    plt.close(fig)

    return


def plot_configuration(u_, v_, icshift, jcshift):
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


def compute_vorticity_xz(u_, w_, kmax):
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

def streamfunction_circular_flow():
    print '--- compute circular stream ---'

    x = np.arange(-10,10,0.1)
    y = np.arange(-10,10,0.1)


    X, Y = np.meshgrid(x, y, indexing='ij')  # if indexing not defined, output will be 'xy'-indexed

    U = Y / (X**2 + Y**2)
    V = - X / (X**2 + Y**2)

    umax = np.amax(U)
    vmax = np.amax(V)
    fig, axes = plt.subplots(2,4, figsize=(20,10))
    ax1 = axes[0,0]
    a = ax1.contourf(x,y,U.T, levels=np.linspace(-umax,umax), cmap=cm_bwr)
    ax1.set_title('contourf(x,y,U.T)')
    ax = axes[1,0]
    ax.contourf(X,Y,U, levels=np.linspace(-umax,umax), cmap=cm_bwr)
    ax.set_title('contoruf(X,Y,U)')
    plt.colorbar(a, ax=ax1)
    plt.colorbar(a, ax = ax)
    ax2 = axes[0,1]
    b = ax2.contourf(x,y,V.T, levels=np.linspace(-10, 10), cmap=cm_bwr)
    ax2.set_title('contourf(x,y,V.T)')
    ax = axes[1,1]
    ax.contourf(X,Y,V, levels=np.linspace(-10,10), cmap=cm_bwr)
    ax.set_title('contourf(X,Y,V)')
    plt.colorbar(b, ax = ax2)
    plt.colorbar(b, ax = ax)

    ax3 = axes[0,2]
    x_ = np.arange(-10,10,0.8)
    y_ = np.arange(-10,10,0.8)
    X_, Y_ = np.meshgrid(x_, y_)
    U_ = Y_ / (X_ ** 2 + Y_ ** 2)
    V_ = - X_ / (X_ ** 2 + Y_ ** 2)
    # q = ax3.quiver(x_, y_, X_, Y_)
    # ax3.quiverkey(q, X=0.3, Y=1.1, U=10,
    #              label='Quiver key, length = 10', labelpos='E')
    q = ax3.quiver(x_, y_, U_, V_)
    ax3.quiverkey(q, X=0.3, Y=1.1, U=3,
                label='Quiver key, length = 3', labelpos='E')
    ax3.set_title('quiver(x,y,U,V)')
    ax = axes[1,2]
    q = ax.quiver(X_, Y_, U_, V_)
    ax.set_title('quiver(X,Y,U,V)')
    # q = ax3.quiver(x, y, U.T, V.T)
    ax.quiverkey(q, X=0.3, Y=1.1, U=3,
                 label='Quiver key, length = 3', labelpos='E')

    # psi = np.log(X**2 + Y**2)
    # a = ax4.contourf(x, y, psi)
    psi_ = np.log(X_**2 + Y_**2)
    ax4 = axes[0,3]
    a = ax4.contourf(x_, y_, psi_)
    plt.colorbar(a, ax=ax4)
    ax4.streamplot(x_, y_, U_, V_)
    # CS = ax4.contour(x, y, psi)
    # ax4.clabel(CS)
    ax4.set_title('streamplot(x,y,U,V)')
    ax = axes[1,3]
    ax.contourf(x_, y_, psi_)
    plt.colorbar(a, ax=ax)
    ax.streamplot(X_,Y_,U_,V_)
    ax.set_title('streamplot(X,Y,U,V)')
    fig.savefig('./circ_flow.png')
    plt.close()

    return


def compute_streamfunction_simps5(v_, w_, rho0, i0, j0):
    # int dy w(y,0)
    yrange = np.arange(w_.shape[0])
    # psi_A_ = integrate.simps(w_[:, 0], yrange[:])
    rhow = np.zeros(shape=w_.shape)
    for k in range(nz):
        rhow[:,k] = rho0[k]*w_[:,k]
    psi_A = integrate.simps(w_[:,:], yrange[:], axis=0)
    print 'shapes: ', psi_A.shape
    # print 'diff: ', (psi_A[0]-psi_A_)
    krange = np.arange(v_.shape[1])
    psi_B = integrate.simps(v_[0,:], krange[:], axis=0)

    return psi_A


def compute_streamfunction_simps4(v_, w_, rho0, t0):
    # Psi(y,z) = \int_y0^y w(y',z)\rho_0(z) dy' + a(z)
    # a(z) = -\int_z0^z u(y0,z) dz'
    if v_.shape != w_.shape:
        print('ERROR: input matrices do not have same shape')
        sys.exit()
    [ly,lz] = v_.shape
    psi_A = np.zeros((ly,lz), dtype=np.double)     # psi_A = int dx v(x,y)
    psi_B = np.zeros((ly,lz), dtype=np.double)     # psi_B = int dy u(x,y)
    # psi_simps = np.zeros((ly,lz))                  # psi = int(v dx) - int(u dy) = psi_A - psi_B
    for j in range(ly):
        for k in range(lz):
            # Simpson Rule
            j0 = np.int(np.round(j/2))
            k0 = np.int(np.round(k/2))
            # (1)  \rho_0(z)*\int_y0^y w(y',z) dy'
            psi_A[j,k] = rho0[k] * (j * dy) / 6 * (w_[0, k] + 4 * w_[j0, k] + w_[j, k])
            # (2) a(z) =  -\int_z0^z u(y0,z) dz'
            psi_B[j,k] = -  j * dy / 6 * (rho0[0]*v_[0, 0] + 4 * rho0[k0]*v_[0, k0] + rho0[k]*v_[0, k])
    psi_simps = psi_A + psi_B

    return psi_simps

def compute_streamfunction_simps3(u_, v_, t0):
    # Psi(x,y) = \int_y0^y w(x,y') dy' + a(x)
    # a(x) = -\int_x0^x v(x',y0) dx'
    if u_.shape != v_.shape:
        print('ERROR: input matrices do not have same shape')
        sys.exit()
    [lx, ly] = u_.shape
    psi_A = np.zeros((lx, ly), dtype=np.double)     # psi_A = int dx v(x,y)
    psi_B = np.zeros((lx, ly), dtype=np.double)     # psi_B = int dy u(x,y)
    psi_simps = np.zeros((lx, ly))                  # psi = int(v dx) - int(u dy) = psi_A - psi_B
    for i in range(lx):
        for j in range(ly):
            # Simpson Rule
            i0 = np.int(np.round(i/2))
            j0 = np.int(np.round(j/2))
            # (1)  \int_y0^y u(x,y') dy'
            psi_B[i, j] = j * dy / 6 * (u_[i, 0] + 4 * u_[i, j0] + u_[i, j])
            # (2) a(x) = -\int_x0^y v(x',y0)
            psi_A[i, j] = - i * dx / 6 * (v_[0, 0] + 4 * v_[i0, 0] + v_[i, 0])
    psi_simps = psi_A + psi_B

    return psi_simps


def compute_streamfunction_simps2(u_, v_, t0):
    # F(x) = int_c^x f(y) dy, with F(c) = 0
    if u_.shape != v_.shape:
        print('ERROR: input matrices do not have same shape')
        sys.exit()
    [lx, ly] = u_.shape
    psi_A = np.zeros((lx, ly), dtype=np.double)     # psi_A = int dx v(x,y)
    psi_B = np.zeros((lx, ly), dtype=np.double)     # psi_B = int dy u(x,y)
    psi_simps = np.zeros((lx, ly))                  # psi = int(v dx) - int(u dy) = psi_A - psi_B
    for i in range(lx):
        for j in range(ly):
            # Simpson Rule
            i0 = np.int(np.round(i/2))
            j0 = np.int(np.round(j/2))
            psi_B[i, j] = j * dy / 6 * (u_[i, 0] + 4 * u_[i, j0] + u_[i, j])
            psi_A[i, j] = i * dx / 6 * (v_[0, j] + 4 * v_[i0, j] + v_[i, j])
    psi_simps = psi_A - psi_B

    return psi_simps


def compute_streamfunction_simps(u_, v_, t0):
    # est, errbound = integrate.quad(func, min, max)
    # funct: callable python object (function, method, class instance); also lambda-function
    # min, max: limits for integration (can us inf, -inf)
    # return:
    #   est: estimated value for integral
    #   errbound: upper bound on the error

    # testing
    # test_num_integration()

    if u_.shape != v_.shape:
        print('ERROR: input matrices do not have same shape')
        sys.exit()

    # phi: vector potential
    # psi: stream function (solenoidal part)
    [lx, ly] = u_.shape
    print(lx, ly)

    # Integrate velocity fields to get potential and streamfunction
    # Use Simpson rule summation(function CUMSIMP)
    # for u = -dpsi/dy, v = dpsi/dx
    # psi = int(v dx) - int(u dy)
    # psi = integrate.quad()

    psi_A = np.zeros((lx, ly),dtype=np.double)
    psi_B = np.zeros((lx, ly),dtype=np.double)
    psi_simps = np.zeros((lx, ly))
    for i in range(1,lx-1):
        for j in range(1,ly-1):
            # Simpson Rule
            psi_A[i,j] = 2*dx/6 * ( v_[i+1,j] + 4*v_[i,j] + v_[i-1,j] )
            psi_B[i,j] = 2*dy/6 * ( u_[i,j+1] + 4*u_[i,j] + u_[i,j-1] )
    psi_simps = psi_A - psi_B

    km = 30
    ym = np.int(ny_/2)
    y_arr = dy*np.arange(0,ny_)
    z_arr = dz*np.arange(0,nz)
    YM_, ZM_ = np.meshgrid(y_arr[:ym], z_arr[:km])
    fig, axes = plt.subplots(1,2, figsize=(10,4))
    ax1 = plt.subplot(121)
    ax1.set_title('u')
    plt.contourf(u_[:ym,:km].T, alpha=0.7)
    plt.colorbar()
    plt.contour(u_[:ym,:km].T, levels=[0], colors='k', linewidth=5)
    CS = ax1.contour(-psi_B[:ym,:km].T, linewidth=3)
    ax1.clabel(CS)
    ax2 = plt.subplot(122)
    ax2.set_title('psi_B = int dy u')
    print('maxs', np.amax(u_[:ym,:km]), np.amax(psi_B[:ym,:km]))
    ax2.contour(-psi_B[:ym,:km].T)
    ax2.contour(u_[:ym,:km].T, linestyles='--')
    fig.savefig(os.path.join(path_out_figs, 'test_u_t'+str(t0)+'s.png'))
    plt.close()

    return psi_simps


def integrand(x, a, b):
    return a*x**2 + b

def test_num_integration():
    x = np.arange(1, 10)
    a = 2
    f = a * x
    I = integrate.simps(f, x)
    print x
    print f
    print I
    x = np.arange(0, 10)
    f = 2 * np.ones(10)
    I = integrate.simps(f, x)
    print f
    print I
    x = np.arange(0, 10)
    y = np.arange(0, 5)
    X, Y = np.meshgrid(x, y, indexing='ij')  # if indexing not defined, output will be 'xy'-indexed
    # 'ij'-indexing: X.shape = Y.shape = [len(x), len(y)]
    # 'xy'-indexing: X.shape = Y.shape = [len(y), len(x)]
    # X[:,j] = x, Y[i,:] = y
    a = 1
    b = 2
    f = a * X + b * Y
    print X
    print Y
    print f

    i = 0
    j = 0
    I = integrate.simps(f[:, j], X[:, j])
    print I
    I = integrate.simps(f[:, j], x, axis=0)
    print I
    I = integrate.simps(f, x, axis=0)
    print I
    I2 = integrate.simps(f, y, axis=1)
    print I2


    a = 2
    b = 1
    # integrand(x, 1, 0)
    psi = integrate.quad(integrand, 0, 1, args=(a,b))
    print('psi test', psi)


    return

# ----------------------------------------------------------------------
# PLOTTING

def plot_comparison_vort_vort_stag(vort_yz, vort_yz_stag, vort_yz_stag_, jcshift, jc, kmax, t0):
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))

    ax = axes[0, 0]
    cf = ax.imshow(vort_yz.T, origin='lower')
    ax.set_title('vort_yz')
    ax.plot([jc, jc], [1, ny_ - 2], 'k')
    plt.colorbar(cf, ax=ax, shrink=0.6)

    ax = axes[0, 1]
    cf = ax.imshow(vort_yz_stag.T, origin='lower')
    ax.set_title('vort_yz staggered')
    plt.colorbar(cf, ax=ax, shrink=0.6)
    ax.plot([jc, jc], [1, ny_ - 2], 'k')

    ax = axes[0, 2]
    cf = ax.imshow(vort_yz_stag_.T, origin='lower')
    ax.set_title('vort_yz_ staggered')
    ax.plot([jcshift, jcshift], [1, ny_ - 2], 'k')
    plt.colorbar(cf, ax=ax, shrink=0.6)

    for i in [0, 1]:
        axes[0, i].set_xlim([jc - np.int(ny_ / 2), jc + np.int(ny_ / 2)])
    for i in [2]:
        axes[0, i].set_xlim([0, ny_])
    for i in range(3):
        axes[0, i].set_ylim([0, kmax])

    for i in range(3):
        ax = axes[1, i]
        ax.plot(np.arange(0, ny_), vort_yz[jc - np.int(ny_ / 2):jc + np.int(ny_ / 2), i], label='vort_yz')
        ax.plot(np.arange(0, ny_), vort_yz_stag[jc - np.int(ny_ / 2):jc + np.int(ny_ / 2), i], label='vort_yz stag')
        ax.plot(np.arange(0, ny_), vort_yz_stag_[:ny_, i], label='vort_yz_ stag')
        ax.legend(loc='best', fontsize=8)
        ax.set_title('k='+str(i))
        ax.set_xlabel('y')


    for i in range(3):
        ax = axes[2, i]
        ax.plot(np.arange(0, ny_), vort_yz_stag[jc - np.int(ny_ / 2):jc + np.int(ny_ / 2), i] - vort_yz_stag_[:ny_, i],
                label='diff')
        ax.set_title('diff k='+str(i))
        ax.set_xlabel('y')

    plt.suptitle('comparison vort_yz unstaggered vs. staggered')
    fig.savefig(os.path.join(path_out_figs, 'test_vort_t' + str(t0) + 's.png'))
    plt.close(fig)
    return


def comparison_vort_yz_vort_xz(vort_xz_, vort_yz_, kmax, t0):
    fig, axes = plt.subplots(1,3, figsize=(14,4), sharey='all')
    y_arr = np.arange(0,ny_)
    z_arr = np.arange(0,kmax)
    Y_, Z_ = np.meshgrid(y_arr, z_arr)
    plt.suptitle('t='+str(t0)+'s')
    ax1 = axes[0]
    # a = ax1.contourf(Y_,Z_,vort_xz_[:,:kmax])
    a = ax1.imshow(vort_xz_[:,:].T, origin='lower')
    ax1.autoscale(False)
    # ax1.set_ylim([0,kmax])
    # ax1.set_xlim([0,ny_])
    ax1.set_title('vort_xz')
    ax2 = axes[1]
    b = ax2.imshow(vort_yz_[:,:].T, origin='lower')
    ax2.set_title('vort_yz')
    ax3 = axes[2]
    c = ax3.imshow(vort_xz_[:,:].T-vort_yz_[:,:].T, origin='lower')
    ax3.set_title('vort_xz - vort_yz')
    plt.colorbar(a, ax=ax1)
    plt.colorbar(b, ax=ax2)
    plt.colorbar(c, ax=ax3)
    plt.savefig(os.path.join(path_out_figs, 'vort_xz_yz_t'+str(t0)+'s.png'))
    plt.close()
    return

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

def set_input_output_parameters(args):
    print('--- set input parameters ---')
    global path, path_in, path_out_figs, path_out_fields, path_stats, path_fields
    path = args.path
    # path = '/Users/bettinameyer/polybox/ClimatePhysics/Copenhagen/Projects/LES_ColdPool/' \
    #        'triple_3D_noise/Out_CPDry_triple_dTh2K/'
    # path = '/nbi/ac/cond1/meyerbe/ColdPools/triple_3D_noise/Out_CPDry_triple_Th10K/'
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

    shift = 60
    nx_half = irstar + shift
    ny_half = irstar + shift
    ishift = np.max(nx_half - ic, 0)
    jshift = np.max(ny_half - jc, 0)
    nx_ = 2 * nx_half
    ny_ = 2 * ny_half


    print('rstar: ' + str(rstar), irstar)
    print('ic,jc,id,jd', ic, jc, nx_half, ny_half)
    print('nx_,ny_', nx_, ny_)
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
    return


def dump_statistics_variable(var_in, var_name, grp_name, file_name, path):
    rootgrp = nc.Dataset(os.path.join(path, file_name), 'r+', format='NETCDF4')
    grp = rootgrp.groups[grp_name]
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