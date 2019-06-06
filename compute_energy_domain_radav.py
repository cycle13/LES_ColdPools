import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import netCDF4 as nc
import argparse
import json as simplejson
import os

# compute potential temperature by integrating over anomaly
#   PE = \int dz g * (th_anomaly(z) - th_env(z)) * z
#   KE ~ v**2 = (v_rad**2 + v_tan**2 + w**2)

# label_size = 8
# plt.rcParams['xtick.labelsize'] = label_size
# plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
# plt.rcParams['legend.fontsize'] = 8
# plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35

# COMPUTATION OF POTENTIAL ENERGY (PE) AND KINETIC ENERGY (KE)
# - compute PE from azimuthally averaged entropy statistics
# - compute PE from 3D field, using buoyancy
# - compute absolute PE from 3D field, using density
# - compute KE from azimuthally averaged velocity statistics
# - compute KE from 3D velocity fields (with and without interpolation)
# - compute SGS KE from Stats-file (ww) and using the viscosity



def main():
    # Parse information from the command line
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

    filename_in = 'stats_radial_averaged.nc'

    nml, times = set_input_output_parameters(args, filename_in)
    define_geometry(case_name, nml)
    path_fields = os.path.join(path, 'fields')
    path_stats = os.path.join(path, 'data_analysis')
    path_figs = os.path.join(path, 'figs_CP_energy')
    if not os.path.exists(path_figs):
        os.mkdir(path_figs)
    if case_name == 'ColdPoolDry_single_3D':
        ic = ic_arr[0]
        jc = jc_arr[0]

    id = os.path.basename(path)
    print('id: ', id)
    filename_out = 'CP_energy_' + id + '_domain_radav.nc'

    ''' create output file '''
    create_output_file(times, filename_in, path_stats, filename_out, path_stats)

    ''' (A) Domain '''
    ''' (A1) Potential Energy (PE) '''
    PE = compute_PE_from_radav(times, id, filename_in, filename_out,
                    path_fields, path_stats, path_figs)
    PE = compute_PE_from_fields(times, id, filename_in, filename_out,
                    path_fields, path_stats, path_figs)
    PE_abs = compute_PE_absolute(times, filename_out, path_fields, path_stats, path_figs)

    ''' (A2) Kinetic Energy (KE) '''
    KE, KEd, KE_r = compute_KE_from_radav(times, id, filename_in, filename_out, path_fields, path_stats, path_figs)
    KE_3D = compute_KE_from_fields(times, id, filename_in, filename_out, path_fields, path_stats, path_figs)
    KE_vel_rad = compute_KE_from_v_rad_2D(times, id, filename_in, filename_out, path_fields, path_stats, path_figs)
    ''' (A3) SGS Kinetic Energy (KE) '''
    # compute SGS momentum fluxes from viscosity profile and 3D velocity fields
    KE_sgs = compute_KE_sgs(times, id, filename_in, filename_out, path_fields, path_stats, path_figs)

    ''' plotting '''
    plot_KE_PE_radav(id, filename_out, path_stats, path_figs)
    plot_KE_PE(id, filename_out, path_stats, path_figs)
    plot_KE(id, filename_out, path_stats, path_figs)
    plot_PE(id, filename_out, path_stats, path_figs)
    return

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def compute_KE_from_v_rad_2D(times, id, filename_in, filename_out, path_fields, path_stats, path_figs):
    print ''
    print('--- compute KE from v_rad ---')
    # read in v_rad 2D fields
    # # read in v_rad-field
    # path_out_data_2D = os.path.join(path, 'fields_v_rad')
    rootgrp = nc.Dataset(os.path.join(path, 'fields_v_rad', 'v_rad.nc'))
    v_rad_2D = rootgrp.variables['v_rad'][:, :, :, :]
    v_tan_2D = rootgrp.variables['v_tan'][:, :, :, :]
    rootgrp.close()

    # ''' output '''
    # rootgrp = nc.Dataset(os.path.join(path_stats, filename_out), 'r+', format='NETCDF4')
    # ts_grp = rootgrp.groups['timeseries']
    # prof_grp = rootgrp.groups['profiles']
    # var = ts_grp.variables['KE']
    # var[:] = KE[:]
    # var = ts_grp.variables['KEd']
    # var[:] = KEd[:]
    # var = ts_grp.variables['KE_sgs_stats']
    # var[:] = KE_sgs[:]
    # var = prof_grp.variables['KE_r']
    # var[:, :] = KE_r[:, :]
    # rootgrp.close()

    return



def compute_KE_from_fields(times, id, filename_in, filename_out, path_fields, path_stats, path_figs):
    print ''
    print('--- compute KE from 3D fields ---')
    nt = len(times)
    print('nt', nt)
    KE = np.zeros((nt), dtype=np.double)
    KE_3D = np.zeros((nt, nx, ny ,kmax), dtype=np.double)
    KE_x = np.zeros((nt, nx), dtype=np.double)  # compute KE[t, r, :] (columnwise integration over z)
    KE_x_int = np.zeros((nt, nx), dtype=np.double)  # compute KE[t, r, :] (columnwise integration over z)
    # interpolated velocity fields
    KE_int = np.zeros((nt), dtype=np.double)
    KE_3D_int = np.zeros((nt, nx, ny ,kmax), dtype=np.double)
    u_int = 0.
    v_int = 0.
    w_int = 0.

    # 1. read in reference density
    rootgrp = nc.Dataset(os.path.join(path, 'stats', 'Stats.' + case_name + '.nc'))
    rho0 = rootgrp.groups['reference'].variables['rho0'][:]
    # rho_unit = rootgrp.groups['reference'].variables['rho0'].units
    rootgrp.close()

    # 2. read in 3D fields & compute KE
    # KE ~ v**2 = (u**2 + v**2 + w**2)
    for it,t0 in enumerate(times):
        print('--t=' + str(t0) + '--')
        u = read_in_netcdf_fields('u', os.path.join(path_fields, str(t0)+'.nc'))[:,:,:kmax]
        v = read_in_netcdf_fields('v', os.path.join(path_fields, str(t0)+'.nc'))[:,:,:kmax]
        w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0)+'.nc'))[:,:,:kmax]
        u2 = u * u
        v2 = v * v
        w2 = w * w
        # del u, v, w

        jc = jc_arr[0]
        for i in range(nx):
            for j in range(ny):
                for k in range(kmax):
                    if i == 0:
                        u_int = 0.5* ( u[i,j,k] + u[-1,j,k] )
                    else:
                        u_int = 0.5* ( u[i,j,k] + u[i-1,j,k] )  # periodic bcs
                    if j == 0:
                        v_int = 0.5* ( v[i,j,k] + v[i,-1,k] )
                    else:
                        v_int = 0.5* ( v[i,j,k] + v[i,j-1,k] )  # periodic bcs
                    if k > 0:
                        w_int = 0.5* ( w[i,j,k] + w[i,j,k-1] )
                    else:
                        w_int = 0.5* w[i,j,k]      # asymmetric bcs: w[i,j,-1]=0

                    KE_3D[it,i,j,k] = 0.5 * dV * rho0[k] * (u2[i,j,k] + v2[i,j,k] + w2[i,j,k])
                    KE_3D_int[it,i,j,k] = 0.5 * dV * rho0[k] * ( u_int*u_int + v_int*v_int + w_int*w_int )
            # KE_x[it, i] = 0.5 * dV * np.sum(rho0[:kmax] * (u2[i, jc, :] + v2[i, jc, :] + w2[i, jc, :]))
            KE_x[it, i] = np.sum(KE_3D[it,i,jc,:])
            KE_x_int[it, i] = np.sum(KE_3D_int[it,i,jc,:])

        # aux = np.sum(np.sum(u2[:,:,:] + v2[:,:,:] + w2[:,:,:], axis=0), axis=0)
        # KE[it] = 0.5 * dV * np.sum(rho0[:kmax] * aux)
        KE[it] = np.sum(KE_3D[it,:,:,:])
        KE_int[it] = np.sum(KE_3D_int[it,:,:,:])

    # compute KE_r from KE[ijk] average_angular
    # compare KE_r to KE_x
    from average_angular import compute_angular_average

    it = 3
    fig, axis = plt.subplots(2,3, figsize=(25,15))
    ax = axis[0,0]
    ax.set_title('diff KE interpolation: max=' + str(np.round(np.amax(KE_3D[it,:,jc,1:]-KE_3D_int[it,:,jc,1:])))
                 + ', min=' + str(np.round(np.amin(KE_3D[it,:,jc,1:]-KE_3D_int[it,:,jc,1:]))) + ')' )
    ax.set_xlabel('z')
    cf = ax.contourf(np.arange(nx)*dx[0], np.arange(kmax)*dx[2], (KE_3D[it,:,jc,:]-KE_3D_int[it,:,jc,:]).T )
    plt.colorbar(cf, ax=ax, shrink=0.6, aspect=12)

    ax = axis[0, 1]
    ax.set_title('KE crosssection')
    ax.set_xlabel('z')
    ax.set_ylabel('x')
    cf = ax.contourf(np.arange(nx) * dx[0], np.arange(kmax) * dx[2], KE_3D[it, :, jc, :].T)
    plt.colorbar(cf, ax=ax, shrink=0.6, aspect=12)

    ax = axis[0, 2]
    ax.set_title('KE')
    ax.plot(times, KE, label='no int')
    ax.plot(times, KE_int, label='int')
    ax.plot(times, np.sum(np.sum(np.sum(KE_3D[:,:,:,1:],axis=1),axis=1),axis=1), label='no int (k>0)')
    ax.plot(times, np.sum(np.sum(np.sum(KE_3D_int[:,:,:,1:],axis=1),axis=1),axis=1), label='int (k>0)')
    ax.legend(loc='best')

    ax = axis[1,0]
    ax.set_title('diff KE interpolation: max=' + str(np.round(np.amax(KE_x-KE_x_int)))
                 + ', min=' + str(np.round(np.amin(KE_x-KE_x_int))) + ')')
    ax.set_xlabel('z')
    cf = ax.contourf(times, np.arange(nx) * dx[0], (KE_x[:, :]-KE_x_int[:, :]).T)
    plt.colorbar(cf, ax=ax, shrink=0.6, aspect=12)

    ax = axis[1,1]
    ax.set_title('KE_x')
    ax.set_xlabel('time')
    ax.set_ylabel('x')
    cf = ax.contourf(times, np.arange(nx)*dx[0], (KE_x[:,:]).T )
    plt.colorbar(cf, ax=ax, shrink=0.6, aspect=12)

    ax = axis[1,2]
    ax.set_title('KE_x int')
    ax.set_xlabel('time')
    ax.set_ylabel('x')
    cf = ax.contourf(times, np.arange(nx)*dx[0], (KE_x_int[:,:]).T )
    plt.colorbar(cf, ax=ax, shrink=0.6, aspect=12)



    plt.suptitle('t='+str(times[it]))
    plt.tight_layout()
    plt.savefig(os.path.join(path_figs, 'KE_3D_' + id + '.png'))
    plt.close()



    ''' output '''
    rootgrp = nc.Dataset(os.path.join(path_stats, filename_out), 'r+', format='NETCDF4')
    ts_grp = rootgrp.groups['timeseries']
    var = ts_grp.variables['KE_3D']
    var[:nt] = KE[:]
    var = ts_grp.variables['KE_3D_int']
    var[:nt] = KE_int[:]
    rootgrp.close()



    return KE



def compute_KE_from_radav(times, id, filename_in, filename_out, path_fields, path_stats, path_figs):
    print ''
    print('--- compute KE ---')
    # 1. read in reference rho
    # 2. read in velocity fields
    # 3. read in SGS momentum fluxes
    # 4. integrate: KE = 0.5*sum_i(rho_i*dV*v_i**2) from center (ic,jc) to rim
    nt = len(times)
    print('nt', nt)

    # 1. read in reference density
    rootgrp = nc.Dataset(os.path.join(path, 'stats', 'Stats.' + case_name + '.nc'))
    rho0 = rootgrp.groups['reference'].variables['rho0'][:]
    # rho_unit = rootgrp.groups['reference'].variables['rho0'].units
    rootgrp.close()


    # 2. (A) read in azimuthally averaged s-field
    file_radav = nc.Dataset(os.path.join(path, path_stats, filename_in))
    time_in = file_radav.groups['timeseries'].variables['time'][:]
    radius = file_radav.groups['stats'].variables['r'][:]
    krange = file_radav.groups['dimensions'].variables['krange'][:]
    v_rad = file_radav.groups['stats'].variables['v_rad'][:, :, :]  # v(nt, nr, nz)
    v_tan = file_radav.groups['stats'].variables['v_tan'][:, :, :]  # v(nt, nr, nz)
    w_rad = file_radav.groups['stats'].variables['w'][:, :, :]      # w(nt, nr, nz)
    file_radav.close()
    nk = np.int(np.minimum(kmax, krange[-1]))
    nr = len(radius)


    # v2_rad = v_rad * v_rad * (2*np.pi*radius**2)
    v2_rad = v_rad * v_rad
    v2_tan = v_tan * v_tan
    w2_av = w_rad * w_rad
    dxi = 1./dx[0]
    rad2 = radius * radius * dxi * dxi
    print 'radius', radius*dxi

    # 3. read in and compute SGS momentum fluxes
    file_stats = nc.Dataset(os.path.join(path, 'stats', 'Stats.'+case_name+'.nc'), 'r')
    # uu_sgs_stats = file_stats.groups['profiles'].variables['u_sgs_flux_x'][:,:]
    # vv_sgs_stats = file_stats.groups['profiles'].variables['v_sgs_flux_y'][:,:]
    ww_sgs_stats = file_stats.groups['profiles'].variables['w_sgs_flux_z'][:,:]
    # visc_mean = file_stats.groups['profiles'].variables['viscosity_mean'][:,:]
    file_stats.close()


    # 4. compute KE
    # KE ~ v**2 = (v_rad**2 + v_tan**2 + w**2)
    KE = np.zeros((nt))
    # KE_aux = np.zeros((nt))
    KEd = np.zeros((nt))
    KE_r = np.zeros((nt, nr))           # compute KE[t, r, :] (columnwise integration over z)
    KE_sgs = np.zeros((nt))
    for it, t0 in enumerate(times):
        print('--t=' + str(t0) + '--')
        # aux = np.sum(v2_rad[it,:,:kmax]+v2_tan[it,:,:kmax]+w2_av[it,:,:kmax], 0)
        aux = 4*np.pi**2 * np.sum(v2_rad[it,:,:kmax].T*rad2[:] + v2_tan[it,:,:kmax].T*rad2[:] + w2_av[it,:,:kmax].T*rad2[:], 1)
        KEd[it] = 0.5*np.sum(aux)
        KE[it] = 0.5 * dV * np.sum(rho0[:kmax] * aux)
        KE_sgs[it] = 0.5 * dV * np.sum(rho0[:kmax] * ww_sgs_stats[it,:kmax])

        for i in range(nr):
            KE_r[it, i] = 0.5 * dV * np.sum(rho0[:kmax] * (v2_rad[it,i,:kmax]*rad2[i]
                                                           + v2_tan[it,i,:kmax]*rad2[i]
                                                           + w2_av[it,i,:kmax]*rad2[i]) )
            # KE_aux[it] += KE_r[it, i]





    fig, axis = plt.subplots(3, 3, figsize=(15, 15), sharey='none')
    rmax = 100
    ax1 = axis[0, 0]
    ax2 = axis[0, 1]
    ax3 = axis[0, 2]
    ax1.set_title('v_rad^2')
    ax2.set_title('v_tan^2')
    ax3.set_title('w^2')
    print(times, len(times), nt, v2_rad.shape)
    cf = ax1.contourf(times, radius[:rmax], v2_rad[:nt, :rmax, 0].T)
    cbar = plt.colorbar(cf, ax=ax1, shrink=0.6, aspect=12)
    cf = ax2.contourf(times, radius[:rmax], v2_tan[:nt, :rmax, 0].T)
    cbar = plt.colorbar(cf, ax=ax2, shrink=0.6, aspect=12)
    cf = ax3.contourf(times, radius[:rmax], w2_av[:nt, :rmax, 0].T)
    cbar = plt.colorbar(cf, ax=ax3, shrink=0.6, aspect=12)
    ax1.set_xlabel('time [s]')
    ax2.set_xlabel('time [s]')
    ax3.set_xlabel('time [s]')
    ax1.set_ylabel('radius [m]')

    ax1 = axis[1, 0]
    ax2 = axis[1, 1]
    ax3 = axis[1, 2]
    ax1.set_title('v_rad^2')
    ax2.set_title('v_tan^2')
    ax3.set_title('w^2')
    print(times, len(times), nt, v2_rad.shape)
    aux = v2_rad[:nt, :rmax, 0]*rad2[:rmax]*2*np.pi
    cf = ax1.contourf(times, radius[:rmax], aux.T)
    cbar = plt.colorbar(cf, ax=ax1, shrink=0.6, aspect=12)
    ax1.contour(times, radius[:rmax], v2_rad[:nt, :rmax, 0].T, colors='w')
    aux = v2_tan[:nt, :rmax, 0]*rad2[:rmax]*2*np.pi
    cf = ax2.contourf(times, radius[:rmax], aux.T)
    cbar = plt.colorbar(cf, ax=ax2, shrink=0.6, aspect=12)
    aux = w2_av[:nt, :rmax, 0]*rad2[:rmax]*2*np.pi
    cf = ax3.contourf(times, radius[:rmax], aux.T)
    cbar = plt.colorbar(cf, ax=ax3, shrink=0.6, aspect=12)
    ax1.set_xlabel('time [s]')
    ax2.set_xlabel('time [s]')
    ax3.set_xlabel('time [s]')
    ax1.set_ylabel('radius [m]')

    ax1 = axis[2, 0]
    ax2 = axis[2, 1]
    ax3 = axis[2, 2]
    ax1.plot(times, KEd)
    ax2.plot(times, KE)
    # ax2.plot(times, KE_aux, 'r')
    # ax2.plot(times, 2*np.pi*KE_aux, 'g')
    cf = ax3.contourf(times, radius[:rmax], KE_r[:nt, :rmax].T)
    cbar = plt.colorbar(cf, ax=ax3, shrink=0.6, aspect=12)
    ax1.set_title('KEd')
    ax2.set_title('KE')
    # # ax2.set_ylabel('radius [m]')
    # # ax3.set_ylabel('radius [m]')

    plt.savefig(os.path.join(path_figs, 'KE_vel_' + id + '.png'))
    plt.close(fig)





    ''' output '''
    rootgrp = nc.Dataset(os.path.join(path_stats, filename_out), 'r+', format='NETCDF4')
    ts_grp = rootgrp.groups['timeseries']
    prof_grp = rootgrp.groups['profiles']
    var = ts_grp.variables['KE']
    var[:] = KE[:]
    var = ts_grp.variables['KEd']
    var[:] = KEd[:]
    var = ts_grp.variables['KE_sgs_stats']
    var[:] = KE_sgs[:]
    var = prof_grp.variables['KE_r']
    var[:, :] = KE_r[:, :]
    rootgrp.close()

    return KE, KEd, KE_r



def compute_KE_sgs(times, id, filename_in, filename_out, path_fields, path_stats, path_figs):
    print ''
    print('--- compute SGS KE ---')
    from compute_energy_domain_sgs_flux import compute_diffusive_flux

    # 0. get dimensions
    file_radav = nc.Dataset(os.path.join(path, path_stats, filename_in))
    # time_in = file_radav.groups['timeseries'].variables['time'][:]
    radius = file_radav.groups['stats'].variables['r'][:]
    krange = file_radav.groups['dimensions'].variables['krange'][:]
    file_radav.close()
    nt = np.int(len(times))
    nk = np.int(np.minimum(kmax, krange[-1]))
    print('nt', nt)


    # 1. read in reference density and mean viscosity
    rootgrp = nc.Dataset(os.path.join(path, 'stats', 'Stats.' + case_name + '.nc'))
    rho0 = rootgrp.groups['reference'].variables['rho0'][:kmax]
    visc_mean = rootgrp.groups['profiles'].variables['viscosity_mean'][:,:]
    rootgrp.close()

    # 2. compute fluxes and KE
    # uu_sgs_flux = np.zeros((nt, nx, ny, nk))
    # vv_sgs_flux = np.zeros((nt, nx, ny, nk))
    # ww_sgs_flux = np.zeros((nt, nx, ny, nk))
    uu_sgs_flux = np.zeros((nt, kmax))
    vv_sgs_flux = np.zeros((nt, kmax))
    ww_sgs_flux = np.zeros((nt, kmax))
    KE_sgs = np.zeros((nt))
    for it, t0 in enumerate(times):
        print('--t=' + str(t0) + '----')
        # strain_rate = compute_strain_rate(d_ing, d_ed, nx, ny, nk, dx, t0, path_fields)
        # uu_sgs_flux[it, :, :, :] = compute_diffusive_flux(visc_mean[it, :nk], rho0[:nk], 0, 0, nx, ny, nk, dx, t0,
        #                                                   path_fields)
        # vv_sgs_flux[it, :, :, :] = compute_diffusive_flux(visc_mean[it, :nk], rho0[:nk], 1, 1, nx, ny, nk, dx, t0,
        #                                                   path_fields)
        # ww_sgs_flux[it, :, :, :] = compute_diffusive_flux(visc_mean[it, :nk], rho0[:nk], 2, 2, nx, ny, nk, dx, t0,
        #                                                   path_fields)

        aux = compute_diffusive_flux(visc_mean[it, :kmax], rho0[:kmax], 0, 0, nx, ny, kmax, dx, t0, path_fields)
        uu_sgs_flux[it,:] = np.sum(np.sum(aux, axis=0), axis=0)
        aux = compute_diffusive_flux(visc_mean[it, :kmax], rho0[:kmax], 1, 1, nx, ny, kmax, dx, t0, path_fields)
        vv_sgs_flux[it,:] = np.sum(np.sum(aux, axis=0), axis=0)
        aux = compute_diffusive_flux(visc_mean[it, :kmax], rho0[:kmax], 2, 2, nx, ny, kmax, dx, t0, path_fields)
        ww_sgs_flux[it,:] = np.sum(np.sum(aux, axis=0), axis=0)

        KE_sgs[it] = np.sum(rho0*uu_sgs_flux[it,:kmax]) \
                     + np.sum(rho0*vv_sgs_flux[it,:kmax])\
                     + np.sum(rho0*ww_sgs_flux[it,:kmax])
    del aux
    KE_sgs = 0.5 * dV * KE_sgs

    ''' output '''
    rootgrp = nc.Dataset(os.path.join(path_stats, filename_out), 'r+', format='NETCDF4')
    ts_grp = rootgrp.groups['timeseries']
    var = ts_grp.variables['KE_sgs']
    var[:] = KE_sgs[:]

    KE = ts_grp.variables['KE'][:]
    KE_sgs_stats = ts_grp.variables['KE_sgs_stats'][:]
    rootgrp.close()




    fig, axis = plt.subplots(1, 3, figsize=(15, 5))
    ax1 = axis[0]
    ax2 = axis[1]
    ax3 = axis[2]
    # ax3 = axis[2]
    # ax4 = axis[3]
    ax1.set_title('KE')
    ax2.set_title('KE sgs')
    ax3.set_title('viscosity')

    ax1.plot(times, KE, label='KE res')
    ax2.plot(times, KE_sgs, label='KE sgs')
    ax2.plot(times, KE_sgs_stats, label='KE sgs (stats)')
    ax1.legend(loc=0)
    ax2.legend(loc=1)
    print  len(times), len(krange[:kmax]), visc_mean.shape, visc_mean[:,:kmax].shape
    # cf = ax3.contourf(times, krange[:kmax], visc_mean[:,:kmax])
    cf = ax3.contourf(visc_mean[:,:kmax])
    # cbar = plt.colorbar(cf, shrink=0.6, ticks=np.arange(min, max + 1, 1), aspect=12)
    cbar = plt.colorbar(cf, ax=ax3, shrink=0.6, aspect=12)
    ax1.set_xlabel('time [s]')
    ax2.set_xlabel('time [s]')
    ax3.set_xlabel('time [s]')
    ax1.set_ylabel('radius [m]')
    ax2.set_ylabel('KE [J]')
    ax3.set_ylabel('KE [J]')

    plt.suptitle('Kinetic Energy')
    plt.tight_layout()
    plt.savefig(os.path.join(path_figs, 'KE_sgs_' + id + '.png'))
    plt.close(fig)


    return KE_sgs



def compute_PE_from_fields(times, id, filename_in, filename_out,
               path_fields, path_stats, path_figs):
    # 3. define background profile (here: take profile at any point outside the anomaly region)
    print ''
    print('--- compute PE from fields---')
    nt = len(times)
    print('nt', nt)

    # 3. define background profile (here: take profile at any point outside the anomaly region)
    s_ = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
    i0 = 10
    j0 = 10
    di = 5
    s_env = np.average(np.average(s_[i0 - di:i0 + di, j0 - di:j0 + di, :kmax], axis=0), axis=0)
    theta_env = theta_s(s_env)
    th_g = theta_env[0]
    del s_, s_env
    rootgrp = nc.Dataset(os.path.join(path, 'stats', 'Stats.' + case_name + '.nc'))
    rho0 = rootgrp.groups['reference'].variables['rho0'][:]
    rootgrp.close()

    # 4. integrate
    g = 9.80665
    # PEd = PE/kg = sum(g*dz*dTh_i) = g*dz*sum(dTh_i)
    # [PE/kg] = m/s^2*m = (m/s)^2
    # PE = m*a*s        >>  [PE] = kg*m/s^2*m = kg*(m/s)^2
    PE_unround = np.zeros(nt, dtype=np.double)
    PE = np.zeros(nt, dtype=np.double)
    PE_t = np.zeros(nt, dtype=np.double)
    theta_env_t = np.zeros((nt,kmax), dtype=np.double)
    theta_diff_unround = np.zeros((nx,ny,kmax), dtype=np.double)
    theta_diff = np.zeros((nx,ny,kmax), dtype=np.double)
    for it, t0 in enumerate(times):
        print('--t=' + str(t0) + '--')
        s_in = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))[:, :, :kmax]
        theta = theta_s(s_in)
        theta_env_t[it] = np.average(np.average(theta[i0 - di:i0 + di, j0 - di:j0 + di, :kmax], axis=0), axis=0)
        for k in range(kmax):
            diff = theta_env[k] - theta[:,:,k]
            diff_t = theta_env_t[it,k] - theta[:,:,k]
            # theta_diff_unround[:, :, k] = theta_env[k] - theta[:,:,k]
            theta_diff_unround[:, :, k] = diff
            theta_diff[:,:,k] = np.ma.masked_where(theta_diff_unround[:,:,k] <= 0.02, theta_diff_unround[:,:,k])
            PE += z_half[k] * np.sum(theta_diff[:,:, k]) * dV * rho0[k]
            PE_unround += z_half[k] * np.sum(theta_diff_unround[:,:,k]) * dV * rho0[k]
            PE_t += z_half[k] * np.sum(diff_t) * dV * rho0[k]
    PE = g/th_g * PE
    PE_unround = g/th_g * PE_unround


    fig, axis = plt.subplots(2, 5, figsize=(20, 12), sharey='row')
    rmax = 60
    for i,it in enumerate([0, 1, 3, 5, 10]):
        if it < nt:
            ax = axis[0,i]
            min = np.amin(theta_diff_unround[it, :,:], axis=1)
            max = np.amax(theta_diff_unround[it,:,:], axis=1)
            # coord = np.argmax(theta_diff_unround[it,:,:], axis=1)
            print('min, max', np.amin(min), np.amax(max))

            ax = axis[1,i]
    #         ax.plot(radius, max, linewidth=3, label='max')
    #         ax.plot(radius, -1e2*min, linewidth=3, label='-min')
    #         ax.plot(radius, -1e2*np.round(min,1), linewidth=3, label='-min rounded')
    #         ax.plot([radius[282],radius[282]], [0,300], 'k-')
    #         ax.legend(loc=0)
    #         if i==0:
    #             ax.set_ylim(0,np.amax(max)+0.1)
    #         ax.set_title('t='+str(times[it]) +
    #                       ' (max(theta_diff): ir='+str(np.argmax(np.amax(theta_diff[it,:,:], axis=1)))+')')
    #
    #         ax = axis[0,i]
    #         min = np.amin(theta[it, :, :kmax])
    #         max = np.amax(theta[it, :, :kmax])
    #         cf = ax.contourf(radius[:nr], np.arange(kmax), theta_diff[it,:nr,:kmax].T,
    #                          levels=np.linspace(-4,4,9))#, norm=LogNorm())
    #         cbar = plt.colorbar(cf, ax=ax, shrink=0.6, aspect=20)
    #         ax.set_title('theta[r,z] - theta_env[z]')
    #         # ax.plot([radius[282],radius[282]], [0,300], 'k-')
    #         # ax.set_title('ir='+str(np.argmax(np.amax(theta_diff[it,:,:], axis=1))))

    plt.tight_layout()
    plt.savefig(os.path.join(path_figs, 'PE_theta_test_3D.png'))
    plt.close(fig)



    fig, axis = plt.subplots(1, 3, figsize=(24, 8))
    ax0 = axis[0]
    ax1 = axis[1]
    ax2 = axis[2]
    ax1.set_title('PE density')
    ax2.set_title('PE')

    for it,t0 in enumerate(times):
        if it == 0 or it == nt-1:
            ax0.plot(theta_env_t[it], np.arange(kmax)*dx[2], color=cm_bwr(np.double(it)/nt), label='theta_env(t='+str(t0)+')')
        else:
            ax0.plot(theta_env_t[it], np.arange(kmax)*dx[2], color=cm_bwr(np.double(it)/nt))
        ax0.plot([theta_env_t[it,0], theta_env_t[it,0]], [0,kmax*dx[2]], '--', color=cm_bwr(np.double(it)/nt))
    ax0.legend(loc=1)
    ax0.set_xlabel('theta  [K]')
    ax0.set_ylabel('z')
    ax0.plot(theta_env, np.arange(kmax) * dx[2], 'k', label='theta_env(z)')
    ax0.plot([th_g, th_g], [0, kmax * dx[2]], 'k', label='theta_env(z=0)')

    ax1.plot(times, PE, label='PE')
    ax1.plot(times, PE_unround, label='PE unround')
    ax1.plot(times, PE_t, label='PE_t')
    # cbar = plt.colorbar(cf, ax=ax1, shrink=0.6, aspect=12)
    ax1.legend()

    ax2.plot(times, PE)

    ax1.set_xlabel('time [s]')
    ax2.set_xlabel('time [s]')
    ax1.set_ylabel('radius [m]')
    ax2.set_ylabel('PE [J]')

    plt.suptitle('Potential Energy')
    plt.tight_layout()
    plt.savefig(os.path.join(path_figs, 'PE_' + id + '.png'))
    plt.close(fig)



    ''' output '''
    rootgrp = nc.Dataset(os.path.join(path_stats, filename_out), 'r+', format='NETCDF4')
    ts_grp = rootgrp.groups['timeseries']
    var = ts_grp.variables['PE_3D']
    var[:nt] = PE[:]
    rootgrp.close()

    return



def compute_PE_from_radav(times, id, filename_in, filename_out,
               path_fields, path_stats, path_figs):
    print ''
    print('--- compute PE from angular average ---')
    # 1. read in initial s-field
    # 2. convert entropy to potential temperature
    # 3. ??? convert potential temperature to density
    # 4. define background profile (here: take profile at any point outside the anomaly region)
    # 5. integrate

    # 1. read in azimuthally averaged s-field
    file_radav = nc.Dataset(os.path.join(path, 'data_analysis', filename_in))
    time_in = file_radav.groups['timeseries'].variables['time'][:]
    radius = file_radav.groups['stats'].variables['r'][:]
    krange = file_radav.groups['dimensions'].variables['krange'][:kmax]
    s_in = file_radav.groups['stats'].variables['s'][:, :, :]       # s(nt, nr, nz)
    file_radav.close()
    nt = len(times)
    print('nt', nt)
    nr = len(radius)

    # 2. convert entropy to potential temperature
    theta = theta_s(s_in)

    # 3. define background profile (here: take profile at any point outside the anomaly region)
    s_ = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
    i0 = 10
    j0 = 10
    di = 5
    s_env = np.average(np.average(s_[i0-di:i0+di,j0-di:j0+di,:kmax], axis=0), axis=0)
    theta_env = theta_s(s_env)
    th_g = theta_env[0]
    del s_, s_env

    rootgrp = nc.Dataset(os.path.join(path, 'stats', 'Stats.'+case_name+'.nc'))
    rho0 = rootgrp.groups['reference'].variables['rho0'][:]
    # rho_unit = rootgrp.groups['reference'].variables['rho0'].units
    # z_half = rootgrp.groups['reference'].variables['z'][:]
    rootgrp.close()


    # 4. integrate
    g = 9.80665
    # PEd = PE/kg = sum(g*dz*dTh_i) = g*dz*sum(dTh_i)
    # [PE/kg] = m/s^2*m = (m/s)^2
    # PE = m*a*s        >>  [PE] = kg*m/s^2*m = kg*(m/s)^2

    # int dz a(z) = sum_i a_i dz_i
    PE = np.zeros(nt)
    PEd = np.zeros(nt)
    PE_r = np.zeros((nt, nr))
    theta_diff = np.zeros((nt, nr, kmax))
    PE_unround = np.zeros(nt)
    theta_diff_unround = np.zeros((nt, nr, kmax))
    dxi = 1./dx[0]

    for ir,r in enumerate(radius[:-2]):
        for k in range(kmax):
            for it, t0 in enumerate(times):
                diff = theta_env[k] - theta[it, ir, k]
                theta_diff_unround[it,ir,k] = theta_env[k] - theta[it, ir, k]
                if np.abs(diff) > 0.02:
                    theta_diff[it,ir,k] = diff
            # PE +=  z_half[k] * theta_diff[:nt, ir, k] * dV*rho0[k]
            PE +=  2*np.pi*radius[ir]*dxi * z_half[k] * theta_diff[:nt, ir, k] * dV*rho0[k]
            PE_unround +=  z_half[k] * theta_diff_unround[:nt, ir, k] * dV*rho0[k]
            PEd += z_half[k] * theta_diff[:nt, ir, k]
            PE_r[:, ir] += z_half[k] * theta_diff[:nt,ir,k] * dV*rho0[k]
            # PE +=  z_half[k]*(theta_env[k] - theta[:nt,ir,k]) * dV*rho0[k]
            # PEd += z_half[k]*(theta_env[k] - theta[:nt,ir,k])
            # PE_r[:, ir] += z_half[k]*(theta_env[k] - theta[:nt,ir,k]) * dV*rho0[k]
    PEd = g/th_g * PEd
    PE = g/th_g * PE
    PE_unround = g / th_g * PE_unround
    # # PE_ = g*dz*PE
    # print('PE', PE, 'PEd', PEd)
    # print('density at 500m: ' + str(rho0[5]) + ' ' + rho_unit)
    # print('mass per grid cell at 500m: ' + str(dV * z_half[5]) + ' kg')



    fig, axis = plt.subplots(3, 5, figsize=(20, 12), sharey='row')
    rmax = 60
    for i,it in enumerate([0, 1, 3, 5, 10]):
        if it < nt:
            ax = axis[0,i]
            min = np.amin(theta[it, :, :kmax])
            max = np.amax(theta[it, :, :kmax])
            cf = ax.contourf(radius[:nr], np.arange(kmax), theta_diff[it,:nr,:kmax].T,
                             levels=np.linspace(-4,4,9))#, norm=LogNorm())
            cbar = plt.colorbar(cf, ax=ax, shrink=0.6, aspect=20)
            ax.set_title('theta[r,z] - theta_env[z]')
            ax.set_xlabel('x')
            ax.set_ylabel('z')

            ax = axis[1, i]
            min = np.amin(theta[it, :, :kmax])
            max = np.amax(theta[it, :, :kmax])
            cf = ax.contourf(radius[:nr], np.arange(kmax), theta_diff[it, :nr, :kmax].T*radius[:nr]*dxi,
                             levels=np.linspace(-4, 40, 9))  # , norm=LogNorm())
            cbar = plt.colorbar(cf, ax=ax, shrink=0.6, aspect=20)
            ax.set_title('theta[r,z] - theta_env[z]')
            ax.set_xlabel('x')
            ax.set_ylabel('z')


            ax = axis[2, i]
            min = np.amin(theta_diff_unround[it, :, :], axis=1)
            max = np.amax(theta_diff_unround[it, :, :], axis=1)
            coord = np.argmax(theta_diff_unround[it, :, :], axis=1)
            print('min, max', np.amin(min), np.amax(max))
            ax.plot(radius, max, linewidth=3, label='max')
            ax.plot(radius, -1e2 * min, linewidth=3, label='-min')
            ax.plot(radius, -1e2 * np.round(min, 1), linewidth=3, label='-min rounded')
            ax.plot([radius[282], radius[282]], [0, 300], 'k-')
            ax.legend(loc=0)
            if i == 0:
                ax.set_ylim(0, np.amax(max) + 0.1)
            ax.set_title('t=' + str(times[it]) +
                         ' (max(theta_diff): ir=' + str(np.argmax(np.amax(theta_diff[it, :, :], axis=1))) + ')')

    plt.tight_layout()
    plt.savefig(os.path.join(path_figs, 'PE_theta_test.png'))
    plt.close(fig)



    fig, axis = plt.subplots(1, 3, figsize=(15, 5))
    ax0 = axis[0]
    ax1 = axis[1]
    ax2 = axis[2]
    ax1.set_title('PE density')
    ax2.set_title('PE')

    ax0.plot(theta_env, krange, label='theta_env(z)')
    ax0.plot([th_g, th_g], [krange[0],krange[-1]], label='theta_env(z=0)')
    ax0.plot([theta_env[10], theta_env[10]], [krange[0],krange[-1]], label='theta_env(z=10)')
    ax0.legend(loc='best')

    cf = ax1.contourf(times, radius, PE_r.T, levels=np.linspace(0,np.amax(PE_r), 10))
    cbar = plt.colorbar(cf, ax=ax1, shrink=0.6, aspect=12)
    ax2.plot(times, PE)

    ax1.set_xlabel('time [s]')
    ax2.set_xlabel('time [s]')
    ax1.set_ylabel('radius [m]')
    ax2.set_ylabel('PE [J]')

    plt.suptitle('Potential Energy')
    plt.tight_layout()
    plt.savefig(os.path.join(path_figs, 'PE_' + id + '_radav.png'))
    plt.close(fig)


    ''' output '''
    rootgrp = nc.Dataset(os.path.join(path_stats, filename_out), 'r+', format='NETCDF4')
    ts_grp = rootgrp.groups['timeseries']
    var = ts_grp.variables['PE']
    var[:] = PE[:]
    var = ts_grp.variables['PEd']
    var[:] = PEd[:]
    rootgrp.close()

    return PE

# ----------------------------------------------------------------------

def compute_PE_absolute(times, filename_out, path_fields, path_stats, path_figs):
    print ''
    print('--- compute PE absolute ---')
    nt = len(times)
    print('kmax', kmax)
    # PE_ref = from s(t=0)
    g = 9.80665
    temp = np.ndarray((nx,ny,kmax+1),dtype=np.double)
    rho_3D = np.ndarray((nx,ny,kmax+1),dtype=np.double)
    PE = np.zeros((nt), dtype=np.double)
    PE_aux = np.zeros((nt,kmax), dtype=np.double)
    # reference slice
    dPdk_ref = np.zeros((nt,kmax), dtype=np.double)

    root_stats = nc.Dataset(os.path.join(path, 'stats', 'Stats.' + case_name + '.nc'))
    p0 = root_stats.groups['reference'].variables['p0'][:]
    z0 = root_stats.groups['reference'].variables['z'][:]
    root_stats.close()

    for it,t0 in enumerate(times):
        print('--t=' + str(t0) + '--')
        root_field = nc.Dataset(os.path.join(path_fields, str(t0)+'.nc'))
        temp[:,:,:] = root_field.groups['fields'].variables['temperature'][:,:,:kmax+1]
        root_field.close()
        for k in range(kmax):
            rho_3D[:, :, k] = compute_rho(p0[k], temp[:, :, k], 0., 0.)
            PE[it] += np.sum(rho_3D[:,:,k]) * z0[k]
            dPdk_ref[it,k] = g * dV * np.sum(rho_3D[:, :, k]) * z0[k]  # added PE for one slice at level k
        for kmax_ in range(kmax+1):
            for k in range(kmax_):
                PE_aux[it,kmax_-1] += np.sum(rho_3D[:,:,k]) * z0[k]

        # # print rho field
        # jc = jc_arr[0]
        # fig, (ax0, ax1, ax2) = plt.subplots(1,3, figsize=(18,5))
        # cf = ax0.contourf((rho_3D[:,jc,:kmax]-rho_3D[20,jc,:kmax]).T, levels=np.linspace(0.0,0.012,1e2))
        # plt.colorbar(cf, ax=ax0)
        # cf = ax1.contourf(rho_3D[:,jc,:kmax].T, levels=np.linspace(0.975,1.175,1e2))
        # plt.colorbar(cf, ax=ax1)
        # cf = ax2.contourf(temp[:,jc,:kmax].T, levels=np.linspace(280,300,1e2))
        # plt.colorbar(cf, ax=ax2)
        # plt.savefig(os.path.join(path_figs, 'PE_test_rho_t'+str(t0)+'.png'))
        # plt.close()

    PE = g * dV * PE
    PE_ref = PE[0]
    # dPE = PE - PE_ref
    PE_aux = g * dV * PE_aux
    PE_ref_aux = PE_aux[0, :]
    # dPE_aux = PE_aux - PE_ref

    # reference slice
    dPdk = np.zeros((nt,kmax), dtype=np.double)
    for k in range(1, kmax):
        dPdk[:, k] = PE_aux[:, k] - PE_aux[:, k-1]

    fig, axis = plt.subplots(1, 4, figsize=(20, 6), sharey='none')
    ax0 = axis[0]
    ax1 = axis[1]
    ax2 = axis[2]
    ax3 = axis[3]
    ax0.plot(times, PE[:], '-x', label='PE')
    ax0.plot(times, PE_aux[:,kmax-1], '--', label='PE aux')
    ax0.legend()
    for k in range(1,kmax):
        ax1.plot(times, PE_aux[:,k] - PE_ref_aux[k], '-x', color=cm_bwr(np.double(k-2)/kmax), label='k='+str(k))
        ax2.plot(times, dPdk[:,k]-dPdk_ref[:,k], '-x', color=cm_bwr(np.double(k-2)/kmax))
    ax1.legend(loc=1, fontsize=8)
    ax3.plot(np.arange(1,kmax), np.amax(np.abs(dPdk-dPdk_ref), axis=0)[1:], label='max over time')
    ax3.plot(np.arange(1,kmax), np.average(np.abs(dPdk-dPdk_ref), axis=0)[1:], label='mean over time')
    ax3.legend(fontsize=8)
    ax3.legend(loc='best', fontsize=8)
    ax3.grid()
    ax0.set_title('PE')
    ax1.set_title('PE_aux - PE_aux(t=0)')
    ax2.set_title('dPdk-dPdk_ref')
    ax3.set_title('max_t(dPdk-dPdk_ref)')
    ax0.set_xlabel('time')
    ax1.set_xlabel('time')
    ax2.set_xlabel('time')
    ax3.set_xlabel('k')
    ax0.set_ylabel('PE [J]')
    ax1.set_ylabel('PE [J]')
    ax2.set_ylabel('PE [J]')
    ax3.set_ylabel('PE [J]')

    plt.tight_layout()
    plt.savefig(os.path.join(path_figs, 'PE_test_abs.png'))
    plt.close(fig)


    ''' output '''
    rootgrp = nc.Dataset(os.path.join(path_stats, filename_out), 'r+', format='NETCDF4')
    ts_grp = rootgrp.groups['timeseries']
    var = ts_grp.variables['PE_abs']
    var[:] = PE[:]
    rootgrp.close()

    return




# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

def plot_KE_PE_radav(id, file_data, path_data, path_figs):

    root = nc.Dataset(os.path.join(path_data, file_data), 'r')
    PE = root.groups['timeseries'].variables['PE'][:]
    KE = root.groups['timeseries'].variables['KE'][:]
    KE_sgs_stats = root.groups['timeseries'].variables['KE_sgs_stats'][:]
    KE_sgs = root.groups['timeseries'].variables['KE_sgs'][:]
    KE_r = root.groups['profiles'].variables['KE_r'][:,:]
    times = root.groups['dimensions'].variables['time'][:]
    radius = root.groups['dimensions'].variables['radius'][:]
    (nt, nr) = KE_r.shape
    root.close()
    print path_figs


    fig, axis = plt.subplots(1, 4, figsize=(15,5))
    ax1 = axis[0]
    ax2 = axis[1]
    ax3 = axis[2]
    ax4 = axis[3]
    ax1.set_title('PE')
    ax2.set_title('KE')
    ax3.set_title('KE, PE')
    ax4.set_title('KE + PE')

    ax2.plot(times, KE, '-o')
    ax1.plot(times, PE, '-o')
    ax3.plot(times, PE, '-o', label='PE')
    ax3.plot(times, KE, '-o', label='KE')
    ax3.legend(loc='best')
    ax4.plot(times, PE+KE, label='KE+PE')
    ax4.plot(times, PE, label='PE')
    ax4.legend(loc='best')

    ax1.set_xlabel('time [s]')
    ax2.set_xlabel('time [s]')
    ax3.set_xlabel('time [s]')
    ax4.set_xlabel('time [s]')
    ax1.set_ylabel('PE [J]')
    ax2.set_ylabel('KE [J]')
    ax4.set_ylabel('KE + PE [J]')

    ax1.grid()
    # ax1.grid(which='minor', linewidth=0.2)#, linestyle='-', linewidth=2)
    ax2.grid()

    plt.suptitle('Kinetic and Potential Energy (kmax='+str(kmax)+')')
    plt.tight_layout()
    plt.savefig(os.path.join(path_figs, 'KE_PE_'+id + '_radav.png'))
    plt.close(fig)




    fig, axis = plt.subplots(1, 4, figsize=(15, 5))
    ax1 = axis[0]
    ax2 = axis[1]
    ax3 = axis[2]
    ax4 = axis[3]
    ax1.set_title('PE')
    ax2.set_title('KE')
    ax3.set_title('KE, PE')
    ax4.set_title('KE + PE')

    ax2.semilogy(times, KE, '-o')
    ax1.semilogy(times, PE, '-o')
    ax3.semilogy(times, PE, '-o', label='PE')
    ax3.semilogy(times, KE, '-o', label='KE')
    ax3.legend(loc='best')
    ax4.semilogy(times, PE + KE, label='KE+PE')
    ax4.semilogy(times, PE, label='PE')
    ax4.legend(loc='best')

    ax1.set_xlabel('time [s]')
    ax2.set_xlabel('time [s]')
    ax3.set_xlabel('time [s]')
    ax4.set_xlabel('time [s]')
    ax1.set_ylabel('PE [J]')
    ax2.set_ylabel('KE [J]')
    ax4.set_ylabel('KE + PE [J]')

    ax1.grid()
    # ax1.grid(which='minor', linewidth=0.2)#, linestyle='-', linewidth=2)
    ax2.grid()

    plt.suptitle('Kinetic and Potential Energy (kmax=' + str(kmax) + ')')
    plt.tight_layout()
    plt.savefig(os.path.join(path_figs, 'KE_PE_' + id + '_radav_log.png'))
    plt.close(fig)

    return




def plot_KE_PE(id, file_data, path_data, path_figs):

    root = nc.Dataset(os.path.join(path_data, file_data), 'r')
    PE = root.groups['timeseries'].variables['PE'][:]
    PE = root.groups['timeseries'].variables['PE_3D'][:]
    KE = root.groups['timeseries'].variables['KE'][:]
    KE = root.groups['timeseries'].variables['KE_3D'][:]
    KE_sgs_stats = root.groups['timeseries'].variables['KE_sgs_stats'][:]
    KE_sgs = root.groups['timeseries'].variables['KE_sgs'][:]
    KE_r = root.groups['profiles'].variables['KE_r'][:,:]
    times = root.groups['dimensions'].variables['time'][:]
    radius = root.groups['dimensions'].variables['radius'][:]
    (nt, nr) = KE_r.shape
    root.close()
    print path_figs


    fig, axis = plt.subplots(1, 4, figsize=(15,5))
    ax1 = axis[0]
    ax2 = axis[1]
    ax3 = axis[2]
    ax4 = axis[3]
    ax1.set_title('PE')
    ax2.set_title('KE')
    ax3.set_title('KE, PE')
    ax4.set_title('KE + PE')

    ax2.plot(times, KE, '-o')
    ax1.plot(times, PE, '-o')
    ax3.plot(times, PE, '-o', label='PE')
    ax3.plot(times, KE, '-o', label='KE')
    ax3.legend(loc='best')
    ax4.plot(times, PE+KE, label='KE+PE')
    ax4.plot(times, PE, label='PE')
    ax4.legend(loc='best')

    ax1.set_xlabel('time [s]')
    ax2.set_xlabel('time [s]')
    ax3.set_xlabel('time [s]')
    ax4.set_xlabel('time [s]')
    ax1.set_ylabel('PE [J]')
    ax2.set_ylabel('KE [J]')
    ax4.set_ylabel('KE + PE [J]')

    ax1.grid()
    # ax1.grid(which='minor', linewidth=0.2)#, linestyle='-', linewidth=2)
    ax2.grid()

    plt.suptitle('Kinetic and Potential Energy (kmax='+str(kmax)+')')
    plt.tight_layout()
    plt.savefig(os.path.join(path_figs, 'KE_PE_'+id + '.png'))
    plt.close(fig)




    fig, axis = plt.subplots(1, 4, figsize=(15, 5))
    ax1 = axis[0]
    ax2 = axis[1]
    ax3 = axis[2]
    ax4 = axis[3]
    ax1.set_title('PE')
    ax2.set_title('KE')
    ax3.set_title('KE, PE')
    ax4.set_title('KE + PE')

    # ax2.semilogy(times, KE, '-o')
    # ax1.semilogy(times, PE, '-o')
    # ax3.semilogy(times, PE, '-o', label='PE')
    # ax3.semilogy(times, KE, '-o', label='KE')
    # ax3.legend(loc='best')
    # ax4.semilogy(times, PE + KE, label='KE+PE')
    # ax4.semilogy(times, PE, label='PE')
    # ax4.legend(loc='best')

    ax1.set_xlabel('time [s]')
    ax2.set_xlabel('time [s]')
    ax3.set_xlabel('time [s]')
    ax4.set_xlabel('time [s]')
    ax1.set_ylabel('PE [J]')
    ax2.set_ylabel('KE [J]')
    ax4.set_ylabel('KE + PE [J]')

    ax1.grid()
    # ax1.grid(which='minor', linewidth=0.2)#, linestyle='-', linewidth=2)
    ax2.grid()

    plt.suptitle('Kinetic and Potential Energy (kmax=' + str(kmax) + ')')
    plt.tight_layout()
    plt.savefig(os.path.join(path_figs, 'KE_PE_' + id + '_log.png'))
    plt.close(fig)

    return



def plot_KE(id, file_data, path_data, path_figs):
    root = nc.Dataset(os.path.join(path_data, file_data), 'r')
    PE = root.groups['timeseries'].variables['PE'][:]
    KE = root.groups['timeseries'].variables['KE'][:]
    KE_3D = root.groups['timeseries'].variables['KE_3D'][:]
    KE_sgs_stats = root.groups['timeseries'].variables['KE_sgs_stats'][:]
    KE_sgs = root.groups['timeseries'].variables['KE_sgs'][:]
    KE_r = root.groups['profiles'].variables['KE_r'][:, :]
    KE_x = root.groups['profiles'].variables['KE_x'][:, :]
    times = root.groups['dimensions'].variables['time'][:]
    radius = root.groups['dimensions'].variables['radius'][:]
    (nt, nr) = KE_r.shape
    root.close()


    fig, axis = plt.subplots(2, 3, figsize=(20, 12))
    ax1 = axis[0,0]
    ax2 = axis[0,1]
    ax3 = axis[0,2]
    ax1.set_title('KE(r)')
    ax2.set_title('KE')
    ax3.set_title('KE_3D / KE')

    cf = ax1.contourf(times, radius, KE_r.T)
    # cbar = plt.colorbar(cf, shrink=0.6, ticks=np.arange(min, max + 1, 1), aspect=12)
    cbar = plt.colorbar(cf, ax=ax1, shrink=0.6, aspect=12)
    ax2.plot(times, KE, label='KE radial')
    ax2.plot(times, KE_3D, label='KE 3D')
    ax3.plot(times[1:], KE_3D[1:]/KE[1:])
    # ax3.plot(times, KE_3D/KE)

    ax2.legend(loc=0)

    ax1.set_xlabel('time [s]')
    ax2.set_xlabel('time [s]')
    ax3.set_xlabel('time [s]')
    ax1.set_ylabel('radius [m]')
    ax3.set_ylabel('KE [J]')

    ax1 = axis[1, 0]
    ax2 = axis[1, 1]
    ax3 = axis[1, 2]
    ax1.set_title('KE sgs')
    ax2.set_title('KE 3D')
    ax1.plot(times, KE_sgs, label='KE sgs')
    ax1.plot(times, KE_sgs_stats, label='KE sgs (stats)')
    ax1.legend(loc=1)
    ax2.plot(times, KE_3D, label='KE 3D')

    ax1.set_xlabel('time [s]')
    ax2.set_xlabel('time [s]')
    ax3.set_xlabel('time [s]')
    ax1.set_ylabel('KE [J]')
    ax3.set_ylabel('KE [J]')

    plt.suptitle('Kinetic Energy')
    plt.tight_layout()
    plt.savefig(os.path.join(path_figs, 'KE_' + id + '.png'))
    plt.close(fig)

    return



def plot_PE(id, filename_out, path_stats, path_figs):

    return


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

def set_input_output_parameters(args, filename_in):
    print('')
    print('--- set input parameters ---')
    global case_name
    global path
    global files

    if args.casename:
        case_name = args.casename
    else:
        case_name = 'ColdPoolDry_triple_3D'
    path = args.path

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
    print('nx, ny, nz', nx, ny, nz)

    ''' determine file range '''
    global kmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = np.int(100)
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = np.int(10000)
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = np.int(20)

    # azimuthally averaged data
    file_radav = nc.Dataset(os.path.join(path, 'data_analysis', filename_in))
    times_ = file_radav.groups['timeseries'].variables['time'][:]
    # radius = file_radav.groups['stats'].variables['r'][:]
    krange = file_radav.groups['dimensions'].variables['krange'][:]
    file_radav.close()
    nk = krange[-1]
    kmax = np.int(np.minimum(kmax, nk))
    del krange

    # path_fields = os.path.join(path, 'fields')
    # times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
    #          and np.int(name[:-3]) >= tmin and np.int(name[:-3]) <= tmax]
    # times.sort()
    print('input radav file: ' + filename_in)
    times = [np.int(t) for t in times_ if np.int(t) >= tmin and np.int(t) <= tmax]
    files = [str(np.int(t)) + '.nc' for t in times]
    print('times', times)
    print('kmax', kmax)
    print ('')
    del times_

    return nml, times




def define_geometry(case_name, nml):
    print('--- define geometry ---')
    global x_half, y_half, z_half
    global ic_arr, jc_arr
    global rstar, irstar, zstar, kstar

    # test file:
    path_fields = os.path.join(path, 'fields')
    var = read_in_netcdf_fields('u', os.path.join(path_fields, files[0]))
    [nx_, ny_, nz_] = var.shape
    del var
    x_half = np.empty((nx_), dtype=np.double, order='c')
    y_half = np.empty((ny_), dtype=np.double, order='c')
    z_half = np.empty((nz_), dtype=np.double, order='c')
    count = 0
    for i in xrange(nx_):
        x_half[count] = (i + 0.5) * dx[0]
        count += 1
    count = 0
    for j in xrange(ny_):
        y_half[count] = (j + 0.5) * dx[1]
        count += 1
    count = 0
    for i in xrange(nz_):
        z_half[count] = (i + 0.5) * dx[2]
        count += 1

    # set coordinates for plots
    if case_name == 'ColdPoolDry_single_3D':
        rstar = nml['init']['r']
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
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
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        kstar = np.int(np.round(zstar / dx[2]))
        isep = 4 * irstar
        jsep = 0
        ic1 = np.int(np.round((nx + 2 * gw) / 3)) - gw
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
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        kstar = np.int(np.round(zstar / dx[2]))
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

    return
# ----------------------------------
def theta_s(s):
    # T_tilde = 298.15
    # thetas_c(s, qt){
    #     return T_tilde * exp((s - (1.0 - qt) * sd_tilde - qt * sv_tilde) / cpm_c(qt));
    # }
    # cpm_c(qt){
    #     return (1.0 - qt) * cpd + qt * cpv;
    # }
    # parameters from pycles
    T_tilde = 298.15
    sd_tilde = 6864.8
    cpd = 1004.0
    th_s = T_tilde * np.exp( (s - sd_tilde)/cpd )
    return th_s


def compute_rho(p0, T, qt, qv):
    Rd = 287.1
    eps_vi = 1.60745384883
    return p0 / (Rd * T) * 1./(1.0 - qt + eps_vi * qv)


def compute_alpha(p0, T, qt, qv):
    Rd = 287.1
    eps_vi = 1.60745384883
    return (Rd * T)/p0 * (1.0 - qt + eps_vi * qv)


# ----------------------------------
def create_output_file(times, filename_in, path_in, filename_out, path_out):
    # output for each CP:
    # - min, max (timeseries)
    # - CP height (field; max=timeseries)
    # - (ok) CP rim (field)
    nt = len(times)
    print ''
    print('--- create output file ---')
    print(path_out)
    print(times, nt)

    root_in = nc.Dataset(os.path.join(path_in, filename_in), 'r')
    radius = root_in.groups['stats'].variables['r'][:]
    root_in.close()
    nr = len(radius)

    rootgrp = nc.Dataset(os.path.join(path_out, filename_out), 'w', format='NETCDF4')

    dim_grp = rootgrp.createGroup('dimensions')
    dim_grp.createDimension('nt', nt)
    dim_grp.createDimension('nr', nr)
    dim_grp.createDimension('nx', nx)
    var = dim_grp.createVariable('time', 'f8', ('nt'))
    var.units = "s"
    var[:] = times
    var = dim_grp.createVariable('radius', 'f8', ('nr'))
    var.units = "m"
    var[:] = radius

    ts_grp = rootgrp.createGroup('timeseries')
    ts_grp.createDimension('nt', nt)
    # KE from radially averaged velocity fields
    var = ts_grp.createVariable('KE', 'f8', ('nt'))
    var.units = "J"
    var.title = "KE from azimuthally averaged velocity fields"
    var = ts_grp.createVariable('KEd', 'f8', ('nt'))
    var.units = "(m/s)^2"
    var = ts_grp.createVariable('KE_3D', 'f8', ('nt'))
    var.units = "(m/s)^2"
    var.title = "KE from 3D velocity fields"
    var = ts_grp.createVariable('KE_3D_int', 'f8', ('nt'))
    var.units = "(m/s)^2"
    var.title = "KE from interpolated 3D velocity fields"
    # SGS KE
    #   KE computed from 3D fields and mean viscosity profile
    var = ts_grp.createVariable('KE_sgs_stats', 'f8', ('nt'))
    var.units = "(m/s)^2"
    #   KE from SGS <ww> profile in stats output
    var = ts_grp.createVariable('KE_sgs', 'f8', ('nt'))
    var.title = "SGS KE"
    var.units = "(m/s)^2"

    # PE from radially averaged velocity fields
    var = ts_grp.createVariable('PE', 'f8', ('nt'))
    var.units = "J"
    var.title = "PE from azimuthally averaged entropy field"
    var = ts_grp.createVariable('PEd', 'f8', ('nt'))
    var.units = "(m/s)^2"
    # PE from 3D fields
    var = ts_grp.createVariable('PE_3D', 'f8', ('nt'))
    var.units = "J"
    var.title = "PE from 3D entropy field"
    # absolute PE from 3D fields
    var = ts_grp.createVariable('PE_abs', 'f8', ('nt'))
    var.units = "J"
    var.title = "absolute PE from 3D temperature field"

    prof_grp = rootgrp.createGroup('profiles')
    prof_grp.createDimension('nt', nt)
    prof_grp.createDimension('nr', nr)
    prof_grp.createDimension('nx', nx)
    var = prof_grp.createVariable('KE_r', 'f8', ('nt', 'nr'))
    var.units = "J"
    var = prof_grp.createVariable('KE_x', 'f8', ('nt', 'nx'))
    var.units = "J"

    # field_grp = rootgrp.createGroup('fields_2D')
    # field_grp.createDimension('nt', nt)
    # field_grp.createDimension('nx', nx)
    # field_grp.createDimension('ny', ny)
    # var = field_grp.createVariable('w_max', 'f8', ('nt', 'nx', 'ny'))
    # var.units = "m/s"
    # var = field_grp.createVariable('w_max_height', 'f8', ('nt', 'nx', 'ny'))
    # var.units = "m"

    rootgrp.close()

    del radius
    del times
    return


# def dump_output_file(filename, KE, KEd, it):
#     rootgrp = nc.Dataset(os.path.join(path_out, filename), 'r+', format='NETCDF4')
#
#     ts_grp = rootgrp.groups['timeseries']
#     var = ts_grp.variables['KE']
#     var[it] = KE[it]
#     var = ts_grp.variables['KEd']
#     var[it] = KEd[it]
#
#     # field_grp = rootgrp.groups['fields_2D']
#     # var = field_grp.variables['w_max']
#     # var[it,:,:] = w_max[0,:,:]
#     # var = field_grp.variables['w_max_height']
#     # var[it,:,:] = w_max[1,:,:]
#
#     rootgrp.close()
#     return



def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:,:,:]
    rootgrp.close()
    return data


if __name__ == '__main__':
    main()

