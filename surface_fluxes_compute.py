import netCDF4 as nc
import argparse
import json as simplejson
import os
import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
import time


from thermodynamic_functions import exner_c, pv_c, sd_c, g, alpha_c

label_size = 10
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.labelsize'] = 15

def main():

    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--path_tracers")
    args = parser.parse_args()



    times, nml = set_input_output_parameters(args)
    ic_arr, jc_arr = define_geometry(nml)
    ic = ic_arr[0]
    jc = jc_arr[0]
    irange = np.minimum(nx - ic, ic)
    jrange = np.minimum(ny - jc, jc)
    rmax = np.int(np.ceil(np.sqrt(irange ** 2 + jrange ** 2)))


    ''' compute surface fluxes from 3D LES output'''
    # # > output 2d-fields and timeseries of mean values
    # filename_2d = 'surface_fluxes.nc'
    # filename_stats = 'surface_fluxes_stats.nc'
    # compute_surface_fluxes(filename_2d, filename_stats, times)
    # print('')

    ''' compute azimuthally averaged fluxes '''
    # # - read in r_field, th_field from fields_v_rad/v_rad.nc
    # # - compute average
    # filename_2d = 'surface_fluxes.nc'
    # filename_stats = 'surface_fluxes_stats.nc'
    # compute_angular_average(filename_2d, filename_stats, rmax)


    ''' compute surface flux from given v_rad profile'''
    # compute_surface_flux_constant()
    # ''' (a) radius of CP from tracers '''
    # (1) read in tracer position = CP rim = position of max(v_rad)
    # (2) assume profile v_rad(r) = m*r, r<=r_tracer
    # (3) compute tracer profile therefrom
    # if args.path_tracers:
    #     path_tracers = os.path.join(path, args.path_tracers, 'output')
    # else:
    #     k0 = 0
    #     path_tracers = os.path.join(path, 'tracer_k' + str(k0), 'output')
    # print('path tracers: ', path_tracers)
    # cp_id = 1  # circle ID that is used for statistics
    # n_tracers = get_number_tracers(path_tracers)
    # n_cps = get_number_cps(path_tracers)
    # print('number of CPs: ', n_cps)
    # print('number of tracers per CP: ', n_tracers)
    # dist_tracers_av = np.zeros((nt))
    # U_tracers_av = np.zeros((nt))
    # for it, t0 in enumerate(times):
    #     dist_tracers_av[it], U_tracers_av[it] = get_radius_vel(path_tracers, it, cp_id, n_tracers, n_cps)
    # r_tracers_av = dist_tracers_av * dx[0]
    # del dist_tracers_av



    ''' (a) compute radius of CP defined as inflection point of radial velocity '''
    iR_CP, R_CP, v_CP = compute_CP_radius(times)

    ''' compute sfc fluxes from linear u_rad model (Romps) '''
    k0 = 0
    u_rad_lin = compute_linear_profiles_Romps(iR_CP, R_CP, v_CP, k0, times)

    shf_rad_lin, shf_mean_lin = compute_surface_flux_1d(u_rad_lin, times, 'test_shf_1d_lin.png')
    root = nc.Dataset(os.path.join(path, 'data_analysis', 'stats_radial_averaged.nc'))
    times_rav = root.groups['timeseries'].variables['time'][:]
    r_range = root.groups['stats'].variables['r'][:]
    u_rad = root.groups['stats'].variables['v_rad'][:,:,0]      # (nt, nr, nz)
    root.close()
    shf_rad, shf_mean = compute_surface_flux_1d(u_rad, times, 'test_shf_1d_normal.png')






    ''' plotting '''

    try:
        try:
            stats_file = nc.Dataset(os.path.join(path, 'stats', 'Stats.' + case_name + '.all.nc'))
        except:
            stats_file = nc.Dataset(os.path.join(path, 'stats', 'Stats.' + case_name + '.nc'))
    except:
        try:
            stats_file = nc.Dataset(os.path.join(path, 'Stats.' + case_name + '.all.nc'))
        except:
            stats_file = nc.Dataset(os.path.join(path, 'Stats.' + case_name + '.nc'))
    filename_fluxes = 'surface_fluxes_stats.nc'
    figname = 'SHF_comparison_stats_vs_computation.png'
    plot_sfc_fluxes(figname, filename_fluxes, stats_file)

    filename_stats = 'surface_fluxes_stats.nc'
    figname = 'fluxes_rav.png'
    plot_sfc_rad_av(figname, filename_fluxes, stats_file)
    stats_file.close()

    figname = 'fluxes_rad_av_comparison.png'
    rmax_plot = 9e3
    fig, (axis) = plt.subplots(1, 3, figsize=(16, 5))
    ax1 = axis[0]
    ax2 = axis[1]
    ax3 = axis[2]
    tstep = 6
    ta = np.where(times_rav == tmin)[0][0]
    tb = np.where(times_rav == tmax)[0][0]
    for it, t0 in enumerate(times_rav[0::tstep]):
        if it >= ta and it <= tb:
            count_color = tstep * np.double(it) / len(times_rav)
            ax1.plot(r_range, shf_rad_lin[it,:],
                     color=mpl.cm.winter(count_color), label='t=' + str(np.int(t0)))
            ax2.plot(r_range, shf_rad[it,:],
                     color=mpl.cm.winter(count_color), label='t=' + str(np.int(t0)))
    ax3.plot(times_rav[:nt], shf_mean_lin, label='SHF lin')
    ax3.plot(times_rav[:nt], shf_mean, label='SHF')
    ax1.set_xlim(0, rmax_plot)
    ax2.set_xlim(0, rmax_plot)
    ax3.set_xlim(0, times_rav[-1])
    ax1.set_xlabel('r  [m]')
    ax2.set_xlabel('r  [m]')
    ax3.set_xlabel('time  [s]')
    ax1.set_title('SHF linear v_rad')
    ax2.set_title('SHF from velocity output')
    ax3.set_title('mean SHF')
    for i,ax in enumerate(axis):
        ax.legend(loc='best')
    plt.subplots_adjust(bottom=0.1, right=.98, left=0.06, top=0.9, wspace=0.2)
    fig.savefig(os.path.join(path, path_out_figs, figname))
    plt.close(fig)

    return

# ----------------------------------------------------------------------

def compute_surface_flux_1d(u_field_1d, times, figname):
    print('--- compute 1D surface flux ---')
    # COEFFICIENTS
    cm = 0.001229  # bulk coefficient for momentum flux (from Rico-case)
    ch = 0.001094  # bulk coefficient for heat flux (from Rico-case)
    # self.cq = 0.001133
    z0 = 0.00015
    z_half0 = 0.5*dx[2]
    CM = cm * (np.log(20.0 / z0) / np.log(z_half0 / z0)) ** 2
    CH = ch * (np.log(20.0 / z0) / np.log(z_half0 / z0)) ** 2

    # REFERENCE VALUES
    Pg = 1.0e5
    Tg = 300.0
    qtg = 0.0
    theta_surface = Tg * exner_c(Pg)
    pv_star = pv_c(Pg, qtg, qtg)
    pd_star = Pg - pv_star
    s_star = sd_c(pd_star, Tg)

    alpha0 = alpha_c(Pg, Tg, qtg, 0.0)
    rho0 = 1. / alpha0

    # LOOP OVER TIME
    # create_output_file_2d(filename_2d, nx, ny, times)
    # create_output_file_stats(filename_stats, times)
    # READ IN FIELDS
    root = nc.Dataset(os.path.join(path, 'data_analysis', 'stats_radial_averaged.nc'))
    nr = root.groups['stats'].dimensions['nr'].size
    s_rad = root.groups['stats'].variables['s'][:,:,0]
    temp_rad = root.groups['stats'].variables['temperature'][:,:,0]
    root.close()

    # OUTPUT FIELDS
    theta_flux = np.zeros((nt,nr), dtype=np.double)
    s_flux = np.zeros((nt, nr), dtype=np.double)
    ur_flux = np.zeros((nt, nr), dtype=np.double)
    shf = np.zeros((nt, nr), dtype=np.double)
    s_flux_surface_mean = np.zeros((nt), dtype=np.double)
    shf_mean = np.zeros((nt), dtype=np.double)
    for it, t0 in enumerate(times):
        # WINDSPEED
        gustiness = 0.0
        u0 = 0.0
        v0 = 0.0
        windspeed = np.zeros((nr), dtype=np.double)
        for i in range(nr):
            windspeed[i] = np.maximum(np.sqrt(u0**2 + v0**2) + u_field_1d[it,i], gustiness)

        # COMPUTE FLUXES
        theta_flux[it,:] = -CH * windspeed[:] * (
                    temp_rad[it,:] * exner_c(Pg) - theta_surface)  # no difference to in-loop calculation
        for i in range(nr-1):
            s_flux[it,i] = -CH * windspeed[i] * (s_rad[it,i] - s_star)
            ur_flux[it,i] = -CM * interp_2(windspeed[i], windspeed[i + 1]) * (u_field_1d[it,i] + u0)

            shf[it,i] = s_flux[it,i] * rho0 * temp_rad[it,i]
            s_flux_surface_mean[it] += s_flux[it,i]
            shf_mean[it] += shf[it,i]

        # DOMAIN MEAN FLUXES
        s_flux_surface_mean[it] /= (nr - 1)
        shf_mean[it] /= (nr - 1)
        # dump_2d(filename_2d, theta_flux, s_flux, u_flux, v_flux, shf, it, times)
        # dump_stats(filename_stats, s_flux_surface_mean, shf_mean, it, times)




    # -------------------------------
    root = nc.Dataset(os.path.join(path, 'data_analysis', 'stats_radial_averaged.nc'))
    nr = root.groups['stats'].dimensions['nr'].size
    r_range = root.groups['stats'].variables['r'][:]
    root.close()
    fig, (axis) = plt.subplots(1, 3, figsize=(16, 5))
    ax1 = axis[0]
    ax2 = axis[1]
    ax3 = axis[2]
    tstep = 6
    # ta = np.where(times_rav == tmin)[0][0]
    # tb = np.where(times_rav == tmax)[0][0]
    rmax_plot = 9e3
    for it in range(nt):
        count_color = np.double(it) / nt
        ax1.plot(r_range, s_flux[it, :],
                 color=mpl.cm.winter(count_color), label='t=' + str(np.int(it)))
        ax2.plot(r_range, shf[it, :],
                 color=mpl.cm.winter(count_color), label='t=' + str(np.int(it)))
    ax3.plot(shf_mean, label='SHF')
    ax1.set_xlim(0, rmax_plot)
    ax2.set_xlim(0, rmax_plot)
    # ax3.set_xlim(0, times_rav[-1])
    ax1.set_title('s-flux')
    ax2.set_title('SHF')
    ax3.set_title('mean SHF')
    for i, ax in enumerate(axis):
        ax.legend(loc='best', fontsize=6, ncol=2)
    plt.subplots_adjust(bottom=0.1, right=.98, left=0.06, top=0.9, wspace=0.2)
    fig.savefig(os.path.join(path, path_out_figs, figname))
    plt.close(fig)


    return shf, shf_mean




def compute_surface_fluxes(filename_2d, filename_stats, times):

    z_half = np.empty((nz), dtype=np.double, order='c')
    count = 0
    for i in xrange(nz):
        z_half[count] = (i + 0.5) * dx[2]
        count += 1
    k0 = 0

    # COEFFICIENTS
    cm = 0.001229  # bulk coefficient for momentum flux (from Rico-case)
    ch = 0.001094  # bulk coefficient for heat flux (from Rico-case)
    # self.cq = 0.001133
    z0 = 0.00015
    CM = cm * (np.log(20.0 / z0) / np.log(z_half[0] / z0)) ** 2
    CH = ch * (np.log(20.0 / z0) / np.log(z_half[0] / z0)) ** 2

    # REFERENCE VALUES
    Pg = 1.0e5
    Tg = 300.0
    qtg = 0.0
    # theta_surface = Ref.Tg * exner_c(Ref.Pg)
    theta_surface = Tg * exner_c(Pg)
    pv_star = pv_c(Pg, qtg, qtg)
    pd_star = Pg - pv_star
    # self.s_star = (1.0-Ref.qtg) * sd_c(pd_star, Ref.Tg) + Ref.qtg * sv_c(pv_star,Ref.Tg)
    s_star = sd_c(pd_star, Tg)

    # temperature_half[k], ql_half[k], qi_half[k] = Thermodynamics.eos(p_half_[k], self.sg, self.qtg)
    # alpha_half[k] = Thermodynamics.alpha(p_half_[k], temperature_half[k], self.qtg, qv_half[k])
    alpha0 = alpha_c(Pg, Tg, qtg, 0.0)
    rho0 = 1. / alpha0


    # LOOP OVER TIME
    create_output_file_2d(filename_2d, nx, ny, times)
    create_output_file_stats(filename_stats, times)
    for it, t0 in enumerate(times):
        print('it, t0', it, t0)
        # READ IN FIELDS
        root = nc.Dataset(os.path.join(path_fields, str(t0)+'.nc'))
        u = root.groups['fields'].variables['u'][:,:,k0]
        v = root.groups['fields'].variables['v'][:,:,k0]
        temperature = root.groups['fields'].variables['temperature'][:,:,k0]
        s = root.groups['fields'].variables['s'][:,:,k0]
        root.close()

        # OUTPUT FIELDS
        theta_flux = np.zeros((nx, ny), dtype=np.double)
        s_flux = np.zeros((nx, ny), dtype=np.double)
        u_flux = np.zeros((nx, ny), dtype=np.double)
        v_flux = np.zeros((nx, ny), dtype=np.double)
        friction_velocity = np.zeros((nx, ny), dtype=np.double)
        shf = np.zeros((nx, ny), dtype=np.double)

        s_flux_surface_mean = 0.0
        shf_mean = 0.0

        # WINDSPEED
        gustiness = 0.0
        u0 = 0.0
        v0 = 0.0
        windspeed = compute_windspeed(u, v, u0, v0, gustiness)

        # COMPUTE FLUXES
        t_ini = time.time()
        theta_flux = -CH * windspeed * (temperature * exner_c(Pg) - theta_surface)  # no difference to in-loop calculation
        buoyancy_flux = g * theta_flux / theta_surface
        for i in range(1,nx-1):
            for j in range(1,ny-1):
                # theta_flux[i,j] = -CH * windspeed[i,j] * (temperature[i,j] * exner_c(Pg) - theta_surface)
                s_flux[i,j] = -CH * windspeed[i,j] * (s[i,j] - s_star)
                # self.qt_flux[ij] = -self.cq * windspeed[ij] * (PV.values[qt_shift + ijk] - Ref.qtg)
                # # buoyancy_flux = g * ((theta_flux + (eps_vi-1.0)*(theta_surface*self.qt_flux[ij] + Ref.qtg * theta_flux))/(theta_surface*(1.0 + (eps_vi-1)*Ref.qtg)))
                # buoyancy_flux = g * theta_flux[i,j] / theta_surface
                u_flux[i,j] = -CM * interp_2(windspeed[i,j], windspeed[i+1,j]) * (u[i,j] + u0)
                v_flux[i,j] = -CM * interp_2(windspeed[i,j], windspeed[i,j + 1]) * (v[i,j] + v0)
                ustar_ = np.sqrt(CM) * windspeed[i,j]
                friction_velocity[i,j] = ustar_

                shf[i,j] = s_flux[i,j] * rho0 * temperature[i,j]
                s_flux_surface_mean += s_flux[i,j]
                shf_mean += shf[i,j]

        # DOMAIN MEAN FLUXES
        s_flux_surface_mean /= (nx-2) * (ny-2)
        shf_mean /= (nx-2) * (ny-2)
        print('--- time: ', time.time()-t_ini)
        dump_2d(filename_2d, theta_flux, s_flux, u_flux, v_flux, shf, it, times)
        dump_stats(filename_stats, s_flux_surface_mean, shf_mean, it, times)





        # # COMPUTE FLUXES without loop
        # filename_new = 'surface_fluxes_new.nc'
        # create_output_file_2d(filename_new, nx, ny, times)
        # filename_stats_new = 'surface_fluxes_stats_new.nc'
        # create_output_file_stats(filename_stats_new, times)
        # t_ini = time.time()
        # theta_flux = -CH * windspeed * (temperature * exner_c(Pg) - theta_surface)
        # s_flux = -CH * windspeed * (s - s_star)
        # shf = s_flux * rho0 * temperature
        # # self.qt_flux[ij] = -self.cq * windspeed[ij] * (PV.values[qt_shift + ijk] - Ref.qtg)
        # # buoyancy_flux = g * ((theta_flux + (eps_vi-1.0)*(theta_surface*self.qt_flux[ij] + Ref.qtg * theta_flux))/(theta_surface*(1.0 + (eps_vi-1)*Ref.qtg)))
        # buoyancy_flux = g * theta_flux / theta_surface
        # for i in range(1, nx - 1):
        #     for j in range(1, ny - 1):
        #         u_flux[i, j] = -CM * interp_2(windspeed[i, j], windspeed[i + 1, j]) * (u[i, j] + u0)
        #         v_flux[i, j] = -CM * interp_2(windspeed[i, j], windspeed[i, j + 1]) * (v[i, j] + v0)
        #         ustar_ = np.sqrt(CM) * windspeed[i, j]
        #         friction_velocity[i, j] = ustar_
        # # compute domain mean fluxes
        # s_flux_surface_mean = np.mean(s_flux)
        # shf_mean = np.mean(shf)
        # print('--- time: ', time.time()-t_ini)
        # dump_2d(filename_new, theta_flux, s_flux, u_flux, v_flux, shf, it, times)
        # dump_stats(filename_stats_new, s_flux_surface_mean, shf_mean, it, times)

    return

# ----------------------------------------------------------------------

def compute_CP_radius(times):
    filename_stats_radav = 'stats_radial_averaged.nc'
    rootgrp = nc.Dataset(os.path.join(path, 'data_analysis', filename_stats_radav), 'r')
    v_rad_av = rootgrp.groups['stats'].variables['v_rad'][:, :, :]  # nt, nr, nz
    r_range = rootgrp.groups['stats'].variables['r'][:]  # nr
    nr = rootgrp.groups['stats'].dimensions['nr'].size
    nk = rootgrp.groups['stats'].dimensions['nz'].size
    times_rav = rootgrp.groups['timeseries'].variables['time'][:]
    rootgrp.close()
    print(nt, nr, nk)
    print('resolution: ', dx, r_range[:5])

    dv_rad_dr = np.zeros((nt, nr, nk), dtype=np.double)
    iR = np.zeros((nt, nk), dtype=np.int)
    v_CP = np.zeros((nt, nk), dtype=np.double)


    for it, t0 in enumerate(times):
        for ir, r0 in enumerate(r_range[1:-1]):
            dv_rad_dr[it, ir, :] = 1. / dx[0] * (v_rad_av[it, ir + 1, :] - v_rad_av[it, ir - 1, :])
        iR[it, :] = np.argmin(dv_rad_dr[it, 10:, :], axis=0) + 10
        for ik in range(nk):
            v_CP[it,ik] = v_rad_av[it,iR[it,ik],ik]

    R = r_range[iR]
    print('--', v_CP.shape)

    figname = 'inflection_point_test.png'
    rmax_plot = 9e3
    fig, (axis) = plt.subplots(1, 2, figsize=(16, 5))
    ax1 = axis[0]
    ax1_ = ax1.twinx()
    ax2 = axis[1]
    tstep = 6
    ta = np.where(times_rav == tmin)[0][0]
    tb = np.where(times_rav == tmax)[0][0]
    for it, t0 in enumerate(times_rav[0::tstep]):
        if it >= ta and it <= tb:
            count_color = tstep * np.double(it) / len(times_rav)
            ax1.plot(r_range, v_rad_av[tstep * it, :, 0],
                     color=mpl.cm.winter(count_color), label='t=' + str(np.int(t0)))
            ax1_.plot(r_range, dv_rad_dr[tstep * it, :, 0],
                      color=mpl.cm.summer(count_color), linewidth=1)
            ax2.plot(r_range, dv_rad_dr[tstep * it, :, 0],
                     color=mpl.cm.summer(count_color), label='t=' + str(np.int(t0)))
            ax2.plot(R[tstep * it,0], dv_rad_dr[tstep * it, iR[tstep * it,0], 0], 'o',
                     color=mpl.cm.summer(count_color))


            ax1.plot(r_range[iR[tstep * it, 0]], v_rad_av[tstep * it, iR[tstep * it, 0], 0], 'o',
                     color=mpl.cm.bone(count_color))
            ax1.plot(R[tstep * it, 0], v_CP[tstep * it, 0], 'xk')

    ax1.set_xlim(0, rmax_plot)
    ax2.set_xlim(0, rmax_plot)
    ax1.set_xlabel('r  [m]')
    ax2.set_xlabel('r  [m]')
    ax1.set_ylabel('v_rad')
    ax1_.set_ylabel('d(v_rad)/dr')
    ax2.set_ylabel('d(v_rad)/dr')
    ax1.set_title('v_rad')
    ax2.set_title('d(v_rad)/dr')
    for i, ax in enumerate(axis):
        ax.legend(loc='best')
    plt.subplots_adjust(bottom=0.1, right=.98, left=0.06, top=0.9, wspace=0.4)
    fig.savefig(os.path.join(path, path_out_figs, figname))
    plt.close(fig)
    return iR, R, v_CP




def compute_linear_profiles_Romps(iR_CP, R_CP, v_CP, k0, times):
    ''' linear model '''
    filename_stats_radav = 'stats_radial_averaged.nc'
    rootgrp = nc.Dataset(os.path.join(path, 'data_analysis', filename_stats_radav), 'r')
    v_rad_av = rootgrp.groups['stats'].variables['v_rad'][:, :, :]  # nt, nr, nz
    r_range = rootgrp.groups['stats'].variables['r'][:]  # nr
    nr = rootgrp.groups['stats'].dimensions['nr'].size
    nk = rootgrp.groups['stats'].dimensions['nz'].size
    times_rav = rootgrp.groups['timeseries'].variables['time'][:]
    rootgrp.close()
    print(nt, nr, nk)
    ur_lin = np.zeros((nt, nr))
    # ir_tracer = np.zeros(nt, dtype=np.int)
    for it, t0 in enumerate(times):
    #     ir_tracer[it] = np.where(r_range == np.int(np.round(R_CP[it], -2)))[0][0]
        ur_lin[it, :] = v_CP[it,k0] / R_CP[it,k0] * r_range[:]

    # print('')
    rmax_plot = 10e3
    figname = 'linear_model_Romps.png'
    fig, axis = plt.subplots(2, 2, figsize=(15, 8))
    ax1 = axis[0, 0]
    ax2 = axis[0, 1]
    ax3 = axis[1, 0]
    ax4 = axis[1, 1]
    ax1.plot(times, R_CP[:,k0], '-o')
    ax2.plot(times, v_CP[:,k0], '-o')
    for it, t0 in enumerate(times):
        count = np.double(it) / len(times)
        ax1.plot(t0, R_CP[it,k0], 'o', color=mpl.cm.bwr(count))
        ax2.plot(t0, v_CP[it,k0], 'o', color=mpl.cm.bwr(count))
        ax3.plot(r_range, ur_lin[it,:], '-', color=mpl.cm.bwr(count))
        ax4.plot(r_range, v_rad_av[it,:,k0], '-', color=mpl.cm.bwr(count))
        ax4.plot(r_range[iR_CP[it,k0]], v_rad_av[it, iR_CP[it,k0],k0], 'o', color=mpl.cm.bwr(count))
        ax4.plot(r_range, ur_lin[it, :], '-', color=mpl.cm.bwr(count), linewidth=1)
        ax4.plot(r_range[iR_CP[it,k0]], v_CP[it,k0], 'd', color=mpl.cm.bwr(count))
    ax1.set_title('R(t)')
    ax2.set_title('U(t)')
    ax3.set_title('u_r=U/R*r (Romps)')
    for i, ax in enumerate(axis[0, :]):
        ax.set_xlabel('time  [s]')
    for i, ax in enumerate(axis[1, :]):
        ax.set_xlabel('r  [m]')
    ax1.set_ylabel('r  [m]')
    ax2.set_ylabel('U  [m/s]')
    ax3.set_ylabel('u_lin  [m/s]')
    ax4.set_ylabel('v_rad  [m/s]')
    ax3.set_xlim(0, rmax_plot)
    ax4.set_xlim(0, rmax_plot)
    ax3.set_ylim(0, np.amax(v_rad_av[:, :, k0]))
    ax4.set_ylim(0, np.amax(v_rad_av[:, :, k0]))
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, wspace=0.3, hspace=0.3)
    fig.savefig(os.path.join(path, path_out_figs, figname))
    plt.close(fig)

    rmax_plot = 10e3
    t_list = [0, 7, 13, 18, 24, 30, 36]
    t_list = [0, 1, 2, 7, 13, 18, 24, 30, 36]
    figname = 'linear_model_Romps_large.png'
    fig, axis = plt.subplots(2, 1, figsize=(15, 8))
    plt.subplots_adjust(bottom=0.25, right=.95, left=0.07, top=0.9, wspace=0.3)
    ax1 = axis[0]
    ax2 = axis[1]
    for it in t_list:
        if it < nt:
            t0 = times[it]
            count = np.double(it) / len(times)
            ax1.plot(r_range, ur_lin[it, :], '-', color=mpl.cm.winter(count))
            ax2.plot(r_range, v_rad_av[it, :, k0], '-', color=mpl.cm.winter(count), label='t=' + str(t0))
            ax2.plot(r_range[iR_CP[it,k0]], v_rad_av[it, iR_CP[it,k0], k0], 'o', color=mpl.cm.winter(count))
            ax2.plot(r_range[:iR_CP[it,k0] + 1], ur_lin[it, :iR_CP[it,k0] + 1], '-', color=mpl.cm.winter(count), linewidth=1)
            ax2.plot(r_range[iR_CP[it,k0]], v_CP[it,k0], 'd', color=mpl.cm.winter(count))
    ax2.set_title('u_r=U/R*r (Romps)')
    ax2.legend(loc='best')
    ax1.set_xlabel('r  [m]')
    ax1.set_xlim(0, rmax_plot)
    ax2.set_xlim(0, rmax_plot)
    ax1.set_ylim(0, np.amax(v_rad_av[:, :, 0]))
    ax2.set_ylim(0, np.amax(v_rad_av[:, :, 0]))
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, wspace=0.3, hspace=0.3)
    fig.savefig(os.path.join(path, path_out_figs, figname))
    plt.close(fig)

    return ur_lin

# ----------------------------------------------------------------------

def compute_angular_average(filename_2d, filename_stats, rmax):
    # (1) read in 2d-surface flux fields
    print(os.path.join(path, 'fields_fluxes', filename_2d))
    sfc_flux_file = nc.Dataset(os.path.join(path, 'fields_fluxes', filename_2d), 'r')

    # (2) read in r_field
    vrad_file = nc.Dataset(os.path.join(path, 'fields_v_rad', 'v_rad.nc'), 'r')
    r_field = vrad_file.variables['r_field'][:,:]
    vrad_file.close()

    # (3) compute azimuthal average
    var_list = ['shf', 'u_flux', 'v_flux', 's_flux', 'th_flux']
    file_out = nc.Dataset(os.path.join(path, 'fields_fluxes', filename_stats), 'r+')
    nt_ = len(file_out.groups['timeseries'].variables['time'][:])
    r_grp = file_out.createGroup('rad_av')
    if not 'nt' in r_grp.dimensions.keys():
        r_grp.createDimension('nt', nt_)
    if not 'nr' in r_grp.dimensions.keys():
        r_grp.createDimension('nr', rmax)
    if not 'r' in r_grp.variables.keys():
        # ri_range = np.arange(0, rmax, dtype=np.int)
        # r_range = np.arange(0, rmax, dtype=np.double) * dx[0]
        var = r_grp.createVariable('r', 'f8', ('nr'))
        var.unit = 'm'
        var[:] = np.arange(0, rmax, dtype=np.double) * dx[0]
        var = r_grp.createVariable('ir', 'f8', ('nr'))
        var.unit = '-'
        var[:] = np.arange(0, rmax, dtype=np.int)

    for var_name in var_list:
        var_2d = sfc_flux_file.variables[var_name][:, :, :]
        count = np.zeros(rmax, dtype=np.int)
        var_av = np.zeros((nt, rmax), dtype=np.double)
        for i in range(nx):
            for j in range(ny):
                r = r_field[i, j]
                count[r] += 1
                var_av[:, r] += var_2d[:nt, i, j]
        for r in range(rmax):
            if count[r] > 0:
                var_av[:, r] /= count[r]

        if var_name in r_grp.variables.keys():
            var = r_grp.variables[var_name]
        else:
            var = r_grp.createVariable(var_name, 'f8', ('nt', 'nr'))
        var[:nt,:] = var_av[:,:]
    file_out.close()
    sfc_flux_file.close()

    return



# ----------------------------------------------------------------------
def plot_sfc_fluxes(figname, filename_fluxes, stats_file):
    flux_file_stats = nc.Dataset(os.path.join(path_out_fields, filename_fluxes))
    shf_mean_stats = stats_file.groups['timeseries'].variables['shf_surface_mean'][:]
    times_stats = stats_file.groups['timeseries'].variables['t'][:]
    s_flux_mean_stats = stats_file.groups['timeseries'].variables['s_flux_surface_mean'][:]
    shf_mean = flux_file_stats.groups['timeseries'].variables['shf_surface_mean'][:]
    s_flux_mean = flux_file_stats.groups['timeseries'].variables['s_flux_surface_mean'][:]
    times_out = flux_file_stats.groups['timeseries'].variables['time'][:]
    flux_file_stats.close()

    nt_ = np.minimum(len(times_out), len(times_stats))
    if times_stats[:nt].any() != times_out[:nt].any():
        print('not same times!! not plotting')
        return
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        ax1.plot(times_stats[:nt], shf_mean_stats[:nt], '-d', label='shf (stats)', linewidth=3)
        ax1.plot(times_stats[:nt], shf_mean[:nt], '--o', label='shf', linewidth=2)
        ax2.plot(times_stats[:nt], s_flux_mean_stats[:nt], '-d', label='s-flux (stats)', linewidth=3)
        ax2.plot(times_stats[:nt], s_flux_mean[:nt], '-o', label='s-flux', linewidth=3)
        ax1.set_xlim(0,times_stats[nt_-1])
        ax2.set_xlim(0,times_stats[nt_-1])
        ax1.set_xlabel('time  [s]')
        ax2.set_xlabel('time  [s]')
        ax1.set_ylabel('mean(SHF)')
        ax2.set_ylabel('mean(entropy flux)')
        ax1.legend()
        plt.suptitle('comparison sfc fluxes computed vs. stats-file')
        plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, wspace=0.3)
        fig.savefig(os.path.join(path, path_out_figs, figname))
        plt.close(fig)
    return



def plot_sfc_rad_av(figname, filename_fluxes, stats_file):
    print('plotting surface fluxes rad av')
    shf_mean_stats = stats_file.groups['timeseries'].variables['shf_surface_mean'][:]
    times_stats = stats_file.groups['timeseries'].variables['t'][:]
    s_flux_mean_stats = stats_file.groups['timeseries'].variables['s_flux_surface_mean'][:]

    flux_file_stats = nc.Dataset(os.path.join(path_out_fields, filename_fluxes))
    r_range = flux_file_stats.groups['rad_av'].variables['r'][:]
    shf_rav = flux_file_stats.groups['rad_av'].variables['shf'][:,:]
    s_flux_rav = flux_file_stats.groups['rad_av'].variables['s_flux'][:,:]
    times_rav = flux_file_stats.groups['timeseries'].variables['time'][:]
    shf_mean = flux_file_stats.groups['timeseries'].variables['shf_surface_mean'][:]
    flux_file_stats.close()

    file_radav = nc.Dataset(os.path.join(path, 'data_analysis', 'stats_radial_averaged.nc'))
    r_rad = file_radav.groups['stats'].variables['r'][:]
    v_rad = file_radav.groups['stats'].variables['v_rad'][:,:,0]
    temp = file_radav.groups['stats'].variables['temperature'][:,:,0]
    file_radav.close()

    nt_ = np.minimum(len(times_rav), len(times_stats))
    ta = np.where(times_rav == tmin)[0][0]
    tb = np.where(times_rav == tmax)[0][0]
    rmax_plot = 10e3
    irmax_plot = np.int(rmax_plot/dx[0])
    tstep = 3
    if times_stats[:nt].any() != times_rav[:nt].any():
        print('not same times!! not plotting')
        return
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        for it, t0 in enumerate(times_rav[0::tstep]):
            if it >= ta and it <= tb:
                count_color = tstep*np.double(it) / len(times_rav)
                ax1.plot(r_range, shf_rav[tstep*it,:], color=mpl.cm.jet(count_color),
                        label='t=' + str(np.int(t0)))
                ax2.plot(r_range, s_flux_rav[tstep*it,:], color=mpl.cm.jet(count_color),
                         label='t=' + str(np.int(t0)))
                if it == ta:
                    ax1.plot(t0, 1e1*shf_mean_stats[tstep*it], 'ok', markersize=8, label='domain mean*1e1')
                    ax2.plot(t0, 1e1*s_flux_mean_stats[tstep*it], 'ok', markersize=8, label='domain mean*1e1')
                else:
                    ax1.plot(t0, 1e1*shf_mean_stats[tstep*it], 'ok', markersize=8)
                    ax2.plot(t0, 1e1*s_flux_mean_stats[tstep*it], 'ok', markersize=8)

        ax1.set_xlim(0,rmax_plot)
        ax2.set_xlim(0,rmax_plot)
        ax1.set_xlabel('r  [m]')
        ax2.set_xlabel('r  [m]')
        ax1.set_ylabel('SHF(r)')
        ax2.set_ylabel('s_flux(r)')
        # ax1.legend()
        ax1.legend(loc='upper left', bbox_to_anchor=(.5, -0.15),fontsize=8, ncol=5)
        plt.suptitle('azimuthally averaged SHF')

        plt.subplots_adjust(bottom=0.25, right=.95, left=0.07, top=0.9, wspace=0.3)
        fig.savefig(os.path.join(path, path_out_figs, figname))
        plt.close(fig)


        fig, (axis) = plt.subplots(1, 2, figsize=(12, 5))
        axis[0].plot(r_range)
        axis[0].plot(r_rad, 'r--')
        plt.subplots_adjust(bottom=0.25, right=.95, left=0.07, top=0.9, wspace=0.3)
        fig.savefig(os.path.join(path, path_out_figs, 'test_radii.png'))
        plt.close(fig)

        figname = 'SHF_rav.png'
        fig, (axis) = plt.subplots(1, 2, figsize=(12, 5))
        ax1 = axis[0]
        ax2 = axis[1]
        tstep = 6
        shf_mean_ = np.zeros((nt_))
        for it, t0 in enumerate(times_rav[0::tstep]):
            if it >= ta and it <= tb:
                count_color = tstep*np.double(it) / len(times_rav)
                # ax3.plot(r_range, shf_rav[tstep*it,:], color=mpl.cm.winter(count_color),
                #          label='t=' + str(np.int(t0)))
                if np.amax(shf_rav[tstep*it,:])>0.:
                    ax2.semilogy(r_range, shf_rav[tstep*it,:], color=mpl.cm.winter(count_color),
                             label='t=' + str(np.int(t0)))
        for it,t0 in enumerate(times_rav):
            shf_mean_[it] = 1./(nx*ny*dx[0])*np.sum(2*np.pi*r_range[:]*shf_rav[it,:])
        ax1.plot(times_stats, shf_mean_stats, label='<SHF> (stats)')
        ax1.plot(times_stats[:tb], shf_mean[:tb], '-o', markersize=4, label='<SHF>')
        ax1.plot(times_stats[:tb], 5e-2*np.amax(v_rad[:tb], axis=1), '-o', markersize=4, label='0.05*max(v_r)')
        ax1.plot(times_stats[:tb], -1e-1*(np.amin(temp[:tb, :irmax_plot], axis=1) - np.amax(temp[:tb, :irmax_plot], axis=1)),
                 '-o', markersize=4, label='-0.1*dT')

        ax1.set_xlim(0, tmax)
        ax1.set_xlabel('time  [s]')
        ax1.set_ylabel('SHF [J/(s*m2)]  /  v [m/s]  /  dT [K]')
        ax2.set_xlim(0, rmax_plot)
        ax2.set_ylim(1e-4, 1e2)
        ax2.set_xlabel('r  [m]')
        ax2.set_ylabel('log[ SHF(r) ]')
        # ax3.set_xlim(0, rmax_plot)
        # ax3.set_xlabel('r  [m]')
        # ax3.set_ylabel('SHF(r)')
        # ax3.set_xlim(0, rmax_plot)
        # ax3.set_xlabel('r  [m]')
        # ax3.set_ylabel('SHF(r)')
        ax1.legend(loc='upper left', bbox_to_anchor=(0., -0.15), fontsize=10, ncol=2)
        ax2.legend(loc='upper left', bbox_to_anchor=(0., -0.15), fontsize=10, ncol=4)
        # plt.suptitle('azimuthally averaged SHF')

        plt.subplots_adjust(bottom=0.25, right=.98, left=0.06, top=0.9, wspace=0.2)
        fig.savefig(os.path.join(path, path_out_figs, figname))
        plt.close(fig)


    return


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
    # count = 0
    # # while CP age is 0 and CP ID is 1
    # while (int(lines[count].split()[0]) == 1):
    #     count += 1
    # cp_number = int(lines[count-1].split()[3])
    cp_number = int(lines[-1].split()[3])
    f.close()

    return cp_number



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

# ----------------------------------------------------------------------


def create_output_file_2d(filename, nx, ny, times):
    # if os.path.exists(os.path.join(path_out_fields, filename)):
    #     print('2d File already existing! Not creating new file')
    # else:
    rootgrp = nc.Dataset(os.path.join(path_out_fields, filename), 'w', format='NETCDF4')
    nt_ = len(times)
    rootgrp.createDimension('nt', nt_)
    var = rootgrp.createVariable('time', 'f8', ('nt'))
    var.units = "s"
    var[:] = times
    rootgrp.createDimension('nx', nx)
    rootgrp.createDimension('ny', ny)

    var = rootgrp.createVariable('th_flux', 'f8', ('nt', 'nx', 'ny'))
    var.units = "K/m2"
    var = rootgrp.createVariable('s_flux', 'f8', ('nt', 'nx', 'ny'))
    var.units = "J/(K m2)"
    var = rootgrp.createVariable('u_flux', 'f8', ('nt', 'nx', 'ny'))
    var.units = "1/(ms)"
    var = rootgrp.createVariable('v_flux', 'f8', ('nt', 'nx', 'ny'))
    var.units = "1/(ms)"
    var = rootgrp.createVariable('shf', 'f8', ('nt', 'nx', 'ny'))
    var.units = "J/m2"
    rootgrp.close()

    return


def create_output_file_stats(filename, times):
    # if os.path.exists(os.path.join(path_out_fields, filename)):
    #     print('Stats file already existing! Not creating new file')
    # else:
    rootgrp = nc.Dataset(os.path.join(path_out_fields, filename), 'w', format='NETCDF4')
    nt_ = len(times)
    ts_grp = rootgrp.createGroup('timeseries')
    ts_grp.createDimension('nt', nt_)
    var = ts_grp.createVariable('time', 'f8', ('nt'))
    var.units = "s"
    var[:] = times

    var = ts_grp.createVariable('s_flux_surface_mean', 'f8', ('nt'))
    var.units = "J/(K m2)"
    var = ts_grp.createVariable('shf_surface_mean', 'f8', ('nt'))
    var.units = "J/m2"
    rootgrp.close()

    return


def dump_2d(filename, th_flux, s_flux, u_flux, v_flux, shf, it, times):
    rootgrp = nc.Dataset(os.path.join(path_out_fields, filename), 'r+', format='NETCDF4')

    t_out = rootgrp.variables['time'][:]
    if times[it] != t_out[it]:
        print('Not correct output time!!')
        sys.exit()

    var = rootgrp.variables['th_flux']
    var[it,:,:] = th_flux[:,:]
    var = rootgrp.variables['s_flux']
    var[it,:,:] = s_flux[:,:]
    var = rootgrp.variables['u_flux']
    var[it,:,:] = u_flux[:,:]
    var = rootgrp.variables['v_flux']
    var[it,:,:] = v_flux[:,:]
    var = rootgrp.variables['shf']
    var[it,:,:] = shf[:,:]

    rootgrp.close()
    return


def dump_stats(filename, s_flux_surface_mean, shf_surface_mean, it, times):
    rootgrp = nc.Dataset(os.path.join(path_out_fields, filename), 'r+', format='NETCDF4')
    ts_grp = rootgrp.groups['timeseries']
    t_out = ts_grp.variables['time'][:]
    if times[it] != t_out[it]:
        print('Not correct output time!!')
        sys.exit()

    var = ts_grp.variables['s_flux_surface_mean']
    var[it] = s_flux_surface_mean
    var = ts_grp.variables['shf_surface_mean']
    var[it] = shf_surface_mean
    rootgrp.close()
    return

# ----------------------------------------------------------------------
def interp_2(phi, phip1):
    return 0.5*(phi + phip1)

def compute_windspeed(u, v, u0, v0, gustiness):
    imin = 1
    jmin = 1
    imax = nx
    jmax = ny

    speed = np.zeros((nx,ny), dtype=np.double)

    for i in range(imin, imax):
        for j in range(jmin, jmax):
            u_interp = interp_2(u[i-1,j],u[i,j]) + u0
            v_interp = interp_2(v[i,j-1],v[i,j]) + v0
            speed[i,j] = np.maximum(np.sqrt(u_interp*u_interp + v_interp*v_interp), gustiness)

    return speed


# ----------------------------------------------------------------------

def set_input_output_parameters(args):
    print('')
    print('--- set input parameters ---')
    global path, path_in, path_out_figs, path_out_fields, path_stats, path_fields
    path = args.path
    path_fields = os.path.join(path, 'fields')
    path_out_figs = os.path.join(path, 'figs_fluxes')
    path_out_fields = os.path.join(path, 'fields_fluxes')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    if not os.path.exists(path_out_fields):
        os.mkdir(path_out_fields)
    print('paths: ')
    print('path_out_fields: ' + path_out_fields)
    print('path_out_figs:   ' + path_out_figs)
    print('')

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
    global nt, tmin, tmax
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

    print('nx, ny, nz', nx, ny, nz)
    print('times', timerange)


    return timerange, nml


def define_geometry(nml):
    '''--- define geometry ---'''
    global rstar
    if case_name == 'ColdPoolDry_double_2D':
        rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        isep = 4 * irstar
        ic1 = np.int(nx / 3)
        ic2 = ic1 + isep
        jc1 = np.int(ny / 2)
        jc2 = jc1
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_single_3D':
        rstar = nml['init']['r']
        try:
            print('(ic,jc) from nml')
            ic = nml['init']['ic']
            jc = nml['init']['jc']
        except:
            print('(ic,jc) NOT from nml')
            ic = np.int(nx / 2)
            jc = np.int(ny / 2)
        ic_arr = [ic]
        jc_arr = [jc]
    elif case_name == 'ColdPoolDry_double_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
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
        d = np.int(np.round(ny / 2))
        a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        ic1 = np.int(np.round(a / 2))  # + gw
        ic2 = ic1
        ic3 = ic1 + np.int(np.round(a))
        jc1 = np.int(np.round(d / 2))  # + gw
        jc2 = jc1 + d
        jc3 = jc1 + np.int(np.round(d / 2))
        ic_arr = [ic1, ic2, ic3]
        jc_arr = [jc1, jc2, jc3]
    print('')
    return ic_arr, jc_arr


# _______________________________________________________
def theta_s(s):
    T_tilde = 298.15
    sd_tilde = 6864.8
    cpd = 1004.0
    th_s = T_tilde * np.exp( (s - sd_tilde)/cpd )
    return th_s
# _______________________________________________________

if __name__ == '__main__':
    main()