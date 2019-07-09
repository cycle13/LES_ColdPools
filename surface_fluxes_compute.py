import netCDF4 as nc
import argparse
import json as simplejson
import os
import numpy as np
import sys
import matplotlib.pyplot as plt

from thermodynamic_functions import exner_c, pv_c, sd_c, g, alpha_c

def main():

    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()


    times, nml = set_input_output_parameters(args)

    compute_surface_fluxes(times)


    print os.path.join(path, 'stats', 'Stats.'+case_name + '.nc')
    stats_file = nc.Dataset(os.path.join(path, 'stats', 'Stats.'+case_name + '.nc'))
    flux_file_stats = nc.Dataset(os.path.join(path_out_fields, 'surface_fluxes_stats.nc'))
    shf_mean_stats = stats_file.groups['timeseries'].variables['shf_surface_mean'][:]
    s_flux_mean_stats = stats_file.groups['timeseries'].variables['s_flux_surface_mean'][:]
    shf_mean = flux_file_stats.variables['shf_surface_mean'][:]
    s_flux_mean = flux_file_stats.variables['s_flux_surface_mean'][:]
    times_out = flux_file_stats.variables['time'][:]
    stats_file.close()
    flux_file_stats.close()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.plot(shf_mean_stats, '-o', label='shf (stats)', linewidth=3)
    ax1.plot(shf_mean, '-d', label='shf', linewidth=2)
    ax2.plot(s_flux_mean_stats, '-o', label='s-flux (stats)', linewidth=3)
    ax2.plot(s_flux_mean, '-d', label='s-flux', linewidth=3)
    ax1.legend()
    plt.show()

    return

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def compute_surface_fluxes(times):

    # x_half = np.empty((nx), dtype=np.double, order='c')
    # y_half = np.empty((ny), dtype=np.double, order='c')
    z_half = np.empty((nz), dtype=np.double, order='c')
    # count = 0
    # for i in xrange(nx):
    #     x_half[count] = (i + 0.5) * dx[0]
    #     count += 1
    # count = 0
    # for j in xrange(ny):
    #     y_half[count] = (j + 0.5) * dx[1]
    #     count += 1
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
    filename = 'surface_fluxes.nc'
    create_output_file_2d(filename, nx, ny, times)
    filename_stats = 'surface_fluxes_stats.nc'
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
        shf_surface_mean = 0.0

        # WINDSPEED
        gustiness = 0.0
        u0 = 0.0
        v0 = 0.0
        windspeed = compute_windspeed(u, v, u0, v0, gustiness)

        # COMPUTE FLUXES
        for i in range(1,nx - 1):
            for j in range(1,ny - 1):
                theta_flux[i,j] = -CH * windspeed[i,j] * (temperature[i,j] * exner_c(Pg - theta_surface))
                s_flux[i,j] = -CH * windspeed[i,j] * (s[i,j] - s_star)
                # self.qt_flux[ij] = -self.cq * windspeed[ij] * (PV.values[qt_shift + ijk] - Ref.qtg)
                # buoyancy_flux = g * ((theta_flux + (eps_vi-1.0)*(theta_surface*self.qt_flux[ij] + Ref.qtg * theta_flux))/(theta_surface*(1.0 + (eps_vi-1)*Ref.qtg)))
                buoyancy_flux = g * theta_flux / theta_surface
                u_flux[i,j] = -CM * interp_2(windspeed[i,j], windspeed[i+1,j]) * (u[i,j] + u0)
                v_flux[i,j] = -CM * interp_2(windspeed[i,j], windspeed[i,j + 1]) * (v[i,j] + v0)
                ustar_ = np.sqrt(CM) * windspeed[i,j]
                friction_velocity[i,j] = ustar_

                shf[i,j] = s_flux[i,j] * rho0 * temperature[i,j]

                s_flux_surface_mean += s_flux[i,j]
                shf_surface_mean += shf[i,j]

        s_flux_surface_mean /= (nx-2)*(ny-2)
        shf_surface_mean /= (nx-2)*(ny-2)
        dump_2d(filename, theta_flux, s_flux, u_flux, v_flux, shf, it, times)
        dump_stats(filename_stats, s_flux_surface_mean, shf_surface_mean, it, times)

    return

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

def create_output_file_2d(filename, nx, ny, times):
    if os.path.exists(os.path.join(path_out_fields, filename)):
        print('2d File already existing! Not creating new file')
    else:
        rootgrp = nc.Dataset(os.path.join(path_out_fields, filename), 'w', format='NETCDF4')
        nt = len(times)
        rootgrp.createDimension('nt', nt)
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
    if os.path.exists(os.path.join(path_out_fields, filename)):
        print('Stats file already existing! Not creating new file')
    else:
        rootgrp = nc.Dataset(os.path.join(path_out_fields, filename), 'w', format='NETCDF4')
        nt = len(times)
        rootgrp.createDimension('nt', nt)
        var = rootgrp.createVariable('time', 'f8', ('nt'))
        var.units = "s"
        var[:] = times

        var = rootgrp.createVariable('s_flux_surface_mean', 'f8', ('nt'))
        var.units = "J/(K m2)"
        var = rootgrp.createVariable('shf_surface_mean', 'f8', ('nt'))
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
    t_out = rootgrp.variables['time'][:]
    if times[it] != t_out[it]:
        print('Not correct output time!!')
        sys.exit()

    var = rootgrp.variables['s_flux_surface_mean']
    var[it] = s_flux_surface_mean
    var = rootgrp.variables['shf_surface_mean']
    var[it] = shf_surface_mean
    rootgrp.close()
    return

# ----------------------------------------------------------------------
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
# ----------------------------------------------------------------------

def set_input_output_parameters(args):
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

    print('nx, ny, nz', nx, ny, nz)
    print('times', timerange)


    return timerange, nml

# _______________________________________________________
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