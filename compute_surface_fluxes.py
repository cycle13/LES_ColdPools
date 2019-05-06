import netCDF4 as nc
import argparse
import json as simplejson
import os
import numpy as np
import sys

from thermodynamic_functions import exner_c, pv_c, sd_c, g

def main():

    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()


    times, nml = set_input_output_parameters(args)

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

    # READ IN FIELDS
    t0 = times[0]
    k0 = 0
    root = nc.Dataset(os.path.join(path_fields, str(t0)+'.nc'))
    u = root.groups['fields'].variables['u'][:,:,k0]
    v = root.groups['fields'].variables['v'][:,:,k0]
    temperature = root.groups['fields'].variables['temperature'][:,:,k0]
    s = root.groups['fields'].variables['s'][:,:,k0]
    root.close()

    # OUTPUT FIELDS
    s_flux = np.zeros((nx, ny), dtype=np.double)
    u_flux = np.zeros((nx, ny), dtype=np.double)
    v_flux = np.zeros((nx, ny), dtype=np.double)
    friction_velocity = np.zeros((nx, ny), dtype=np.double)

    # COEFFICIENTS
    cm = 0.001229  # bulk coefficient for momentum flux (from Rico-case)
    ch = 0.001094  # bulk coefficient for heat flux (from Rico-case)
    # self.cq = 0.001133
    z0 = 0.00015
    CM = cm * (np.log(20.0 / z0) / np.log(z_half[0] / z0)) ** 2
    CH = ch * (np.log(20.0 / z0) / np.log(z_half[0] / z0)) ** 2

    # WINDSPEED
    gustiness = 0.0
    u0 = 0.0
    v0 = 0.0
    windspeed = compute_windspeed(u, v, u0, v0, gustiness)

    # COMPUTE FLUXES
    # theta_surface = 298.0
    Pg = 1.0e5
    Tg = 300.0
    qtg = 0.0
    # theta_surface = Ref.Tg * exner_c(Ref.Pg)
    theta_surface = Tg * exner_c(Pg)
    pv_star = pv_c(Pg, qtg, qtg)
    pd_star = Pg - pv_star
    # self.s_star = (1.0-Ref.qtg) * sd_c(pd_star, Ref.Tg) + Ref.qtg * sv_c(pv_star,Ref.Tg)
    s_star = sd_c(pd_star, Tg)
    for i in range(nx - 1):
        for j in range(ny - 1):
            # ijk = i * istride + j * jstride + gw
            # ij = i * istride_2d + j
            theta_flux = -CH * windspeed[i,j] * (temperature[i,j] * exner_c(Pg - theta_surface))

            s_flux[i,j] = -CH * windspeed[i,j] * (s[i,j] - s_star)
            # self.qt_flux[ij] = -self.cq * windspeed[ij] * (PV.values[qt_shift + ijk] - Ref.qtg)
            # buoyancy_flux = g * ((theta_flux + (eps_vi-1.0)*(theta_surface*self.qt_flux[ij] + Ref.qtg * theta_flux))/(theta_surface*(1.0 + (eps_vi-1)*Ref.qtg)))
            buoyancy_flux = g * theta_flux / theta_surface
            u_flux[i,j] = -CM * interp_2(windspeed[i,j], windspeed[i+1,j]) * (u[i,j] + u0)
            v_flux[i,j] = -CM * interp_2(windspeed[i,j], windspeed[i,j + 1]) * (v[i,j] + v0)
            ustar_ = np.sqrt(CM) * windspeed[i,j]
            friction_velocity[i,j] = ustar_

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