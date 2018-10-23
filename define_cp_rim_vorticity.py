import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os

from define_cp_rim_plottingfct import set_colorbars

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--casename")
    parser.add_argument("--path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    parser.add_argument("--perc")
    args = parser.parse_args()

    global path_fields, path_out, path_stats
    if args.path:
        path = args.path
    else:
        # path = '/Users/bettinameyer/polybox/ClimatePhysics/Copenhagen/Projects/LES_ColdPool/' \
        #        'triple_3D_noise/Out_CPDry_triple_dTh2K/'
        path = '/nbi/ac/cond1/meyerbe/ColdPools/triple_3D_noise/Out_CPDry_triple_Th10K/'
    if os.path.exists(os.path.join(path, 'fields')):
        path_fields = os.path.join(path, 'fields')
    elif os.path.exists(os.path.join(path, 'fields_k120')):
        path_fields = os.path.join(path, 'fields_k120')
    path_out = os.path.join(path, 'figs_cp_rim')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    path_stats = os.path.join(path, 'fields_cp_rim')
    if not os.path.exists(path_stats):
        os.mkdir(path_stats)

    global case_name
    if args.casename:
        case_name = args.casename
    else:
        case_name = 'ColdPoolDry_triple_3D'

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

    if args.kmin:
        kmin = np.int(args.kmin)
    else:
        kmin = 5
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = 5
    krange = np.arange(kmin, kmax + 1, 1)
    nk = len(krange)

    # percentile for threshold
    if args.perc:
        perc = args.perc
    else:
        # perc = 95     # tested for triple 3D, dTh=3K, t=400s
        perc = 98       # tested for triple 3D, dTh=10K, t=100-400s

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dy, dz, gw
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    gw = nml['grid']['gw']

    global cm_bwr, cm_grey, cm_vir
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_vir = plt.cm.get_cmap('jet')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    set_colorbars(cm_bwr, cm_vir, cm_grey)      # to set colorbars as global functions in define_cp_rim_plottingfct.py

    # define subdomain to scan
    nx_half, ny_half = define_geometry(case_name, nml)
    icshift = nx_half
    jcshift = ny_half

    # define general arrays
    dphi = 6  # angular resolution for averaging of radius
    n_phi = 360 / dphi
    angular_range = np.arange(0, 360, dphi)
    # - rim_intp_all = (phi(t,k,i_phi)[deg], phi(t,k,i_phi)[rad], r_out(t,k,i_phi)[m], r_int(t,k,i_phi)[m], D(t,k,i_phi)[m])
    #   (phi: angles at interval of 6 deg; r_out,int: outer,inner boundary of convergence zone; D: thickness of convergence zone)
    # - rim_vel = (phi(t,k,i_phi)[deg], phi(t,k,i_phi)[rad], r_out(t,k,i_phi)[m], U(t,k,i_phi)[m/s], dU(t,k,i_phi)[m/s**2])
    # - rim_vel_av = (r_av(t,k), U_av(t,k), dU_av/dt(t,k))
    rim_intp_all = np.zeros(shape=(5, nt, nk, n_phi), dtype=np.double)
    rim_vel = np.zeros(shape=(5, nt, nk, n_phi), dtype=np.double)
    rim_vel_av = np.zeros(shape=(3, nt, nk))
    for it,t0 in enumerate(timerange):
        for ik,k0 in enumerate(krange):
            rim_intp_all[0, it, ik, :] = angular_range
            rim_intp_all[1, it, ik, :] = np.pi * angular_range / 180

    # create statistics file
    stats_file_name = 'rimstats_vort_perc' + str(perc) + 'th.nc'
    create_statistics_file(stats_file_name, path_stats, angular_range, n_phi, nt, timerange, nk, krange)

    for it, t0 in enumerate(timerange):
        if it > 0:
            dt = t0-timerange[it-1]
        else:
            dt = t0
        print('----- time: ' + str(t0) + ' (dt=' + str(dt) + ') -----')

        # create mask file for time step t0
        mask_file_name = 'rimmask_vort_perc' + str(perc) + 'th' + '_t' + str(t0) + '.nc'
        create_mask_file(mask_file_name, path_stats, nk, krange, ic, jc, ishift, jshift)


        # reading in file; selecting subdomain
        w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0) + '.nc'))
        w_roll = np.roll(np.roll(w[:, :, :], ishift, axis=0), jshift, axis=1)
        w_ = w_roll[ic - nx_half + ishift:ic + nx_half + ishift, jc - ny_half + jshift:jc + ny_half + jshift, :]
        del w

        for ik,k0 in enumerate(krange):
            print('level: k=' + str(k0), '(z=' + str(k0 * dz) + 'm)')

            # ''' (b) find inner & outer rim '''
            rim_int = np.zeros((nx_, ny_), dtype=np.int)
            rim_out = np.zeros((nx_, ny_), dtype=np.int)

            # ''' dump mask and rim '''
            '''_____ '''
            mask = np.zeros((nx_, ny_))
            dump_mask(mask, rim_int, rim_out, mask_file_name, path_stats, k0, ik)
            '''_____ '''

            # dump statistics
            dump_statistics_file(rim_intp_all[:, it, ik, :], rim_vel[:, it, ik, :], rim_vel_av[:, it, ik],
                                 stats_file_name, path_stats, k0, ik, t0, it)

            print('')


    del w_roll


    return


# ----------------------------------------------------------------------

def create_statistics_file(file_name, path, angles, n_phi, nt, timerange, nk, krange):
    print('-------- create statistics file -------- ')
    print(path + ', ' + file_name)
    print('')

    rootgrp = nc.Dataset(os.path.join(path, file_name), 'w', format='NETCDF4')

    ts_grp = rootgrp.createGroup('timeseries')
    ts_grp.createDimension('nt', nt)
    ts_grp.createDimension('nz', nk)
    var = ts_grp.createVariable('time', 'f8', ('nt'))
    var.units = "s"
    var[:] = timerange
    var = ts_grp.createVariable('r_av', 'f8', ('nt', 'nz'))
    var.units = "m"
    var = ts_grp.createVariable('U_av', 'f8', ('nt', 'nz'))
    var.units = "m/s"
    var = ts_grp.createVariable('dU_av', 'f8', ('nt', 'nz'))
    var.units = "m/s^2"

    stats_grp = rootgrp.createGroup('stats')
    stats_grp.createDimension('nt', nt)
    stats_grp.createDimension('nz', nk)
    stats_grp.createDimension('nphi', n_phi)  # number of segments
    var = stats_grp.createVariable('r_out', 'f8', ('nt', 'nz', 'nphi'))
    var.units = "m"
    var = stats_grp.createVariable('r_int', 'f8', ('nt', 'nz', 'nphi'))
    var.units = "m"
    var = stats_grp.createVariable('D', 'f8', ('nt', 'nz', 'nphi'))  # thickness of rim
    var.units = "m"
    var = stats_grp.createVariable('U', 'f8', ('nt', 'nz', 'nphi'))  # velocity of outer rim
    var.units = "m/s"
    var = stats_grp.createVariable('dU', 'f8',('nt', 'nz', 'nphi'))  # change of rim velocity with time ('acceleration of rim')
    var.units = "m/s^2"

    pol_grp = rootgrp.createGroup('angles')
    pol_grp.createDimension('nphi', n_phi)
    var = pol_grp.createVariable('phi_deg', 'f8', ('nphi'))
    var.units = "degrees"
    var[:] = angles[:]
    var = pol_grp.createVariable('phi_rad', 'f8', ('nphi'))
    var.units = "radian"
    var[:] = np.pi / 180 * angles[:]

    prof_grp = rootgrp.createGroup('profiles')
    prof_grp.createDimension('nz', nk)
    prof_grp.createDimension('nt', nt)
    var = prof_grp.createVariable('k_dumped', 'f8', ('nt', 'nz'))
    var[:] = np.zeros((nt,nk))
    var.description = "0: level not dumped, 1: level dumped"
    var = prof_grp.createVariable('krange', 'f8', ('nz'))
    var[:] = krange[:]
    var.units = "-"
    var = prof_grp.createVariable('zrange', 'f8', ('nz'))
    var[:] = krange[:] * dz
    var.units = "m"

    rootgrp.close()

    return



def dump_statistics_file(rim_intp_all, rim_vel, rim_vel_av, file_name, path, k0, ik, t0, it):
    # statistics are dumped after each (t,k)-tuple (i.e. for each level)

    # - rim_intp_all = (phi(i_phi)[deg], phi(i_phi)[rad], r_out(i_phi)[m], r_int(i_phi)[m], D(i_phi)[m])
    #   (phi: angles at interval of 6 deg; r_out,int: outer,inner boundary of convergence zone; D: thickness of convergence zone)
    # - rim_vel = (phi(i_phi)[deg], phi(i_phi)[rad], r_out(i_phi)[m], U(i_phi)[m/s], dU(i_phi)[m/s**2])
    # - rim_vel_av = (r_av, U_av, dU_av/dt)
    rootgrp = nc.Dataset(os.path.join(path, file_name), 'r+', format='NETCDF4')

    ts_grp = rootgrp.groups['timeseries']
    var = ts_grp.variables['r_av']
    var[it,ik] = rim_vel_av[0]
    var = ts_grp.variables['U_av']
    var[it,ik] = rim_vel_av[1]
    var = ts_grp.variables['dU_av']
    var[it,ik] = rim_vel_av[2]

    stats_grp = rootgrp.groups['stats']
    var = stats_grp.variables['r_out']
    var[it,ik,:] = rim_intp_all[2,:]
    var = stats_grp.variables['r_int']
    var[it,ik,:] = rim_intp_all[3,:]
    var = stats_grp.variables['D']
    var[it,ik,:] = rim_intp_all[4,:]
    var = stats_grp.variables['U']
    var[it,ik,:] = rim_vel[3,:]
    var = stats_grp.variables['dU']
    var[it,ik,:] = rim_vel[4,:]

    prof_grp = rootgrp.groups['profiles']
    var = prof_grp.variables['k_dumped']
    var[ik] = 1

    rootgrp.close()
    return





def create_mask_file(file_name, path, nk, krange, ic, jc, ishift, jshift):
    print('-------- create mask file --------')
    rootgrp = nc.Dataset(os.path.join(path, file_name), 'w', format='NETCDF4')

    descr_grp = rootgrp.createGroup('description')
    var = descr_grp.createVariable('ic', 'f8', )
    var[:] = ic
    var = descr_grp.createVariable('jc', 'f8', )
    var[:] = jc
    var = descr_grp.createVariable('ishift', 'f8', )
    var[:] = ishift
    var = descr_grp.createVariable('jshift', 'f8', )
    var[:] = jshift

    mask_grp = rootgrp.createGroup('fields')
    mask_grp.createDimension('nx', nx_)
    mask_grp.createDimension('ny', ny_)
    mask_grp.createDimension('nz', nk)
    mask_grp.createVariable('mask', 'f8', ('nx', 'ny', 'nz'))
    mask_grp.createVariable('rim_inner', 'f8', ('nx', 'ny', 'nz'))
    mask_grp.createVariable('rim_outer', 'f8', ('nx', 'ny', 'nz'))

    prof_grp = rootgrp.createGroup('profiles')
    prof_grp.createDimension('nz', nk)
    var = prof_grp.createVariable('k_dumped', 'f8', ('nz'))
    var[:] = np.zeros(nk)
    var.description = "0: level not dumped, 1: level dumped"
    var = prof_grp.createVariable('krange', 'f8', ('nz'))
    var[:] = krange[:]
    var.units = "-"
    var = prof_grp.createVariable('zrange', 'f8', ('nz'))
    var[:] = krange[:] * dz
    var.units = "m"
    rootgrp.close()
    return


def dump_mask(mask, rim_int, rim_out, file_name, path, k0, ik):
    # mask = (nx_, ny_)
    # rim_int = (nx_, ny_)
    # rim_out = (nx_, ny_)
    # k0: level
    print('dumping mask, k0=', k0)
    rootgrp = nc.Dataset(os.path.join(path, file_name), 'r+', format='NETCDF4')

    mask_grp = rootgrp.groups['fields']
    var = mask_grp.variables['mask']
    var[:,:,ik] = mask[:,:]
    var = mask_grp.variables['rim_inner']
    var[:,:,ik] = rim_int[:,:]
    var = mask_grp.variables['rim_outer']
    var[:,:,ik] = rim_out[:,:]

    levels_grp = rootgrp.groups['profiles']
    var = levels_grp.variables['k_dumped']
    var[ik] = np.int(k0)

    rootgrp.close()
    print('')
    return


# ----------------------------------------------------------------------
def define_geometry(case_name, nml):
    global nx_, ny_
    global ic, jc, shift, ishift, jshift
    if case_name == 'ColdPoolDry_triple_3D':
        flag = 'triple'
        # d = np.int(np.round(ny / 2))
        d = np.int(np.round((ny + gw) / 2))
        a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        ic = np.int(np.round(a / 2))
        jc = np.int(np.round(d / 2))
        shift = 60
        nx_half = irstar + shift
        ny_half = irstar + shift
        ishift = np.max(nx_half - ic, 0)
        jshift = np.max(ny_half - jc, 0)
        nx_ = 2 * nx_half
        ny_ = 2 * ny_half
    elif case_name == 'ColdPoolDry_double_3D':
        flag = 'double'
        rstar = 5000.0
        irstar = np.int(np.round(rstar / dx))
        isep = 4 * irstar
        jsep = 0
        ic1 = np.int(nx / 3)
        jc1 = np.int(ny / 2)
        # ic2 = ic1 + isep
        # jc2 = jc1 + jsep
        ic = ic1
        jc = jc1
        shift = 40
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
    return nx_half, ny_half

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