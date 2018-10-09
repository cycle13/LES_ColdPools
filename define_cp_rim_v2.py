import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os

from define_cp_rim_plottingfct import set_colorbars
from define_cp_rim_plottingfct import plot_yz_crosssection, plot_w_field, plot_s, \
    plot_outlines, plot_rim_mask, plot_angles, plot_cp_outline_alltimes, \
    plot_cp_rim_velocity, plot_cp_rim_averages, plot_rim_thickness

def main():
    """
    Find inner and outer rim of mask based on a threshold (usually 95th percentile) of the vertical velocity.
    The rim is found by number of neighbours and filling interior of mask

    :param --path: The full path to the files
    :param --casename: casename
    :param --tmin: minimum time taken into account
    :param --tmax: maximum time taken into account
    :param --k0: level, at which the mask is computed

    :return: figures


    Details:

    1. read in w-field, shift field (roll) and define partial domain where to look for cold pool
    2. mask 2D field and turn mask from boolean (True: w>w_c) into integer (1: w>w_c)
    3. Define rim of cold pool as the outline of the mask; based on number of neighbours
    """


    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--casename")
    parser.add_argument("--path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    # parser.add_argument("--k0", nargs = '+', type = int)
    parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    global path_fields, path_out
    if args.path:
        path = args.path
    else:
        path = '/Users/bettinameyer/polybox/ClimatePhysics/Copenhagen/Projects/LES_ColdPool/' \
                'triple_3D_noise/Out_CPDry_triple_dTh2K/'
        # path = '/nbi/ac/cond1/meyerbe/ColdPools/triple_3D_noise/Out_CPDry_triple_Th3K/'
    if os.path.exists(os.path.join(path, 'fields')):
        path_fields = os.path.join(path, 'fields')
    elif os.path.exists(os.path.join(path, 'fields_k120')):
        path_fields = os.path.join(path, 'fields_k120')
    path_out = os.path.join(path, 'figs_cp_rim')
    if not os.path.exists(path_out):
        os.mkdir(path_out)

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
    timerange = np.arange(tmin,tmax + 100,100)
    nt = len(timerange)

    # if args.k0:
    #     k0 = np.int(args.k0)
    # else:
    #     k0 = 5  # level
    if args.kmin:
        kmin = np.int(args.kmin)
    else:
        kmin = 5
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = 5
    krange = np.arange(kmin, kmax+1, 1)
    nk = len(krange)

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dy, dz
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    gw = nml['grid']['gw']

    global cm_bwr, cm_grey, cm_vir
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_vir = plt.cm.get_cmap('viridis')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    set_colorbars(cm_bwr, cm_vir, cm_grey)      # to set colorbars as global functions in define_cp_rim_plottingfct.py

    # define subdomain to scan
    global nx_, ny_
    if case_name == 'ColdPoolDry_triple_3D':
        flag = 'triple'
        # d = np.int(np.round(ny / 2))      # nbi
        d = np.int(np.round(ny / 2)) + gw
        a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        ic = np.int(np.round(a / 2))
        jc = np.int(np.round(d / 2))
        shift = 40
        id = irstar + shift
        jd = irstar + shift
        ishift = np.max(id - ic, 0)
        jshift = np.max(jd - jc, 0)
        nx_ = 2 * id
        ny_ = 2 * jd
    elif case_name == 'ColdPoolDry_double_3D':
        flag = 'double'
        rstar = 5000.0
        irstar = np.int(np.round(rstar / dx))
        isep = 4*irstar
        jsep = 0
        ic1 = np.int(nx / 3)
        jc1 = np.int(ny / 2)
        # ic2 = ic1 + isep
        # jc2 = jc1 + jsep
        ic = ic1
        jc = jc1
        shift = 40
        id = irstar + shift
        jd = irstar + shift
        ishift = np.max(id - ic, 0)
        jshift = np.max(jd - jc, 0)
        nx_ = 2 * id
        ny_ = 2 * jd

    print('ic,jc,id,jc,nx_,ny_', ic, jc, id, jd, nx_, ny_)
    # percentile for threshold
    perc = 95  # tested for triple 3D, t=400s

    # define general arrays
    dphi = 6        # angular resolution for averaging of radius
    n_phi = 360 / dphi
    # - rim_intp_all = (phi(t,i_phi)[deg], phi(t,i_phi)[rad], r_out(t,i_phi)[m], r_int(t,i_phi)[m], D(t,i_phi)[m])
    #   (phi: angles at interval of 6 deg; r_out,int: outer,inner boundary of convergence zone; D: thickness of convergence zone)
    # - rim_vel = (phi(t,i_phi)[deg], phi(t,i_phi)[rad], r_out(t,i_phi)[m], U(t,i_phi)[m/s], dU(t, i_phi)[m/s**2])
    # - rim_vel_av = (r_av(t), U_av(t), dU_av/dt(t))
    rim_intp_all = np.zeros(shape=(5, nt, n_phi), dtype=np.double)
    rim_vel = np.zeros(shape=(5, nt, n_phi), dtype=np.double)
    rim_vel_av = np.zeros(shape=(3, nt))

    # create statistics file
    stats_file_name = 'rimstats_perc' + str(perc) + 'th.nc'
    create_statistics_file(stats_file_name, path_out, n_phi, nt, timerange, nk, krange)

    for it,t0 in enumerate(timerange):
        if it > 0:
            dt = t0-timerange[it-1]
        else:
            dt = t0
        print('time: '+ str(t0), '(dt='+str(dt)+')')

        ''' create mask file for time step t0'''
        mask_file_name = 'rimmask_perc' + str(perc) + 'th' + '_t' + str(t0) + '.nc'
        create_mask_file(mask_file_name, path_out, nk, krange)
        # stats_file_name = 'rimstats_perc' + str(perc) + 'th' + '_t' + str(t0) + '.nc'
        # create_statistics_file(stats_file_name, path_out, nk, krange)


        '''(A) read in w-field, shift domain and define partial domain '''
        w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0) + '.nc'))

        for ik,k0 in enumerate(krange):
            # w_roll = np.roll(np.roll(w[:, :, k0], ishift, axis=0), jshift, axis=1)    # nbi
            w_roll = np.roll(w[:, :, k0], [ishift, jshift], [0, 1])
            w_ = w_roll[ic - id + ishift:ic + id + ishift, jc - jd + jshift:jc + jd + jshift]
            icshift = id
            jcshift = jd

            # plot_yz_crosssection(w, ic, path_out, t0)


            ''' (B) mask 2D field and turn mask from boolean (True: w>w_c) into integer (1: w>w_c)'''
            # Note:  no difference if percentile of total field w or subdomain w_
            w_c = np.percentile(w_, perc)
            # w_mask = True, if w<w_c
            # w_mask_r = True, if w>w_c
            w_mask = np.ma.masked_less(w_, w_c)
            w_mask_r = np.ma.masked_where(w_ > w_c, w_)
            if not w_mask_r.mask.any():
                print('STOP (t='+str(t0)+')' )
                continue
            else:
                w_bin_r = np.asarray(
                    [np.int(w_mask_r.mask.reshape(nx_ * ny_)[i]) for i in range(nx_ * ny_)]).reshape(nx_, ny_)

            # plot_s(w, w_c, t0, k0, path_fields, path_out)
            plot_w_field(w_c, perc, w, w_roll, w_, w_mask,
                         ishift, jshift, id, jd, ic, jc, icshift, jcshift,
                         k0, t0, dz, gw, nx_, ny_, ny, ny, path_out)
            del w_roll

            ''' (C) define outline of cold pool '''

            mask_aux = np.array(w_bin_r, copy=True)

            ''' (a) fill interior of mask '''
            imin = icshift
            imax = icshift
            jmin = jcshift
            jmax = jcshift
            di = 0
            dj = 0
            while (w_mask.mask[icshift+di, jcshift] or w_mask.mask[icshift-di, jcshift]):
                imin = np.minimum(icshift - di, imin)-1
                imax = np.maximum(icshift + di, imax)+1
                di += 1
            while (w_mask.mask[icshift, jcshift+dj] or w_mask.mask[icshift, jcshift-dj]):
                jmin = np.minimum(jcshift - dj, jmin)-1
                jmax = np.maximum(jcshift + dj, jmax)+1
                dj += 1
            rmax2 = np.maximum(np.maximum(imax-icshift,icshift-imin),np.maximum(jmax-jcshift,jcshift-jmin))**2

            di = 0
            while (icshift - di > imin or icshift + di < imax):
                dj = 0
                r2 = di ** 2 + dj ** 2
                while (r2 <= rmax2):
                    for si in [-1, 1]:
                        for sj in [-1, 1]:
                            r2 = di ** 2 + dj ** 2
                            if w_mask.mask[icshift+si*di,jcshift+sj*dj]:
                                mask_aux[icshift + si * di, jcshift + sj * dj] = 2
                    dj += 1
                di += 1

            # plt.figure()
            # plt.subplot(131)
            # plt.contourf(w_mask.mask.T, origin='lower')
            # plt.colorbar()
            # plt.title('w_mask')
            # ax = plt.subplot(132)
            # ax.imshow(w_mask.mask.T, origin='lower')
            # plt.plot([imin, imin], [0, ny_ - 1], 'w', linewidth=1)
            # plt.plot([imax, imax], [0, ny_ - 1], 'w', linewidth=1)
            # plt.title('w mask')
            # ax = plt.subplot(133)
            # ax.imshow(mask_aux.T, origin='lower')
            # circle1 = plt.Circle((icshift, jcshift), np.sqrt(rmax2), fill=False, color='w')
            # ax.add_artist(circle1)
            # plt.title('mask_aux')
            # plt.savefig('./test_mask_aux.png')

            ''' (b) find inner & outer rim '''
            rim_int = np.zeros((nx_, ny_), dtype=np.int)
            rim_out = np.zeros((nx_, ny_), dtype=np.int)
            rim_aux = np.zeros((nx_, ny_), dtype=np.int)
            rim_list_int = []
            rim_list_out = []

            di = 0
            dj = 0
            imin = icshift
            jmin = jcshift
            imax = icshift
            jmax = jcshift
            while (mask_aux[icshift+di, jcshift]>0 and icshift+di<nx_):
                imin = np.minimum(icshift - di, imin)-1
                di += 1
            while (mask_aux[icshift - di, jcshift] > 0 and icshift - di >= 0):
                imax = np.maximum(icshift + di, imax)+1
                di += 1
            while (mask_aux[jcshift, jcshift+dj]>0 and jcshift+dj<ny_):
                jmin = np.minimum(jcshift - dj, jmin)-1
                dj += 1
            while (mask_aux[jcshift, jcshift - dj] > 0 and jcshift - dj >= 0):
                jmax = np.maximum(jcshift + dj, jmax)+1
                dj += 1
            rmax2 = np.maximum(np.maximum(imax-icshift,icshift-imin),np.maximum(jmax-jcshift,jcshift-jmin))**2
            # plot_outlines(perc, w_mask, rim_int, rim_out, rim_list_out, rim_aux, rmax2, icshift, jcshift, imin, imax, jmin, jmax,
            #               nx_, ny_, k0, t0, path_out)
            for si in [-1, 1]:
                for sj in [-1, 1]:
                    for di in range(imax):
                        i = icshift + si*di
                        for dj in range(jmax):
                            j = jcshift + sj*dj
                            r2 = di ** 2 + dj ** 2
                            if r2 <= rmax2:
                                rim_aux[i,j] = 1
                                if w_mask_r.mask[i,j]:
                                    a = np.count_nonzero(w_bin_r[i - 1:i + 2, j - 1:j + 2])
                                    if a > 5 and a < 9:
                                        if np.sum(mask_aux[i-1:i+2,j-1:j+2]) > 9:
                                            rim_int[i, j] = 1
                                            rim_list_int.append((i, j))
                                        else:
                                            rim_out[i,j] = 1
                                            rim_list_out.append((i, j))
                                            # a = np.count_nonzero(w_bin_r[i - 1:i + 2, j - 1 + sj:j + 2 + sj])
                                            # if a <= 5 or a >= 9:
                                            #     print('breaking')
                                            #     break

            plot_outlines(perc, w_mask, rim_int, rim_out, rim_list_out, rim_aux, rmax2, icshift, jcshift, imin, imax, jmin, jmax,
                          nx_, ny_, k0, t0, path_out)
            del mask_aux
            # ''' dump mask and rim '''
            dump_mask(w_mask_r.mask, rim_int, rim_out, mask_file_name, path_out, k0, ik)

            ''' (D) Polar Coordinates & sort according to angle '''
            # (1) find/define center of mass (here = (ic/jc))
            # (2)
            # Once you create a tuple, you cannot edit it, it is immutable. Lists on the other hand are mutable,
            #   you can edit them, they work like the array object in JavaScript or PHP. You can add items,
            #   delete items from a list; but you can't do that to a tuple, tuples have a fixed size.
            nrim_out = len(rim_list_out)
            nrim_int = len(rim_list_int)
            for i, coord in enumerate(rim_list_out):
                rim_list_out[i] = (coord, (polar(coord[0] - icshift, coord[1] - jcshift)))
            for i, coord in enumerate(rim_list_int):
                rim_list_int[i] = (coord, (polar(coord[0] - icshift, coord[1] - jcshift)))
            # if rim already very close to subdomain (nx_,ny_), make domain larger
            if coord[0] >= nx_ - 3 or coord[1] >= ny_ - 3:
                print('!!! changing domain size', nx_, nx_ + 4)
                shift += 10
                id = irstar + shift
                jd = irstar + shift
                ishift = np.max(id - ic, 0)
                jshift = np.max(jd - jc, 0)
                nx_ = 2 * id
                ny_ = 2 * jd

            # sort list according to angle
            rim_list_out.sort(key=lambda tup: tup[1][1])
            rim_list_int.sort(key=lambda tup: tup[1][1])
            plot_rim_mask(w_, w_mask, rim_out, rim_int, rim_list_out, rim_list_int, icshift, jcshift, nx_, ny_,
                          t0, k0, path_out)

            del w_mask
            del rim_out, rim_int

            # average and interpolate for bins of 6 degrees
            angular_range = np.arange(0, 361, dphi)
            # - rim_intp_all = (phi(t,i_phi)[deg], phi(t,i_phi)[rad], r_out(t,i_phi)[m], r_int(t,i_phi)[m], D(t,i_phi)[m])
            #   (phi: angles at interval of 6 deg; r_out,int: outer,inner boundary of convergence zone; D: thickness of convergence zone)
            # - rim_vel = (phi(t,i_phi)[deg], phi(t,i_phi)[rad], r_out(t,i_phi)[m], U(t,i_phi)[m/s], dU(t, i_phi)[m/s**2])
            # - rim_vel_av = (r_av(t), U_av(t), dU_av/dt(t))
            rim_intp_all[0, it, :] = angular_range[:-1]
            rim_intp_all[1, it, :] = np.pi * rim_intp_all[0, it, :] / 180
            print('')
            i = 0
            for n, phi in enumerate(rim_intp_all[0, it, :]):
                phi_ = rim_list_out[i][1][1]
                r_aux = 0.0
                count = 0
                while (phi_ >= phi and phi_ < angular_range[n + 1]):
                    r_aux += rim_list_out[i][1][0]
                    count += 1
                    i += 1
                    if i < nrim_out:
                        phi_ = rim_list_out[i][1][1]
                    else:
                        # phi_ = angular_range[n+1]
                        i = 0  # causes the rest of n-, phi-loop to run through without entering the while-loop
                        # >> could probably be done more efficiently
                        break
                if count > 0:
                    rim_intp_all[2, it, n] = dx * r_aux / count
            i = 0
            for n, phi in enumerate(rim_intp_all[0, it, :]):
                phi_ = rim_list_int[i][1][1]
                r_aux = 0.0
                count = 0
                while (phi_ >= phi and phi_ < angular_range[n + 1]):
                    r_aux += rim_list_int[i][1][0]
                    count += 1
                    i += 1
                    if i < nrim_int:
                        phi_ = rim_list_int[i][1][1]
                    else:
                        # phi_ = angular_range[n+1]
                        i = 0  # causes the rest of n-, phi-loop to run through without entering the while-loop
                        # >> could probably be done more efficiently
                        break
                if count > 0:
                    rim_intp_all[3, it, n] = dx * r_aux / count
            print('')

            # dump statistics
            # dump_statistics_file(rim_intp_all, stats_file_name, path_out)
            dump_statistics_file(rim_intp_all[:,it,:], rim_vel[:,it,:], rim_vel_av[:,it], angular_range[:-1], stats_file_name, path_out, k0, ik, t0, it)


            # plot outline in polar coordinates r(theta)
            plot_angles(rim_list_out, rim_list_int, rim_intp_all[:,it,:], t0, path_out)
            plot_cp_outline_alltimes(rim_intp_all[:,0:it+1,:], timerange, dx, k0, path_out)

            rim_intp_all[4,:,:] = rim_intp_all[2, :, :] - rim_intp_all[3, :, :]     # rim thickness

            plot_rim_thickness(rim_intp_all[:,0:it+1,:], timerange[:it+1], dx, k0, path_out)
            del rim_list_out, rim_list_int


            ''' Compute radial velocity of rim '''
            rim_vel[0:3, it, :] = rim_intp_all[0:3, it, :]  # copy phi [deg + rad], r(phi)


            if it == 0:
                rim_vel_av[0, it] = np.average(np.ma.masked_less(rim_intp_all[2, it, :], 1.))
                rim_vel_av[1, it] = 0.0
            elif it > 0:
                # for n, phi in enumerate(rim_intp_all[0,it,:]):
                rim_vel[3, it, :] = (rim_intp_all[2, it, :] - rim_intp_all[2, it-1, :]) / dt
                rim_vel[4, it, :] = (rim_vel[3, it, :] - rim_vel[3, it-1, :]) / dt
                rim_vel_av[0, it] = np.average(np.ma.masked_less(rim_intp_all[2,it,:],1.))
                rim_vel_av[1, it] = np.average(np.ma.masked_where(rim_intp_all[2,it,:]>1., rim_vel[3,it,:]).data)
                rim_vel_av[2, it] = np.average(np.ma.masked_where(rim_intp_all[2,it,:]>1., rim_vel[4,it,:]).data)

                plot_cp_rim_averages(rim_vel[:, 0:it+1, :], rim_vel_av[:, :it+1], timerange[:it+1], k0, path_out)

            plot_cp_rim_velocity(rim_vel[:, 0:it + 1, :], rim_vel_av, k0, timerange, path_out)

    return

# ----------------------------------
import math
def polar(x, y):
    """
    :param x: x-coordinate
    :param y: y-coordinate
    :return: polar coordinates r, theta(degrees)
    """
    r = (x ** 2 + y ** 2) ** .5
    if y == 0:
        theta = 180 if x < 0 else 0
    elif x == 0:
        theta = 90 if y > 0 else 270
    elif x > 0:
        theta = math.degrees(math.atan(float(y) / x)) if y > 0 \
            else 360 + math.degrees(math.atan(float(y) / x))
    elif x < 0:
        theta = 180 + math.degrees(math.atan(float(y) / x))
    return r, theta
# ----------------------------------

# def create_statistics_file(file_name, path, n_phi, nt, timerange, nk, krange):
#     # output list with coordinates of rim points; mean velocity, radius etc.
#
#     print('-------- create statistics file -------- ' )
#     print(path + ', ' + file_name)
#     print('')
#     rootgrp = nc.Dataset(os.path.join(path, file_name), 'w', format='NETCDF4')
#     # # dimgrp = rootgrp.createGroup('dims')
#     ts_grp = rootgrp.createGroup('time')
#     ts_grp.createDimension('nt', nt)
#     var = ts_grp.createVariable('time', 'f8', ('nt'))
#     var[:] = timerange
#
#     # fields_grp = rootgrp.createGroup('fields')
#     # fields_grp.createDimension('nx', nx)
#     # fields_grp.createDimension('ny', ny)
#     # fields_grp.createDimension('nz', nz)
#
#     stats_grp = rootgrp.createGroup('stats')
#     stats_grp.createDimension('nt', nt)
#     stats_grp.createDimension('nphi', n_phi)    # number of segments
#     stats_grp.createVariable('r_out', 'f8', ('nt', 'nphi'))
#     stats_grp.createVariable('r_int', 'f8', ('nt', 'nphi'))
#     stats_grp.createVariable('D', 'f8', ('nt', 'nphi'))     # thickness of rim
#
#     pol_grp = rootgrp.createGroup('angles')
#     pol_grp.createDimension('nphi', n_phi)
#     var = pol_grp.createVariable('degree', 'f8', ('nphi'))
#     var = pol_grp.createVariable('radian', 'f8', ('nphi'))
#
#     prof_grp = rootgrp.createGroup('profiles')
#     prof_grp.createDimension('nz', nk)
#     var = prof_grp.createVariable('k_dumped', 'f8', ('nz'))
#     var[:] = np.zeros(nk)
#     var = prof_grp.createVariable('krange', 'f8', ('nz'))
#     var[:] = krange[:]
#     var = prof_grp.createVariable('zrange', 'f8', ('nz'))
#     var[:] = krange[:] * dz
#
#     rootgrp.close()
#
#     return

# def dump_statistics_file(rim_intp_all, file_name, path):
#     # - rim_intp_all = (phi(t,i_phi)[deg], phi(t,i_phi)[rad], r_out(t,i_phi)[m], r_int(t,i_phi)[m], D(t,i_phi)[m])
#     #   (phi: angles at interval of 6 deg; r_out,int: outer,inner boundary of convergence zone; D: thickness of convergence zone)
#     # - rim_vel = (phi(t,i_phi)[deg], phi(t,i_phi)[rad], r_out(t,i_phi)[m], U(t,i_phi)[m/s], dU(t, i_phi)[m/s**2])
#     # - rim_vel_av = (r_av(t), U_av(t), dU_av/dt(t))
#     rootgrp = nc.Dataset(os.path.join(path, file_name), 'r+', format='NETCDF4')
#
#     prof_grp = rootgrp.createGroup('profiles')
#     var = levels_grp.variables['k_dumped']
#     var[ik] = 1
#
#     rootgrp.close()
#     return

def create_statistics_file(file_name, path, n_phi, nt, timerange, nk, krange):
    # output list with coordinates of rim points; mean velocity, radius etc.

    print('-------- create statistics file -------- ' )
    print(path + ', ' + file_name)
    print('')
    rootgrp = nc.Dataset(os.path.join(path, file_name), 'w', format='NETCDF4')
    # # dimgrp = rootgrp.createGroup('dims')
    ts_grp = rootgrp.createGroup('time')
    ts_grp.createDimension('nt', nt)
    var = ts_grp.createVariable('time', 'f8', ('nt'))
    var[:] = timerange

    ts_grp = rootgrp.createGroup('timeseries')
    ts_grp.createDimension('nt', nt)
    ts_grp.createDimension('nz', nk)
    ts_grp.createVariable('r_av', 'f8', ('nt', 'nz'))
    ts_grp.createVariable('U_av', 'f8', ('nt', 'nz'))
    ts_grp.createVariable('dU_av', 'f8', ('nt', 'nz'))

    stats_grp = rootgrp.createGroup('stats')
    stats_grp.createDimension('nt', nt)
    stats_grp.createDimension('nz', nk)
    stats_grp.createDimension('nphi', n_phi)    # number of segments
    stats_grp.createVariable('r_out', 'f8', ('nt', 'nz', 'nphi'))
    stats_grp.createVariable('r_int', 'f8', ('nt', 'nz', 'nphi'))
    stats_grp.createVariable('D', 'f8', ('nt', 'nz', 'nphi'))     # thickness of rim
    stats_grp.createVariable('U', 'f8', ('nt', 'nz', 'nphi'))     # velocity of outer rim
    stats_grp.createVariable('dU', 'f8', ('nt', 'nz', 'nphi'))    # change of rim velocity with time ('acceleration of rim')

    pol_grp = rootgrp.createGroup('angles')
    pol_grp.createDimension('nphi', n_phi)
    var = pol_grp.createVariable('degree', 'f8', ('nphi'))
    var = pol_grp.createVariable('radian', 'f8', ('nphi'))

    prof_grp = rootgrp.createGroup('profiles')
    prof_grp.createDimension('nz', nk)
    var = prof_grp.createVariable('k_dumped', 'f8', ('nz'))
    var[:] = np.zeros(nk)
    var = prof_grp.createVariable('krange', 'f8', ('nz'))
    var[:] = krange[:]
    var = prof_grp.createVariable('zrange', 'f8', ('nz'))
    var[:] = krange[:] * dz

    rootgrp.close()

    return


def dump_statistics_file(rim_intp_all, rim_vel, rim_vel_av, angles, file_name, path, k0, ik, t0, it):
    # - rim_intp_all = (phi(i_phi)[deg], phi(i_phi)[rad], r_out(i_phi)[m], r_int(i_phi)[m], D(i_phi)[m])
    #   (phi: angles at interval of 6 deg; r_out,int: outer,inner boundary of convergence zone; D: thickness of convergence zone)
    # - rim_vel = (phi(i_phi)[deg], phi(i_phi)[rad], r_out(i_phi)[m], U(i_phi)[m/s], dU(i_phi)[m/s**2])
    # - rim_vel_av = (r_av(t), U_av(t), dU_av/dt(t))
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


    pol_grp = rootgrp.groups['angles']
    var = pol_grp.variables['degree']
    var[:] = angles[:]
    var = pol_grp.variables['radian']
    var[:] = np.pi / 180 * angles[:]

    prof_grp = rootgrp.groups['profiles']
    var = prof_grp.variables['k_dumped']
    var[ik] = 1

    rootgrp.close()
    return

def create_mask_file(file_name, path, nk, krange):
    print('-------- create mask file --------', nx_, ny_, nk)
    # file_name = 'rimmask_perc' + str(perc) + 'th' + '_t' + str(time) + '.nc'
    rootgrp = nc.Dataset(os.path.join(path, file_name), 'w', format='NETCDF4')

    mask_grp = rootgrp.createGroup('fields')
    mask_grp.createDimension('nx', nx_)
    mask_grp.createDimension('ny', ny_)
    mask_grp.createDimension('nz', nk)
    mask_grp.createVariable('mask', 'f8', ('nx', 'ny', 'nz'))
    mask_grp.createVariable('rim_inner', 'f8', ('nx', 'ny', 'nz'))
    mask_grp.createVariable('rim_outer', 'f8', ('nx', 'ny', 'nz'))

    levels_grp = rootgrp.createGroup('profiles')
    levels_grp.createDimension('nz', nk)
    var = levels_grp.createVariable('k_dumped', 'f8', ('nz'))
    var[:] = np.zeros(nk)
    var = levels_grp.createVariable('krange', 'f8', ('nz'))
    var[:] = krange[:]
    var = levels_grp.createVariable('zrange', 'f8', ('nz'))
    var[:] = krange[:] * dz
    rootgrp.close()
    print('')
    return

def dump_mask(mask, rim_int, rim_out, file_name, path, k0, ik):
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
    var[ik] = 1
    rootgrp.close()
    print('')
    return


#     # ----------------------------------------------------------------------
#
#     def create_statistics_file(self, path, file_name, time, ncomp, nvar, nz_):
#         print('create statistics file: '+ path+', '+ file_name)
#         # ncomp: number of Gaussian components in EM
#         # nvar: number of variables of multi-variate Gaussian components
#         rootgrp = nc.Dataset(os.path.join(path,file_name), 'w', format='NETCDF4')
#         dimgrp = rootgrp.createGroup('dims')
#         ts_grp = rootgrp.createGroup('time')
#         ts_grp.createDimension('nt',len(time)-1)
#         means_grp = rootgrp.createGroup('means')
#         means_grp.createDimension('nz', nz_)
#         means_grp.createDimension('ncomp', ncomp)
#         means_grp.createDimension('nvar', nvar)
#         cov_grp = rootgrp.createGroup('covariances')
#         cov_grp.createDimension('nz', nz_)
#         cov_grp.createDimension('ncomp', ncomp)
#         cov_grp.createDimension('nvar', nvar)
#         weights_grp = rootgrp.createGroup('weights')
#         weights_grp.createDimension('nz', nz_)
#         weights_grp.createDimension('ncomp', ncomp)
#         error_grp = rootgrp.createGroup('error')
#         error_grp.createDimension('nz', nz_)
#
#         var = ts_grp.createVariable('t','f8',('nt'))
#         for i in range(len(time)-1):
#             var[i] = time[i+1]
#         z_grp = rootgrp.createGroup('profiles')
#         z_grp.createDimension('nz', nz_)
#         var = z_grp.createVariable('z', 'f8', ('nz'))
#         for i in range(nz_):
#             var[i] = self.zrange[i]
#         var = z_grp.createVariable('k', 'f8', ('nz'))
#         for i in range(nz_):
#             var[i] = self.krange[i]
#         rootgrp.close()
#         return
#
#
# def dump_variable(path, group_name, data_, var_name, ncomp, nvar, nz_):
#     print('-------- dump variable --------', var_name, group_name, path)
#     # print('dump variable', path, group_name, var_name, data_.shape, ncomp, nvar)
#     rootgrp = nc.Dataset(path, 'r+')
#     if group_name == 'means':
#         # rootgrp = nc.Dataset(path, 'r+')
#         var = rootgrp.groups['means'].createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar'))
#         # var = nc.Dataset(path, 'r+').groups['means'].createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar'))
#         # var = nc.Dataset(path, 'r+').groups['means'].createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar'))[:,:,:]
#         var[:,:,:] = data_[:,:,:]
#
#     elif group_name == 'covariances':
#         var = rootgrp.groups['covariances'].createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar', 'nvar'))
#         var[:,:,:,:] = data_[:,:,:,:]
#
#     elif group_name == 'weights':
#         var = rootgrp.groups['weights'].createVariable(var_name, 'f8', ('nz', 'ncomp'))
#         var[:,:] = data_[:,:]
#
#     elif group_name == 'error':
#         var = rootgrp.groups['error'].createVariable(var_name, 'f8', ('nz'))
#         var[:] = data_[:]
#
#     elif group_name == 'profiles':
#         var = rootgrp.groups['profiles'].createVariable(var_name, 'f8', ('nz'))
#         var[:] = data_[:]
#
#     # # write_field(path, group_name, data, var_name)
#     # # print('--------')
#     rootgrp.close()
#     print('')
#     return

#----------------------------------------------------------------------

def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:, :, :]
    rootgrp.close()
    return data

if __name__ == '__main__':
    main()