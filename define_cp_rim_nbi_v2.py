import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os

# todo
# - problems for mask output, if domain size changes

from define_cp_rim_plottingfct import set_colorbars
from define_cp_rim_plottingfct import plot_yz_crosssection, plot_w_field, plot_s, \
    plot_outlines, plot_rim_mask, plot_angles, plot_cp_outline_alltimes, \
    plot_cp_rim_velocity, plot_cp_rim_averages, plot_cp_rim_thickness

def main():
    """
    Find inner and outer rim of mask based on a threshold (usually 95th percentile) of the vertical velocity.
    The rim is found by number of neighbours and filling interior of mask
    :param --path: The full path to the files
    :param --casename: casename
    :param --tmin: minimum time taken into account
    :param --tmax: maximum time taken into account
    :param --k0: level, at which the mask is computed

    :return: figures in repository 'figs_cp_rim'; 2D field output with mask, inner and outer rim at levels k=[kmin..kmax]; statistics with cold pool radius, rim velocity

    Details:

    1. read in w-field, shift field (roll) and define partial domain where to look for cold pool
    2. mask 2D field and turn mask from boolean (True: w>w_c) into integer (1: w>w_c)
    3. Define rim of cold pool as the outline of the mask; based on number of neighbours

    Output 3D fields:

    1. mask(x,y,k)
    2. rim_inner(x,y,k)
    3. rim_outer(x,y,k)
     (x=[0..nx_-1], y=[0..ny_-1], k=[kmin..kmax])
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

    global path_fields, path_out, path_stats
    if args.path:
        path = args.path
    else:
        # path = '/Users/bettinameyer/polybox/ClimatePhysics/Copenhagen/Projects/LES_ColdPool/' \
        #        'triple_3D_noise/Out_CPDry_triple_dTh2K/'
        path = '/nbi/ac/cond1/meyerbe/ColdPools/triple_3D_noise/Out_CPDry_triple_Th3K/'
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
    cm_vir = plt.cm.get_cmap('jet')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    set_colorbars(cm_bwr, cm_vir, cm_grey)      # to set colorbars as global functions in define_cp_rim_plottingfct.py

    # define subdomain to scan
    global nx_, ny_
    if case_name == 'ColdPoolDry_triple_3D':
        flag = 'triple'
        d = np.int(np.round( (ny+gw) / 2))
        a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        ic = np.int(np.round(a / 2))
        jc = np.int(np.round(d / 2))
        shift = 60
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

    print('ic,jc,id,jc', ic, jc, id, jd)
    print('nx_,ny_', nx_, ny_)
    print('shift, ishift, jshift', shift, ishift, jshift)
    # percentile for threshold
    # perc = 95     # tested for triple 3D, dTh=3K, t=400s
    perc = 98       # tested for triple 3D, dTh=10K, t=100-400s


    # (A) read in w-field
    #       - shift field (roll) and define partial domain where to look for cold pool
    # (B) mask 2D field and turn mask from boolean (True: w>w_c) into integer (1: w>w_c)
    # (C) Define rim of cold pool as the outline of the mask; based on number of neighbours

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
    create_statistics_file(stats_file_name, path_stats, n_phi, nt, timerange, nk, krange)

    for it,t0 in enumerate(timerange):
        if it > 0:
            dt = t0-timerange[it-1]
        else:
            dt = t0
        print('--- time: '+ str(t0), '(dt='+str(dt)+') ---')
        print('nx_, ny_', nx_, ny_, 'id, jd', id, jd, shift)

        ''' create mask file for time step t0'''
        mask_file_name = 'rimmask_perc' + str(perc) + 'th' + '_t' + str(t0) + '.nc'
        create_mask_file(mask_file_name, path_stats, nk, krange, ic, jc, ishift, jshift)
        # stats_file_name = 'rimstats_perc' + str(perc) + 'th' + '_t' + str(t0) + '.nc'
        # create_statistics_file(stats_file_name, path_stats, nk, krange)


        '''(A) read in w-field, shift domain and define partial domain '''
        w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0) + '.nc'))

        for ik,k0 in enumerate(krange):
            print('level: k=' + str(k0), '(z=' + str(k0 * dz) + 'm)')
            w_roll = np.roll(np.roll(w[:, :, k0], ishift, axis=0), jshift, axis=1)
            w_ = w_roll[ic - id + ishift:ic + id + ishift, jc - jd + jshift:jc + jd + jshift]
            icshift = id -1
            jcshift = jd
            print('icshift', icshift, ic + ishift)

            if flag == 'triple':
                plot_yz_crosssection(w, ic, path_out, t0)


            ''' (B) mask 2D field and turn mask from boolean (True: w>w_c) into integer (1: w>w_c)'''
            # Note:  no difference if percentile of total field w or subdomain w_
            # w_c = np.percentile(w_, perc)
            w_c = np.percentile(w, perc)
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
            # (imin, imax, jmin, jmax) enclose inner circle of cold pool
            imin1 = icshift
            imax1 = icshift
            jmin1 = jcshift
            jmax1 = jcshift
            i = 0
            j = 0
            while (w_mask.mask[icshift+i, jcshift] or w_mask.mask[icshift-i, jcshift]):
                imin1 = np.minimum(icshift - i, imin1)-1
                imax1 = np.maximum(icshift + i, imax1)+1
                i += 1
            while (w_mask.mask[icshift, jcshift+j] or w_mask.mask[icshift, jcshift-j]):
                jmin1 = np.minimum(jcshift - j, jmin1)-1
                jmax1 = np.maximum(jcshift + j, jmax1)+1
                j += 1
            rmax2_int = (3+ np.maximum(np.maximum(imax1-icshift,icshift-imin1),
                                   np.maximum(jmax1-jcshift,jcshift-jmin1)) )**2
            print('rmax2', rmax2_int, np.sqrt(rmax2_int))
            print('imin, imax, jmin, jmax', imin1, imax1, jmin1, jmax1)
            print('icshift, jcshift', icshift, jcshift)

            # fill interior of cold pool (within inner rim)
            i = 0
            while (icshift - i >= imin1-1 or icshift + i < imax1):
                # while (icshift - i >= 0 or icshift + i < nx_):
                # while (icshift - i > imin or icshift + i < imax): # older version
                j = 0
                r2 = i ** 2 + j ** 2
                while (r2 <= rmax2_int):
                    for si in [-1, 1]:
                        for sj in [-1, 1]:
                            r2 = i ** 2 + j ** 2
                            if w_mask.mask[icshift+si*i,jcshift+sj*j]:
                                mask_aux[icshift + si * i, jcshift + sj * j] = 2
                    j += 1
                i += 1



            ''' (b) find inner & outer rim '''
            rim_int = np.zeros((nx_, ny_), dtype=np.int)
            rim_out = np.zeros((nx_, ny_), dtype=np.int)
            rim_aux = np.zeros((nx_, ny_), dtype=np.int)
            rim_list_int = []
            rim_list_out = []



            # find radius of circle that encloses whole cp rim (rmax2 >= rim)
            i = 0
            j = 0
            imin = icshift
            jmin = jcshift
            imax = icshift
            jmax = jcshift
            while (mask_aux[icshift+i, jcshift]>0 and icshift+i<nx_):
                imin = np.minimum(icshift - i, imin)-1
                i += 1
            while (mask_aux[icshift - i, jcshift] > 0 and icshift - i >= 0):
                imax = np.maximum(icshift + i, imax)+1
                i += 1
            while (mask_aux[jcshift, jcshift+j]>0 and jcshift+j<ny_):
                jmin = np.minimum(jcshift - j, jmin)-1
                j += 1
            while (mask_aux[jcshift, jcshift - j] > 0 and jcshift - j >= 0):
                jmax = np.maximum(jcshift + j, jmax)+1
                j += 1
            rmax2 = np.maximum(np.maximum(imax-icshift,icshift-imin),np.maximum(jmax-jcshift,jcshift-jmin))**2
            # plot_outlines(perc, w_mask, rim_int, rim_out, rim_list_out, rim_aux, rmax2, icshift, jcshift, imin, imax, jmin, jmax,
            #               nx_, ny_, k0, t0, path_out)

            plt.figure()
            plt.subplot(131)
            plt.imshow(w_mask.mask.T, origin='lower')
            plt.colorbar(shrink=0.3)
            plt.title('w_mask.mask')
            ax = plt.subplot(132)
            ax.imshow(mask_aux.T, origin='lower')
            plt.colorbar(shrink=0.3)
            plt.plot([imin1, imin1], [0, ny_ - 1], 'w', linewidth=1)
            plt.plot([imax1, imax1], [0, ny_ - 1], 'w', linewidth=1)
            circle1 = plt.Circle((icshift, jcshift), np.sqrt(rmax2_int), fill=False, color='w')
            ax.add_artist(circle1)
            plt.tight_layout()
            plt.title('mask_aux, (before)')
            ax = plt.subplot(133)
            plt.imshow(mask_aux.T, origin='lower')
            plt.colorbar(shrink=0.3)
            plt.plot([imin, imin], [0, ny_ - 1], 'w', linewidth=1)
            plt.plot([imax, imax], [0, ny_ - 1], 'w', linewidth=1)
            circle1 = plt.Circle((icshift, jcshift), np.sqrt(rmax2), fill=False, color='w')
            ax.add_artist(circle1)
            plt.plot(icshift, jcshift, 'ow', markersize=5)
            plt.tight_layout()
            plt.title('mask_aux (after)')
            plt.savefig(os.path.join(path_out, 'test_mask_aux_t' + str(t0) + '.png'))
            plt.close()

            for si in [-1, 1]:
                for sj in [-1, 1]:
                    for di in range(imax):
                        i = icshift + si*di
                        for dj in range(jmax):
                            j = jcshift + sj*dj
                            r2 = di ** 2 + dj ** 2
                            if r2 <= rmax2:
                                rim_aux[i,j] = 1        # just to see which pixels where tested
                                if w_mask_r.mask[i,j] and r2 > rmax2_int:       # outer rim
                                    a = np.count_nonzero(w_bin_r[i - 1:i + 2, j - 1:j + 2])
                                    if a > 5 and a < 9:     # btw 4-8 neighbouring pixels belong to rim
                                        rim_out[i,j] = 1
                                        rim_list_out.append((i,j))
                                elif w_mask_r.mask[i,j]:
                                    a = np.count_nonzero(w_bin_r[i - 1:i + 2, j - 1:j + 2])
                                    if a > 5 and a < 9:  # btw 4-8 neighbouring pixels belong to rim
                                        rim_int[i, j] = 1
                                        rim_list_int.append((i, j))
                                        # if np.sum(mask_aux[i-1:i+2,j-1:j+2]) > 9:
                                        #     rim_int[i, j] = 1
                                        #     rim_list_int.append((i, j))
                                        # else:
                                        #     rim_out[i,j] = 1
                                        #     rim_list_out.append((i, j))
                                        #     # a = np.count_nonzero(w_bin_r[i - 1:i + 2, j - 1 + sj:j + 2 + sj])
                                        #     # if a <= 5 or a >= 9:
                                        #     #     print('breaking')
                                        #     #     break

            plot_outlines(perc, w_mask, rim_int, rim_out, rim_list_out, rim_aux, rmax2_int, icshift, jcshift, imin, imax, jmin, jmax,
                          nx_, ny_, k0, t0, path_out)
            del mask_aux
            # ''' dump mask and rim '''
            dump_mask(w_mask_r.mask, rim_int, rim_out, mask_file_name, path_stats, k0, ik)


            ''' (D) Polar Coordinates & sort according to angle '''
            # (1) find/define center of mass (here = (ic/jc))
            # (2)
            # Once you create a tuple, you cannot edit it, it is immutable. Lists on the other hand are mutable,
            #   you can edit them, they work like the array object in JavaScript or PHP. You can add items,
            #   delete items from a list; but you can't do that to a tuple, tuples have a fixed size.
            nrim_out = len(rim_list_out)
            nrim_int = len(rim_list_int)
            for i, coord in enumerate(rim_list_int):
                rim_list_int[i] = (coord, (polar(coord[0] - icshift, coord[1] - jcshift)))
            for i, coord in enumerate(rim_list_out):
                rim_list_out[i] = (coord, (polar(coord[0] - icshift, coord[1] - jcshift)))
            # if rim already very close to subdomain (nx_,ny_), make domain larger
            print('rim size: ', coord[0], coord[1], 'domain size nx_, ny_', nx_, ny_)
            if coord[0] >= nx_ - 5 or coord[1] >= ny_ - 5:
                print('!!! changing domain size', nx_, nx_ +10)
                shift += 10
                id = irstar + shift
                jd = irstar + shift
                ishift = np.max(id - ic, 0)
                jshift = np.max(jd - jc, 0)
                if 2 * id <= nx and 2 * jd <= ny:
                    nx_ = 2 * id
                    ny_ = 2 * jd
                else:
                    print('!!! reached domain size')

            # sort list according to angle
            rim_list_out.sort(key=lambda tup: tup[1][1])
            rim_list_int.sort(key=lambda tup: tup[1][1])
            plot_rim_mask(w_, w_mask, rim_out, rim_int, rim_list_out, rim_list_int,
                          icshift, jcshift, nx_, ny_,
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
            rim_intp_all[4, :, :] = rim_intp_all[2, :, :] - rim_intp_all[3, :, :]    # rim thickness
            print('')

            # dump statistics
            # dump_statistics_file(rim_intp_all, stats_file_name, path_stats)
            dump_statistics_file(rim_intp_all[:,it,:], rim_vel[:,it,:], rim_vel_av[:,it], angular_range[:-1],
                                 stats_file_name, path_stats, k0, ik, t0, it)


            # plot outline in polar coordinates r(theta)
            plot_angles(rim_list_out, rim_list_int, rim_intp_all[:,it,:], t0, path_out)
            plot_cp_outline_alltimes(rim_intp_all[:,0:it+1,:], perc, k0, timerange, dx, path_out)


            plot_cp_rim_thickness(rim_intp_all[:,0:it+1,:], perc, k0, timerange[:it+1], dx, path_out)
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

                plot_cp_rim_averages(rim_vel[:, 0:it+1, :], rim_vel_av[:, :it+1], perc, k0, timerange[:it+1], path_out)

            plot_cp_rim_velocity(rim_vel[:, 0:it + 1, :], rim_vel_av, perc, k0, timerange, path_out)
            print('')

    return

# ----------------------------------------------------------------------
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

# ----------------------------------------------------------------------
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
    var.units = "s"
    var[:] = timerange

    ts_grp = rootgrp.createGroup('timeseries')
    ts_grp.createDimension('nt', nt)
    ts_grp.createDimension('nz', nk)
    var = ts_grp.createVariable('r_av', 'f8', ('nt', 'nz'))
    var.units = "m"
    var = ts_grp.createVariable('U_av', 'f8', ('nt', 'nz'))
    var.units = "m/s"
    var = ts_grp.createVariable('dU_av', 'f8', ('nt', 'nz'))
    var.units = "m/s^2"

    stats_grp = rootgrp.createGroup('stats')
    stats_grp.createDimension('nt', nt)
    stats_grp.createDimension('nz', nk)
    stats_grp.createDimension('nphi', n_phi)    # number of segments
    stats_grp.createVariable('r_out', 'f8', ('nt', 'nz', 'nphi'))
    stats_grp.createVariable('r_int', 'f8', ('nt', 'nz', 'nphi'))
    var = stats_grp.createVariable('D', 'f8', ('nt', 'nz', 'nphi'))     # thickness of rim
    var.units = "m"
    var = stats_grp.createVariable('U', 'f8', ('nt', 'nz', 'nphi'))     # velocity of outer rim
    var.units = "m/s"
    var = stats_grp.createVariable('dU', 'f8', ('nt', 'nz', 'nphi'))    # change of rim velocity with time ('acceleration of rim')
    var.units = "m/s^2"

    pol_grp = rootgrp.createGroup('angles')
    pol_grp.createDimension('nphi', n_phi)
    var = pol_grp.createVariable('phi_deg', 'f8', ('nphi'))
    var.units = "degrees"
    var = pol_grp.createVariable('phi_rad', 'f8', ('nphi'))
    var.units = "radian"

    prof_grp = rootgrp.createGroup('profiles')
    prof_grp.createDimension('nz', nk)
    var = prof_grp.createVariable('k_dumped', 'f8', ('nz'))
    var[:] = np.zeros(nk)
    var.description = "0:level not dumped, 1:level dumped"
    var = prof_grp.createVariable('krange', 'f8', ('nz'))
    var[:] = krange[:]
    var.units = "-"
    var = prof_grp.createVariable('zrange', 'f8', ('nz'))
    var[:] = krange[:] * dz
    var.units = "m"

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
    var = pol_grp.variables['phi_deg']
    var[:] = angles[:]
    var = pol_grp.variables['phi_rad']
    var[:] = np.pi / 180 * angles[:]

    prof_grp = rootgrp.groups['profiles']
    var = prof_grp.variables['k_dumped']
    var[ik] = 1

    rootgrp.close()
    return


def create_mask_file(file_name, path, nk, krange, ic, jc, ishift, jshift):
    print('-------- create mask file --------')
    # file_name = 'rimmask_perc' + str(perc) + 'th' + '_t' + str(time) + '.nc'
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

    levels_grp = rootgrp.createGroup('profiles')
    levels_grp.createDimension('nz', nk)
    var = levels_grp.createVariable('k_dumped', 'f8', ('nz'))
    var[:] = np.zeros(nk)
    var.description = "0:level not dumped, 1:level dumped"
    var = levels_grp.createVariable('krange', 'f8', ('nz'))
    var[:] = krange[:]
    var.units = "-"
    var = levels_grp.createVariable('zrange', 'f8', ('nz'))
    var[:] = krange[:] * dz
    var.units = "m"
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