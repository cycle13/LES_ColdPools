import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os

from define_cp_rim_plottingfct import set_colorbars
from define_cp_rim_plottingfct import plot_yz_crosssection, plot_w_field, plot_s, \
    plot_outlines, plot_rim_mask, plot_angles, plot_cp_outline_alltimes, \
    plot_cp_rim_velocity, plot_cp_rim_averages

def main():
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_vir = plt.cm.get_cmap('viridis')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    set_colorbars(cm_bwr, cm_vir, cm_grey)

    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("--casename")
    parser.add_argument("--path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--k0")
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
        tmax = tmin + 100
    timerange = np.arange(tmin,tmax,100)
    nt = len(timerange)

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

    # define subdomain to scan
    # --- for triple coldpool ---
    # d = np.int(np.round(ny / 2))      # nbi
    d = np.int(np.round(ny / 2)) + gw
    a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
    rstar = 5000.0  # half of the width of initial cold-pools [m]
    irstar = np.int(np.round(rstar / dx))
    ic = np.int(np.round(a / 2))
    jc = np.int(np.round(d / 2))
    shift = 40
    id = irstar + shift
    jd = irstar + shift
    ishift = np.max(id-ic,0)
    jshift = np.max(jd-jc,0)
    nx_ = 2 * id
    ny_ = 2 * jd

    print('ic,jc,id,jc,nx_,ny_',ic,jc,id,jd,nx_,ny_)

    # (A) read in w-field
    #       - shift field (roll) and define partial domain where to look for cold pool
    # (B) mask 2D field and turn mask from boolean (True: w>w_c) into integer (1: w>w_c)
    # (C) Define rim of cold pool as the outline of the mask; based on number of neighbours

    # define general arrays
    if args.k0:
        k0 = np.int(args.k0)
    else:
        k0 = 5      # level
    dphi = 6        # angular resolution for averaging of radius
    n_phi = 360 / dphi
    # rim_intp_int: inner rim of mask
    # - rim_intp_int = (phi(t,i_phi)[deg], phi(t,i_phi)[rad], r(t,i_phi))
    # rim_intp_out: outer rim of mask
    # - rim_intp_out = (phi(t,i_phi)[deg], phi(t,i_phi)[rad], r(t,i_phi))
    # - rim_intp_all = (phi(t,i_phi)[deg], phi(t,i_phi)[rad], r_out(t,i_phi), r_int(t,i_phi)
    # - rim_vel = (phi(t,i_phi)[deg], phi(t,i_phi)[rad], r(t,i_phi), U(t,i_phi), dU(t, i_phi))
    # - rim_vel_av = (r_av(t), U_av(t), dU_av/dt(t))
    rim_intp_all = np.zeros(shape=(4, nt, n_phi), dtype=np.double)
    rim_vel = np.zeros(shape=(5, nt, n_phi), dtype=np.double)
    rim_vel_av = np.zeros(shape=(3, nt))

    for it,t0 in enumerate(timerange):
        if it > 0:
            dt = t0-timerange[it-1]
        else:
            dt = t0
        print('time: '+ str(t0), '(dt='+str(dt)+')')

        '''(A) read in w-field, shift domain and define partial domain '''
        w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0) + '.nc'))
        # w_roll = np.roll(np.roll(w[:, :, k0], ishift, axis=0), jshift, axis=1)    # nbi
        w_roll = np.roll(w[:, :, k0], [ishift, jshift], [0, 1])
        w_ = w_roll[ic - id + ishift:ic + id + ishift, jc - jd + jshift:jc + jd + jshift]
        icshift = id
        jcshift = jd

        plot_yz_crosssection(w, ic, path_out, t0)


        ''' (B) mask 2D field and turn mask from boolean (True: w>w_c) into integer (1: w>w_c)'''
        perc = 90
        # ??? or use percentile of total field w: np.percentile(w, perc)
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
        del w, w_roll

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
                            mask_aux[icshift+si*di,jcshift+sj*dj] = 2
                dj += 1
            di += 1

        plt.figure()
        plt.subplot(131)
        plt.contourf(w_mask.mask.T, origin='lower')
        plt.colorbar()
        plt.title('w_mask')
        ax = plt.subplot(132)
        ax.imshow(w_mask.mask.T, origin='lower')
        plt.plot([imin, imin], [0, ny_ - 1], 'w', linewidth=1)
        plt.plot([imax, imax], [0, ny_ - 1], 'w', linewidth=1)
        plt.title('w mask')
        ax = plt.subplot(133)
        ax.imshow(mask_aux.T, origin='lower')
        circle1 = plt.Circle((icshift,jcshift), np.sqrt(rmax2), fill=False, color='w')
        ax.add_artist(circle1)
        plt.title('mask_aux')
        plt.savefig('./test_mask_aux.png')

        ''' (b) find inner&outer rim '''
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
        print('imin,etc', imin, imax, jmin, jmax)
        plot_outlines(perc, w_mask, rim_int, rim_out, rim_list_out, rim_aux, rmax2, icshift, jcshift, imin, imax, jmin, jmax,
                      nx_, ny_, t0, path_out)
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
                      nx_, ny_, t0, path_out)
        del mask_aux

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
        plot_rim_mask(w_, w_mask, rim_out, rim_int, rim_list_out, rim_list_int, icshift, jcshift, nx_, ny_, t0, path_out)
        del w_mask

        # average and interpolate for bins of 6 degrees
        angular_range = np.arange(0, 361, dphi)
        # - rim_intp_all = (phi[t,deg], phi[t,rad], r_out(t,phi))
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
                rim_intp_all[2, it, n] = dx*r_aux / count
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
                rim_intp_all[3, it, n] = dx*r_aux / count
        print('')


        # plot outline in polar coordinates r(theta)
        plot_angles(rim_list_out, rim_list_int, rim_intp_all[:,it,:], t0, path_out)
        plot_cp_outline_alltimes(rim_intp_all[:,0:it+1,:], timerange, dx, path_out)
        del rim_list_out, rim_list_int


        ''' Compute radial velocity of rim '''
        rim_vel[0:3, it, :] = rim_intp_all[0:3, it, :]  # copy phi [deg + rad], r(phi)


        if it == 0:
            rim_vel_av[0, it] = np.average(np.ma.masked_less(rim_intp_all[2, it, :], 1.))
            rim_vel_av[1, it] = 0.0
        elif it > 0:
            # for n, phi in enumerate(rim_intp_all[0,it,:]):
            rim_vel[3, it, :] = (rim_intp_all[2, it, :] - rim_intp_all[2, it-1, :])*dx / dt
            rim_vel_av[0, it] = np.average(np.ma.masked_less(rim_intp_all[2,it,:],1.))
            rim_vel_av[1, it] = np.average(np.ma.masked_where(rim_intp_all[2,it,:]>1., rim_vel[3,it,:]).data)

            plot_cp_rim_averages(rim_vel[:, 0:it+1, :], rim_vel_av[:, :it+1], timerange[:it+1], path_out)

        plot_cp_rim_velocity(rim_vel[:, 0:it + 1, :], rim_vel_av, timerange, path_out)

    return

# ----------------------------------
import math
def polar(x, y):
    """returns r, theta(degrees)
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

def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:, :, :]
    rootgrp.close()
    return data

if __name__ == '__main__':
    main()