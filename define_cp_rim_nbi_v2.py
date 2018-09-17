import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os

from define_cp_rim_plottingfct import plot_yz_crosssection, plot_w_field, plot_s, plot_outlines

def main():
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
        # path = '/Users/bettinameyer/polybox/ClimatePhysics/Copenhagen/Projects/LES_ColdPool/' \
        #        'triple_3D_noise/Out_CPDry_triple_dTh2K/'
        path = '/nbi/ac/cond1/meyerbe/ColdPools/triple_3D_noise/Out_CPDry_triple_Th3K/'
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
        tmax = 200
    timerange = np.arange(tmin, tmax, 100)
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
    cm_vir = plt.cm.get_cmap('jet')
    cm_grey = plt.cm.get_cmap('gist_gray_r')

    # define subdomain to scan
    # --- for triple coldpool ---
    d = np.int(np.round(ny / 2))
    a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
    rstar = 5000.0  # half of the width of initial cold-pools [m]
    irstar = np.int(np.round(rstar / dx))
    ic = np.int(np.round(a / 2))
    jc = np.int(np.round(d / 2))
    shift = 20
    id = irstar + shift
    jd = irstar + shift
    ishift = np.max(id - ic, 0)
    jshift = np.max(jd - jc, 0)
    nx_ = 2 * id
    ny_ = 2 * jd

    print('ic,jc,id,jc,nx_,ny_', ic, jc, id, jd, nx_, ny_)

    # (A) read in w-field
    #       - shift field (roll) and define partial domain where to look for cold pool
    # (B) mask 2D field and turn mask from boolean (True: w>w_c) into integer (1: w>w_c)
    # (C) Define rim of cold pool as the outline of the mask; based on number of neighbours

    # define general arrays
    if args.k0:
        k0 = np.int(args.k0)
    else:
        k0 = 5  # level
    dphi = 6  # angular resolution for averaging of radius
    n_phi = 360 / dphi
    # rim_intp_int: inner rim of mask
    # - rim_intp_int = (phi(t,i_phi)[deg], phi(t,i_phi)[rad], r(t,i_phi))
    # rim_intp_out: outer rim of mask
    # - rim_intp_out = (phi(t,i_phi)[deg], phi(t,i_phi)[rad], r(t,i_phi))
    # - rim_vel = (phi(t,i_phi)[deg], phi(t,i_phi)[rad], r(t,i_phi), U(t,i_phi), dU(t, i_phi))
    # - rim_vel_av 0 (r_av(t), U_av(t), dU_av/dt(t))
    rim_intp_out = np.ndarray(shape=(3, nt, n_phi), dtype=np.double)
    rim_vel = np.ndarray(shape=(5, nt, n_phi), dtype=np.double)
    rim_vel_av = np.ndarray(shape=(2, nt))

    for it,t0 in enumerate(timerange):
        if it > 0:
            dt = t0-timerange[it-1]
        else:
            dt = t0
        print('time: '+ str(t0), '(dt='+str(dt)+')')

        '''(A) read in w-field, shift domain and define partial domain '''
        w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0) + '.nc'))
        w_roll = np.roll(np.roll(w[:, :, k0], ishift, axis=0), jshift, axis=1)
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
        del w, w_roll, w_mask

        ''' (C) define outline of cold pool '''
        rim = np.zeros((nx_, ny_), dtype=np.int)
        rim_aux = np.zeros((nx_, ny_), dtype=np.int)
        nx_2 = np.int(nx_/2)
        ny_2 = np.int(ny_/2)
        rim_list = []
        # plt.figure()
        for si in [-1,1]:
            for sj in [-1,1]:
                stop_flag = False
                for i in range(icshift,nx_2 + si*nx_2,si):
                    if not stop_flag:
                        for j in range(jcshift,ny_2 + sj*ny_2, sj):
                            if i==70:
                                print(i,j)
                            if i == 70:
                                print('yes')
                            rim_aux[i,j] = 1
                            if w_mask_r.mask[i,j]:
                                a = np.count_nonzero(w_bin_r[i - 1:i + 2, j - 1:j + 2])
                                if a > 5 and a < 9:
                                    rim[i, j] = 1
                                    rim_list.append((i, j))
                                    if j == jcshift:
                                        stop_flag = True
                                    a = np.count_nonzero(w_bin_r[i - 1:i + 2, j-1+sj:j+2+sj])
                                    if i==70:
                                        print(i,j,'a(i=70,j+1)',a)
                                    if a <= 5 or a >= 9:
                                        print('breaking')
                                        break

        print('')
        print('icshift, jcshicft', icshift, jcshift)
        print(rim_list)


        # plt.xlim([0,nx])
        # plt.ylim([0,ny])
        # plt.savefig('./test_coverage.png')
        # plt.close()

        plot_outlines(perc, w_mask_r, rim, rim_list, rim_aux, icshift, jcshift, nx_, ny_, t0, path_out)

    return


# ----------------------------------

def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:,:,:]
    rootgrp.close()
    return data

if __name__ == '__main__':
    main()