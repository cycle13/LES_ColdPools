import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os

label_size = 10
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.labelsize'] = 15

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    # parser.add_argument("casename")
    # parser.add_argument("path")
    # parser.add_argument("--k0")
    # parser.add_argument("--tmin")
    # parser.add_argument("--tmax")
    # parser.add_argument("--path_tracers")
    parser.add_argument("--t_shift")
    parser.add_argument("--x_shift")
    args = parser.parse_args()
    # set_input_parameters(args)

    path_olga = '/nbi/ac/conv1/henneb/results/coldpool/lindp2K_13/output/cp/'
    path_out_figs = os.path.join('/nbi/home/meyerbe/paper_olga/figs_2D_real/')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    file_olga = 'coldpool_tracer_out_all.txt'
    CP_id_olga = 3
    col_id = 4  # column in textfile of CP-ID
    # tracer file:
    # - col=0:      timestep of simulation
    # - col=1:      age of CP (timestep since beginning of first CP)
    # - col=2:      tracer id
    # - col=3:      CP id
    # - col=4,5:
    # - col=8:      tracer radius
    # - col=14,15:  CP center (xc,yc)


    path_tracer_file = os.path.join(path_olga, file_olga)

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')

    print('???? dt_fields, dx')
    dt_fields = 100
    dx = 200
    dy = dx

    if args.t_shift:
        t_shift = np.int(args.t_shift)
    else:
        t_shift = 0
    if args.x_shift:
        x_shift = np.int(args.x_shift)
    else:
        x_shift = 0
    print('t_shift: '+str(t_shift))
    print('x_shift: '+str(x_shift))


    print('')
    print('path figs:    ' + path_out_figs)
    print('path tracers: ' + path_tracer_file)
    print('')




    cp_id = CP_id_olga
    print('CP ID: ' + str(cp_id))
    # n_cps = get_number_cps(path_tracer_file)
    n_cps = 2265
    print('number of CPs: ' + str(n_cps))
    # n_tracers = get_number_tracers(path_tracer_file)
    n_tracers = 999
    print('number of tracers: ' + str(n_tracers))

    # tau, t_ini, t_end = get_cp_lifetime(cp_id, path_tracer_file)
    tau = 5
    t_ini = 47
    t_end = 51
    lifetime = np.arange(t_ini, t_end+1)
    print('CP lifetime: ' + str(tau), lifetime)
    xc, yc = get_cp_center(cp_id, tau, n_tracers, path_tracer_file)
    ic = np.asarray([np.int(i) for i in xc])
    jc = np.asarray([np.int(i) for i in yc])
    print('CP center: ')
    print('xc: ', xc)
    print('yc: ', yc)
    # print('           '+ str(ic)+', '+str(jc))
    print('')


    # ''' tracer coordinates '''
    # coordinates = get_tracer_coords(cp_id, n_cps, n_tracers, lifetime, path_tracer_file)
    # print('')


    # ''' averaging tracer radius'''
    # aux = np.zeros(shape=coordinates.shape)
    # for it in range(tau):
    #     aux[it,:,0] = coordinates[it,:,0] - xc[it]
    #     aux[it,:,1] = coordinates[it,:,0] - yc[it]
    # rad_ = np.sqrt(aux[:,:,0]**2 + aux[:,:,1]**2)
    # rad = np.average(rad_, axis=1)
    # del aux
    # print('averaging: ', rad_.shape, rad.shape, coordinates.shape)
    # # fig_name = 'dist_tracers.png'
    # # fig, axis = plt.subplots(1, 2, figsize=(11, 6))
    # # ax0 = axis[0]
    # # ax0.plot(rad[:])
    # # plt.tight_layout()
    # # fig.savefig(os.path.join(path_out_figs, fig_name))
    # # plt.close(fig)


    '''compute & average radial velocity field'''
    # read in fields
    k0 = 0
    path_fields = path_out_figs
    rootgrp = nc.Dataset(os.path.join(path_fields, 'input_u.nc'))
    u_in = rootgrp.variables['var1'][:, k0, :, :]
    rootgrp.close()
    rootgrp = nc.Dataset(os.path.join(path_fields, 'input_v.nc'))
    v_in = rootgrp.variables['var1'][:, k0, :, :]
    rootgrp.close()

    # define subdomain for which to compute radial veloctiy wrt CP center (xc,yc)
    [nx,ny] = u_in.shape[1:3]
    lx = 10e3
    irange = np.int(np.minimum(np.int(lx / dx), np.minimum(nx - xc[0], xc[0])))
    jrange = np.int(np.minimum(np.int(lx / dx), np.minimum(ny - yc[0], yc[0])))
    nx_ = 2 * irange
    ny_ = 2 * jrange
    print('subdomain: ')
    print(nx, xc, xc / dx, lx, lx / dx, irange)
    print(ny, yc, yc / dx)

    # compute theta, r
    # th_field, r_field = compute_radius(50, 50, 50, 50, 100, 100, path_out_figs)     # ic at center; whole domain
    # th_field, r_field = compute_radius(50, 50, 40, 40, 100, 100, path_out_figs)   # ic at cetner; part domain
    # n = 100
    # xc = 30
    # nxh_ = np.int(np.minimum(np.int(n / 2), np.minimum(nx - xc, xc)))
    # print(n, xc, nxh_)
    # th_field, r_field = compute_radius(xc, xc, nxh_, nxh_, n, n, path_out_figs)   # ic not at center; part of domain
    print(nx, xc[0], irange)
    th_field, r_field = compute_radius(ic[0], jc[0], irange, jrange, nx, nx, path_out_figs)

    # interpolate horizontal velocity field to cell centre in subdomain
    v_hor_int = compute_v_hor_int(u_in, v_in, ic, jc, irange, jrange, nx, ny, t_ini, t_end, tau, path_out_figs)


    # v_rad_int = np.zeros((tau, nx, ny))
    # v_tan_int = np.zeros((tau, nx, ny))
    v_rad_int = np.zeros((tau, nx_, ny_))
    v_tan_int = np.zeros((tau, nx_, ny_))
    it = 0
    for it, t0 in enumerate(lifetime):
        v_rad_int[it,:,:], v_tan_int[it,:,:] = compute_radial_vel(v_hor_int[:,it,:,:], th_field,
                                                                  irange, jrange, nx_, ny_, ic[it], jc[it],
                                                                  t0, path_out_figs)






    # '''' plot input fields '''
    # # plot_tracers(coordinates, cp_id, n_tracers, x_shift, t_shift, lifetime, path_out_figs)
    # path_fields = path_out_figs
    # plot_input_fields(coordinates, xc, yc, cp_id, n_tracers, x_shift, t_shift, lifetime,
    #                   'u', path_fields, path_out_figs)





    # var_list = ['s', 'w']
    # for it,t0 in enumerate(times):
    #     print('-plot time: '+str(t0))
    #     fig_name = 's_w' + '_t' + str(t0) + '_tracers.png'
    #     fig, axis = plt.subplots(1, 2, figsize=(11, 6), sharey='all')
    #     rootgrp = nc.Dataset(os.path.join(path_fields, str(t0)+'.nc'))
    #     for j, var_name in enumerate(var_list):
    #         print(var_name, j)
    #         var = rootgrp.groups['fields'].variables[var_name][:,:,k0]
    #         max = np.amax(var)
    #         if var_name in ['w', 'v_rad', 'v_tan']:
    #             min = -max
    #             cm_ = cm_bwr
    #         else:
    #             min = np.amin(var)
    #             cm_ = cm_hsv
    #         # axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_)
    #         axis[j].imshow(var.T, vmin=min, vmax=max, cmap=cm_, origin='lower')
    #         axis[j].set_title(var_name)
    #         axis[j].set_aspect('equal')
    #     rootgrp.close()
    #     for i in range(n_tracers):
    #         for j in range(len(var_list)):
    #             axis[j].plot(coordinates[it,i,0]+shift, coordinates[it,i,1]+shift, 'ok', markersize=2)
    #     plt.tight_layout()
    #     fig.savefig(os.path.join(path_out_figs, fig_name))
    #     plt.close(fig)
    #
    #
    # var_list = ['v_rad', 'v_tan']
    # rootgrp = nc.Dataset(os.path.join(path, 'fields_v_rad', 'v_rad.nc'))
    # for it, t0 in enumerate(times):
    #     print('-plot time: ' + str(t0))
    #     fig_name = 'v_rad_tan' + '_t' + str(t0) + '_tracers.png'
    #     fig, axis = plt.subplots(1, 2, figsize=(11, 6), sharey='all')
    #     for j, var_name in enumerate(var_list):
    #         print(var_name, j)
    #         var = rootgrp.variables[var_name][it, :, :, k0]
    #         max = np.amax(var)
    #         if var_name in ['w', 'v_rad', 'v_tan']:
    #             min = -max
    #         else:
    #             min = np.amin(var)
    #         axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_bwr)
    #         axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_bwr)
    #         axis[j].set_title(var_name)
    #         axis[j].set_aspect('equal')
    #     for i in range(n_tracers):
    #         for j in range(len(var_list)):
    #             axis[j].plot(coordinates[it,i,0]+shift, coordinates[it,i,1]+shift, 'ok', markersize=2)
    #     plt.tight_layout()
    #     fig.savefig(os.path.join(path_out_figs, fig_name))
    #     plt.close(fig)
    # rootgrp.close()
    #
    #
    #
    # var_list = ['u', 'v']
    # rootgrp = nc.Dataset(os.path.join(path_tracers, 'input', 'uv_alltimes.nc'))
    # for it, t0 in enumerate(times):
    #     print('-plot time: ' + str(t0))
    #     fig_name = 'uv_alltimes' + '_t' + str(t0) + '_tracers.png'
    #     fig, axis = plt.subplots(1, 2, figsize=(11, 6), sharey='all')
    #     for j, var_name in enumerate(var_list):
    #         print(var_name, j)
    #         var = rootgrp.variables[var_name][it, :, :]
    #         max = np.amax(var)
    #         if var_name in ['w', 'v_rad', 'v_tan']:
    #             min = -max
    #         else:
    #             min = np.amin(var)
    #         axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_bwr)
    #         axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_bwr)
    #         axis[j].set_title(var_name)
    #         axis[j].set_aspect('equal')
    #     for i in range(n_tracers):
    #         for j in range(len(var_list)):
    #             axis[j].plot(coordinates[it,i,0]+shift, coordinates[it,i,1]+shift, 'ok', markersize=2)
    #     plt.tight_layout()
    #     fig.savefig(os.path.join(path_out_figs, fig_name))
    #     plt.close(fig)
    # rootgrp.close()


    return
# ----------------------------------------------------------------------

def plot_tracers(coordinates, cp_id, n_tracers, x_shift, t_shift, lifetime, path_out_figs):
    print('plot tracers: ', n_tracers)
    lvls = np.linspace(0,1,n_tracers)
    for it, t0 in enumerate(lifetime):
        print('-plot time: ' + str(t0))
        fig_name = 'cp'+str(cp_id)+'_tracers_t' + str(t0) + '_tracers.png'
        fig, axis = plt.subplots(2, 2, figsize=(7, 6), sharey='none')
        ax0 = axis[0,0]
        ax1 = axis[0,1]
        ax2 = axis[1,0]
        ax0.set_aspect('equal')
        for i in range(n_tracers):
            ax0.plot(coordinates[it, i, 0] + shift, coordinates[it, i, 1] + shift,
                    'o', color=cm_hsv(lvls[i]), markeredgecolor='w', markersize=5)
            ax1.plot(i, coordinates[it, i, 1],
                     'o', color=cm_hsv(lvls[i]), markeredgecolor='w', markersize=5)
            ax2.plot(coordinates[it,i,0], i,
                     'o', color=cm_hsv(lvls[i]), markeredgecolor='w', markersize=5)
        plt.tight_layout()
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)
    return




def plot_input_fields(coordinates, xc, yc, cp_id, n_tracers, x_shift, t_shift, lifetime,
                      field_name, path_fields, path_out_figs):
    k0 = 0
    rootgrp = nc.Dataset(os.path.join(path_fields, 'input_u.nc'))
    var1 = rootgrp.variables['var1'][:, k0, :, :]
    rootgrp.close()
    rootgrp = nc.Dataset(os.path.join(path_fields, 'input_v.nc'))
    var2 = rootgrp.variables['var1'][:, k0, :, :]
    rootgrp.close()
    print('plot input fields: ')
    print(lifetime)

    for it, t0 in enumerate(lifetime):
        print('-plot time: ' + str(t0), it)
        fig_name = 'uv_input_fields' + '_t' + str(t0) + '_tracers.png'
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 5), sharey='all')
        max = np.amax(var1[t0,:,:])
        min = -max
        lvls = np.linspace(min, max, 10)
        cf = ax0.contourf(var1[t0,:,:], levels=lvls, cmap=cm_bwr)
        ax0.set_title('u')
        plt.colorbar(cf, ax=ax0)
        max = np.amax(var2[t0,:,:])
        min = -max
        lvls = np.linspace(min, max, 10)
        cf = ax1.contourf(var2[t0,:,:], levels=lvls, cmap=cm_bwr)
        ax1.set_title('v')
        plt.colorbar(cf, ax=ax1)
        # # ax.contour(var2.T, levels=np.linspace(min, max, 2e1), colors='k', linewidth=0.5)
        # ax.contour(var2.T, colors='k', linewidths=0.5)
        # ax.set_title('u and v')

        # print(coordinates[it,:,0])
        print(np.amin(coordinates[it,:,0]), np.amax(coordinates[it,:,0]))
        print(np.amin(coordinates[it,:,1]), np.amax(coordinates[it,:,1]))
        for i in range(n_tracers):
            ax0.plot(coordinates[it, i, 0] + x_shift, coordinates[it, i, 1] + x_shift, 'ok', markersize=2)
            ax1.plot(coordinates[it, i, 0] + x_shift, coordinates[it, i, 1] + x_shift, 'ok', markersize=2)
        for ax in [ax0, ax1]:
            ax.set_aspect('equal')
            # ax.set_xlim(np.amin(coordinates[it, :, 0] - 1), np.amax(coordinates[it, :, 0] + 1))
            # ax.set_ylim(np.amin(coordinates[it, :, 1] - 1), np.amax(coordinates[it, :, 1] + 1))
            ax.set_xlim(630, 670)
            ax.set_ylim(580, 620)
            ax.plot(xc[it],yc[it], 'xk', markersize=20)
        plt.suptitle('t='+str(t0)+'s, t-shift '+str(t_shift)+', x-shift '+str(x_shift))
        plt.tight_layout()
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)

    return



# ------------------------- CP STATISTICS ---------------------------------------------


def get_tracer_coords(cp_id, n_cps, n_tracers, times, fullpath_in):
    print('get tracer coordinates')
    f = open(fullpath_in, 'r')
    lines = f.readlines()

    nt = len(times)
    coords = np.zeros((nt, n_tracers, 2))

    count_start = 0
    for it,t0 in enumerate(times):
        print('----t0='+str(t0), it, '----', count_start)
        timestep = int(lines[count_start].split()[0])
        age = int(lines[count_start].split()[1])
        cp_ID = int(lines[count_start].split()[3])
        while (cp_ID < cp_id):
            count_start += 1
            cp_ID = int(lines[count_start].split()[3])
        # print('t0 - cp_ID:', t0, timestep, cp_ID, count_start)
        while (timestep < t0):
            count_start += 1
            timestep = int(lines[count_start].split()[0])
        age = int(lines[count_start].split()[1])
        cp_ID = int(lines[count_start].split()[3])
        # print('t_end - timestep: ', t_ini, timestep, cp_ID, count_start)

        count = count_start
        i = 0
        while (timestep == t0 and age == it and cp_ID == cp_id):
            columns = lines[count].split()
            coords[it,i,0] = float(columns[4])
            coords[it,i,1] = float(columns[5])
            i += 1
            count += 1
            timestep = int(lines[count].split()[0])
            age = int(lines[count].split()[1])
            cp_ID = int(lines[count].split()[3])

    f.close()
    print('')
    return coords


def get_number_cps(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    cp_number = int(lines[-1].split()[3])
    f.close()

    return cp_number


def get_number_tracers(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    count = 0
    cp_age = int(lines[count].split()[0])
    cp_ID = int(lines[count].split()[3])
    # print('cp_age', cp_age)
    # print('cp_ID', cp_ID)
    age = cp_age
    ID = cp_ID
    while (age == cp_age and ID == cp_ID):
        count += 1
        age = int(lines[count].split()[0])
        ID = int(lines[count].split()[3])
    n_tracers = count
    f.close()

    return n_tracers


def get_cp_lifetime(cp_ID, fullpath_in):
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    count = 0
    ID = int(lines[count].split()[3])
    while (ID < cp_ID):
        count += 1
        ID = int(lines[count].split()[3])
    t0 = int(lines[count].split()[0])
    ID = int(lines[count].split()[3])

    while (ID == cp_ID):
        count += 1
        ID = int(lines[count].split()[3])
    t1 = int(lines[count-1].split()[0])
    tau = t1-t0+1


    f.close()
    return tau, t0, t1


def get_cp_center(cp_ID, tau, n_tracers, fullpath_in):
    xc = np.zeros(tau)
    yc = np.zeros(tau)

    f = open(fullpath_in, 'r')
    lines = f.readlines()
    count = 0
    ID = int(lines[count].split()[3])
    while (ID < cp_ID):
        count += 1
        ID = int(lines[count].split()[3])
    count_start = count

    for it in range(tau):
        xc[it] = float(lines[count_start+it*n_tracers].split()[14])
        yc[it] = float(lines[count_start+it*n_tracers].split()[15])

    f.close()
    # print('xc: ', xc)
    # print('yc: ', yc)
    return xc, yc



# def get_radius_vel(fullpath_in, t0, cp_id, n_tracers, n_cps):
#     f = open(fullpath_in, 'r')
#     # f = open(DIR+EXPID+'/'+child+'/output/irt_tracks_output_pure_sort.txt', 'r')
#     lines = f.readlines()
#     count = 0
#     dist = []
#     vel = []
#
#     count = t0 * n_cps * n_tracers + (cp_id - 1)*n_tracers
#     # while CP age is 0 and CP ID is cp_id
#     timestep = int(lines[count].split()[0])
#     cp_ID = int(lines[count].split()[3])
#     while (timestep-1 == t0 and int(lines[count].split()[3])==cp_id):
#         columns = lines[count].split()
#         dist.append(float(columns[8]))
#         # vel.append(np.sqrt(float(columns[10])**2 + float(columns[11])**2))
#         vel.append(float(columns[12]))
#         count += 1
#         timestep = int(lines[count].split()[0])
#     f.close()
#     r_av = np.average(dist)
#     vel_av = np.average(vel)
#
#     return r_av, vel_av


# ----------------------------------------------------------------------
def compute_radius(ic, jc, irange, jrange, nx, ny, path_out_figs):
    r_field = np.zeros((nx, ny), dtype=np.int)          # radius
    th_field = np.zeros((nx, ny), dtype=np.double)      # angle
    for i in range(irange):
        for j in range(jrange):
            r_field[ic+i, jc+j] = np.round(np.sqrt(i**2+j**2))
            r_field[ic-i, jc+j] = r_field[ic+i,jc+j]
            r_field[ic-i, jc-j] = r_field[ic+i,jc+j]
            r_field[ic+i, jc-j] = r_field[ic+i,jc+j]
            if i == 0:
                # i = 1e-9
                aux = np.arctan(np.double(j)/1e-9)
            else:
                aux = np.arctan(np.double(j)/i)
            # print('i,j', i, j, aux, np.pi-aux, np.pi + aux, 2*np.pi - aux)
            th_field[ic+i, jc+j] = aux
            th_field[ic-i, jc+j] = np.pi - aux
            th_field[ic-i, jc-j] = np.pi + aux
            th_field[ic+i, jc-j] = 2*np.pi - aux

    fig, axis = plt.subplots(3, 2)
    # cf = ax1.imshow(r_field[ic - irange:ic + irange + 1, jc - jrange:jc + jrange + 1].T, origin='lower')
    for ax in axis[:,0].flat:
        cf = ax.imshow(r_field.T, origin='lower')
        plt.colorbar(cf, ax=ax)
        ax.set_title('r(x,y)')
    for ax in axis[:,1].flat:
        cf = ax.imshow(th_field.T, origin='lower')
        # cf = ax2.imshow(th_field[ic - irange:ic + irange + 1, jc - jrange:jc + jrange + 1].T, origin='lower')
        plt.colorbar(cf, ax=ax)
        ax.set_title('th(x,y)')

    for ax in axis[1:,:].flat:
        ax.plot([ic,ic], [0, ny], 'k', linewidth=1)
        ax.plot([0, nx], [jc,jc], 'k', linewidth=1)
    # ax4.plot([ic,ic], [0, ny], 'k')
    # ax4.plot([0, nx], [jc,jc], 'k')

    for ax in axis[:2,:].flat:
        ax.set_xlim(0, nx)
        ax.set_ylim(0, ny)
    for ax in axis[2,:].flat:
        ax.set_xlim(ic-irange-5, ic+irange+5)
        ax.set_ylim(jc-jrange-5, jc+jrange+5)

    plt.suptitle('ic,jc='+str(ic) + ', '+str(jc) + ', (nxh = '+str(irange)+')')
    plt.savefig(os.path.join(path_out_figs, 'r_th_field.png'))
    plt.close()
    return th_field, r_field


def compute_v_hor_int(u_in, v_in, xc, yc, irange, jrange, nx, ny, t0, t1, tau, path_out_figs):
    v_hor_int = np.zeros((2, tau, nx, ny))
    print('')
    print('compute v hor')
    print('v_hor_int: ', v_hor_int.shape)
    imin = xc[0] - irange
    imax = xc[0] + irange
    jmin = yc[0] - jrange
    jmax = yc[0] + jrange
    for i in range(irange):
        v_hor_int[0, :, xc[0] + i, jmin:jmax] = 0.5 * (
                    u_in[t0:t1 + 1, xc[0] + i, jmin:jmax] + u_in[t0:t1 + 1, xc[0] + i - 1, jmin:jmax])
        v_hor_int[0, :, xc[0] - i, jmin:jmax] = 0.5 * (
                    u_in[t0:t1 + 1, xc[0] - i, jmin:jmax] + u_in[t0:t1 + 1, xc[0] - i - 1, jmin:jmax])
    for j in range(jrange):
        v_hor_int[1, :, imin:imax, yc[0] + j] = 0.5 * (
                    v_in[t0:t1 + 1, imin:imax, yc[0] + j] + v_in[t0:t1 + 1, imin:imax, yc[0] + j - 1])
        v_hor_int[1, :, imin:imax, yc[0] - j] = 0.5 * (
                    v_in[t0:t1 + 1, imin:imax, yc[0] - j] + v_in[t0:t1 + 1, imin:imax, yc[0] - j - 1])

    # for i in range(1, nx - 1):
    #     v_hor_int[0,:,i,:] = 0.5 * (u_in[t0:t1+1,i,:] + u_in[t0:t1+1,i-1,:])
    # for j in range(1, ny - 1):
    #     v_hor_int[1,:,:,j] = 0.5 * (v_in[t0:t1+1,:,j] + v_in[t0:t1+1,:,j-1])

    for it, t0_ in enumerate(np.arange(t0, t1+1)):
        fig, axis = plt.subplots(2, 2, figsize=(10, 10))
        ax11 = axis[0, 0]
        ax12 = axis[0, 1]
        ax21 = axis[1, 0]
        ax22 = axis[1, 1]
        cf = ax11.imshow(v_hor_int[0, 0, :, :].T, origin='lower')
        plt.colorbar(cf, ax=ax11, shrink=0.8)
        ax11.plot([xc[0], xc[0]], [0, ny], 'k')
        ax11.plot([xc[0]+irange, xc[0]+irange], [0, ny], 'k')
        ax11.plot([xc[0]-irange, xc[0]-irange], [0, ny], 'k')
        ax11.set_title('u hor int')
        cf = ax12.imshow(v_hor_int[1, 0, :, :].T, origin='lower')
        plt.colorbar(cf, ax=ax12, shrink=0.8)
        ax12.plot([0, nx], [yc[0], yc[0]], 'k')
        ax12.plot([0, nx], [yc[0]+jrange, yc[0]+jrange], 'k')
        ax12.plot([0, nx], [yc[0]-jrange, yc[0]-jrange], 'k')
        ax12.set_title('v hor int')
        cf = ax21.contourf(v_hor_int[0, 0, :, :].T)
        plt.colorbar(cf, ax=ax21, shrink=0.8)
        ax21.set_title('u hor int')
        cf = ax22.contourf(v_hor_int[1, 0, :, :].T)
        plt.colorbar(cf, ax=ax22, shrink=0.8)
        ax22.set_title('v hor int')
        for ax in axis[0,:].flat:
            ax.set_xlim(0, nx)
            ax.set_ylim(0, ny)
        for ax in axis[1,:].flat:
            ax.set_xlim(xc[0]-irange-10, xc[0]+irange+10)
            ax.set_ylim(yc[0]-jrange-10, yc[0]+jrange+10)
        plt.tight_layout()
        plt.savefig(os.path.join(path_out_figs, 'test_field_vhor_int_t' + str(t0_) + '.png'))
        plt.close()
    return v_hor_int


def compute_radial_vel(uv, th_field, irange, jrange, nx, ny, ic, jc, t0, path_out_figs):
    ur = np.zeros((nx,ny), dtype=np.double)
    utan = np.zeros((nx,ny), dtype=np.double)
    print('compute radial velocity: ')
    print(uv.shape, th_field.shape)
    for i in range(nx):
        for j in range(ny):
            ii = ic-np.int(nx/2)+i
            jj = jc-np.int(ny/2)+j
            th = th_field[ii,jj]
            # # # clockwise rotation
            # # ur[i,j] = uv[0,i,j]*np.cos(th) + uv[1,i,j]*np.sin(th)
            # counter-clockwise rotation
            ur[i,j] = uv[0,ii,jj]*np.cos(th) + uv[1,ii,jj]*np.sin(th)
            utan[i,j] = -uv[0,ii,jj]*np.sin(th) + uv[1,ii,jj]*np.cos(th)


    print('ij', i, ii, ic-np.int(nx/2), j, jj, jc-np.int(ny/2))
    print(nx, np.int(nx/2))
    print('ur, utan ', np.amin(ur), np.amax(ur), np.amin(utan), np.amax(utan))
    # print(np.amax(uv[0,:,:]), np.amax(uv[1,:,:]))
    # print()




    fig, axis = plt.subplots(2, 2, figsize=(10, 15))
    ax11 = axis[0, 0]
    ax12 = axis[0, 1]
    ax21 = axis[1, 0]
    ax22 = axis[1, 1]
    cf = ax11.imshow(uv[0, :, :].T, origin='lower')
    plt.colorbar(cf, ax=ax11, shrink=0.5)
    ax11.set_title('u')
    cf = ax12.imshow(uv[1, :, :].T, origin='lower')
    plt.colorbar(cf, ax=ax12, shrink=0.5)
    ax12.set_title('v')
    # circle1 = plt.Circle((xc[0], yc[0]), rmax / 2, fill=False, color='r', linewidth=2)
    cf = ax21.imshow(ur[:, :].T, origin='lower')
    plt.colorbar(cf, ax=ax21, shrink=0.5)
    ax21.set_title('radial velocity')
    # ax11.add_artist(circle1)
    # ax11.plot(xc[0] - 0.5, jc - 0.5, 'ow', markersize=7)
    cf = ax22.imshow(utan[:, :].T, origin='lower')
    plt.colorbar(cf, ax=ax22, shrink=0.5)
    ax22.set_title('tangential velocity')
    for ax in axis[0, :].flat:
        ax.set_xlim(ic - irange, ic + irange)
        ax.set_ylim(jc - jrange, jc + jrange)
    plt.tight_layout()
    plt.savefig(os.path.join(path_out_figs, 'test_field_vrad_vtan_t' + str(t0) + '.png'))
    plt.close()

    return ur, utan
# ----------------------------------------------------------------------

def set_input_parameters(args):
    print('--- set input parameters ---')
    global path, path_fields, case_name
    # path = args.path
    # path_fields = os.path.join(path, 'fields')
    # case_name = args.casename

    # nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    # global nx, ny, nz, dx, dV, gw
    # nx = nml['grid']['nx']
    # ny = nml['grid']['ny']
    # nz = nml['grid']['nz']
    # dx = np.zeros(3, dtype=np.int)
    # dx[0] = nml['grid']['dx']
    # dx[1] = nml['grid']['dy']
    # dx[2] = nml['grid']['dz']
    # gw = nml['grid']['gw']
    # dV = dx[0] * dx[1] * dx[2]
    #
    # global tmin, tmax
    # if args.tmin:
    #     tmin = np.int(args.tmin)
    # else:
    #     tmin = 100
    # if args.tmax:
    #     tmax = np.int(args.tmax)
    # else:
    #     tmax = 100
    #
    # ''' time range '''
    # global times, files
    # times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
    #          and tmin <= np.int(name[:-3]) <= tmax]
    # times.sort()
    # files = [str(t) + '.nc' for t in times]
    # print('tmin, tmax: ', tmin, tmax)
    # print('times: ', times)
    # print('')

    return

# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()
