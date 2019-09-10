import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os

import plot_CP_singleCP
from plot_CP_singleCP import get_radius_vel

label_size = 12
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['text.usetex'] = 'true'

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    # parser.add_argument("--k0")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--timerange", nargs='+', type=int)
    parser.add_argument("--dx")
    parser.add_argument("--imin")
    parser.add_argument("--shift_x")
    parser.add_argument("--shift_t")
    args = parser.parse_args()

    k0 = 0
    if args.dx:
        dx = np.int(args.dx)
    else:
        dx = 50

    case_name = 'ColdPoolDry_single_3D'
    if dx == 25:
        path_data = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run4_dx25m/dTh3_z1000_r1000'
        path_tracer_file = os.path.join(path_data, 'tracer_k' + str(k0) + '/output/', 'coldpool_tracer_out.txt')
    elif dx == 50:
        path_data = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run3_dx50m/dTh3_z1000_r1000'
        path_tracer_file = os.path.join(path_data, 'tracer_k' + str(k0) + '_nointerpol/output/', 'coldpool_tracer_out.txt')
    path_fields = os.path.join(path_data, 'fields')
    path_fields_merged = os.path.join(path_data, 'fields_merged')
    # path_out_figs = '/Users/bettinameyer/Dropbox/ClimatePhysics/Copenhagen/Projects/ColdPool_LES/data_olga/figs_2D_1CP'
    path_out_figs = os.path.join('/nbi/home/meyerbe/paper_olga/figs_2D_1CP_dx'+str(dx)+'m_nointerpol/')
    # path_out_figs = os.path.join('/nbi/home/meyerbe/paper_olga/figs_2D_1CP_dx'+str(dx)+'m_interpol/')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print('path figures: ' + path_out_figs)
    print('path tracers: ' + path_tracer_file)
    print('')
    set_input_parameters(case_name, path_data)



    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gray')
    cm_grey_r = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')
    cm_winter = plt.cm.get_cmap('winter')
    cm_summer = plt.cm.get_cmap('summer')
    # cm_fall = plt.cm.get_cmap('fall')

    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 0
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 2500
    times = np.arange(tmin, tmax, 100)
    if args.timerange:
        timerange = [np.int(t) for t in args.timerange]
    else:
        timerange = [0, 900, 2400]
    print('tmin: ', tmin)
    print('tmax: ', tmax)
    print('times: ' + str(times))
    print('timerange: ' + str(timerange))
    nml = simplejson.loads(open(os.path.join(path_data, case_name + '.in')).read())
    global dt_fields
    dt_fields = nml['fields_io']['frequency']

    ''' get tracer coordinates '''
    cp_id = 1
    n_cps = get_number_cps(path_tracer_file)
    n_tracers = get_number_tracers(path_tracer_file)
    print('CP ID: ' + str(cp_id))
    print('number of CPs: ' + str(n_cps))
    print('number of tracers: ' + str(n_tracers))
    print('')
    coordinates = get_tracer_coords(cp_id, n_cps, n_tracers, times, dt_fields, path_tracer_file)

    ''' get 2D field '''
    var_name = 'w'
    # print('')
    # print('!!!! old merged fields file !!!')
    # print('')
    root = nc.Dataset(os.path.join(path_fields_merged, 'fields_allt_xy_k'+str(k0)+'.nc'))
    # root = nc.Dataset(os.path.join(path_fields_merged, 'fields_allt_xy_k'+str(k0)+'_v1.nc'))
    var = root.variables[var_name][:,:,:]
    root.close()
    root = nc.Dataset(os.path.join(path_data, 'fields_v_rad', 'v_rad.nc'))
    v_rad = root.variables['v_rad'][:, :, :, k0]
    root.close()

    # parameters plotting
    itmax = np.int(tmax/dt_fields)
    print('itmax: ', itmax)
    max_var = np.amax(var[:itmax,:,:])
    if var_name == 'w':
        min_var = -max_var
    else:
        min_var = np.amin(var[:itmax,:,:])
    lvls_var = np.linspace(min_var, max_var, 25)
    # lvls_var = np.linspace(min_var, max_var, 10)
    max = np.amax(np.abs(v_rad[:,:,:]))
    lvls_v_rad = np.linspace(-max, max, 25)
    # lvls_v_rad = np.linspace(-max, max, 10)
    colmap = cm_hsv
    colmap = cm_bwr
    print('min max: ', min_var, max_var)

    if args.imin:
        imin = np.int(args.imin)
    else:
        if dx == 25:
            imin = 100
        elif dx == 50:
            imin = 0
    imax = nx-imin
    print('imin, imax', imin, imax)

    if args.shift_x:
        shift = np.double(args.shift_x)
    elif dx == 50:
        shift = -0.5
    else:
        shift = 0
    if args.shift_t:
        shift_t = np.int(args.shift_t)
    elif dx == 50:
        shift_t = 1
    else:
        shift_t = 1
    print('x-shift: '+str(shift))
    print('t-shift: '+str(shift_t))
    print('')


    '''------- plot timerange -------'''
    fig_name = 'timerange_t' + str(timerange) + '_k' + str(k0) + '_tshift'+str(shift_t)+'.png'
    textprops = dict(facecolor='white', alpha=0.5, linewidth=0.)
    label_list = ['a)', 'b)', 'c)']
    fig, axes = plt.subplots(nrows=2, ncols=len(timerange), figsize=(12,7), sharey='row', sharex='col')
    for i,t0 in enumerate(timerange):
        it = np.int(t0 / dt_fields)
        ax = axes[0,i]
        im0 = ax.contourf(var[it, :, :].T, cmap=colmap, levels=lvls_var)
        textstr = label_list[i] + ' t='+str(t0)+'s'
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=textprops)
        ax = axes[1,i]
        im1 = ax.contourf(v_rad[it, :, :].T, cmap=colmap, levels=lvls_v_rad)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=textprops)
        for ax in axes[:,i].flat:
            for i in range(n_tracers):
                ax.plot(coordinates[it + shift_t, i, 0] + shift, coordinates[it + shift_t, i, 1] + shift, 'ok', markersize=1)
    for ax in axes.flat:
        # ax.set_aspect('equal')
        ax.set_xlim(imin,imax)
        ax.set_ylim(imin,imax)
    for ax in axes[1,:].flat:
        ax.set_xlabel('Distance / km')
        print('xticks', ax.get_xticks())
        x_ticks = [np.int((ti-imin) * dx * 1e-3) for ti in ax.get_xticks()]
        # x_ticks = ax.get_xticks()
        # x_ticks = x_ticks * dx *1e-3
        print('      ', x_ticks)
        ax.set_xticklabels(x_ticks)
        for label in ax.xaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
    for ax in axes[:,0].flat:
        ax.set_ylabel('Distance / km')
        y_ticks = [np.int((ti-imin) * dx * 1e-3) for ti in ax.get_yticks()]
        ax.set_yticklabels(y_ticks)
        for label in ax.yaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
    # for ax in axes.flat:
    #     ax.set_aspect('equal')
    cbar_ax0 = fig.add_axes([0.9, 0.58, 0.015, 0.3])
    cbar_ax1 = fig.add_axes([0.9, 0.14, 0.015, 0.3])
    cb0 = fig.colorbar(im0, cax=cbar_ax0, ticks=np.arange(np.floor(min_var), np.floor(max_var) + 1, 1))
    cb1 = fig.colorbar(im1, cax=cbar_ax1, ticks=[-10,-5,0,5,10])
    cb0.set_label('vertical velocity / ms' + r'$^{-1}$')
    cb1.set_label('radial velocity / ms'+r'$^{-1}$')
    fig.subplots_adjust(top=0.95, bottom=0.07, left=0.07, right=0.87, hspace=0.1, wspace=0.1)
    plt.savefig(os.path.join(path_out_figs, fig_name))


    colmap = cm_grey_r
    # tracercolor = 'royalblue'
    # tracercolor = 'blue'
    tracercolor = 'mediumblue'
    fig_name = 'timerange2_t' + str(timerange) + '_k' + str(k0) + '_tshift'+str(shift_t)+ '.png'
    textprops = dict(boxstyle='round', facecolor='white', alpha=0.5, linewidth=0.)
    fig, axes = plt.subplots(nrows=2, ncols=len(timerange), figsize=(11, 7), sharey='row', sharex='col')
    for i, t0 in enumerate(timerange):
        it = np.int(t0 / dt_fields)
        ax = axes[0, i]
        im0 = ax.contourf(var[it, :, :].T, cmap=colmap, levels=lvls_var)
        textstr = 't=' + str(t0) + 's'
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=textprops)
        ax = axes[1, i]
        im1 = ax.contourf(v_rad[it, :, :].T, cmap=colmap, levels=lvls_v_rad)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=textprops)
        for ax in axes[:, i].flat:
            for i in range(n_tracers):
                ax.plot(coordinates[it + shift_t, i, 0] + shift, coordinates[it + shift_t, i, 1] + shift, 'o',
                        markerfacecolor=tracercolor, markeredgecolor=tracercolor, markersize=1)
    for ax in axes.flat:
        # ax.set_aspect('equal')
        ax.set_xlim(imin, imax)
        ax.set_ylim(imin, imax)
    for ax in axes[1,:].flat:
        ax.set_xlabel('Distance / km')
        print('xticks', ax.get_xticks())
        x_ticks = [np.int((ti-imin) * dx * 1e-3) for ti in ax.get_xticks()]
        print('      ', x_ticks)
        ax.set_xticklabels(x_ticks)
        for label in ax.xaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
    for ax in axes[:,0].flat:
        ax.set_ylabel('Distance / km')
        y_ticks = [np.int((ti-imin) * dx * 1e-3) for ti in ax.get_yticks()]
        ax.set_yticklabels(y_ticks)
        for label in ax.yaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
    # for ax in axes.flat:
    #     # ax.set_aspect('equal')
    #     ax.set_xlim(imin,imax)
    #     ax.set_ylim(imin,imax)
    cbar_ax0 = fig.add_axes([0.92, 0.60, 0.015, 0.3])
    cbar_ax1 = fig.add_axes([0.92, 0.14, 0.015, 0.3])
    cb0 = fig.colorbar(im0, cax=cbar_ax0, ticks=np.arange(np.floor(min_var), np.floor(max_var) + 1, 1))
    cb1 = fig.colorbar(im1, cax=cbar_ax1, ticks=[-10, -5, 0, 5, 10])
    cb0.set_label('vertical velocity / ms' + r'$^{-1}$')
    cb1.set_label('radial velocity / ms' + r'$^{-1}$')
    # fig.subplots_adjust(top=0.95, bottom=0.07, left=0.07, right=0.87, hspace=0.1, wspace=0.1)
    fig.subplots_adjust(top=0.95, bottom=0.1, left=0.07, right=0.90, hspace=0.1, wspace=0.1)
    plt.savefig(os.path.join(path_out_figs, fig_name))




    # '''------- plot each timestep -------'''
    # plot_each_timestep(coordinates, n_tracers, var_name, var, v_rad, lvls_var, lvls_v_rad,
    #                    k0, shift_t, shift, times, path_out_figs, colmap)
    #
    #
    #
    # '''------- plot each timestep incl. circle of average radius -------'''
    # plot_each_timestep_circles(coordinates, cp_id, n_tracers, n_cps,
    #                            var_name, var, v_rad, lvls_var, lvls_v_rad,
    #                            k0, shift_t, shift, dx, times, path_data, path_tracer_file, path_out_figs, colmap)


    return
# ----------------------------------------------------------------------

def plot_each_timestep(coordinates, n_tracers, var_name, var, v_rad, lvls_var, lvls_v_rad,
                       k0, shift_t, shift, times, path_out_figs, colmap):
    '''------- plot each timestep -------'''
    for t0 in times[:-1]:
        it = np.int(t0 / dt_fields)
        print('plot: t='+str(t0) + ', '+str(it))
        fig_name = var_name + '_t' + str(np.int(t0)) + '_k'+str(k0) + '.png'
        fig, (ax0, ax1) = plt.subplots(1,2, figsize=(20,10))
        cf = ax0.contourf(var[it,:,:].T, cmap=colmap, levels=lvls_var)
        plt.colorbar(cf, ax=ax0, shrink=0.5)
        cf = ax1.contourf(v_rad[it,:,:].T, cmap=colmap, levels=lvls_v_rad)
        plt.colorbar(cf, ax=ax1, shrink=0.5)
        for ax in [ax0, ax1]:
            ax.set_aspect('equal')
            for i in range(n_tracers):
                ax.plot(coordinates[it+shift_t, i, 0] + shift, coordinates[it+shift_t, i, 1] + shift, 'ok', markersize=2)
        ax0.set_title('w')
        ax1.set_title('v_{rad}')
        plt.suptitle('t-shift '+str(shift_t)+', x-shift '+str(shift))
        plt.tight_layout()
        plt.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)



def plot_each_timestep_circles(coordinates, cp_id, n_tracers, n_cps,
                               var_name, var, v_rad, lvls_var, lvls_v_rad,
                               k0, shift_t, shift, dx, times, path_data, path_tracer_file, path_out_figs, colmap):
    '''------- plot each timestep incl. circle of average radius -------'''
    ''' a) read in CP radius r_torus(t)=r_av and CP spreading velocity u_torus(t)=U_rad_av from tracer statistics '''
    print('path tracers: ', path_tracer_file)
    nt = len(times)
    nk = 1
    dist_av = np.zeros((nt, nk))
    U_rad_av = np.zeros((nt, nk))
    for it, t0 in enumerate(times):
        # print('---t0: ' + str(t0) + '---', it)
        dist_av[it, k0], U_rad_av[it, k0] = get_radius_vel(path_tracer_file, it, cp_id, n_tracers, n_cps)
    r_tracers_av = dist_av * dx
    del dist_av

    ''' b) get azimuthally averaged fields '''
    rad_stats_file = nc.Dataset(os.path.join(path_data, 'data_analysis', 'stats_radial_averaged.nc'))
    v_rad_av = rad_stats_file.groups['stats'].variables['v_rad'][:, :, :]
    r_av = rad_stats_file.groups['stats'].variables['r'][:]  # nr
    rad_stats_file.close()
    v_rad_gradient = np.zeros(shape=v_rad_av.shape)
    for ir, r in enumerate(r_av[1:-1]):
        v_rad_gradient[:, ir, :] = (v_rad_av[:, ir + 1, :] - v_rad_av[:, ir - 1, :]) / (r_av[ir + 1] - r_av[ir - 1])
    ir_vrad_grad_max = np.argmin(v_rad_gradient, axis=1)  # nt, nk

    print ''
    for t0 in times[:-1]:
        if t0 < 1000:
            imin = 100
        elif t0 < 2000:
            imin = 50
        else:
            imin = 0
        it = np.int(t0 / dt_fields)
        print('plot: t=' + str(t0) + ', ' + str(it))
        fig_name = 'circle_' + var_name + '_t' + str(np.int(t0)) + '_k' + str(k0) + '.png'
        fig, axis = plt.subplots(2, 2, figsize=(12, 10))
        ax11 = axis[0, 0]
        ax12 = axis[0, 1]
        ax21 = axis[1, 0]
        ax22 = axis[1, 1]
        radius_tracers = r_tracers_av[it+shift_t, k0] / dx
        radius_gradient  = ir_vrad_grad_max[it,k0]
        circle1_grad = plt.Circle((ic+0.5, jc+0.5), radius_gradient, fill=False, color='b', linewidth=2)
        circle2_grad = plt.Circle((ic+0.5, jc+0.5), radius_gradient, fill=False, color='b', linewidth=2)
        circle1_tr = plt.Circle((ic+0.5, jc+0.5), radius_tracers, fill=False, color='g', linewidth=3)
        circle2_tr = plt.Circle((ic+0.5, jc+0.5), radius_tracers, fill=False, color='g', linewidth=3)
        circle3 = plt.Circle((ic, jc), radius_tracers, fill=False, color='g', linewidth=3)
        circle4 = plt.Circle((ic, jc), radius_tracers, fill=False, color='g', linewidth=3)
        cf = ax11.contourf(var[it, :, :].T, cmap=colmap, levels=lvls_var)
        plt.colorbar(cf, ax=ax11, shrink=0.5)
        cf = ax12.contourf(v_rad[it, :, :].T, cmap=colmap, levels=lvls_v_rad)
        plt.colorbar(cf, ax=ax12, shrink=0.5)
        ax11.set_title('vertical velocity')
        ax11.add_artist(circle1_grad)
        ax11.add_artist(circle1_tr)
        ax12.add_artist(circle2_grad)
        ax12.add_artist(circle2_tr)
        ax11.plot(ic - 0.5, jc - 0.5, 'ow', markersize=7)
        ax12.set_title('radial velocity')
        cf = ax21.contourf(var[it, :, :].T, cmap=colmap, levels=lvls_var)
        plt.colorbar(cf, ax=ax21, shrink=0.5)
        cf = ax22.contourf(v_rad[it, :, :].T, cmap=colmap, levels=lvls_v_rad)
        plt.colorbar(cf, ax=ax22, shrink=0.5)

        for ax in [ax11, ax12, ax21, ax22]:
            # ax.set_aspect('equal')
            for i in range(n_tracers):
                ax.plot(coordinates[it + shift_t, i, 0] + shift, coordinates[it + shift_t, i, 1] + shift, 'ok',
                        markersize=1)

        ax21.add_artist(circle3)
        ax22.add_artist(circle4)
        ax22.set_title('v')
        for ax in axis.flat:
            ax.set_xlim(imin, nx-imin)
            ax.set_ylim(imin, nx-imin)
        plt.tight_layout()
        plt.suptitle('t-shift ' + str(shift_t) + ', x-shift ' + str(shift))
        plt.savefig(os.path.join(path_out_figs, fig_name))
        plt.close()

# ----------------------------------------------------------------------

def plot_tracer_input_fields(coordinates, n_tracers, shift, path_tracers, path_out_figs):
    print(os.path.join(path_tracers, 'input', 'uv_alltimes.nc'))
    rootgrp = nc.Dataset(os.path.join(path_tracers, 'input', 'uv_alltimes.nc'))
    for it, t0 in enumerate(times):
        if it > 0:
            print('-plot time: ' + str(t0))
            fig_name = 'uv_input_fields' + '_t' + str(t0) + '_tracers.png'
            fig, ax = plt.subplots(1, 1, figsize=(7, 6), sharey='all')
            var1 = rootgrp.variables['u'][it, :, :]
            var2 = rootgrp.variables['v'][it, :, :]
            max = np.amax(var1)
            min = -max
            # min = np.amin(var)
            ax.contourf(var1.T, levels=np.linspace(min, max, 1e2), cmap=cm_bwr)
            # ax.contour(var2.T, levels=np.linspace(min, max, 2e1), colors='k', linewidth=0.5)
            ax.contour(var2.T, colors='k', linewidths=0.5)
            ax.set_title('u and v')
            ax.set_aspect('equal')
            # for i in range(n_tracers):
            #     ax.plot(coordinates[it, i, 0] + shift, coordinates[it, i, 1] + shift, 'ok', markersize=2)
            plt.tight_layout()
            fig.savefig(os.path.join(path_out_figs, fig_name))
            plt.close(fig)
    rootgrp.close()

    return



# ----------------------------------------------------------------------


def get_tracer_coords(cp_id, n_cps, n_tracers, times, dt_fields, fullpath_in):
    print('get tracer coordinates')
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    column = lines[0].split()

    nt = len(times)
    coords = np.zeros((nt, n_tracers, 2))

    for it, t0 in enumerate(times):
        print('----t0='+str(t0), it, '----')
        i = 0
        # count = t0 * n_cps * n_tracers + (cp_id - 1) * n_tracers
        # count = it * n_cps * n_tracers + (cp_id - 1) * n_tracers
        count = np.int(t0/dt_fields) * n_cps * n_tracers + (cp_id - 1) * n_tracers
        # while CP age is 0 and CP ID is cp_id
        timestep = int(lines[count].split()[0])
        cp_ID = int(lines[count].split()[3])
        # while (timestep - 1 == it and cp_ID == cp_id):
        while (timestep - 1 == t0/dt_fields and cp_ID == cp_id):
            columns = lines[count].split()
            coords[it,i,0] = float(columns[4])
            coords[it,i,1] = float(columns[5])
            i += 1
            count += 1
            cp_ID = int(lines[count].split()[3])
            timestep = int(lines[count].split()[0])

    f.close()
    # print ''
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
# ----------------------------------------------------------------------

def set_input_parameters(case_name, path_data):
    nml = simplejson.loads(open(os.path.join(path_data, case_name + '.in')).read())
    global nx, ny
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    # nz = nml['grid']['nz']
    # dx = np.zeros(3, dtype=np.int)
    # dx[0] = nml['grid']['dx']
    # dx[1] = nml['grid']['dy']
    # dx[2] = nml['grid']['dz']
    # gw = nml['grid']['gw']
    # dV = dx[0] * dx[1] * dx[2]

    global ic, jc
    try:
        ic = nml['init']['ic']
        jc = nml['init']['jc']
    except:
        ic = np.int(nx/2)
        jc = np.int(ny/2)
    print('ic, jc', ic, jc)

    return

# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()
