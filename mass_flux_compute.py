import netCDF4 as nc
import argparse
import json as simplejson
import os
import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
import matplotlib.patches as patch
import time



# from convert_fields_smaller_k import convert_file_for_varlist_horsection

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
    parser.add_argument("--level")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()


    path, nml, times, files = set_input_output_parameters(args)
    tmin = times[0]
    tmax = times[-1]
    itmin = np.int(tmin / dt)
    itmax = np.int(tmax / dt)
    nt = len(times)
    ic_arr, jc_arr, ic, jc, ncp = define_geometry(nml)
    path_out_figs = os.path.join(path, 'figs_massflux')
    path_out_data = os.path.join(path, 'data_analysis')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    if not os.path.exists(path_out_data):
        os.mkdir(path_out_data)

    if args.level:
        zlevel = np.int(args.level)
    else:
        zlevel = 1000
    k0 = np.int(np.double(zlevel) / dx[2])
    print('-----------------------------------')
    print('considering mass flux at level: z=' + str(zlevel) + ', k=' + str(k0))
    print('')


    path_field_k0 = os.path.join(path, 'fields_merged', 'fields_allt_xy_k'+str(k0)+'.nc')
    if not os.path.exists(path_field_k0):
        print('merging file for k='+str(k0), times)
        convert_file_for_varlist_horsection(['w'], times, files, os.path.join(path, 'fields'), os.path.join(path, 'fields_merged'), k0)
    print('path: ', path_field_k0)
    field_k0 = nc.Dataset(path_field_k0, 'r')
    t_merged = field_k0.variables['time'][:]
    print('time merged field: ', t_merged, t_merged[0], t_merged[-1])
    if t_merged[0] > tmin or t_merged[-1] < tmax:
        print('removing file: ', os.path.join(path, 'fields_merged', 'fields_allt_xy_k'+str(k0)+'.nc'))
        os.remove(os.path.join(path, 'fields_merged', 'fields_allt_xy_k'+str(k0)+'.nc'))
        convert_file_for_varlist_horsection(['w'], times, files, os.path.join(path, 'fields'),
                                            os.path.join(path, 'fields_merged'), k0)
    else:
        print('merged file at z='+str(zlevel)+' already existing')
        print(path_field_k0)
    print('')



    ''' computing mass flux '''
    filename_out = 'mass_flux_z'+str(zlevel)+'.nc'
    create_output_file(filename_out, times, path_out_data)
    rho = np.zeros((nt), dtype=np.double)
    # read in density
    print('stats-file: ', os.path.join(path, 'stats', 'Stats.' + case_name + '.nc'))
    rootgrp = nc.Dataset(os.path.join(path, 'stats', 'Stats.' + case_name + '.nc'))
    time_stats = rootgrp.groups['profiles'].variables['t'][:]
    # rho0_stat = rootgrp.groups['reference'].variables['rho0_full'][:]
    rho_mean = 1./rootgrp.groups['profiles'].variables['alpha_mean'][:,k0]
    rootgrp.close()
    print('time stats: ', time_stats)
    if time_stats[-1] < tmax:
        print('Using restart stats-file')
        nt_ = len(time_stats)
        rho[:nt_] = rho_mean[:]
        try:
            rootgrp = nc.Dataset(os.path.join(path, 'stats', 'Stats.' + case_name + '.Restart_0.nc'))
            rho_mean_restart = 1. / rootgrp.groups['profiles'].variables['alpha_mean'][:, k0]
            rootgrp.close()
            rho[nt_:] = rho_mean_restart[:(nt-nt_)]
        except:
            tmax = np.int(time_stats[-1])
            print('no restart-file, changing tmax to '+str(tmax))
    else:
        rho[:] = rho_mean[:nt]
    # compute mass flux
    w_2D = field_k0.variables['w'][:,:,:]
    print('shapes: rho=', rho.shape, 'w_2D=',w_2D.shape, tmin, tmax)
    mf_w = np.sum(w_2D[:itmax+1,:,:], axis=0)
    flux = np.zeros((w_2D.shape), dtype=np.double)
    for it,t0 in enumerate(range(tmin, tmax+dt, dt)):
        # it_fields = t0/dt_fields
        # it_stats = t0/dt_stats
        flux[it,:,:] = rho[it]*w_2D[it,:,:]
    flux_acc = np.sum(flux, axis=0)
    flux_mean = np.mean(flux, axis=0)
    # only account for positive (i.e., upward) mass flux
    mf_w_pos = np.sum(np.where(w_2D[:itmax + 1, :, :] > 0, w_2D[:itmax + 1, :, :], 0), axis=0)
    flux_pos = np.where(w_2D>0, flux, 0)
    flux_pos_mean = np.mean(flux_pos, axis=0)
    # del flux
    field_k0.close()

    ''' output flux '''
    dump_output_file('mass_flux_2D', 'fields_2D', flux, filename_out, path_out_data)
    dump_output_file('mass_flux_2D_positive', 'fields_2D', flux_pos, filename_out, path_out_data)
    # dump_output_file('mass_flux_accumulated', 'fields_2D', flux_acc, filename_out, path_out_data)
    # dump_output_file('mass_flux_mean', 'fields_2D', flux_mean, filename_out, path_out_data)
    # dump_output_file('mass_flux_positive_accumulated', 'fields_2D', np.sum(flux_pos, axis=0), filename_out, path_out_data)
    # dump_output_file('mass_flux_positive_mean', 'fields_2D', flux_pos_mean, filename_out, path_out_data)


    ''' compute mean flux in box - single, double, triple '''
    dx_box = 20.
    dy_box = dx_box
    # d = jc_arr[1] - jc_arr[0]
    # mf_mean_single:       MF in box   x0=ic1-(r0 + itmax*v0)/dx, ..., x1=x0+20*dx; v0=1m/s
    #                                   y0=jc1-(r0 + itmax*v0)/dx,  ..., y1=y0+20*dx
    # mf_mean_double:       MF in box   x0=ic1-20*dx, ..., x1=ic1
    #                                   y0=jc1-20*dx, ..., y1=jc1
    # mf_mean_triple:       MF in box   x0=icoll-10*dx, ..., x1=icoll+10*dx
    #                                   y0=jcoll-10*dx, ..., y1=jcoll+10*dx
    x_single = ic_arr[0] - (3*rstar) / dx[0]
    y_single = jc_arr[0] - (3*rstar) / dx[1]
    mf_mean_single = np.mean(np.mean(flux[:, x_single:x_single + dx_box, y_single:y_single + dy_box], axis=1), axis=1)
    mf_mean_single_pos = np.mean(np.mean(flux_pos[:, x_single:x_single + dx_box, y_single:y_single + dy_box], axis=1),
                                 axis=1)
    mf_mean_single_pos_acc = [np.sum(mf_mean_single_pos[:it + 1]) for it in range(itmin, itmax + 1)]
    if ncp > 1:
        x_double = ic_arr[0] - dx_box / 2
        y_double = 0.5 * (jc_arr[0] + jc_arr[1]) - dy_box / 2
        x_triple = ic - dx_box / 2
        y_triple = jc - dy_box / 2

        mf_mean_double = np.mean(np.mean(flux[:,x_double:x_double+dx_box, y_double:y_double+dy_box], axis=1), axis=1)
        mf_mean_triple = np.mean(np.mean(flux[:,x_triple:x_triple+dx_box, y_triple:y_triple+dy_box], axis=1), axis=1)
        mf_mean_double_pos = np.mean(np.mean(flux_pos[:,x_double:x_double+dx_box, y_double:y_double+dy_box], axis=1), axis=1)
        mf_mean_triple_pos = np.mean(np.mean(flux_pos[:,x_triple:x_triple+dx_box, y_triple:y_triple+dy_box], axis=1), axis=1)
    else:
        mf_mean_double = np.zeros(itmax)
        mf_mean_triple = np.zeros(itmax)
        mf_mean_double_pos = np.zeros(itmax)
        mf_mean_triple_pos = np.zeros(itmax)


    # accumulated flux
    mf_mean_single_acc = [np.sum(mf_mean_single[:it+1]) for it in range(itmin,itmax+1)]
    if ncp > 1:
        mf_mean_double_acc = [np.sum(mf_mean_double[:it+1]) for it in range(itmin,itmax+1)]
        mf_mean_triple_acc = [np.sum(mf_mean_triple[:it+1]) for it in range(itmin,itmax+1)]
        mf_mean_double_pos_acc = [np.sum(mf_mean_double_pos[:it+1]) for it in range(itmin,itmax+1)]
        mf_mean_triple_pos_acc = [np.sum(mf_mean_triple_pos[:it+1]) for it in range(itmin,itmax+1)]
    else:
        mf_mean_double_acc = 0
        mf_mean_triple_acc = 0
        mf_mean_double_pos_acc = 0
        mf_mean_triple_pos_acc = 0

    # ''' output files '''
    # # create file
    #
    #
    ''' plotting '''
    # fig_name = 'massflux_k' + str(k0) + '_tmax' + str(tmax) + '.png'
    # plot_massflux_allversions(flux_mean, flux_pos_mean, mf_w, mf_w_pos, k0, path_out_figs, fig_name)
    fig_name = 'massflux_locations_k' + str(k0) + '_tmax' + str(tmax) + '.png'
    if ncp == 3:
        plot_massflux_locations(flux_mean, flux_pos_mean,
                            mf_mean_single, mf_mean_double, mf_mean_triple,
                            mf_mean_single_pos, mf_mean_double_pos, mf_mean_triple_pos,
                            mf_mean_single_acc, mf_mean_double_acc, mf_mean_triple_acc,
                            mf_mean_single_pos_acc, mf_mean_double_pos_acc, mf_mean_triple_pos_acc,
                            ic_arr, jc_arr, ic, jc,
                            x_single, y_single, x_double, y_double, x_triple, y_triple, dx_box, dy_box,
                            k0, t_merged, path_out_figs, fig_name)


    return

# ---------------------------------------------------------------------
def plot_massflux_allversions(mf_mean, mf_mean_pos, mf_w, mf_w_pos, k0, path_out_figs, fig_name):

    ncol = 4
    fig, axis = plt.subplots(2, ncol, figsize=(ncol * 5, 10), sharey='all')
    ax = axis[0, 0]
    cf = ax.contourf(mf_mean)
    plt.colorbar(cf, ax=ax)
    ax.set_title('(A) sum(w[z]*rho[z]) (z=' + str(k0 * dx[2]) + ')')
    ax = axis[0, 1]
    cf = ax.contourf(mf_mean_pos)
    plt.colorbar(cf, ax=ax)
    ax.set_title('(B) sum(w[z]*rho[z], w>0) (z=' + str(k0 * dx[2]) + ')')
    ax = axis[1, 1]
    cf = ax.contourf(mf_mean - mf_mean_pos)
    plt.colorbar(cf, ax=ax)
    ax.set_title('diff=A-B')

    ax = axis[0, 2]
    cf = ax.contourf(mf_w)
    plt.colorbar(cf, ax=ax)
    ax.set_title('(C) sum(w[z=' + str(k0 * dx[2]) + '])')
    ax = axis[1, 2]
    cf = ax.contourf(mf_mean - mf_w)
    plt.colorbar(cf, ax=ax)
    ax.set_title('diff=A-C')
    ax = axis[0, 3]
    cf = ax.contourf(mf_w_pos)
    plt.colorbar(cf, ax=ax)
    ax.set_title('(D) sum(w[z]>=0) (z=' + str(k0 * dx[2]) + ')')
    ax = axis[1, 3]
    cf = ax.contourf(mf_mean - mf_w_pos)
    plt.colorbar(cf, ax=ax)
    ax.set_title('diff=B-D')

    # for ax in axis:
    #     ax.set_aspect('equal')
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, hspace=0.2, wspace=0.1)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return




def plot_massflux_locations(mf_mean, mf_mean_pos,
                            mf_mean_single, mf_mean_double, mf_mean_triple,
                            mf_mean_single_pos, mf_mean_double_pos, mf_mean_triple_pos,
                            mf_mean_single_acc, mf_mean_double_acc, mf_mean_triple_acc,
                            mf_mean_single_pos_acc, mf_mean_double_pos_acc, mf_mean_triple_pos_acc,
                            ic_arr, jc_arr, ic, jc,
                            x_single, y_single, x_double, y_double, x_triple, y_triple, dx_box, dy_box,
                            k0, time_range, path_out_figs, fig_name):


    mini = np.maximum(-20, np.int(np.amin(mf_mean)))
    maxi = np.maximum(np.amax(mf_mean), np.amax(mf_mean_pos))
    print('---', mini, maxi)
    lvls = np.arange(mini, maxi+1, dtype=np.int)
    print('---', mini, maxi, lvls[0], lvls[-1])
    del mini, maxi

    ncol = 4
    fig, axis = plt.subplots(2, ncol, figsize=(ncol * 5, 10), sharex='none')
    ax = axis[0, 0]
    cf = ax.contourf(mf_mean.T, levels=lvls, extend='min')
    plt.colorbar(cf, ax=ax)
    ax.set_title('(A) sum(w[z]*rho[z]) (z=' + str(k0 * dx[2]) + ')')
    ax = axis[1, 0]
    cf = ax.contourf(mf_mean_pos.T, levels=lvls, extend='min')
    plt.colorbar(cf, ax=ax)
    ax.set_title('(B) sum(w[z]*rho[z], w>0) (z=' + str(k0 * dx[2]) + ')')

    # rect_single = patch.Rectangle((x_single,y_single), dx_box, dy_box, color='w', fill=False, linewidth=2)
    # rect_double = patch.Rectangle((x_double, y_double), dx_box, dy_box, fill=False, linewidth=2)
    # rect_triple = patch.Rectangle((x_triple, y_triple), dx_box, dy_box, fill=False, linewidth=2)
    # ax = axis[0, 0]
    # print('--- x single, double: ', x_single, x_double)
    # print('--- y single, double: ', y_single, y_double)
    # ax.plot(x_single, y_single, 'o', color='w', markersize=12)
    # ax.plot(x_double, y_double, 'o', color='r', markersize=12)
    # ax.add_artist(rect_single)
    # ax.add_artist(rect_double)
    # ax.add_artist(rect_triple)
    # rect_single = patch.Rectangle((x_single, y_single), dx_box, dy_box, fill=False, linewidth=2)
    # rect_double = patch.Rectangle((x_double, y_double), dx_box, dy_box, fill=False, linewidth=2)
    # rect_triple = patch.Rectangle((x_triple, y_triple), dx_box, dy_box, fill=False, linewidth=2)
    # ax = axis[1, 0]
    # ax.add_artist(rect_single)
    # ax.add_artist(rect_double)
    # ax.add_artist(rect_triple)
    # for i in range(2):
    #     axis[i, 0].plot(ic_arr[:2], jc_arr[:2], 'ko', markersize=7)
    #
    # ax = axis[0, 1]
    # ax.plot(time_range, mf_mean_single, label='flux')
    # ax.plot(time_range, mf_mean_single_pos, '--', label='flux w>0')
    # ax = axis[0, 2]
    # ax.plot(time_range, mf_mean_double, label='flux')
    # ax.plot(time_range, mf_mean_double_pos, '--', label='flux w>0')
    # ax = axis[0, 3]
    # ax.plot(time_range, mf_mean_triple, label='flux')
    # ax.plot(time_range, mf_mean_triple_pos, '--', label='flux w>0')
    # y0 = np.amax(mf_mean_single)
    # ax.plot([tmin,tmax+dt_fields], [y0,y0], 'k', linewidth=0.5)
    # y0 = np.amax(mf_mean_single_pos)
    # ax.plot([tmin,tmax+dt_fields], [y0,y0], 'k--', linewidth=0.5)
    # y0 = np.amax(mf_mean_double)
    # ax.plot([tmin,tmax+dt_fields], [y0,y0], 'k', linewidth=0.5)
    # y0 = np.amax(mf_mean_double_pos)
    # ax.plot([tmin,tmax+dt_fields], [y0,y0], 'k--', linewidth=0.5)
    #
    # ax = axis[1, 1]
    # ax.plot(np.arange(tmin, tmax+dt, dt), mf_mean_single_acc, label='accumulated flux')
    # ax = axis[1, 2]
    # ax.plot(np.arange(tmin, tmax+dt, dt), mf_mean_double_acc, label='accumulated flux')
    # ax.plot(np.arange(tmin, tmax+dt, dt), mf_mean_double_pos_acc, '--', label='accumulated positive flux')
    # ax = axis[1, 3]
    # ax.plot(np.arange(tmin, tmax+dt, dt), mf_mean_triple_acc, label='accumulated flux')
    # ax.plot(np.arange(tmin, tmax+dt, dt), mf_mean_triple_pos_acc, '--', label='accumulated positive flux')
    # y0 = np.amax(mf_mean_single_acc)
    # ax.plot([tmin,tmax+dt], [y0,y0], 'k', linewidth=0.5)
    # y0 = np.amax(mf_mean_single_pos_acc)
    # ax.plot([tmin,tmax+dt], [y0,y0], 'k--', linewidth=0.5)
    # y0 = np.amax(mf_mean_double_acc)
    # ax.plot([tmin,tmax+dt], [y0,y0], 'k', linewidth=0.5)
    # y0 = np.amax(mf_mean_double_pos_acc)
    # ax.plot([tmin,tmax+dt], [y0,y0], 'k--', linewidth=0.5)
    #
    # axis[1,1].legend(loc='best')
    # axis[0,1].set_title('flux single')
    # axis[0,2].set_title('flux double')
    # axis[0,3].set_title('flux triple')
    # axis[1,1].set_title('accumulated flux single')
    # axis[1,2].set_title('accumulated flux double')
    # axis[1,3].set_title('accumulated flux triple')
    # for i in range(1,ncol):
    #     axis[0,i].set_xlim(tmin,tmax)
    #     axis[1,i].set_xlim(tmin,tmax)
    #     axis[1,i].set_xlabel('time [s]')
    # delta = (10.*rstar)/dx[0]
    # for i in range(2):
    #     axis[i,0].set_aspect('equal')
    #     axis[i,0].set_xlim(ic-delta, ic+delta)
    #     axis[i,0].set_ylim(jc-delta, jc+delta)
    #     axis[i,ncol-1].legend(loc='best')

    plt.suptitle('mass flux (t='+str(tmin)+',..'+str(tmax)+', z='+str(k0*dx[2])+'m)', fontsize=21)
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, hspace=0.2, wspace=0.1)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)




    # ncol=3
    # fig, axis = plt.subplots(2, ncol, figsize=(ncol * 5, 10), sharex='none')
    # ax = axis[0, 0]
    # cf = ax.contourf(mf_mean.T, levels=lvls)  # , extend='max')
    # plt.colorbar(cf, ax=ax)
    # ax = axis[0, 1]
    # cf = ax.contourf(mf_mean.T, vmin=mini, vmax=maxi)
    # plt.colorbar(cf, ax=ax)
    # ax = axis[0, 2]
    # cf = ax.contourf(mf_mean.T, levels=lvls, extend='min')
    # plt.colorbar(cf, ax=ax)
    # ax.set_title('(A) sum(w[z]*rho[z]) (z=' + str(k0 * dx[2]) + ')')
    # ax = axis[1, 0]
    # cf = ax.contourf(mf_mean_pos.T, levels=lvls)
    # plt.colorbar(cf, ax=ax)
    # ax = axis[1, 1]
    # cf = ax.contourf(mf_mean_pos.T, vmin=mini, vmax=maxi)
    # plt.colorbar(cf, ax=ax)
    # ax = axis[1, 2]
    # cf = ax.contourf(mf_mean_pos.T, levels=lvls, extend='min')
    # plt.colorbar(cf, ax=ax)
    # ax.set_title('(B) sum(w[z]*rho[z], w>0) (z=' + str(k0 * dx[2]) + ')')
    # plt.suptitle('mass flux (t=' + str(tmin) + ',..' + str(tmax) + ', z=' + str(k0 * dx[2]) + 'm)', fontsize=21)
    # plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, hspace=0.2, wspace=0.1)
    # fig.savefig(os.path.join(path_out_figs, 'test_'+fig_name))
    # plt.close(fig)

    return
# ----------------------------------------------------------------------

def set_input_output_parameters(args):
    print('')
    print('--- set input parameters ---')
    path = args.path

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

    global dt_fields, dt_stats, dt
    dt_fields = np.int(nml['fields_io']['frequency'])
    dt_stats = np.int(nml['stats_io']['frequency'])
    dt = np.int(np.maximum(dt_stats, dt_fields))

    global tmin, tmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100
    print(tmax)

    # times = [np.int(name[:-3]) for name in os.listdir(os.path.join(path, 'fields')) if name[-2:] == 'nc'
    #          and tmin <= np.int(name[:-3]) <= tmax]
    times = [np.int(name[:-3]) for name in os.listdir(os.path.join(path, 'fields')) if name[-2:] == 'nc'
             and tmin <= np.int(name[:-3]) <= tmax and np.mod(np.int(name[:-3]), dt) == 0]
    times.sort()
    nt = len(times)
    print('times: ' + str(times))
    print('tmin, tmax: ' + str(tmin) + ', ' + str(tmax))
    print('dt: ' + str(dt))
    print('nt:', nt)
    files = [str(t) + '.nc' for t in times]
    # print(files)
    print('')

    return path, nml, times, files
# _______________________________________________________

def define_geometry(nml):
    global rstar

    '''--- define geometry ---'''
    if case_name[:21] == 'ColdPoolDry_single_3D':
        rstar = nml['init']['r']
        try:
            print('(ic,jc) from nml')
            ic = np.int(nml['init']['ic'])
            jc = np.int(nml['init']['jc'])
        except:
            print('(ic,jc) NOT from nml')
            ic = np.int(nx / 2)
            jc = np.int(ny / 2)
        ic_arr = [ic]
        jc_arr = [jc]
        ncp = 1
    elif case_name[:21] == 'ColdPoolDry_double_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        sep = nml['init']['sep']
        isep = np.int(np.round(sep / dx[0]))
        jsep = 0
        ic = np.int(np.round(nx / 2))
        jc = np.int(np.round(ny / 2))
        ic1 = ic - np.int(np.round(isep / 2))
        jc1 = jc
        ic2 = ic1 + isep
        jc2 = jc1 + jsep
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
        ncp = 2
    elif case_name[:21] == 'ColdPoolDry_triple_3D':
        print(case_name)
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        try:
            d = nml['init']['d']
        except:
            d = 1e3
        i_d = np.int(np.round(d / dx[0]))
        idhalf = np.int(np.round(i_d / 2))
        a = np.int(np.round(i_d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        r_int = np.int(np.sqrt(3.) / 6 * i_d)  # radius of inscribed circle
        # point of 3-CP collision (ic, jc)
        ic = np.int(np.round(nx / 2))
        jc = np.int(np.round(ny / 2))
        ic1 = ic - r_int
        ic2 = ic1
        ic3 = ic + (a - r_int)
        jc1 = jc - idhalf
        jc2 = jc + idhalf
        jc3 = jc
        ic_arr = np.asarray([ic1, ic2, ic3])
        jc_arr = np.asarray([jc1, jc2, jc3])
        ncp = 3
        print('ic, ic1, ic2, ic3: ', ic, ic1, ic2, ic3)
        print('nx,ny: ', nx, ny)
        print('distance btw. CPs: ', d, i_d, idhalf)
        print('r*: ', rstar, r_int)
        print('')

    return ic_arr, jc_arr, ic, jc, ncp

# _______________________________________________________
def convert_file_for_varlist_horsection(var_list, times, files, path_fields, path_out, level):
    print('................................')
    print('converting: ', times, files)

    # read in test fields file
    fullpath_in = os.path.join(path_fields, files[0])
    rootgrp_in = nc.Dataset(fullpath_in, 'r')
    # field_keys = rootgrp_in.groups['fields'].variables.keys()
    # dims_keys = rootgrp_in.groups['fields'].dimensions.keys()
    dims = rootgrp_in.groups['fields'].dimensions
    nx_ = dims['nx'].size
    ny_ = dims['ny'].size
    rootgrp_in.close()

    k0 = level
    file_name = 'fields_allt_xy_k' + str(np.int(k0)) + '.nc'
    fullpath_out = os.path.join(path_out, file_name)
    print('filename', file_name)

    if os.path.exists(fullpath_out):
        print('')
        print('file ' + fullpath_out + ' already exists! ')
        print('')
    else:
        rootgrp_out = nc.Dataset(fullpath_out, 'w', format='NETCDF4')
        rootgrp_out.createDimension('time', None)
        rootgrp_out.createDimension('nx', nx_)
        rootgrp_out.createDimension('ny', ny_)
        descr_grp = rootgrp_out.createGroup('description')
        var = descr_grp.createVariable('k0', 'f8', )
        var[:] = k0

        time_out = rootgrp_out.createVariable('time', 'f8', ('time',))
        time_out.long_name = 'Time'
        time_out.units = 's'
        time_out[:] = times

        # create variables
        var_list_all = np.append(var_list, 'theta')
        for var in np.append(var_list, 'theta'):
            rootgrp_out.createVariable(var, 'f8', ('time', 'nx', 'ny'))

        # fill variables
        for it, file in enumerate(files):
            print('file: ', file)
            fullpath_in = os.path.join(path_fields, file)
            rootgrp_in = nc.Dataset(fullpath_in, 'r')
            for var in var_list:
                print('var', var)
                var_out = rootgrp_out.variables[var]
                data = rootgrp_in.groups['fields'].variables[var][:, :, k0]
                var_out[it, :, :] = data[:, :]
            var = 'theta'
            data = rootgrp_in.groups['fields'].variables['s'][:, :, k0]
            data_th = theta_s(data)
            del data
            var_out = rootgrp_out.variables[var]
            var_out[it, :,:] = data_th

        rootgrp_out.close()
    return
# _______________________________________________________

def create_output_file(filename, times, path_out):
    nt = len(times)
    print('create output file in: ', path_out)
    print('size: ', nz, nt)

    rootgrp = nc.Dataset(os.path.join(path_out, filename), 'w', format='NETCDF4')

    ts_grp = rootgrp.createGroup('timeseries')
    ts_grp.createDimension('nt', nt)
    var = ts_grp.createVariable('time', 'f8', ('nt'))
    var.units = "s"
    var[:] = times


    field_grp = rootgrp.createGroup('fields_2D')
    field_grp.createDimension('nt', nt)
    field_grp.createDimension('nx', nx)
    field_grp.createDimension('ny', ny)
    var = field_grp.createVariable('mass_flux_2D', 'f8', ('nt', 'nx', 'ny'))
    var.long_name = 'mass flux'
    var.units = "m^(-2)"
    var = field_grp.createVariable('mass_flux_2D_positive', 'f8', ('nt', 'nx', 'ny'))
    var.long_name = 'upward mass flux'
    var.units = "m^(-2)"
    # var = field_grp.createVariable('mass_flux_accumulated', 'f8', ('nx', 'ny'))
    # var.long_name = 'mass flux time accumulated'
    # var = field_grp.createVariable('mass_flux_mean', 'f8', ('nx', 'ny'))
    # var.long_name = 'mass flux time mean'
    # var = field_grp.createVariable('mass_flux_positive_accumulated', 'f8', ('nx', 'ny'))
    # var.long_name = 'upward mass flux time accumulated'
    # var = field_grp.createVariable('mass_flux_positive_mean', 'f8', ('nx', 'ny'))
    # var.long_name = 'upward mass flux time mean'

    rootgrp.close()
    print ''
    return


def dump_output_file(var_name, group_name, var_in, filename, path_out):
    print('dump ' + var_name )
    rootgrp = nc.Dataset(os.path.join(path_out, filename), 'r+', format='NETCDF4')
    grp = rootgrp.groups[group_name]
    var = grp.variables[var_name]
    if group_name == 'timeseries':
        var[:] = var_in
    elif group_name == 'fields_2D':
        var[:,:,:] = var_in
    rootgrp.close()
    return
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