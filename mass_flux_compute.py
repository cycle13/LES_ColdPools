import netCDF4 as nc
import argparse
import json as simplejson
import os
import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
import time



from convert_fields_smaller_k import convert_file_for_varlist_horsection

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


    path, nml = set_input_output_parameters(args)
    ic_arr, jc_arr = define_geometry(nml)
    path_out_figs = os.path.join(path, 'figs_massflux')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)

    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100
    itmin = np.int(tmin / dt_fields)
    itmax = np.int(tmax / dt_fields)
    times = [np.int(name[:-3]) for name in os.listdir(os.path.join(path, 'fields')) if name[-2:] == 'nc'
             and tmin <= np.int(name[:-3]) <= tmax]
    times.sort()
    nt = len(times)
    print('times: ' + str(times))
    print('nt:', nt)
    files = [str(t) + '.nc' for t in times]
    print(files)
    print('')

    if args.level:
        zlevel = np.double(args.level)
    else:
        zlevel = 1000.
    k0 = np.int(zlevel / dx[2])
    print('-----------------------------------')
    print('considering mass flux at level: z=' + str(zlevel) + ', k=' + str(k0))
    print('')


    path_field_k0 = os.path.join(path, 'fields_merged', 'fields_allt_xy_k'+str(k0)+'.nc')
    if not os.path.exists(path_field_k0):
        print('merging file for k='+str(k0))
        convert_file_for_varlist_horsection(['w'], times, files, os.path.join(path, 'fields'), os.path.join(path, 'fields_merged'), k0)

    field_k0 = nc.Dataset(path_field_k0, 'r')
    t_merged = field_k0.variables['time'][:]
    print('time merged field: ', t_merged, t_merged[0], t_merged[-1])
    if t_merged[0] > tmin or t_merged[-1] < tmax:
        print('removing file: ', os.path.join(path, 'fields_merged', 'fields_allt_xy_k'+str(k0)+'.nc'))
        os.remove(os.path.join(path, 'fields_merged', 'fields_allt_xy_k'+str(k0)+'.nc'))
        convert_file_for_varlist_horsection(['w'], times, files, os.path.join(path, 'fields'),
                                            os.path.join(path, 'fields_merged'), k0)
    print('')



    ''' computing mass flux '''
    rho = np.zeros((nt), dtype=np.double)
    # read in density
    rootgrp = nc.Dataset(os.path.join(path, 'stats', 'Stats.' + case_name + '.nc'))
    time_stats = rootgrp.groups['profiles'].variables['t'][:]
    # rho0_stat = rootgrp.groups['reference'].variables['rho0_full'][:]
    rho_mean = 1./rootgrp.groups['profiles'].variables['alpha_mean'][:,k0]
    rootgrp.close()
    if time_stats[-1] < tmax:
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
        print(rho.shape, rho_mean.shape, nt)
        rho[:] = rho_mean[:nt]
    print('rho: ', rho)
    # compute mass flux
    w = field_k0.variables['w'][:,:,:]
    mf_w = np.sum(w[:itmax+1,:,:], axis=0)
    # mf_ref = rho0_stats[k0]*np.sum(field_k0.variables['w'], axis=0)
    aux = np.zeros((w.shape), dtype=np.double)
    for it in range(itmin, itmax+1):
        aux[it,:,:] = rho[k0]*w[it,:,:]
    mf_mean = np.sum(aux, axis=0)
    # only account for positive (i.e., upward) mass flux
    mf_w_pos = np.sum(np.where(w[:itmax + 1, :, :] > 0, w[:itmax + 1, :, :], 0), axis=0)
    mf_mean_pos = np.sum(np.where(w>0, aux, 0), axis=0)
    del aux
    field_k0.close()





    # ''' plotting '''
    # fig_name = 'massflux_k' + str(k0) + '_tmax' + str(tmax) + '.png'
    # plot_massflux_allversions(mf_mean, mf_mean_pos, mf_w, mf_w_pos, k0, path_out_figs, fig_name)
    fig_name = 'massflux_locations_k' + str(k0) + '_tmax' + str(tmax) + '.png'
    plot_massflux_locations(mf_mean, mf_mean_pos, mf_w, mf_w_pos,
                            ic_arr, jc_arr, k0, path_out_figs, fig_name)


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




def plot_massflux_locations(mf_mean, mf_mean_pos, mf_w, mf_w_pos,
                            ic_arr, jc_arr, k0, path_out_figs, fig_name):
    dx_box = 10
    dy_box = 10
    x_single = ic_arr[0]

    ncol = 2
    fig, axis = plt.subplots(2, ncol, figsize=(ncol * 5, 10), sharey='all')
    ax = axis[0, 0]
    cf = ax.contourf(mf_mean.T)
    for i in range(3):
        ax.plot(ic_arr[i], jc_arr[i], 'ko', markersize=10)
    plt.colorbar(cf, ax=ax)
    ax.set_title('(A) sum(w[z]*rho[z]) (z=' + str(k0 * dx[2]) + ')')
    ax = axis[0, 1]
    cf = ax.contourf(mf_mean.T)
    plt.colorbar(cf, ax=ax)
    ax.set_title('(A) sum(w[z]*rho[z]) (z=' + str(k0 * dx[2]) + ')')


    ax = axis[1, 0]
    cf = ax.contourf(mf_mean_pos.T)
    plt.colorbar(cf, ax=ax)
    ax.set_title('(B) sum(w[z]*rho[z], w>0) (z=' + str(k0 * dx[2]) + ')')

    # for ax in axis:
    #     ax.set_aspect('equal')
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.07, top=0.9, hspace=0.2, wspace=0.1)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

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

    global dt_fields
    dt_fields = nml['fields_io']['frequency']

    return path, nml
# _______________________________________________________

def define_geometry(nml):

    '''--- define geometry ---'''
    global rstar
    if case_name == 'ColdPoolDry_double_2D':
        rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        isep = 4 * irstar
        ic1 = np.int(nx / 3)
        ic2 = ic1 + isep
        jc1 = np.int(ny / 2)
        jc2 = jc1
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name[:21] == 'ColdPoolDry_single_3D':
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
    elif case_name == 'ColdPoolDry_double_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        isep = 4 * irstar
        jsep = 0
        ic1 = np.int(np.round((nx + 2 * gw) / 3)) - gw
        jc1 = np.int(np.round((ny + 2 * gw) / 2)) - gw
        ic2 = ic1 + isep
        jc2 = jc1 + jsep
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_triple_3D':
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
        print(ic1, ic2, ic3)
        print(nx, ny, id, idhalf)
        print(rstar, r_int, ic)
        print('')



    return ic_arr, jc_arr

# _______________________________________________________

if __name__ == '__main__':
    main()