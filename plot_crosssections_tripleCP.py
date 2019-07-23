import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35

def main():

    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    case, times, files = set_input_parameters(args)
    define_geometry(case_name, files)
    kmax = np.int(2000. / dx[2])
    krange = np.arange(kmax)
    print('krange: ', krange)
    print('kmax: ', kmax)
    nt = len(times)
    icoll, jcoll, ic_arr, jc_arr, x_half, y_half, z_half = define_geometry(case_name, files)

    filename_out = os.path.join('data_analysis', 'stats_crosssections.nc')

    ''' --- plot configuration ---'''
    rootgrp = nc.Dataset(os.path.join(path_fields, '0.nc'), 'r')
    var = rootgrp.groups['fields'].variables['s'][:, :, 0]
    rootgrp.close()
    plot_configuration(var, icoll, jcoll, ic_arr, jc_arr)

    #
    # ''' (a) compute horizontal convergence''''
    #


    ''' (a) compute and plot maxima in single, double and triple time window '''
    # time_single: time, when w(icoll,jcoll,k=0)<1.m/s for the last time (single for t<=time_single)
    # time_double: timestep when w(icol,jcoll,k=0) decays again for the first time
    var_list = ['w', 's']
    for var_name in var_list:
        print('')
        print('------ var: '+var_name)
        min_single, min_double, min_triple, max_single, max_double, max_triple, \
            time_single, time_double, it_single, it_double = compute_minima_maxima(var_name, icoll, jcoll, ic_arr, jc_arr, times, nt, kmax)
        output_minima_maxima(var_name, min_single, min_double, min_triple, max_single, max_double, max_triple,
                         time_single, time_double, times, kmax, filename_out)
        plot_xz_crosssection_single_double_triple(var_name, min_single, min_double, min_triple,
                                              max_single, max_double, max_triple,
                                              icoll, jcoll, ic_arr, jc_arr,
                                              time_single, time_double, times, x_half)
    print('')


    # ''' compute mass-flux '''
    # # mass_flux[i,j,k]= dm*w[i,j,k]*dt_fields
    # # mass_flux[i,j,k]= rho*dA*w[i,j,k]*dt_fields = rho[i,j,k]*dx*dy*w[i,j,k]*dt_fields
    # # define area:
    # # - band along xz-crosssection of width 300m
    # # - single: between x=r*+delta and x=icoll+10; time btw. t<=time_single
    # # - double: between x=ic_arr[0]-r_int=2*ic_arr[0]-icoll and x=ic_arr[0]+r_int=icoll; time time_single<t<=time_double
    # # - triple: between x=icoll-10 and x=icoll+10; time_double<t
    # width=400.
    # delta_y = np.int(width/(2*dx[1]))
    # rstar = 2000.
    # ir = rstar / dx[0]
    # # path = '/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2/dTh3_z1000_r1000/stats/'
    # # case_name = 'ColdPoolDry_single_3D'
    # try:
    #     rootgrp = nc.Dataset(os.path.join(path_in, 'stats', 'Stats.' + case_name + '.nc'))
    # except:
    #     rootgrp = nc.Dataset(os.path.join(path_in, 'Stats.' + case_name + '.nc'))
    # rho0_stats = rootgrp.groups['reference'].variables['rho0'][:]
    # rootgrp.close()
    # mass_flux_ = np.ndarray((3, nt, kmax), dtype=np.double)
    # mass_flux = np.ndarray((3, kmax), dtype=np.double)
    # for it, t0 in enumerate(times):
    #     rootgrp = nc.Dataset(os.path.join(path_fields, str(t0) + '.nc'), 'r')
    #     var = rootgrp.groups['fields'].variables['w'][:, :, :kmax]
    #     rootgrp.close()
    #     mass_flux_[0,it,:] = np.sum(np.sum(var[icoll+10:ic_arr[2]-ir,jcoll-delta_y:jcoll+delta_y,:], axis=0), axis=0)
    #     mass_flux_[1,it,:] = np.sum(np.sum(var[2*ic_arr[0]-icoll:icoll-3,jcoll-delta_y:jcoll+delta_y,:], axis=0), axis=0)
    #     mass_flux_[2,it,:] = np.sum(np.sum(var[icoll-3:icoll+10,jcoll-delta_y:jcoll+delta_y,:], axis=0), axis=0)
    # mass_flux_ = dx[0]*dx[1]*mass_flux_
    # for k in range(kmax):
    #     mass_flux_[:,:,k] = rho0_stats[k]*mass_flux_[:,:,k]
    # mass_flux[0,:] = dt_fields*np.sum(mass_flux_[0,:time_single,:], axis=0)
    # mass_flux[1,:] = dt_fields*np.sum(mass_flux_[1,time_single:time_double+1,:], axis=0)
    # mass_flux[2,:] = dt_fields*np.sum(mass_flux_[2,:time_single+1,:], axis=0)




    ''' (b) plot crosssections at level k0 for w, s for all time-steps '''
    var_list = ['w', 's']
    # var_list = ['s']
    for k0 in krange:
        print('---- k0: '+str(k0) + '----')
        if case == 'triple':
            j0 = jcoll
            imin = (icoll - 200)
            imax = (icoll + 200)
            delta = r_int
            plot_xz_crosssections(var_list,
                                  min_single, min_double, min_triple,
                                  max_single, max_double, max_triple,
                                  j0, k0, ic_arr, jc_arr, icoll, jcoll,
                                  x_half, y_half, imin, imax, delta,
                                  files, path_out, path_fields)

    ''' (c) plot crosssections at different levels for collision time '''
    krange_ = [1, 4, 10, 15]
    if case == 'triple':
        j0 = jcoll
        imin = (icoll - 200)
        imax = (icoll + 200)
        plot_xz_crosssections_multilevel('w', j0, krange_, imin, imax,
                                         icoll, jcoll, ic_arr, jc_arr,
                                         files, path_out, path_fields)
    # print ''

    return


# --------------------------------------------------------------------

def compute_minima_maxima(var_name, icoll, jcoll, ic_arr, jc_arr, times, nt, kmax):
    print('compute and plot maxima / minima (kmax='+str(kmax)+')')

    max_single_ = np.zeros((nt, kmax), dtype=np.double)
    max_double_ = np.zeros((nt, kmax), dtype=np.double)
    max_triple_ = np.zeros((nt, kmax), dtype=np.double)
    min_single_ = np.zeros((nt, kmax), dtype=np.double)
    min_double_ = np.zeros((nt, kmax), dtype=np.double)
    min_triple_ = np.zeros((nt, kmax), dtype=np.double)
    time_single = 0

    for it, t0 in enumerate(times):
        rootgrp = nc.Dataset(os.path.join(path_fields, str(t0) + '.nc'), 'r')
        var = rootgrp.groups['fields'].variables[var_name][:, :, :kmax]
        rootgrp.close()
        if var_name == 's':
            for k in range(kmax):
                if var_name == 's':
                    min_single_[it, k] = np.maximum(6000, np.amin(var[icoll+10:, jcoll, k]))
                    min_double_[it, k] = np.maximum(6000, np.amin(var[:icoll-3, jcoll, k]))
                    min_triple_[it, k] = np.maximum(6000, np.amin(var[icoll-3:icoll+3, jcoll, k]))
        else:
            min_single_[it, :] = np.amin(var[icoll + 10:, jcoll, :], axis=0)
            min_double_[it, :] = np.amin(var[:icoll - 3, jcoll, :], axis=0)
            min_triple_[it, :] = np.amin(var[icoll - 3:, jcoll, :], axis=0)
        max_single_[it, :] = np.amax(var[icoll + 10:, jcoll, :], axis=0)
        max_double_[it, :] = np.amax(var[:icoll - 3, jcoll, :], axis=0)
        max_triple_[it, :] = np.amax(var[icoll - 3:, jcoll, :], axis=0)
        if var_name == 'w':
            if time_single == 0 and var[ic_arr[0], jcoll, 0] > 1.:
                time_single = t0 - dt_fields
            time_double = (np.argmax(np.amax(max_double_[:, :], axis=1))) * dt_fields + tmin
        elif var_name == 's':
            if time_single == 0 and var[ic_arr[0], jcoll, 0] < (np.mean(var[:10, jcoll, 0])-5):
                time_single = t0 - dt_fields
            time_double = (np.argmin(np.amin(min_triple_[:, :], axis=1))) * dt_fields + tmin
    it_single = (time_single-tmin) / dt_fields
    it_double = (time_double-tmin) / dt_fields
    print('it single, double: ', it_single, it_double)
    print('time single, double: ', time_single, time_double, it_single, it_double)

    max_single = np.amax(max_single_[:it_single + 1, :], axis=0)
    max_double = np.amax(max_double_[it_single:it_double+1,:], axis=0)
    max_triple = np.amax(max_triple_[it_double:,:], axis=0)
    min_single = np.amin(min_single_[:it_single+1,:], axis=0)
    min_double = np.amin(min_double_[it_single:it_double+1,:], axis=0)
    min_triple = np.amin(min_triple_[it_double:,:], axis=0)
    # if var_name == 'w':
    #     max_single = np.amax(max_single_[:it_single+1,:], axis=0)
    #     max_double = np.amax(max_double_[it_single:it_double+1,:], axis=0)
    #     max_triple = np.amax(max_triple_[it_double:,:], axis=0)
    #     min_single = np.amin(min_single_[:it_single+1,:], axis=0)
    #     min_double = np.amin(min_double_[it_single:it_double+1,:], axis=0)
    #     min_triple = np.amin(min_triple_[it_double:,:], axis=0)
    # elif var_name == 's':
    #     max_single = np.amax(max_single_[:it_single+1, :], axis=0)
    #     max_double = np.amax(max_double_[it_single:it_double+1, :], axis=0)
    #     max_triple = np.amax(max_triple_[it_double:, :], axis=0)
    #     min_single = np.amin(min_single_[:it_single+1, :], axis=0)
    #     min_double = np.amin(min_double_[it_single:it_double+1, :], axis=0)
    #     min_triple = np.amin(min_triple_[it_double:, :], axis=0)



    # ---- plotting -----
    cm = plt.cm.get_cmap('coolwarm')

    fig_name = 'xz_plane_' + var_name + '_single_double_triple.png'
    fig, axis = plt.subplots(3, 1, figsize=(12, 12))
    max = np.maximum(np.amax(max_double), np.amax(max_triple)) + 0.5
    min = np.amin(min_double) - 0.2
    for i, ax in enumerate(axis):
        ax.set_xlabel('x')
        ax.set_ylabel(var_name)
        # ax.set_xlim(imin * dx[0], imax * dx[0])
        ax.set_ylim(min, max)
    #     ax.grid()
    # plot min/max
    for it,t0 in enumerate([time_double-dt_fields, time_double, time_double+dt_fields]):
        rootgrp = nc.Dataset(os.path.join(path_fields, str(t0) + '.nc'), 'r')
        var = rootgrp.groups['fields'].variables[var_name][:, jcoll, :]
        rootgrp.close()
        for k0 in range(10):
            count_color = np.double(k0)/10
            axis[it].plot(var[:, k0], color=cm(count_color), label='k=' + str(k0))
    axis[0].set_title('t=time_double-dt='+str(time_double-dt_fields))
    axis[1].set_title('t=time_double='+str(time_double))
    axis[2].set_title('t=time_double+dt='+str(time_double+dt_fields))
    del var
    axis[-1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
               fancybox=True, shadow=False, ncol=6, fontsize=10)
    # plt.suptitle(var_name + ' (z=' + str(k0 * dx[2]) + 'm, k=' + str(k0) + ')', fontsize=21)
    plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.9, hspace=0.26)
    plt.savefig(os.path.join(path_out, fig_name))
    plt.close()




    print('plot minmax levels')
    fig_name = 'minmax_levels_' + var_name + '_' + case_name + '.png'
    zrange = np.arange(kmax) * dx[2]
    fig, axis = plt.subplots(1, 2, figsize=(6, 8))
    ax0 = axis[0]
    ax1 = axis[1]
    ax0.plot(max_single, zrange, 'o-', label='single')
    ax0.plot(max_double, zrange, 'o--', label='double')
    ax0.plot(max_triple, zrange, 'o:', linewidth=2, label='triple')
    ax1.plot(min_single, zrange, 'o-')
    ax1.plot(min_double, zrange, 'o--')
    ax1.plot(min_triple, zrange, 'o:', linewidth=2)
    for ax in axis:
        ax.set_xlabel(var_name)
        ax.set_ylabel('height  [m]')
        ax.grid()
    ax0.set_title('maxima')
    ax1.set_title('minima')
    ax0.legend(loc='upper center', bbox_to_anchor=(0.75, -0.05),
               fancybox=True, shadow=False, ncol=6, fontsize=10)
    plt.suptitle('min/max for ' + var_name, fontsize=21)
    plt.subplots_adjust(bottom=0.1, right=.95, left=0.1, top=0.9, hspace=0.4)
    plt.savefig(os.path.join(path_out, fig_name))
    plt.close()

    return min_single, min_double, min_triple, max_single, max_double, max_triple, \
           time_single, time_double, it_single, it_double

# --------------------------------------------------------------------

def plot_configuration(var, icoll, jcoll, ic_arr, jc_arr):
    fig, (ax1, ax2)= plt.subplots(1,2)
    ax1.imshow(var[:, :].T, origin='lower')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    plt.tight_layout()
    ax1.set_xlim(0,nx)
    ax1.set_ylim(0,ny)
    ax2.contourf(var[:, :].T, origin='lower')
    ax2.set_xlabel('x')

    for i,ic in enumerate(ic_arr):
        col = str(np.double(i)/len(ic_arr))
        ax1.plot([ic, ic], [0, ny-1], color=col, linewidth=3, label='ic'+str(i))
        ax2.plot([ic, ic], [0, ny-1], color=col, linewidth=3)
    for j,jc in enumerate(jc_arr):
        col = str(np.double(j) / len(jc_arr))
        ax1.plot([0, nx], [jc,jc], color=col, linewidth=3, label='jc'+str(j))
        ax2.plot([0, nx], [jc,jc], color=col, linewidth=3)
    ax1.plot([icoll, icoll], [0, ny], 'g--', linewidth=2, label='icoll')
    ax2.plot([icoll, icoll], [0, ny], 'g--', linewidth=2)
    ax1.plot([0, nx], [jcoll, jcoll], 'g--', linewidth=2, label='jcoll')
    ax2.plot([0, nx], [jcoll, jcoll], 'g--', linewidth=2)

    ax1.plot([ic_arr[0]+r_int, ic_arr[0]+r_int], [0, ny], 'b:', linewidth=2, label='ic0 + r_int')
    # ax2.plot([icoll, icoll], [0, ny], 'g--', linewidth=2)
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
               fancybox=True, shadow=False, ncol=3, fontsize=10)

    plt.tight_layout()
    ax2.set_xlim(icoll-200,icoll+200)
    ax2.set_ylim(jcoll-200,jcoll+200)
    ax2.set_aspect('equal')
    plt.savefig(os.path.join(path_out, 'initial.png'))
    plt.close()

    return

# ----------------------------------

def plot_xz_crosssection_single_double_triple(var_name,
                                              min_single, min_double, min_triple,
                                              max_single, max_double, max_triple,
                                              icoll, jcoll, ic_arr, jc_arr,
                                              time_single, time_double, times, x_half):
    print('plot xz-crosssection single, double, triple')

    cm = plt.cm.get_cmap('coolwarm')
    imin = (icoll - 200)
    imax = (icoll + 200)
    for k0 in [0, 2, 4, 6, 8, 10, 12, 14, 16]:
        fig_name = 'xz_plane_' + var_name + '_single_double_triple_k' + str(k0) + '.png'
        fig, axis = plt.subplots(3, 1, figsize=(12, 12))
        ax1 = axis[0]
        ax2 = axis[1]
        ax3 = axis[2]
        max = np.maximum(max_double[k0], max_triple[k0]) + 0.5
        min = np.minimum(min_single[k0], min_double[k0]) - 0.2
        for i, ax in enumerate(axis):
            # plot boundary of area considered for double, triple collision
            ax.plot([dx[0] * (icoll - 3), dx[0] * (icoll - 3)], [min, max], 'b--', linewidth=0.5)
            ax.plot([dx[0] * icoll, dx[0] * icoll], [min, max], 'b--', linewidth=1)
            ax.plot([dx[0] * (icoll + 10), dx[0] * (icoll + 10)], [min, max], 'b--', linewidth=0.5)
            ax.plot([dx[0] * ic_arr[0], dx[0] * ic_arr[0]], [min, max], 'k--', linewidth=1, label='ic1')
            ax.set_xlabel('y')
            ax.set_ylabel(var_name)
            ax.set_xlim(imin * dx[0], imax * dx[0])
            ax.set_ylim(min, max)
            ax.grid()
        # plot min/max
        ax1.plot([0,nx*dx[0]],[max_single[k0], max_single[k0]], 'k', linewidth=0.5)
        ax2.plot([0,nx*dx[0]],[max_double[k0], max_double[k0]], 'k', linewidth=0.5)
        ax3.plot([0,nx*dx[0]],[max_double[k0], max_double[k0]], 'k', linewidth=0.5)
        ax3.plot([0,nx*dx[0]],[max_triple[k0], max_triple[k0]], 'k--', linewidth=0.5)
        for it, t0 in enumerate(times):
            rootgrp = nc.Dataset(os.path.join(path_fields, str(t0) + '.nc'), 'r')
            var = rootgrp.groups['fields'].variables[var_name][:, :, k0]
            rootgrp.close()
            count_color = np.double(it) / len(times)
            if t0 <= time_single:
                ax1.plot(x_half, var[:, jcoll], color=cm(count_color), label='t=' + str(t0))
            elif t0 <= time_double:
                ax2.plot(x_half, var[:, jcoll], color=cm(count_color), label='t=' + str(t0))
            else:
                ax3.plot(x_half, var[:, jcoll], color=cm(count_color), label='t=' + str(t0))
        ax1.set_title('single: t<=' + str(time_single))
        ax2.set_title('double: t<=' + str(time_double))
        ax3.set_title('triple: t>' + str(time_double))
        # ax3.set_title('through collision point (max(2CP)=' + str(np.round(var_max2, 1))
        #               + ', max(3CP)=' + str(np.round(var_max3, 1)) + ')')
        ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
                   fancybox=True, shadow=False, ncol=6, fontsize=10)
        plt.suptitle(var_name + ' (z=' + str(k0 * dx[2]) + 'm, k=' + str(k0) + ')', fontsize=21)
        plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.9, hspace=0.26)
        plt.savefig(os.path.join(path_out, fig_name))
        plt.close()
    return

# ----------------------------------

def plot_xz_crosssections(var_list,
                          min_single, min_double, min_triple,
                          max_single, max_double, max_triple,
                          j0, k0, ic_arr, jc_arr,
                          icoll, jcoll, x_arr, y_arr, imin, imax, delta, files,
                          path_out, path_fields):
    # cm = plt.cm.get_cmap('bwr')
    cm = plt.cm.get_cmap('coolwarm')
    # cm = plt.cm.get_cmap('winter')
    # cm = plt.cm.get_cmap('viridis')
    # cm = plt.cm.get_cmap('viridis_r')

    print('xz-crosssection at j0='+str(j0))
    ic1 = ic_arr[0]
    ic2 = ic_arr[2]
    jc1 = jc_arr[0]
    jc2 = jc_arr[2]

    for var_name in var_list:
        fig, axis = plt.subplots(3, 1, figsize=(12, 12))
        ax1 = axis[0]
        ax2 = axis[1]
        ax3 = axis[2]
        s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
        ax1.imshow(s0[:, :, k0].T, origin="lower")
        ax1.plot([0, nx-1], [jc1, jc1], 'k', linewidth=2)
        ax1.plot([0, nx-1], [j0, j0], 'k', linewidth=2)
        ax1.plot([ic1, ic1], [0, ny-1], 'k--', linewidth=1)
        ax1.plot([ic2, ic2], [0, ny-1], 'k--', linewidth=1)
        ax1.plot([icoll, icoll], [0, ny-1], 'k:', linewidth=2)
        ax1.set_xlim([imin,imax])
        ax1.set_ylim([imin,imax])
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        count = 0
        for file in files:
            rootgrp = nc.Dataset(os.path.join(path_fields, file), 'r')
            var = rootgrp.groups['fields'].variables[var_name][:, :, k0]
            rootgrp.close()
            count_color = np.double(count) / len(files)
            ax2.plot(x_arr,var[:,jc1], color=cm(count_color), label='t='+str(file[:-3]))
            ax3.plot(x_arr,var[:,j0], color=cm(count_color), label='t='+str(file[:-3]))
            ax2.set_xlim(imin*dx[0], imax*dx[0])
            ax3.set_xlim(imin*dx[0], imax*dx[0])
            count += 1

        max = np.maximum(max_double[k0], max_triple[k0]) + 0.1
        min = np.minimum(min_single[k0], min_double[k0]) - 0.1
        ax3.plot([dx[0]*icoll, dx[0]*icoll], [min, max], ':', linewidth=2)
        ax3.fill_betweenx(np.linspace(min, max, 2), dx[0]*(ic1-delta)*np.ones(2), dx[0]*(ic1+delta-2)*np.ones(2),
                          facecolor='green', alpha=0.5, linewidth=0)
        ax3.fill_betweenx(np.linspace(min, max, 2), dx[0]*(ic1+delta-2)*np.ones(2), dx[0]*nx*np.ones(2),
                          facecolor='gray', linewidth=0)
        for i,ax in enumerate(axis[1:]):
            ax.plot([dx[0]*ic1, dx[0]*ic1], [min, max], 'k--', linewidth=1.5)
            ax.plot([dx[0]*ic2, dx[0]*ic2], [min, max], 'k--', linewidth=1.5)
            ax.set_xlabel('y')
            ax.set_ylabel(var_name)
            ax.grid()
        ax2.set_title('through center of coldpool #1 (max='+str(np.round(max_single[k0],1))+')')
        ax3.set_title('through collision point (max(2CP)='+str(np.round(max_double[k0],1))
                      + ', max(3CP)='+str(np.round(max_triple[k0],1))+')')
        ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
                   fancybox=True, shadow=False, ncol=6, fontsize=10)
        plt.suptitle(var_name + ' (z='+str(k0*dx[2])+'m, k='+str(k0)+')', fontsize=21  )
        plt.subplots_adjust(bottom=0.18, right=.95, left=0.1, top=0.95, hspace=0.26)
        plt.savefig(os.path.join(path_out, 'xz_plane_' + var_name + '_z' + str(np.int(k0*dx[2])) + 'm.png'))
        plt.close()
    return



def plot_y_crosssections(var_list, k0, i0, ic1, jc1, ic2, jc2, files, path_out, path_fields):
    cm = plt.cm.get_cmap('coolwarm')
    # cm = plt.cm.get_cmap('viridis')

    var_min1 = 9999.9
    var_max1 = -9999.9
    for var_name in var_list:
        f1 = plt.figure(figsize=(6, 12))
        ax1 = plt.subplot(3, 1, 1)
        ax1.set_title('entropy')
        ax2 = plt.subplot(3, 1, 2)
        ax2.set_title('y-crossection through center of coldpool #1')
        ax3 = plt.subplot(3, 1, 3)
        ax3.set_title('y-crossection through collision point')
        s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
        ax1.imshow(s0[:,:,k0].T)
        ax1.plot([ic1,ic1],[0,2*jc1], 'k', linewidth=2)
        ax1.plot([i0,i0],[0,2*jc1], 'k', linewidth=2)

        for file in files:
            var = read_in_netcdf_fields(var_name, os.path.join(path_fields, file))
            var_min1 = np.minimum(var_min1, np.amin(var[:, jc1, k0]))
            var_max1 = np.maximum(var_max1, np.amax(var[:, jc1, k0]))
            count_color = np.double(file[:-3]) / np.double(files[-1][:-3])
            ax2.plot(var[ic1, :, k0], color=cm(count_color), label='t=' + str(file[:-3]))
            ax3.plot(var[i0, :, k0], color=cm(count_color), label='t=' + str(file[:-3]))

        plt.subplot(3, 1, 1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.subplot(3, 1, 2)
        plt.xlabel('y')
        plt.ylabel(var_name)
        plt.subplot(3, 1, 3)
        plt.xlabel('y')
        plt.ylabel(var_name)
        #
        # ax3.legend( )
        # Put a legend below current axis
        ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                   fancybox=True, shadow=True, ncol=5)
        plt.suptitle(var_name + '  (k=' + str(k0) + ')', fontsize=21)
        plt.savefig(os.path.join(path_out, 'yz_plane_' + var_name + '_k'+str(k0)+'.png'))
        plt.close()
    return



def plot_xz_crosssections_multilevel(var_name, j0, krange, imin, imax,
                                     icoll, jcoll, ic_arr, jc_arr,
                                     files, path_out, path_fields):
    cm = plt.cm.get_cmap('coolwarm')
    cm_cw = plt.cm.get_cmap('coolwarm')

    nline = 5

    fig, axis = plt.subplots(nline, 1, figsize=(8, 18))
    ax1 = axis[0]
    s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
    ax1.imshow(s0[:, :, 1].T, origin='lower', cmap=cm_cw)
    ax1.plot([0, nx], [j0, j0], 'k', linewidth=2)
    ax1.plot([icoll, icoll], [0, ny], 'k:', linewidth=2)
    ax1.set_xlim(imin, imax)
    ax1.set_ylim(imin, imax)

    max = np.zeros(len(krange))
    var_min = 9999.9
    var_max = -9999.9
    for file in files:
        var = read_in_netcdf_fields(var_name, os.path.join(path_fields, file))
        var_min = np.minimum(var_min, np.amin(var))
        var_max = np.maximum(var_max, np.amax(var))
        count_color = np.double(file[:-3]) / np.double(files[-1][:-3])
        for i,ax in enumerate(axis[1:]):
            ax.plot(var[:, j0, krange[i]], color=cm(count_color), label='t=' + str(file[:-3]))
            ax.plot([icoll, icoll], [var_min, var_max], ':', linewidth=2)
            k = krange[i]
            max[i] = np.maximum(max[i], np.amax(var[:,j0,k]))
    for i,ax in enumerate(axis[1:]):
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_xlim(imin, imax)
        k = krange[i]
        ax.set_title('k=' + str(k) + ', z=' + str(dx[2] * k) + 'm (max=' + str(np.round(max[i], 1)) + ')')
        ax.set_xlabel('y')
        ax.set_ylabel(var_name)
    axis[-1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
                   fancybox=True, shadow=True, ncol=5)
    plt.suptitle(var_name + ' multilevel')
    plt.subplots_adjust(bottom=0.11, right=.95, left=0.1, top=0.95, hspace=0.3)
    plt.savefig(os.path.join(path_out, 'xz_multilevel_'+var_name+'.png'))
    plt.close(fig)
    return



def plot_yz_crosssections_multilevel(var_name, i0, krange, jmin, jmax,
                                     jc_arr, files, path_out, path_fields):
    cm = plt.cm.get_cmap('coolwarm')
    jc1 = jc_arr[0]
    nline = 5

    fig, axis = plt.subplots(nline, 1, figsize=(6, 18))
    ax1 = axis[0]
    s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
    ax1.imshow(s0[:, :, 1].T, origin='lower')
    ax1.plot([i0, i0], [0, 4 * jc1], 'k', linewidth=2)
    ax1.set_xlim(jmin,jmax)
    ax1.set_ylim(jmin,jmax)

    max = np.zeros(len(krange))
    for file in files:
        var = read_in_netcdf_fields(var_name, os.path.join(path_fields, file))
        count_color = np.double(file[:-3]) / np.double(files[-1][:-3])
        for i,ax in enumerate(axis[1:]):
            ax.plot(var[i0, :, krange[i]], color=cm(count_color), label='t=' + str(file[:-3]))
            k = krange[i]
            max[i] = np.maximum(max[i], np.amax(var[i0,:,k]))

    for i,ax in enumerate(axis[1:]):
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_ylim(jmin, jmax)
        k = krange[i]
        ax.set_title('k=' + str(k) + ', z=' + str(dx[2] * k) + 'm (max=' + str(np.round(max[i], 1)) + ')')
        ax.set_xlabel('y')
        ax.set_ylabel(var_name)
    axis[-1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
                    fancybox=True, shadow=True, ncol=5)
    plt.suptitle(var_name + ' multilevel')
    plt.subplots_adjust(bottom=0.12, right=.95, left=0.1, top=0.95, hspace=0.3)
    plt.savefig(os.path.join(path_out, 'yz_multilevel_'+var_name+'.png'))
    plt.close()
    return




# ----------------------------------

def output_minima_maxima(var_name, min_single, min_double, min_triple, max_single, max_double, max_triple,
                         time_single, time_double, times, kmax, filename_out):
    print('')
    print(os.path.join(path_in, filename_out))
    print(path_in)
    if not os.path.exists(os.path.join(path_in, 'data_analysis')):
        os.mkdir(os.path.join(path_in, 'data_analysis'))
    root_out = nc.Dataset(os.path.join(path_in, filename_out), 'w', format='NETCDF4')

    rootgrp = root_out.createGroup(var_name, )
    rootgrp.createDimension('time', None)
    rootgrp.createDimension('nz', kmax)

    time_out = rootgrp.createVariable('time', 'f8', ('time',))
    time_out.long_name = 'Time'
    time_out.units = 's'
    time_out[:] = times
    var = rootgrp.createVariable('time_single', 'f8', )
    var.long_name = 'single for t<=time_single'
    var[:] = time_single
    var = rootgrp.createVariable('time_double', 'f8', )
    var.long_name = 'double for time_single<t<=time_double'
    var[:] = time_double

    var = rootgrp.createVariable('min_single', 'f8', ('nz'))
    var[:] = min_single[:]
    var = rootgrp.createVariable('min_double', 'f8', ('nz'))
    var[:] = min_double[:]
    var = rootgrp.createVariable('min_triple', 'f8', ('nz'))
    var[:] = min_triple[:]
    var = rootgrp.createVariable('max_single', 'f8', ('nz'))
    var[:] = max_single[:]
    var = rootgrp.createVariable('max_double', 'f8', ('nz'))
    var[:] = max_double[:]
    var = rootgrp.createVariable('max_triple', 'f8', ('nz'))
    var[:] = max_triple[:]

    root_out.close()

    return

# ----------------------------------


def set_input_parameters(args):
    print ''' setting parameters '''
    global path_in, path_fields, path_out
    path_in = args.path
    if os.path.exists(os.path.join(path_in, 'fields')):
        path_fields = os.path.join(path_in, 'fields')
    path_out = os.path.join(path_in, 'figs_crosssections')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    print('path: ', path_in)
    print('path out: ', path_out)
    print('')

    global case_name
    case_name = args.casename
    case = case_name[12:-3]
    print('')
    print('casename: ' + case_name)
    print('case: ' + case)
    print('')


    ''' determine time range '''
    global tmin, tmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = np.int(0)
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = np.int(10000)
    times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
             and np.int(name[:-3]) >= tmin and np.int(name[:-3]) <= tmax]
    times.sort()
    print('times', times)
    files = [str(t) + '.nc' for t in times]
    print('')

    # if args.kmax:
    #     kmax = np.int(args.kmax)
    # else:
    #     kmax = 0
    print('')

    return case, times, files



def define_geometry(case_name, files):
    print 'define geometry'
    global nx, ny, nz, dx, gw
    nml = simplejson.loads(open(os.path.join(path_in, case_name + '.in')).read())
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = np.ndarray(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']
    global dt_fields
    dt_fields = np.int(nml['fields_io']['frequency'])

    # set coordinates for plots
    if case_name[:21] == 'ColdPoolDry_triple_3D':
        global d, r_int
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        d = nml['init']['d']
        i_d = np.int(np.round(d / dx[0]))
        a = np.int(np.round(i_d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        idhalf = np.int(np.round(i_d / 2))
        r_int = np.int(np.round(np.sqrt(3.) / 6 * i_d))  # radius of inscribed circle
        # point of 3-CP collision (ic, jc)
        ic = np.int(np.round(nx / 2))
        jc = np.int(np.round(ny / 2))
        ic1 = ic - r_int
        ic2 = ic1
        ic3 = ic + (a - r_int)
        jc1 = jc - idhalf
        jc2 = jc + idhalf
        jc3 = jc
        ic_arr = [ic1, ic2, ic3]
        jc_arr = [jc1, jc2, jc3]
        icoll = ic
        jcoll = jc

    # print('ic, jc: ', ic_arr, jc_arr)

    ''' --- auxiliary arrays (since no Grid.pyx) ---'''
    global nx_, ny_, nz_
    # test file:
    var = read_in_netcdf_fields('u', os.path.join(path_fields, files[0]))
    [nx_, ny_, nz_] = var.shape
    x_half = np.empty((nx_), dtype=np.double, order='c')
    y_half = np.empty((ny_), dtype=np.double, order='c')
    z_half = np.empty((nz_), dtype=np.double, order='c')
    count = 0
    for i in xrange(nx_):
        x_half[count] = (i + 0.5) * dx[0]
        count += 1
    count = 0
    for j in xrange(ny_):
        y_half[count] = (j + 0.5) * dx[1]
        count += 1
    count = 0
    for i in xrange(nz_):
        z_half[count] = (i + 0.5) * dx[2]
        count += 1

    return icoll, jcoll, ic_arr, jc_arr, x_half, y_half, z_half

# ----------------------------------
def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    # shape = var.shape
    # data = np.ndarray(shape = var.shape)
    data = var[:,:,:]
    rootgrp.close()
    return data



if __name__ == '__main__':
    main()

