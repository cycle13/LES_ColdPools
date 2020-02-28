import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import netCDF4 as nc
import argparse
import json as simplejson
import os
import time
import scipy
import scipy.ndimage



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
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--s_crit")
    args = parser.parse_args()

    global cm_bwr, cm_grey, cm_vir, cm_hsv, cm_cw
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')
    cm_fall = plt.cm.get_cmap('winter')
    cm_summer = plt.cm.get_cmap('spring')
    cm_cw = plt.cm.get_cmap('coolwarm')

    nml, times = set_input_output_parameters(args)
    i0_center, j0_center, i0_coll, j0_coll, xmin_plt, xmax_plt, ymin_plt, ymax_plt = define_geometry(case_name, nml)
    kmax = kstar + 20

    ''' threshold for entropy '''
    if args.s_crit:
        s_crit = args.scrit
    else:
        s_crit = 5e-1
    print('threshold for ds=s-s_bg: ' + str(s_crit) + 'J/K')
    print('')

    ID = os.path.basename(path)
    if ID == '':
        ID = os.path.basename(path[:-1])
    print('id: ', ID)


    ''' background entropy '''
    try:
        s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
    except:
        s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '100.nc'))[ic_arr[2],jc_arr[0],5]
    s_bg = np.average(np.average(s0[:10,:10,:kmax], axis=0), axis=0)
    smax = np.amax(s0)
    smin = np.amin(s0)
    plot_geometry(s0, i0_coll, j0_coll)
    del s0
    print('sminmax', smin, smax)
    print ''

    # ''' create output file '''
    filename = 'CP_height_' + ID + '_sth' + str(s_crit) + '.nc'
    create_output_file(filename, s_crit, nx, ny, times)
    ''' define CP height & maximal updraft velocity'''
    # Output:
    #       CP_top[it, x, y]: 2D-field with CP-height for each xy-value
    #       w_max[0,it,x,y]:  2D-field with maximum value of w for each column
    #       w_max[1,it,x,y]:  height where maximum value of w for each column
    w_max, w_max_height, w_max_2d = compute_w_max(kmax, times)
    dump_output_file('w_max', 'timeseries', w_max, filename)
    dump_output_file('w_max_2d', 'fields_2D', w_max_2d[0,:,:,:], filename)
    dump_output_file('w_max_height', 'timeseries', w_max_height, filename)
    dump_output_file('w_max_height_2d', 'fields_2D', w_max_2d[1,:,:,:], filename)
    CP_top, CP_top_max = compute_CP_height_threshold(s_bg, s_crit, kmax, times)
    dump_output_file('CP_height_max', 'timeseries', CP_top_max, filename)
    dump_output_file('CP_height_2d', 'fields_2D', CP_top[:,:,:], filename)
    CP_top_grad, CP_top_grad_max, CP_top_gradgrad = compute_CP_height_gradient(i0_coll, j0_coll, kmax, times)
    dump_output_file('CP_height_gradient_max', 'timeseries', CP_top_grad_max, filename)
    dump_output_file('CP_height_gradient_2d', 'fields_2D', CP_top_grad[:,:,:], filename)
    print('')



    '''plot contour-figure of CP_top, w_max, height of w_max (xy-plane)'''
    print('plotting')
    filename = 'CP_height_' + ID + '_sth' + str(s_crit) + '.nc'
    root = nc.Dataset(os.path.join(path_out, filename), 'r')
    time_range = root.groups['timeseries'].variables['time'][:]
    CP_top = root.groups['timeseries'].variables['CP_height_max'][:]
    CP_top_2D = root.groups['fields_2D'].variables['CP_height_2d'][:,:,:]
    CP_top_gradient = root.groups['timeseries'].variables['CP_height_gradient_max'][:]
    CP_top_gradient_2D = root.groups['fields_2D'].variables['CP_height_gradient_2d'][:,:,:]
    w_max_height = root.groups['timeseries'].variables['w_max_height'][:]
    w_max_2D = root.groups['fields_2D'].variables['w_max_2d'][:,:,:]
    w_max_height_2D = root.groups['fields_2D'].variables['w_max_height_2d'][:,:,:]
    root.close()



    for it, t0 in enumerate(times):
        print('plot --- t: ', it, t0)
        print(xmin_plt, xmax_plt, ymin_plt, ymax_plt)
        figname = 'CP_height_t'+str(t0)+'.png'
        plot_contourf_xy(CP_top_2D[it, xmin_plt:xmax_plt, ymin_plt:ymax_plt],
                         w_max_2D[it, xmin_plt:xmax_plt, ymin_plt:ymax_plt],
                         w_max_height_2D[it, xmin_plt:xmax_plt, ymin_plt:ymax_plt], t0,
                         figname)
        figname = 'CP_height_gradient_t'+str(t0)+'.png'
        plot_contourf_xy(CP_top_gradient_2D[it, xmin_plt:xmax_plt, ymin_plt:ymax_plt],
                         w_max_2D[it, xmin_plt:xmax_plt, ymin_plt:ymax_plt],
                         w_max_height_2D[it, xmin_plt:xmax_plt, ymin_plt:ymax_plt], t0,
                         figname)
        plot_contourf_test_yz(smin, smax, CP_top_2D[it, :, :], CP_top_gradient_2D[it, :, :],
                              CP_top_gradgrad[it,:,:], kmax, j0_coll, t0)



    ''' plotting crosssections '''
    i1 = i0_center
    i2 = i0_coll
    j1 = i0_center
    j2 = j0_coll
    s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
    # plot_x_crosssections(s0, CP_top_2D, w_max_2D, j1, j2, times)
    # plot_x_crosssections(s0, CP_top_2D, w_max_2D, j1, j2, times)
    # plot_x_crosssections_CPtop_w(s0, CP_top_2D, w_max_2D, j1, j2, times)
    # plot_y_crosssections(s0, CP_top_2D, w_max_2D, i1, i2, times)
    # plot_y_crosssections_CPtop_w(s0, CP_top_2D, w_max_2D, i1, i2, times)
    #
    ''' plotting timeseries '''
    figname = 'CP_height_timeseries.png'
    plot_timeseries(CP_top_2D, CP_top, w_max_height, s0, i0_coll, j0_coll, time_range, figname)
    figname = 'CP_height_timeseries_gradient.png'
    plot_timeseries(CP_top_gradient_2D, CP_top_gradient, w_max_height, s0, i0_coll, j0_coll, time_range, figname)



    return


# ----------------------------------------------------------------------
def compute_w_max(kmax, times):
    print('--- compute w_max ---')
    nt = len(times)
    xmin = 0
    xmax = nx

    w_max_2d = np.zeros((2, nt, nx, ny))
    w_max = np.zeros((nt))
    w_max_height = np.zeros((nt))
    for it, t0 in enumerate(times):
        print('--- t: ', it, t0)
        w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0) + '.nc'))[xmin:xmax, xmin:xmax, :kmax+1]
        w_ = np.array(w, copy=True)
        w_[w_ < 0.5] = 0
        w_max_2d[0, it, :, :] = np.amax(w, axis=2)              # maximum value
        w_max_2d[1, it, :, :] = np.argmax(w_, axis=2)*dx[2]     # height of maximum value

        w_max[it] = np.average(w_max_2d[0, it, :, :])
        w_max_height[it] = np.average(w_max_2d[1, it, :, :])*dx[2]
    del w_, w

    return w_max, w_max_height, w_max_2d


def compute_CP_height_threshold(s_bg, s_crit, kmax, times):
    print('--- compute CP height by threshold ---')
    nt = len(times)
    xmin = 0
    xmax = nx
    # define CP height by threshold in entropy directly (s > s_crit)
    CP_top = np.zeros((nt, nx, ny), dtype=np.int)
    CP_top_max = np.zeros((nt))

    for it, t0 in enumerate(times):
        print('--- t: ', it, t0)
        s = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))[xmin:xmax, xmin:xmax, :kmax]
        s_diff = s - s_bg
        for i in range(nx):
            for j in range(ny):
                # for k in range(kmax-1, 0, -1):
                k = kmax
                while k >= 0:
                    k -= 1
                    # if np.abs(s_diff[i, j, k]) > s_crit and CP_top[it, i, j] == 0:
                    if np.abs(s[i, j, k]-s_bg[k]) > s_crit and CP_top[it, i, j] == 0:
                        # print('if-in', k)
                        CP_top[it, i, j] = (k+0.5)*dx[2]    # staggered grid for entropy
                        k = 0
        print(kmax, (kmax+0.5)*dx[2], np.amax(CP_top[it,:,:]), np.amin(CP_top[it,:,:]))
        CP_top_max[it] = np.amax(CP_top[it, :, :])



    return CP_top, CP_top_max


def compute_CP_height_gradient(i0_coll, j0_coll, kmax, times):
    print('--- compute CP height by gradient ---')
    nt = len(times)
    xmin = 0
    xmax = nx
    dzi = 1. / dx[2]

    # define CP height by inflection point (i.e., the zero-point the vertical gradient) of entropy
    CP_top_grad = np.zeros((nt, nx, ny), dtype=np.int)
    CP_top_gradgrad = np.zeros((nt, nx, ny), dtype=np.int)

    # irstar used for determining the averaging domain
    if case_name[:21] == 'ColdPoolDry_single_3D':
        i_eps_min = ic_arr[0] - irstar
        j_eps_min = jc_arr[0] - 2 * irstar
        deltax = 3 * irstar
        deltay = irstar
    elif case_name[:21] == 'ColdPoolDry_double_3D':
        i_eps_min = i0_coll - irstar
        j_eps_min = jc_arr[0] - 2 * irstar
        deltax = 3 * irstar
        deltay = irstar
    elif case_name[:21] == 'ColdPoolDry_triple_3D':
        i_eps_min = i0_coll - irstar
        j_eps_min = j0_coll - 2 * irstar
        deltax = 3 * irstar
        deltay = irstar

    for it, t0 in enumerate(times):
        fig, axes = plt.subplots(3, 3, figsize=(15, 12))
        jmin = 100
        jmax = 300

        print('--- t: ', it, t0)
        s = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))[xmin:xmax, xmin:xmax, :kmax+1]
        s_grad = np.zeros((nx, ny, kmax))
        s_gradgrad = np.zeros((nx, ny, kmax), dtype=np.double)

        # # s_grad[:, :, 0] = dzi * (s[:, :, 1] - 2 * s[:, :, 0] + s[:, :, 0])
        s_grad[:,:,0] = dzi * (s[:,:,1] - s[:,:,0])
        s_gradgrad[:,:,0] = dzi ** 2 * (s[:, :, 1] - 2 * s[:, :, 0] + s[:, :, 0])  # symmetric bcs: s[k=-1]=s[k=0]
        for k in range(1, kmax):
            s_grad[:,:,k] = dzi * (s[:,:,k + 1] - s[:,:,k])
            s_gradgrad[:,:,k] = dzi ** 2 * (s[:, :, k+1] - 2 * s[:, :, k] + s[:, :, k - 1])



        ax = axes[0, 0]
        levels = np.linspace(np.amin(s), np.amax(s), 1e2)
        a = ax.contourf(s[ic_arr[0], jmin:jmax, :kstar + 10].T, levels=levels)
        aux = np.ndarray((2,ny))
        for j in range(xmin,xmax):
            aux[0,j] = np.argmax(s_grad[ic_arr[0],j,:kmax])
            aux[1,j] = np.argmin(s_gradgrad[ic_arr[0],j,:kmax])
        ax.plot(aux[0,jmin:jmax], '-w')
        ax.plot(aux[1,jmin:jmax], '-r')
        plt.colorbar(a, ax=ax)
        ax.set_title('s')
        ax = axes[0, 1]
        cf = ax.contourf(s_grad[ic_arr[0], jmin:jmax, :kmax].T, levels=np.linspace(0, 0.09, 10))
        plt.colorbar(cf, ax=ax)
        ax.set_title('s_grad')
        ax = axes[0, 2]
        cf = ax.contourf(np.abs(s_gradgrad[ic_arr[0], jmin:jmax, :kmax]).T,
                         levels=np.linspace(-4e-4, 4e-4, 17), cmap=cm_bwr, extend="both")

        plt.colorbar(cf, ax=ax)
        ax.set_title('s_gradgrad')




        s_grad[s_grad <= 0.02] = 0.
        # s_gradgrad[s_grad <= 0.01] = 0.     # no good result in interior of CP
        k_init = np.argmax(s_gradgrad[:,:,:kmax], axis=2)
        k_end = np.argmin(s_gradgrad[:,:,:kmax], axis=2)
        # k_init_smooth = scipy.ndimage.uniform_filter(k_init, size=3, mode='constant')
        # k_end_smooth = scipy.ndimage.uniform_filter(k_end, size=3, mode='constant')
        activity = np.zeros((nx, ny), dtype=np.int)
        for i in range(xmin,xmax):
            for j in range(xmin,xmax):
                CP_top_grad[it,i,j] = np.argmax(s_grad[i,j,:kmax])
                min_aux = 9999.9
                # test continuity of k_init
                if np.abs(np.mean(k_init[i-1:i+2,j-1:j+2])-k_init[i,j]) <= 8.*2/9:
                    if np.abs(np.mean(k_end[i-1:i+2, j-1:j + 2]) - k_end[i, j]) <= 8.*2/9:
                        if k_init[i,j] < k_end[i,j]:
                            k = k_init[i,j]
                            while k <= k_end[i,j]:
                                if np.abs(s_gradgrad[i,j,k]) <= min_aux:
                                    min_aux = np.abs(s_gradgrad[i,j,k])
                                    CP_top_gradgrad[it,i,j] = k
                                k += 1
                # CP_top_gradgrad[it,i,j] = np.argmin(s_gradgrad[i,j,:kmax])





        ax = axes[1, 0]
        levels = np.linspace(np.amin(s), np.amax(s), 1e2)
        a = ax.contourf(s[ic_arr[0], jmin:jmax, :kstar + 10].T, levels=levels)
        ax.plot(CP_top_grad[it, ic_arr[0], jmin:jmax], '-w')
        ax.plot(CP_top_gradgrad[it, ic_arr[0], jmin:jmax], '-r')
        plt.colorbar(a, ax=ax)
        ax = axes[1,1]
        cf = ax.contourf(s_grad[ic_arr[0], jmin:jmax, :kmax].T, levels=np.linspace(0,0.09,10))
        ax.plot(CP_top_grad[it, ic_arr[0], jmin:jmax], '-w')
        plt.colorbar(cf, ax=ax)
        ax.set_title('s_grad (threshold=0.02)')
        ax = axes[1,2]
        cf = ax.contourf(s_gradgrad[ic_arr[0], jmin:jmax, :kmax].T,
                         levels=np.linspace(-4e-4, 4e-4, 17), cmap=cm_bwr, extend="both")
        cf.cmap.set_under('black')
        cf.cmap.set_over('black')
        plt.colorbar(cf, ax=ax)
        ax.plot(k_init[ic_arr[0], jmin:jmax], 'r', linewidth=1.)
        ax.plot(k_end[ic_arr[0], jmin:jmax], 'b', linewidth= 1.)
        ax.plot(CP_top_gradgrad[it, ic_arr[0], jmin:jmax], '-g', linewidth=2)
        ax.set_title('s_gradgrad (threshold=s_grad>0.02)')
        ax = axes[2,0]
        cf = ax.contourf(np.amin(np.abs(s_gradgrad[jmin:jmax, jmin:jmax,:]), axis=2))
        plt.colorbar(cf, ax=ax)
        ax.set_title('min(s_gradgrad)')
        ax = axes[2,1]
        cf = ax.contourf(CP_top_grad[it, jmin:jmax, jmin:jmax].T)
        plt.colorbar(cf, ax=ax)
        ax.set_title('CP top = argmax(gradient)')
        ax = axes[2,2]
        cf = ax.contourf(CP_top_gradgrad[it, jmin:jmax, jmin:jmax].T)
        plt.colorbar(cf, ax=ax)
        ax.set_title('CP top = argmin(gradgrad)')

        for ax in axes[0, :]:
            ax.set_ylim(0, kstar + 8)
            ax.set_xlim(100,200)
        for ax in axes[1, :]:
            ax.set_ylim(0, kstar + 8)
            ax.set_xlim(100,200)
        for ax in axes[2, :]:
            ax.set_aspect('equal')
        fig.tight_layout()
        fig.savefig(os.path.join(path_out_figs, 'gradient_yz_t' + str(t0) + '.png'))
        plt.close(fig)



    CP_top_grad = (CP_top_grad + 0.5) * dx[2]
    CP_top_gradgrad = (CP_top_gradgrad + 0.5) * dx[2]
    CP_top_max = np.amax(CP_top_grad[it, :, :])
    return CP_top_grad, CP_top_max, CP_top_gradgrad

# ----------------------------------------------------------------------
def plot_timeseries(CP_height_2d, CP_height, w_max_height, s0, i0_coll, j0_coll, time_range, figname):

    # rootgrp = nc.Dataset(os.path.join(path_out, filename), 'r')
    # ts_grp = rootgrp.groups['timeseries']
    # field_grp = rootgrp.groups['fields_2D']
    # time_range = ts_grp.variables['time'][:]
    # CP_height = ts_grp.variables['CP_height_max'][:]
    # w_max = ts_grp.variables['w_max'][:]
    # w_max_height = ts_grp.variables['w_max_height'][:]
    # CP_height_2d = field_grp.variables['CP_height_2d'][:, :]
    # rootgrp.close()

    cm = plt.cm.get_cmap('coolwarm')

    f1, axes = plt.subplots(3, 1, figsize=(6, 12))
    ax1 = axes[0]
    ax1.set_title('entropy')
    ax1.set_xlabel('x  (dx=' + str(dx[0]) + 'm)')
    ax1.set_ylabel('y  (dy=' + str(dx[1]) + 'm)')
    ax1.imshow(s0[:, :, 1].T)
    ax1.plot([0, nx], [j0_coll, j0_coll], 'k-', linewidth=2, label='center CP #1')
    ax1.plot([i0_coll, i0_coll], [0, ny], 'k-', linewidth=2, label='center CP #1')
    ax1.set_xlim(0, nx)
    ax1.set_ylim(0, ny)

    ax2 = axes[1]
    ax2.plot(time_range, CP_height, 'k', linewidth=2, label='max CP height')
    # ax2.plot(time_range, CP_height_2d[:, i0_coll, j0_coll], 'g', linewidth=2, label='CP height collision')
    ax2.legend(loc='best', fontsize=12)
    ax2.set_title('CP height')
    ax2.set_ylabel('CP height [m]')

    ax3 = axes[2]
    ax3.plot(time_range, w_max_height, 'k', linewidth=2, label='height max(w)')
    ax3.set_title('domain-max of w')
    ax3.legend(loc='best', fontsize=12)
    ax3.set_xlabel('y  (dy=' + str(dx[1]) + 'm)')
    ax3.set_ylabel('CP height [m]')


    plt.suptitle('CP height', fontsize=21)
    plt.subplots_adjust(bottom=0.05, right=.95, left=0.1, top=0.85, hspace=0.25)
    plt.savefig(os.path.join(path_out_figs, figname))
    plt.close()

    # del CP_height, CP_height_2d, w_max, w_max_height
    return

# ----------------------------------------------------------------------
def plot_contourf_test_yz(smin, smax, CP_top, CP_top_grad, CP_top_gradgrad, kmax, j0_coll, t0):
    ic1 = ic_arr[0]
    levels = np.arange(smin, smax, 0.1)
    s = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0)+'.nc'))

    fig, axes = plt.subplots(1, 4, figsize=(16, 5))
    ax1 = axes[0]
    a = ax1.contourf(s[ic1, :, :kmax].T, levels=levels)
    ax1.plot(CP_top[ic1,:], '-w')
    ax1.plot(CP_top_grad[ic1,:], '-r')
    ax1.plot(CP_top_grad[ic1,:], '-y')
    plt.colorbar(a, ax=ax1)
    ax2 = axes[1]
    cf = ax2.contourf(CP_top[:,:].T)
    plt.colorbar(cf, ax=ax2)
    ax2.set_title('CP top (threshold)')
    # ax2.contourf(s[ic1, :, :kmax].T, levels=levels)
    # ax2.plot(CP_top[ic1, :], 'w-')
    ax3 = axes[2]
    cf = ax3.contourf(CP_top_grad[:,:].T)
    plt.colorbar(cf, ax=ax3)
    ax3.set_title('CP top (gradient)')
    # ax3.contourf(s[ic1, :, :kmax].T, levels=levels)
    # ax3.plot(CP_top_grad[ic1, :], 'w-')
    ax4 = axes[3]
    cf = ax4.contourf(CP_top_gradgrad[:, :].T)
    plt.colorbar(cf, ax=ax4)
    ax4.set_title('CP top (gradgrad)')

    # fig.suptitle('t=' + str(t0) + 's')
    ax2.set_aspect('equal')
    ax3.set_aspect('equal')
    ax4.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, 'CP_height_yz_t' + str(t0) + '.png'))
    plt.close(fig)
    return


def plot_contourf_xy(CP_top, w_max, w_max_height, t0, fig_name):
    wmax = np.ceil(np.maximum(np.amax(w_max[:,:]), np.abs(np.amin(w_max[:,:]))))
    if wmax > 0.0:
        levels = np.arange(-wmax, wmax+0.1, 0.1)
    else:
        levels = np.arange(-1,2,1)

    fig, axes = plt.subplots(1,3, figsize=(16,5), sharey=True)
    ax1 = axes[0]
    # if CP_height not from file; it's not in metres but in number of levels
    # ax1.set_title('CP height  (max='+str(dx[2]*np.amax(CP_top))+'m)')
    # a = ax1.contourf(dx[2]*CP_top.T)
    a = ax1.contourf(CP_top.T)
    plt.colorbar(a, ax=ax1)
    ax1.set_title('CP height  (max='+str(np.amax(CP_top))+'m)')
    ax1.set_ylabel('y  (dy=' + str(dx[1]) + 'm)')

    ax2 = axes[1]
    b = ax2.contourf(w_max[:,:].T, levels=levels, cmap=cm_bwr)
    plt.colorbar(b, ax=ax2)
    ax2.set_title('max(w), (max='+str(np.round(np.amax(w_max[:,:]),2))+'m/s)')


    ax3 = axes[2]
    b = ax3.contourf(w_max_height[:, :].T)
    plt.colorbar(b, ax=ax3)
    ax3.set_title('height of max(w), (max='+str(dx[2]*np.int(np.amax(w_max_height[:,:])))+'m)')

    for ax in axes:
        ax.set_aspect('equal')
        ax.set_xlabel('x  (dx='+str(dx[0])+'m)')

    fig.suptitle('t='+str(t0)+'s', fontsize=15)
    fig.tight_layout()
    plt.subplots_adjust(bottom=0.1, right=.98, left=0.04, top=0.85, wspace=0.15)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return


def plot_geometry(s, i0, j0):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    ax1 = axes[0]
    ax2 = axes[1]
    cf1 = ax1.contourf(s[:,:,0].T)
    cf2 = ax2.contourf(s[:,:,kstar].T)
    eps = 100
    plt.colorbar(cf1, ax=ax1)
    plt.colorbar(cf2, ax=ax2)
    ax1.plot(ic_arr, jc_arr, 'yo')
    ax1.plot(i0, j0, 'ko')
    ax2.plot(i0, j0, 'ko')

    textprops = dict(facecolor='white', alpha=0.9, linewidth=0.)
    for i in range(len(ic_arr)):
        txt = 'CP'+str(i)
        ax1.text(ic_arr[i]+10, jc_arr[i]+10, txt, fontsize=18, bbox=textprops)

    ax1.set_xlim(eps, nx-eps)
    ax1.set_ylim(eps, ny-eps)
    ax2.set_xlim(eps, nx-eps)
    ax2.set_ylim(eps, ny-eps)
    ax1.plot([i0, i0],[0, ny], 'k-')
    ax2.plot([i0, i0],[0, ny], 'k-')
    ax1.plot([0, ny],[j0, j0], 'k-')
    ax2.plot([0, ny],[j0, j0], 'k-')
    ax1.set_aspect('equal')  # ax.set_aspect(1.0)
    ax2.set_aspect('equal')  # ax.set_aspect(1.0)
    ax1.set_xlabel('x  (dx=' + str(dx[0]) + 'm)')
    ax1.set_ylabel('y  (dy=' + str(dx[1]) + 'm)')
    ax1.set_title('k=0')
    ax2.set_title('k='+str(kstar))

    fig.tight_layout()
    fig.savefig(os.path.join(path_out_figs, 'CP_geometry.png'))
    plt.close(fig)

    return



def plot_x_crosssections(s0, CP_top, w_max, jp1, jp2, times):
    cm = plt.cm.get_cmap('coolwarm')

    f1, axes = plt.subplots(3, 1, figsize=(6, 12))
    ax1 = axes[0]
    ax1.set_title('entropy')
    ax1.set_xlabel('x  (dx=' + str(dx[0]) + 'm)')
    ax1.set_ylabel('y  (dy=' + str(dx[1]) + 'm)')
    ax2 = axes[1]
    ax2.set_title('x-crossection through collision point')
    ax2.set_ylabel('CP height [m]')
    ax3 = axes[2]
    ax3.set_title('x-crossection through center of coldpool #1')
    ax3.set_xlabel('y  (dy='+str(dx[1])+'m)')
    ax3.set_ylabel('CP height [m]')
    ax1.imshow(s0[:,:,1].T)
    ax1.set_xlim(0,nx)
    ax1.set_ylim(0,ny)
    ax1.plot([0,nx],[jp1,jp1], 'k', linewidth=2, label='collision center')
    ax1.plot([0,nx],[jp2,jp2], 'k--', linewidth=2, label='center CP #1')
    ax1.legend()
    for it,t0 in enumerate(times):
        count_color = np.double(it) / len(times)
        ax2.plot(CP_top[it, :, jp1], color=cm(count_color), label='t=' + str(t0))
        ax3.plot(CP_top[it, :, jp2], color=cm(count_color), label='t=' + str(t0))
    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
               fancybox=True, shadow=True, ncol=5)
    plt.suptitle('CP height', fontsize=21)
    plt.savefig(os.path.join(path_out_figs, 'CP_height_xz_plane.png'))
    plt.close()
    return



def plot_y_crosssections(s0, CP_top, w_max, ip1, ip2, times):
    cm = plt.cm.get_cmap('coolwarm')

    f1, axes = plt.subplots(3, 1, figsize=(6, 12))
    ax1 = axes[0]
    ax1.set_title('entropy')
    ax1.set_xlabel('x  (dx=' + str(dx[0]) + 'm)')
    ax1.set_ylabel('y  (dy=' + str(dx[1]) + 'm)')
    ax2 = axes[1]
    ax2.set_title('y-crossection through 2 CP-collision point')
    ax2.set_ylabel('CP height [m]')
    ax3 = axes[2]
    ax3.set_title('y-crossection through center of coldpool #3')
    ax3.set_xlabel('y  (dy='+str(dx[1])+'m)')
    ax3.set_ylabel('CP height [m]')
    ax1.imshow(s0[:,:,1].T)
    ax1.set_xlim(0,nx)
    ax1.set_ylim(0,ny)
    ax1.plot([ip1,ip1],[0,ny], 'k', linewidth=2, label='2-CP collision')
    ax1.plot([ip2,ip2],[0,ny], 'k--', linewidth=2, label='center of CP #3')
    ax1.legend()
    for it,t0 in enumerate(times):
        count_color = np.double(it) / len(times)
        ax2.plot(CP_top[it, ip1, :], color=cm(count_color), label='t=' + str(t0))
        ax3.plot(CP_top[it, ip2, :], color=cm(count_color), label='t=' + str(t0))
    # # ax3.legend( )
    # # Put a legend below current axis
    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
               fancybox=True, shadow=True, ncol=5)
    plt.suptitle('CP height', fontsize=21)
    plt.savefig(os.path.join(path_out_figs, 'CP_height_yz_plane.png'))
    plt.close()
    return




def plot_x_crosssections_CPtop_w(s0, CP_top, w_max, jp1, jp2, times):
    cm = plt.cm.get_cmap('coolwarm')

    f1, axes = plt.subplots(3, 1, figsize=(6, 12))
    ax1 = axes[0]
    ax1.set_title('entropy')
    ax1.set_xlabel('x  (dx='+str(dx[0])+'m)')
    ax1.set_ylabel('y  (dy='+str(dx[1])+'m)')
    ax2 = axes[1]
    ax2.set_title('CP height')
    ax2.set_ylabel('CP height [m]')
    ax3 = axes[2]
    ax3.set_title('max w')
    ax3.set_xlabel('y  (dy='+str(dx[1])+'m)')
    ax3.set_ylabel('max(w) [m/s]')
    ax1.imshow(s0[:,:,1].T)
    ax1.set_xlim(0,nx)
    ax1.set_ylim(0,ny)
    ax1.plot([0,nx],[jp1,jp1], 'k', linewidth=2, label='jp1')
    for it,t0 in enumerate(times):
        count_color = np.double(it) / len(times)
        ax2.plot(CP_top[it, :, jp1], color=cm(count_color), label='t=' + str(t0))
        ax3.plot(w_max[it, :, jp1], color=cm(count_color), label='t=' + str(t0))
    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
               fancybox=True, shadow=True, ncol=5)
    plt.suptitle('x-crossection through collision point', fontsize=21)
    plt.savefig(os.path.join(path_out_figs, 'CP_height_wmax_xz_plane.png'))
    plt.close()
    return



def plot_y_crosssections_CPtop_w(s0, CP_top, w_max, ip1, ip2, times):
    cm = plt.cm.get_cmap('coolwarm')

    f1, axes = plt.subplots(3, 1, figsize=(6, 12))
    ax1 = axes[0]
    ax1.set_title('entropy')
    ax1.set_xlabel('x  (dx='+str(dx[0])+'m)')
    ax1.set_ylabel('y  (dy='+str(dx[1])+'m)')
    ax2 = axes[1]
    ax2.set_title('CP height')
    ax2.set_ylabel('CP height [m]')
    ax3 = axes[2]
    ax3.set_title('max w')
    ax3.set_xlabel('y  (dy='+str(dx[1])+'m)')
    ax3.set_ylabel('max(w) [m/s]')
    ax1.imshow(s0[:,:,1].T)
    ax1.set_xlim(0,nx)
    ax1.set_ylim(0,ny)
    ax1.plot([ip1,ip1],[0,ny], 'k', linewidth=2, label='ip1')
    for it,t0 in enumerate(times):
        count_color = np.double(it) / len(times)
        ax2.plot(CP_top[it, ip1, :], color=cm(count_color), label='t=' + str(t0))
        ax3.plot(w_max[it, ip1, :], color=cm(count_color), label='t=' + str(t0))
    # # ax3.legend( )
    # # Put a legend below current axis
    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
               fancybox=True, shadow=True, ncol=5)
    plt.suptitle('y-crossection through 2 CP-collision point', fontsize=21)
    plt.savefig(os.path.join(path_out_figs, 'CP_height_wmax_yz_plane.png'))
    plt.close()

    return



# ----------------------------------------------------------------------

def set_input_output_parameters(args):
    print('--- set input parameters ---')
    global case_name
    global path, path_fields, path_out, path_out_figs
    global files

    if args.casename:
        case_name = args.casename
    else:
        case_name = 'ColdPoolDry_triple_3D'

    path = args.path
    if os.path.exists(os.path.join(path, 'fields')):
        path_fields = os.path.join(path, 'fields')
    path_out = os.path.join(path, 'data_analysis')
    print(path_out)
    path_out_figs = os.path.join(path, 'figs_CP_height')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print ''
    print('path data: ' + path_out)
    print('path figs: ' + path_out_figs)
    print ''

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dV, gw
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = np.zeros(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx[0] * dx[1] * dx[2]
    print('nx, ny, nz', nx, ny, nz)


    ''' determine file range '''
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = np.int(100)
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = np.int(10000)
    times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
             and np.int(name[:-3]) >= tmin and np.int(name[:-3]) <= tmax]
    # times = [np.int(name[:-3]) for name in files]
    times.sort()
    print('times', times)
    files = [str(t) + '.nc' for t in times]
    print('files', files)

    print ''

    return nml, times




def define_geometry(case_name, nml):
    print('--- define geometry ---')
    global x_half, y_half, z_half
    global ic_arr, jc_arr
    global i0_center, j0_center
    global rstar, irstar, zstar, kstar

    # test file:
    var = read_in_netcdf_fields('u', os.path.join(path_fields, files[0]))
    [nx_, ny_, nz_] = var.shape
    del var
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

    # set coordinates for plots
    if case_name[:21] == 'ColdPoolDry_single_3D':
        rstar = nml['init']['r']
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        kstar = np.int(np.round(zstar / dx[2]))
        try:
            ic = nml['init']['ic']
            jc = nml['init']['jc']
            # print('(ic,jc) from nml')
        except:
            ic = np.int(nx/2)
            jc = np.int(ny/2)
            # print('(ic,jc) NOT from nml')
        ic_arr = np.zeros(1)
        jc_arr = np.zeros(1)
        ic_arr[0] = np.int(ic)
        jc_arr[0] = np.int(jc)
    elif case_name[:21] == 'ColdPoolDry_double_2D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        isep = 4 * irstar
        ic1 = np.int(nx / 3)  # np.int(Gr.dims.ng[0] / 3)
        ic2 = ic1 + isep
        jc1 = np.int(ny / 2)
        jc2 = jc1
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name[:21] == 'ColdPoolDry_double_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        kstar = np.int(np.round(zstar / dx[2]))
        isep = 4 * irstar
        jsep = 0
        ic1 = np.int(np.round((nx + 2 * gw) / 3)) - gw
        jc1 = np.int(np.round((ny + 2 * gw) / 2)) - gw
        ic2 = ic1 + isep
        jc2 = jc1 + jsep
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name[:21] == 'ColdPoolDry_triple_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        kstar = np.int(np.round(zstar / dx[2]))
        try:
            d = nml['init']['d']
        except:
            d = 10000
        dhalf = np.int(np.round(ny / 4))
        id = np.round(d / dx[0])
        idhalf = np.int(np.round(id / 2))
        a = np.int(np.round(id * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        r_int = np.int(np.sqrt(3.) / 6 * id)  # radius of inscribed circle
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




    ''' plotting parameters '''
    if case_name[:21] == 'ColdPoolDry_single_3D':
        i0_coll = 1
        j0_coll = 1
        i0_center = np.int(ic_arr[0])
        j0_center = np.int(jc_arr[0])
        if nx > 200:
            xmin_plt = 100
        else:
            xmin_plt = 0
        xmax_plt = nx - xmin_plt
        ymin_plt = xmin_plt
        ymax_plt = xmax_plt
    elif case_name[:21] == 'ColdPoolDry_double_3D':
        i0_coll = 0.5 * (ic_arr[0] + ic_arr[1])
        i0_center = ic_arr[0]
        j0_coll = jc_arr[0]
        j0_center = jc_arr[0]
        # domain boundaries for plotting
        xmin_plt = 30
        xmax_plt = 230
        ymin_plt = xmin_plt
        ymax_plt = xmax_plt
    elif case_name[:21] == 'ColdPoolDry_triple_3D':
        i0_coll = ic
        i0_center = ic_arr[0]
        j0_coll = jc
        j0_center = jc_arr[0]
        # domain boundaries for plotting
        xmin_plt = 0
        xmax_plt = nx
        ymin_plt = xmin_plt
        ymax_plt = xmax_plt

    print ''
    return i0_center, j0_center, i0_coll, j0_coll, xmin_plt, xmax_plt, ymin_plt, ymax_plt


# ----------------------------------
def theta_s(s):
    # T_tilde = 298.15
    # thetas_c(s, qt){
    #     return T_tilde * exp((s - (1.0 - qt) * sd_tilde - qt * sv_tilde) / cpm_c(qt));
    # }
    # cpm_c(qt){
    #     return (1.0 - qt) * cpd + qt * cpv;
    # }
    # parameters from pycles
    T_tilde = 298.15
    sd_tilde = 6864.8
    cpd = 1004.0
    th_s = T_tilde * np.exp( (s - sd_tilde)/cpd )
    return th_s

# ----------------------------------
def create_output_file(filename, s_crit, nx, ny, times):
    # output for each CP:
    # - min, max (timeseries)
    # - CP height (field; max=timeseries)
    # - (ok) CP rim (field)
    nt = len(times)
    print('create output file: ', os.path.join(path_out, filename))
    print('size: ', nz, nt)

    rootgrp = nc.Dataset(os.path.join(path_out, filename), 'w', format='NETCDF4')

    ts_grp = rootgrp.createGroup('timeseries')
    ts_grp.createDimension('nt', nt)
    var = ts_grp.createVariable('time', 'f8', ('nt'))
    var.units = "s"
    var[:] = times

    var = ts_grp.createVariable('CP_height_max', 'f8', ('nt'))
    var.long_name = 'domain max CP height (from threshold)'
    var.units = "m"
    var = ts_grp.createVariable('CP_height_gradient_max', 'f8', ('nt'))
    var.long_name = 'domain max CP height (from gradient)'
    var.units = "m"

    var = ts_grp.createVariable('w_max', 'f8', ('nt'))
    var.units = "m/s"
    var = ts_grp.createVariable('w_max_height', 'f8', ('nt'))
    var.units = "m/s"

    field_grp = rootgrp.createGroup('fields_2D')
    field_grp.createDimension('nt', nt)
    field_grp.createDimension('nx', nx)
    field_grp.createDimension('ny', ny)
    var = field_grp.createVariable('CP_height_2d', 'f8', ('nt', 'nx', 'ny'))
    var.units = "m"
    var = field_grp.createVariable('CP_height_gradient_2d', 'f8', ('nt', 'nx', 'ny'))
    var.units = "m"
    var = field_grp.createVariable('w_max_2d', 'f8', ('nt', 'nx', 'ny'))
    var.units = "m/s"
    var = field_grp.createVariable('w_max_height_2d', 'f8', ('nt', 'nx', 'ny'))
    var.units = "m"

    rootgrp.close()
    print ''
    return



def dump_output_file(var_name, group_name, var_in, filename):
    # print('dumping', os.path.join(path_out, filename))
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



def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:,:,:]
    rootgrp.close()
    return data


if __name__ == '__main__':
    main()
