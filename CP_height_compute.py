import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import netCDF4 as nc
import argparse
import json as simplejson
import os
import time


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

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')
    cm_fall = plt.cm.get_cmap('winter')
    cm_summer = plt.cm.get_cmap('spring')

    nml, times = set_input_output_parameters(args)
    i0_center, j0_center, i0_coll, j0_coll, xmin_plt, xmax_plt, ymin_plt, ymax_plt = define_geometry(case_name, nml)
    kmax = kstar + 10

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
    plot_geometry(s0, i0_center, j0_center)
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
    CP_top_grad, CP_top_grad_max = compute_CP_height_gradient(i0_coll, j0_coll, kmax, times)
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
        # plot_contourf_test_yz(s, smin, smax, CP_top_2D[it, :, :], CP_top_grad[it, :, :], kmax, t0)



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
    print('s_bg', s_bg)

    for it, t0 in enumerate(times):
        print('--- t: ', it, t0)
        s = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))[xmin:xmax, xmin:xmax, :kmax]
        s_diff = s - s_bg
        print('s_diff', s_diff.shape)
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
        # print(CP_top[it,:,:])
        print(kmax, (kmax+0.5)*dx[2], np.amax(CP_top[it,:,:]), np.amin(CP_top[it,:,:]))
        CP_top_max[it] = np.amax(CP_top[it, :, :])



    return CP_top, CP_top_max




def compute_CP_height_gradient(i0_coll, j0_coll, kmax, times):
    print('--- compute CP height by threshold ---')
    nt = len(times)
    xmin = 0
    xmax = nx
    dzi = 1. / dx[2]

    # define CP height by inflection point (i.e., the zero-point the vertical gradient) of entropy
    CP_top_grad = np.zeros((nt, nx, ny))
    CP_top_max = np.zeros((nt))

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
        print('--- t: ', it, t0)
        s = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))[xmin:xmax, xmin:xmax, :kmax+1]
        s_grad = np.zeros((nx, ny, kmax))
        s_gradgrad = np.zeros((nx, ny, kmax))

        # # s_grad[:, :, 0] = dzi * (s[:, :, 1] - 2 * s[:, :, 0] + s[:, :, 0])
        # s_grad[:,:,0] = dzi * (s[:,:,1] - s[:,:,0])
        s_gradgrad[:,:,0] = dzi ** 2 * (s[:, :, 1] - 2 * s[:, :, 0] + s[:, :, 0])  # symmetric bcs: s[k=-1]=s[k=0]
        for k in range(1, kmax):
            # s_grad[:,:,k] = dzi * (s[:,:,k + 1] - s[:,:,k])
            s_gradgrad[:,:,k] = dzi ** 2 * (s[:, :, k + 1] - 2 * s[:, :, k] + s[:, :, k - 1])

        # epsilon = np.average(np.average(s_grad[i_eps_min:i_eps_min + deltax, j_eps_min:j_eps_min + deltay, :], axis=0),
        #                      axis=0)
        for i in range(xmin,xmax):
            for j in range(xmin,xmax):
                CP_top_grad[it,i,j] = (np.argmin(np.abs(s_gradgrad[i,j,:]))+0.5)*dx[2]
                # for k in range(kmax, 0, -1):
                #     if s_grad[i, j, k] > epsilon[k] + 1e-3 and CP_top_grad[it, i, j] == 0:
                #         CP_top_grad[it, i, j] = k
        CP_top_max[it] = np.amax(CP_top_grad[it, :, :])

    return CP_top_grad, CP_top_max



# ----------------------------------------------------------------------
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
# ----------------------------------------------------------------------
def plot_contourf_test_yz(s, smin, smax, CP_top, CP_top_grad, kmax, t0):
    ic1 = ic_arr[0]
    levels = np.arange(smin, smax, 0.1)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    ax1 = axes[0]
    a = ax1.contourf(s[:, :, 1].T, levels=levels)
    ax1.plot([ic1, ic1], [0, ny], 'k-')
    plt.colorbar(a, ax=ax1)
    ax2 = axes[1]
    ax2.contourf(s[ic1, :, :kmax].T, levels=levels)
    ax2.plot(CP_top[ic1, :], 'w-')
    ax2.set_title('CP top (threshold)')
    ax3 = axes[2]
    ax3.contourf(s[ic1, :, :kmax].T, levels=levels)
    ax3.plot(CP_top_grad[ic1, :], 'w-')
    ax3.set_title('CP top (gradient)')
    fig.suptitle('t=' + str(t0) + 's')
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
    a1 = ax1.contourf(s[:,:,kstar].T)
    a2 = ax2.contourf(s[:,:,0].T)
    eps = 20
    plt.colorbar(a1, ax=ax1)
    plt.colorbar(a2, ax=ax2)
    ax1.set_xlim(i0-eps, i0+eps)
    ax1.set_ylim(j0-eps, j0+eps)
    ax2.set_xlim(i0-eps, i0+eps)
    ax2.set_ylim(j0-eps, j0+eps)
    ax1.plot(i0, j0, 'ko')
    ax2.plot(i0, j0, 'ko')
    ax1.plot([i0, i0],[0, ny], 'k-')
    ax2.plot([i0, i0],[0, ny], 'k-')
    ax1.plot([0, ny],[j0, j0], 'k-')
    ax2.plot([0, ny],[j0, j0], 'k-')
    ax1.set_aspect('equal')  # ax.set_aspect(1.0)
    ax2.set_aspect('equal')  # ax.set_aspect(1.0)
    ax1.set_xlabel('x  (dx=' + str(dx[0]) + 'm)')
    ax1.set_ylabel('y  (dy=' + str(dx[1]) + 'm)')

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
        d = np.int(np.round(ny / 2))
        dhalf = np.int(np.round(ny / 4))
        a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        ic1 = np.int(np.round(a / 2))  # + gw
        ic2 = ic1
        ic3 = ic1 + np.int(np.round(a))
        jc1 = np.int(np.round(d / 2))  # + gw
        jc2 = jc1 + d
        jc3 = jc1 + np.int(np.round(d / 2))
        ic_arr = [ic1, ic2, ic3]
        jc_arr = [jc1, jc2, jc3]
        i0_coll = 0.5 * (ic_arr[0] + ic_arr[2])
        i0_center = ic_arr[0]
        j0_coll = jc_arr[2]
        j0_center = jc_arr[0]

        isep = dhalf

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
    print('create output file in: ', path_out)
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
    var.units = "m"
    var = ts_grp.createVariable('w_max_height', 'f8', ('nt'))
    var.units = "m"

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
