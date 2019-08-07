import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import netCDF4 as nc
import argparse
import json as simplejson
import os
import sys


label_size = 12
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 15
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 8
# plt.rcParams['lines.linewidth'] = 3
# plt.rcParams['grid.linewidth'] = 20
# plt.rcParams['xtick.major.size'] = 8.5
# plt.rcParams['xtick.minor.size'] = 5
# plt.rcParams['ytick.major.size'] = 8.5
# plt.rcParams['ytick.minor.size'] = 5
# plt.rcParams['xtick.major.width'] = 2
# plt.rcParams['xtick.minor.width'] = 1.5
# plt.rcParams['ytick.major.width'] = 2
# plt.rcParams['ytick.minor.width'] = 1.5
# plt.rcParams['pdf.fonttype'] = 42         # Output Type 3 (Type3) or Type 42 (TrueType)

def main():

    ''' set paths & parameters '''
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='LES_CP')
    # parser.add_argument("casename")
    # parser.add_argument("path_root")
    # parser.add_argument("path1")
    # parser.add_argument("path2")
    # parser.add_argument("path3")
    parser.add_argument("--path_out")
    parser.add_argument("--dTh")
    parser.add_argument("--zstar")
    parser.add_argument("--rstar")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--kmax")
    parser.add_argument("--everysecond")
    args = parser.parse_args()

    path_root = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise'
    dTh, zstar, rstar = set_input_parameters(args, path_root)
    ID = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
    path1 = os.path.join(path_root, 'run2_dx100m', ID)
    path2 = os.path.join(path_root, 'run3_dx50m', ID)
    path3 = os.path.join(path_root, 'run4_dx25m', ID)
    print ''
    print 'paths:'
    print 'path 1: ' + path1
    print 'path 2: ' + path2
    print 'path 3: ' + path3
    print ''
    print 'path figs: ' + path_out_figs
    print ''

    define_geometry(args, path1, path2, path3)
    # ic_arr, jc_arr = define_geometry(nml)
    if args.everysecond:
        every_second = np.int(args.everysecond)
    else:
        every_second = 1
    print('every_second', every_second)

    file_name = 'stats_radial_averaged.nc'
    file1 = nc.Dataset(os.path.join(path1, 'data_analysis', file_name), 'r')
    file2 = nc.Dataset(os.path.join(path2, 'data_analysis', file_name), 'r')
    file3 = nc.Dataset(os.path.join(path3, 'data_analysis', file_name), 'r')

    dx1 = file1.groups['dimensions'].dimensions['dx'].size
    dx2 = file2.groups['dimensions'].dimensions['dx'].size
    dx3 = file3.groups['dimensions'].dimensions['dx'].size
    dz1 = file1.groups['dimensions'].dimensions['dz'].size
    dz2 = file2.groups['dimensions'].dimensions['dz'].size
    dz3 = file3.groups['dimensions'].dimensions['dz'].size
    # if dx1 < dx2:   # want to have file2 higher resolution
    #     file_aux = file2
    #     file2 = file1
    #     file1 = file_aux
    #     del file_aux
    print('dx: ', dx1_nml[0], dx2_nml[0], dx3_nml[0])
    print('dz: ', dx1_nml[2], dx2_nml[2], dx3_nml[2])
    print('nx: ', nx1[0], nx2[0], nx3[0])
    print('nz: ', nx1[2], nx2[2], nx3[2])
    print('dt fields: ', dt_fields1, dt_fields2, dt_fields3)


    time1 = file1.groups['timeseries'].variables['time'][:]
    time2 = file2.groups['timeseries'].variables['time'][:]
    time3 = file3.groups['timeseries'].variables['time'][:]
    print('times: ')
    print time1
    print time2
    print time3
    nt = np.minimum(np.minimum(len(time1), len(time2)), len(time3))
    time1 = time1[:nt]
    time2 = time2[:nt]
    time3 = time3[:nt]
    if time1.any() != time2.any() or time1.any() != time3.any() or time2.any() != time3.any():
        print'problem with times'
        sys.exit()

    r1 = file1.groups['stats'].variables['r'][:]
    r2 = file2.groups['stats'].variables['r'][:]
    r3 = file3.groups['stats'].variables['r'][:]
    ri1 = file1.groups['stats'].variables['ri'][:]
    ri2 = file2.groups['stats'].variables['ri'][:]

    try:
        krange1 = file1.groups['dimensions'].variables['krange'][:]
        krange2 = file2.groups['dimensions'].variables['krange'][:]
        krange3 = file3.groups['dimensions'].variables['krange'][:]
    except:
        krange1 = np.arange(kmax)
        krange2 = np.arange(kmax)
        krange3 = np.arange(kmax)
    zrange1 = krange1*dz1
    zrange2 = krange2*dz2
    zrange3 = krange3*dz3
    print ''
    # print 'kranges:'
    # print krange1
    # print krange2
    # print krange3
    zmax = np.minimum(np.minimum(zrange1[-1], zrange2[-1]), zrange3[-1])

    global dz12, dz13, dz23
    dz12 = np.double(dz1) / dz2
    dz13 = np.double(dz1) / dz3
    dz23 = np.double(dz2) / dz3
    print 'dz: ', dz1, dz2, dz3
    print 'dz12, dz13, dz23: ', dz12, dz13, dz23
    if dz2 > dz1 or dz3 > dz2:
        print('Stopping bcs. opposite order of resultion required')
        sys.exit()
    kmax2_ = np.int(np.double(kmax) / dz23)
    kmax1_ = np.int(np.double(kmax) / dz13)
    krange1 = krange1[:kmax1_]
    krange2 = krange2[0::np.int(dz23)][:kmax2_]
    krange3 = krange3[0::np.int(dz13)]
    kmax_ = np.minimum(np.minimum(len(krange1), len(krange2)), len(krange3))
    krange1 = krange1[:kmax_]
    krange2 = krange2[:kmax_]
    krange3 = krange3[:kmax_]
    # if dz1 > dz2:
    #     print 'dz1 > dz2'
    #     dz12 = np.double(dz1) / dz2
    #     kmax_ = np.int(np.minimum(np.int(np.round(kmax / dz12)), len(krange2)/dz12))
    #     krange1 = krange1[:kmax_]
    #     krange2 = krange2[0::np.int(dz12)]
    # elif dz1 < dz2:
    #     print 'dz1 < dz2'
    #     dz12 = np.double(dz2) / dz1
    #     kmax_ = np.int(np.round(kmax / dz12))
    #     krange1 = krange1[0::2]
    #     krange2 = krange2[:kmax_]
    # else:
    #     dz12 = 1.
    #     kmax_ = kmax


    print 'kmax: ', kmax, kmax1_, kmax2_
    print 'zranges:'
    print krange1*dz1
    print krange2*dz2
    print krange3*dz3
    print('')



    ''' plotting rim width '''
    file_name = 'stats_radial_averaged_rimwidth.nc'
    fig_name = 'rim_geometry'
    plot_rim_width(path1, path2, path3, time1, file_name, fig_name, every_second, kmax)


    # ''' plotting profiles (radial dependence) '''
    # var_list = ['w', 'v_rad', 's']
    # ncol = len(var_list)
    # grp_stats1 = file1.groups['stats'].variables
    # grp_stats2 = file2.groups['stats'].variables
    # grp_stats3 = file3.groups['stats'].variables
    # rmax = 9e3
    # irmax1 = np.where(r1 == rmax)[0]
    # irmax2 = np.where(r2 == rmax)[0]
    # irmax3 = np.where(r3 == rmax)[0]
    # # irmax1 = -1
    # # irmax2 = -1
    # trange = np.arange(0, 2500, 600)
    # trange[0] = 100
    # trange = [100, 600, 1200, 2400]
    # print ''
    # print r1[-1], r2[-1]
    # print 'irmax: ', irmax1, irmax2
    # print 'rmax', r1[irmax1], r2[irmax2]
    # print('trange: ', trange)
    # print ''
    #
    # cm1 = cm.copper
    # cm2 = cm.winter
    # cm3 = cm.spring
    # ls_list = ['-', '--', '-', '--', '-']
    # ls_list = ['-', 'densely dashed', 'dashdot', 'densely dotted', '-']
    # dash_list = [(1,0,0,0), (7,1,7,1), (5,3,5,3), (5,2,1,2), (4,2,1,2)]
    # cm0 = plt.cm.get_cmap('coolwarm')
    # c1 = cm0(1.)
    # c2 = cm0(.75)
    # c3 = cm0(.01)
    # colorlist = [c1, c2, c3]
    #
    #
    # ''' at same k '''
    # for k0 in range(kmax_):
    #     k1 = np.int(krange1[k0])
    #     k2 = np.int(krange2[k0])
    #     k3 = np.int(krange3[k0])
    #     print('z1, z2, z3: ', (k1+1)*dz1, (k2+1)*dz2, (k3+1)*dz3)
    #
    #     fig_name = 'radial_average_k' + str(np.int(k1)) + '.png'
    #     fig, axes = plt.subplots(1, ncol, sharey='none', figsize=(5 * ncol, 5))
    #     for i, ax in enumerate(axes):
    #         var1 = grp_stats1[var_list[i]][:, :, k1]
    #         var2 = grp_stats2[var_list[i]][:, :, k2]
    #         var3 = grp_stats3[var_list[i]][:, :, k3]
    #         for c,t0 in enumerate(trange):
    #             it = np.int(t0/dt_fields1)
    #             # print('- t='+str(t0), it)
    #             # count_color = np.double(c) / len(trange)
    #             # ax.plot(r1[:irmax1], var1[it, :irmax1], color=cm1(count_color), linewidth=3,
    #             #             label='dx='+str(dx1) + ', t=' + str(t0))
    #             ax.plot(r1[:irmax1], var1[it, :irmax1], color=colorlist[0],
    #                     linestyle=ls_list[c], dashes=dash_list[c], linewidth=2,
    #                     label='dx=' + str(dx1) + ', t=' + str(t0))
    #
    #         for c, t0 in enumerate(trange):
    #             it = np.int(t0/dt_fields1)
    #             # print('- t='+str(t0), it)
    #             # count_color = np.double(c) / len(trange)
    #             # ax.plot(r2[:irmax2], var2[it, :irmax2], '-', color=cm2(count_color), linewidth=3,
    #             #             label='dx='+str(dx2))
    #             ax.plot(r2[:irmax2], var2[it, :irmax2], '-', color=colorlist[1],
    #                     linestyle=ls_list[c],  dashes=dash_list[c], linewidth=2,
    #                         label='dx='+str(dx2))
    #         for c, t0 in enumerate(trange):
    #             it = np.int(t0/dt_fields1)
    #             # print('- t='+str(t0), it)
    #             # count_color = np.double(c) / len(trange)
    #             # ax.plot(r3[:irmax3], var3[it, :irmax3], '-', color=cm3(count_color), linewidth=3,
    #             #             label='dx='+str(dx2))
    #             ax.plot(r3[:irmax3], var3[it, :irmax3], '-', color=colorlist[2],
    #                     linestyle=ls_list[c],  dashes=dash_list[c], linewidth=2,
    #                         label='dx='+str(dx2))
    #
    #         ax.set_title(var_list[i])
    #         ax.set_xlabel('radius r  [m]')
    #         ax.set_ylabel(var_list[i])
    #     axes[1].legend(loc='upper center', bbox_to_anchor=(0, -0.2),
    #                    fancybox=False, ncol=3, fontsize=8)
    #     plt.subplots_adjust(bottom=0.3, right=.98, left=0.07 , top=0.9, wspace=0.25)
    #     fig.suptitle('radially averaged variables   (k='+str(k0)+')')
    #     fig.savefig(os.path.join(path_out_figs, fig_name))
    #     plt.close(fig)
    #
    #
    #
    #
    #
    # print('')
    # print('')
    # ''' at same height '''
    # for k1 in range(kmax_):
    #     k2 = (k1+1)*dz12-1
    #     k3 = (k1+1)*dz13-1
    #     z1 = (k1+1)*dz1
    #     z2 = (k2+1)*dz2
    #     z3 = (k3+1)*dz3
    #     print(k1, k2, k3)
    #     print('z1, z2, z3: ', z1, z2, z3)
    #
    #     fig_name = 'radial_average_z' + str(z1) + '.png'
    #     fig, axes = plt.subplots(1, ncol, sharey='none', figsize=(6 * ncol, 5))
    #     for i, ax in enumerate(axes):
    #         var1 = grp_stats1[var_list[i]][:, :, k1]
    #         var2 = grp_stats2[var_list[i]][:, :, k2]
    #         var3 = grp_stats3[var_list[i]][:, :, k3]
    #         for c,t0 in enumerate(trange):
    #             it = np.int(t0/dt_fields1)
    #             ax.plot(r1[:irmax1], var1[it, :irmax1], color=colorlist[0],
    #                     linestyle=ls_list[c], dashes=dash_list[c], linewidth=2,
    #                     label='dx=' + str(dx1) + ', t=' + str(t0))
    #
    #         for c, t0 in enumerate(trange):
    #             it = np.int(t0/dt_fields1)
    #             ax.plot(r2[:irmax2], var2[it, :irmax2], '-', color=colorlist[1],
    #                     linestyle=ls_list[c],  dashes=dash_list[c], linewidth=2,
    #                         label='dx='+str(dx2))
    #         for c, t0 in enumerate(trange):
    #             it = np.int(t0/dt_fields1)
    #             ax.plot(r3[:irmax3], var3[it, :irmax3], '-', color=colorlist[2],
    #                     linestyle=ls_list[c],  dashes=dash_list[c], linewidth=2,
    #                         label='dx='+str(dx2))
    #
    #         ax.set_title(var_list[i])
    #         ax.set_xlabel('radius r  [m]')
    #         ax.set_ylabel(var_list[i])
    #     axes[1].legend(loc='upper center', bbox_to_anchor=(0, -0.2),
    #                    fancybox=False, ncol=3, fontsize=8)
    #     plt.subplots_adjust(bottom=0.3, right=.98, left=0.07 , top=0.9, wspace=0.25)
    #     fig.suptitle('radially averaged variables   (z='+str(z1)+')')
    #     fig.savefig(os.path.join(path_out_figs, fig_name))
    #     plt.close(fig)
    #
    # file1.close()
    # file2.close()
    # file3.close()

    return


# _______________________________

def plot_rim_width(path1, path2, path3, time1, file_name, fig_name_root, everysecond, kmax):

    file1_rw = nc.Dataset(os.path.join(path1, 'data_analysis', file_name), 'r')
    file2_rw = nc.Dataset(os.path.join(path2, 'data_analysis', file_name), 'r')
    file3_rw = nc.Dataset(os.path.join(path3, 'data_analysis', file_name), 'r')
    dx1 = file1_rw.groups['dimensions'].dimensions['dx'].size
    dx2 = file2_rw.groups['dimensions'].dimensions['dx'].size
    dx3 = file3_rw.groups['dimensions'].dimensions['dx'].size
    zrange1_rw = dx1_nml[2]*file1_rw.groups['dimensions'].variables['krange'][:]
    zrange2_rw = dx2_nml[2]*file2_rw.groups['dimensions'].variables['krange'][:]
    zrange3_rw = dx3_nml[2]*file3_rw.groups['dimensions'].variables['krange'][:]
    wmin1 = file1_rw.groups['rim_width'].variables['wmin'][:,:]
    wmin2 = file2_rw.groups['rim_width'].variables['wmin'][:,:]
    wmin3 = file3_rw.groups['rim_width'].variables['wmin'][:,:]
    wmax1 = file1_rw.groups['rim_width'].variables['wmax'][:,:]
    wmax2 = file2_rw.groups['rim_width'].variables['wmax'][:,:]
    wmax3 = file3_rw.groups['rim_width'].variables['wmax'][:,:]
    rmin1 = file1_rw.groups['rim_width'].variables['r_wmin'][:,:]
    rmin2 = file2_rw.groups['rim_width'].variables['r_wmin'][:,:]
    rmin3 = file3_rw.groups['rim_width'].variables['r_wmin'][:,:]
    rmax1 = file1_rw.groups['rim_width'].variables['r_wmax'][:,:]
    rmax2 = file2_rw.groups['rim_width'].variables['r_wmax'][:,:]
    rmax3 = file3_rw.groups['rim_width'].variables['r_wmax'][:,:]
    rc1 = file1_rw.groups['rim_width'].variables['r_wcenter'][:,:]
    rc2 = file2_rw.groups['rim_width'].variables['r_wcenter'][:,:]
    rc3 = file3_rw.groups['rim_width'].variables['r_wcenter'][:,:]

    print('')
    kmax1_ = np.int(np.double(kmax) / dz13 + 1)
    kmax2_ = np.int(np.double(kmax) / dz12 + 1)
    kmax3_ = kmax + 1
    zrange1_rw = zrange1_rw[:kmax1_]
    zrange2_rw = zrange2_rw[:kmax2_]
    zrange3_rw = zrange3_rw[:kmax3_]
    print('kmax: ', kmax, kmax1_, kmax2_, kmax3_)
    print(len(zrange1_rw), len(zrange2_rw), len(zrange3_rw))
    print zrange1_rw
    print zrange2_rw
    print zrange3_rw

    if everysecond==0:
        delta = 1
    else:
        delta = 2
    for it, t0 in enumerate(time1[1::delta]):
        print('plotting rim geometry', t0, delta)
        cm = plt.cm.get_cmap('coolwarm')
        c1 = cm(1.)
        c2 = cm(.75)
        c3 = cm(.01)
        fig_name = fig_name_root +'_t' + str(np.int(t0)) + '.png'
        # fig_name = 'rim_geometry_t' + str(np.int(t0)) + '.png'
        ncol = 3
        fig, axes = plt.subplots(1, ncol, sharey='none', figsize=(5 * ncol, 5))
        count_color = 2 * np.double(it) / len(time1)
        ax = axes[0]
        # ax.plot(rcenter[1, it, :], krange, '-o', label='r center')
        print(rmin1.shape, zrange1_rw.shape)
        ax.plot(rmin1[it, :kmax1_], zrange1_rw, '-o', label='dx='+str(dx1), color=c1)
        ax.plot(rmin2[it, :kmax2_], zrange2_rw, '-o', label='dx='+str(dx2), color=c2)
        ax.plot(rmin3[it, :kmax3_], zrange3_rw, '-o', label='dx='+str(dx3), color=c3)
        ax.plot(rmax1[it, :kmax1_], zrange1_rw, '-o', label='dx='+str(dx1), color=c1)
        ax.plot(rmax2[it, :kmax2_], zrange2_rw, '-o', label='dx='+str(dx2), color=c2)
        ax.plot(rmax3[it, :kmax3_], zrange3_rw, '-o', label='dx='+str(dx3), color=c3)
        ax.legend(loc='best', fontsize=10)
        ax.set_xlabel('radius r  [m]')
        ax.set_ylabel('height k  (dz=' + str(dx1) + ', ' + str(dx2) + ', ' + str(dx1) + 'm)')
        ax.set_title('r(min(w)), r(max(w))')

        ax = axes[1]
        ax.plot(rmin1[it, :kmax1_]-rc1[it, :kmax1_], zrange1_rw, '-o', label='dx='+str(dx1), color=c1)
        ax.plot(rmin2[it, :kmax2_]-rc2[it, :kmax2_], zrange2_rw, '-o', label='dx='+str(dx2), color=c2)
        ax.plot(rmin3[it, :kmax3_]-rc3[it, :kmax3_], zrange3_rw, '-o', label='dx='+str(dx3), color=c3)
        ax.plot(rmax1[it, :kmax1_]-rc1[it, :kmax1_], zrange1_rw, '-o', label='dx='+str(dx1), color=c1)
        ax.plot(rmax2[it, :kmax2_]-rc2[it, :kmax2_], zrange2_rw, '-o', label='dx='+str(dx2), color=c2)
        ax.plot(rmax3[it, :kmax3_]-rc3[it, :kmax3_], zrange3_rw, '-o', label='dx='+str(dx3), color=c3)
        ax.set_xlim(-1.5e3, 1e3)

        ax = axes[2]
        ax.set_title('min(w), max(w)')
        ax.plot([0,0], [0, zrange1_rw[-1]], 'k', linewidth=0.5)
        ax.plot(wmin1[it, :kmax1_], zrange1_rw, '-o', label='dx='+str(dx1), color=c1)
        ax.plot(wmin2[it, :kmax2_], zrange2_rw, '-o', label='dx='+str(dx2), color=c2)
        ax.plot(wmin3[it, :kmax3_], zrange3_rw, '-o', label='dx='+str(dx3), color=c3)
        ax.plot(wmax1[it, :kmax1_], zrange1_rw, '-o', label='dx='+str(dx1), color=c1)
        ax.plot(wmax2[it, :kmax2_], zrange2_rw, '-o', label='dx='+str(dx2), color=c2)
        ax.plot(wmax3[it, :kmax3_], zrange3_rw, '-o', label='dx='+str(dx3), color=c3)
        # ax.plot([0, 0], [krange[0], krange[-1]], 'k-')
        # ax.legend(fontsize=10)
        ax.legend(loc='lower right', bbox_to_anchor=(1.35, 0.3),
                       fancybox=False, ncol=1, fontsize=10)
        ax.set_xlabel('w [m/s]')
        ax.set_ylabel('height z  (dz=' + str(dx1_nml[2]) + 'm)')

        plt.suptitle('t=' + str(t0) + 's')
        plt.subplots_adjust(bottom=0.1, right=.9, left=0.07 , top=0.9, wspace=0.25)
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)
    file1_rw.close()
    file2_rw.close()
    file3_rw.close()
    return
# _______________________________
def set_input_parameters(args, path_root):
    global case_name
    case_name = 'ColdPoolDry_single_3D'


    print ''' setting parameters '''
    global path_out_data, path_out_figs
    if args.path_out:
        path_out_data = os.path.join(path_root, args.path_out)
    else:
        path_out_data = os.path.join(path_root, 'figs_comp_dx25_dx50_dx100')
    if not os.path.exists(path_out_data):
        os.mkdir(path_out_data)
    path_out_figs = path_out_data
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)


    if args.dTh:
        dTh = args.dTh
    else:
        dTh = 3
    if args.zstar:
        zstar = np.int(args.zstar)
    else:
        zstar = 1000
    if args.rstar:
        rstar = np.int(args.rstar)
    else:
        rstar = 1000

    return dTh, zstar, rstar



def define_geometry(args, path1, path2, path3):
    # global case_name
    # case_name = args.casename
    nml1 = simplejson.loads(open(os.path.join(path1, case_name + '.in')).read())
    nml2 = simplejson.loads(open(os.path.join(path2, case_name + '.in')).read())
    nml3 = simplejson.loads(open(os.path.join(path3, case_name + '.in')).read())
    global nx1, nx2, nx3, dx1_nml, dx2_nml, dx3_nml
    nx1 = np.ndarray(3, dtype=np.int)
    nx2 = np.ndarray(3, dtype=np.int)
    nx3 = np.ndarray(3, dtype=np.int)
    dx1_nml = np.ndarray(3, dtype=np.int)
    dx2_nml = np.ndarray(3, dtype=np.int)
    dx3_nml = np.ndarray(3, dtype=np.int)
    nx1[0] = nml1['grid']['nx']
    nx1[1] = nml1['grid']['ny']
    nx1[2] = nml1['grid']['nz']
    nx2[0] = nml2['grid']['nx']
    nx2[1] = nml2['grid']['ny']
    nx2[2] = nml2['grid']['nz']
    nx3[0] = nml2['grid']['nx']
    nx3[1] = nml2['grid']['ny']
    nx3[2] = nml2['grid']['nz']
    dx1_nml[0] = nml1['grid']['dx']
    dx1_nml[1] = nml1['grid']['dy']
    dx1_nml[2] = nml1['grid']['dz']
    dx2_nml[0] = nml2['grid']['dx']
    dx2_nml[1] = nml2['grid']['dy']
    dx2_nml[2] = nml2['grid']['dz']
    dx3_nml[0] = nml3['grid']['dx']
    dx3_nml[1] = nml3['grid']['dy']
    dx3_nml[2] = nml3['grid']['dz']
    gw = nml1['grid']['gw']

    global dt_fields1, dt_fields2, dt_fields3
    dt_fields1 = nml1['fields_io']['frequency']
    dt_fields2 = nml2['fields_io']['frequency']
    dt_fields3 = nml3['fields_io']['frequency']


    global ic_arr1, jc_arr1, ic_arr2, jc_arr2
    try:
        print('(ic,jc) from nml')
        ic1 = nml1['init']['ic']
        jc1 = nml1['init']['jc']
        ic2 = nml2['init']['ic']
        jc2 = nml2['init']['jc']
        if ic1 != ic2 or jc1 != jc2:
            print('Problem: different geometries')
            sys.exit()
        else:
            ic = ic1
            jc = jc1
            ic_arr = [ic]
            jc_arr = [jc]
    except:
        print('(ic,jc) NOT from nml')
        if case_name == 'ColdPoolDry_single_3D':
            ic = np.int(nx1/2)
            jc = np.int(ny1/2)
            ic_arr1 = [ic]
            jc_arr1 = [jc]
            ic = np.int(nx2/2)
            jc = np.int(ny2/2)
            ic_arr2 = [ic]
            jc_arr2 = [jc]
        else:
            print('ic, jc not defined')
    global kmax
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = nx1[2]

    global tmin, tmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100

    print('')
    print('kmax ', kmax, 'nz1 ', nx1[2])
    print('')

    return

# _______________________________

# def define_geometry(nml):
#     a = nml['grid']['nx']
#     '''--- define geometry ---'''
#     global rstar
#     if case_name == 'ColdPoolDry_double_2D':
#         rstar = 5000.0  # half of the width of initial cold-pools [m]
#         irstar = np.int(np.round(rstar / dx))
#         # zstar = nml['init']['h']
#         isep = 4 * irstar
#         ic1 = np.int(nx / 3)
#         ic2 = ic1 + isep
#         jc1 = np.int(ny / 2)
#         jc2 = jc1
#         ic_arr = [ic1, ic2]
#         jc_arr = [jc1, jc2]
#     elif case_name == 'ColdPoolDry_single_3D':
#         rstar = nml['init']['r']
#         # irstar = np.int(np.round(rstar / dx))
#         # zstar = nml['init']['h']
#         dTh = nml['init']['dTh']
#         ic = np.int(nx / 2)
#         jc = np.int(ny / 2)
#         # xc = Gr.x_half[ic + Gr.dims.gw]  # center of cold-pool
#         # yc = Gr.y_half[jc + Gr.dims.gw]  # center of cold-pool
#         ic_arr = [ic]
#         jc_arr = [jc]
#     elif case_name == 'ColdPoolDry_double_3D':
#         try:
#             rstar = nml['init']['r']
#         except:
#             rstar = 5000.0  # half of the width of initial cold-pools [m]
#         irstar = np.int(np.round(rstar / dx))
#         # zstar = nml['init']['h']
#         isep = 4 * irstar
#         jsep = 0
#         ic1 = np.int(np.round((nx + 2 * gw) / 3)) - gw
#         jc1 = np.int(np.round((ny + 2 * gw) / 2)) - gw
#         ic2 = ic1 + isep
#         jc2 = jc1 + jsep
#         ic_arr = [ic1, ic2]
#         jc_arr = [jc1, jc2]
#     elif case_name == 'ColdPoolDry_triple_3D':
#         try:
#             rstar = nml['init']['r']
#         except:
#             rstar = 5000.0  # half of the width of initial cold-pools [m]
#         irstar = np.int(np.round(rstar / dx))
#         d = np.int(np.round(ny / 2))
#         dhalf = np.int(np.round(ny / 4))
#         a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
#         ic1 = np.int(np.round(a / 2))  # + gw
#         ic2 = ic1
#         ic3 = ic1 + np.int(np.round(a))
#         jc1 = np.int(np.round(d / 2))  # + gw
#         jc2 = jc1 + d
#         jc3 = jc1 + np.int(np.round(d / 2))
#         ic_arr = [ic1, ic2, ic3]
#         jc_arr = [jc1, jc2, jc3]
#
#         isep = dhalf
#
#     return ic_arr, jc_arr

# _______________________________

if __name__ == '__main__':
    main()