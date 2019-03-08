import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import netCDF4 as nc
import argparse
import json as simplejson
import os
import sys

def main():
    ''' set paths & parameters '''
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    # parser.add_argument("--t0")
    parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    # times, nml = set_input_parameters(args)
    nml = set_input_parameters(args)



    # (1) read in averaged w
    #       path_data = os.path.join(path, 'data_analysis')
    #       file_name = 'stats_radial_averaged.nc'
    # (2) start at rmax, find irmax
    # - outer radius r0: find first point where w > w0
    #     >> compare to tracer radius
    #     ??? test choice for w0; maybe depending on k
    # - center radius r1: find first zero-crossing for r<r0
    # - innter radius r1: find second zero-crossing for r<r1 (?? no sure, if given)
    # (3) test
    # - rcenter indpendent of height?
    # - linear increase of w from rmin throught rcetner to rmax >> solid body
    # - function of decay of w from rmax to r..

    # read in file
    # file_name = 'stats_radial_averaged.nc'
    print('!!!!!!!!!!!!!!!!! filename')
    file_name = 'stats_radial_averaged_.nc'
    data_stats = nc.Dataset(os.path.join(path_data, file_name), 'r')
    nz_stats = data_stats.groups['dimensions'].dimensions['nz']
    krange_ = data_stats.groups['dimensions'].variables['krange'][:]
    stats_grp = data_stats.groups['stats'].variables
    times_ = data_stats.groups['timeseries'].variables['time'][:]
    nt = len(times_)
    r_range = stats_grp['r'][:]
    ''' ----- w ----- '''
    var_name = 'w'
    var = stats_grp[var_name][:, :, :]  # var(t, r, k)
    data_stats.close()


    # if args.t0:
    #     t0 = np.int(args.t0)
    # else:
    #     t0 = 800
    # it = np.where(times_ == t0)[0][0]

    # k0 = np.int(args.k0)
    ''' krange '''
    kmin = np.int(args.kmin)
    kmax = np.int(args.kmax)
    # krange = np.arange(kmin, kmax+1)
    krange = [np.int(k) for k in krange_ if k>=kmin and k<=kmax]
    print('krange: ', krange_)
    print('krange: ', krange)
    print ''

    ''' time range '''
    times = [np.int(t) for t in times_ if t>=tmin and t<=tmax]
    itmin = np.where(times_ == tmin)
    itmax = np.where(times_ == tmax)
    print times
    print times_
    print('tmin etc', tmin, tmax, itmin, itmax)

    # declare arrays
    wmin = np.zeros((nt, kmax - kmin+1))
    wmax = np.zeros((nt, kmax - kmin+1))
    rcenter = np.zeros((2, nt, kmax - kmin + 1))
    rmin = np.zeros((nt, kmax - kmin + 1))
    rmax = np.zeros((nt, kmax - kmin + 1))



    rmax_plot = 6e3
    irmax = np.where(r_range == rmax_plot)[0][0]

    for it, t0 in enumerate(times):
        print('---- t0='+str(t0)+' ----')


        omega_plus = np.zeros((kmax-kmin+1))
        omega_minus = np.zeros((kmax-kmin+1))

        for k0 in range(kmin, kmax+1):
            # indices of max / min
            wmin[it, k0] = np.amin(var[it, 3:, k0])
            imin = np.argmin(var[it, 3:, k0])+3
            rmin[it, k0] = r_range[imin]
            wmax[it, k0] = np.amax(var[it, :, k0])
            imax = np.argmax(var[it, :, k0])
            rmax[it, k0] = r_range[imax]
            print ''

            # find r1 (zero point)
            i = imax
            while (i>imin and var[it, i-1, k0]>0):
                i -= 1
            if np.abs(var[it, i, k0]) < np.abs(var[it, i-1, k0]):
                icenter = i
            else:
                icenter = i-1
            wcenter = var[it, icenter, k0]

            rcenter[0, it, k0] = icenter
            rcenter[1, it, k0] = r_range[icenter]
            print('wmin', wmin[it, k0])
            print('wcenter', wcenter)
            print('wmax', wmax[it, k0])
            # print 'imax', imax
            # print 'icenter', icenter
            # print 'imin', imin
            print ''




            # compute linear fitting functions
            # solid body rotation: u = r*omega
            omega_plus[k0] = (wmax[it, k0] - wcenter) / (rmax[it, k0] - rcenter[1, it, k0])
            omega_minus[k0] = (wcenter - wmin[it, k0]) / (rcenter[1, it, k0] - rmin[it, k0])
            dyplus = -omega_plus[k0]*rcenter[1,it, k0]
            dyminus = -omega_minus[k0]*rcenter[1,it, k0]
            lin_plus = omega_plus[k0]*r_range + dyplus
            lin_minus = omega_minus[k0]*r_range + dyminus
            print('mplus', omega_plus[k0])
            print('mminus', omega_minus[k0])
            print('dyplus', dyplus)
            print('dyminus', dyminus)
            print ''


            # plotting test_fig
            # fig_name = 'test_fig_t'+str(t0)+'_z'+str(k0*dx[2])+'m.png'
            # ncol = 3
            # fig, axes = plt.subplots(1, ncol, sharey='none', figsize=(5 * ncol, 5))
            # count_color = 2 * np.double(it) / len(times)
            # ax = axes[0]
            # ax.plot([0, r_range[irmax]], [0., 0.], 'k')
            # ax.plot(rmax[k0], wmax[it, k0], 'kx')
            # ax.plot(rmin[it, k0], wmin[it, k0], 'kx')
            # ax.plot([rmin[it, k0], rmin[it, k0]], [wmin[it, k0], wmax[it, k0]], '0.5')
            # ax.plot([rmax[it, k0], rmax[it, k0]], [wmin[it, k0], wmax[it, k0]], '0.5')
            # ax.plot([rcenter[1,it, k0], rcenter[1,it, k0]], [wmin[it, k0], wmax[it, k0]], 'k')
            # ax.plot(r_range[:irmax], var[it, :irmax, k0], color=cm.jet(count_color), label='t=' + str(np.int(t0)))
            # ax = axes[1]
            # ax.plot([0, r_range[irmax]], [0., 0.], 'k')
            # ax.plot([0, r_range[irmax]], [wmax[it, k0], wmax[it, k0]], '0.5')
            # ax.plot([0, r_range[irmax]], [wmin[it, k0], wmin[it, k0]], '0.5')
            # ax.plot([rmin[it, k0], rmin[it, k0]], [wmin[it, k0], wmax[it, k0]], '0.5')
            # ax.plot([rmax[it, k0], rmax[it, k0]], [wmin[it, k0], wmax[it, k0]], '0.5')
            # ax.plot([rcenter[1,it, k0], rcenter[1,it, k0]], [wmin[it, k0], wmax[it, k0]], 'k')
            # ax.plot(r_range[:irmax], var[it, :irmax, k0], color=cm.jet(count_color), label='t=' + str(np.int(t0)))
            # ax.plot(r_range[:irmax], lin_plus[:irmax], label='om+='+str(np.round(omega_plus[k0],3)))
            # ax.plot(r_range[:irmax], lin_minus[:irmax], label='om-='+str(np.round(omega_minus[k0],3)))
            # ax.legend(fontsize=10)
            # ax.set_xlim(r_range[imin-10], r_range[imax+10])
            # ax.set_ylim([wmin[it, k0]-1e-1, wmax[it, k0]+1e-1])
            # ax.legend(loc=4, fontsize=10)
            # ax = axes[2]
            # min = imin -30
            # max = imax +30
            # ax.plot([r_range[min], r_range[max]], [0., 0.], 'k')
            # ax.plot(rmax[it, k0], wmax[it, k0], 'kx')
            # ax.plot(rmin[it, k0], wmin[it, k0], 'kx')
            # ax.plot([rmin[it, k0], rmin[it, k0]], [wmin[it, k0], wmax[k0]], '0.5')
            # ax.plot([rmax[it, k0], rmax[it, k0]], [wmin[it, k0], wmax[k0]], '0.5')
            # ax.plot([rcenter[1,it, k0], rcenter[1,it, k0]], [wmin[it, k0], wmax[k0]], 'k')
            # ax.plot(r_range, var[it, :, k0], 'o-',
            #         color=cm.jet(count_color), label='t=' + str(np.int(t0)))
            # ax.set_xlim(r_range[imin - 10], r_range[imax + 10])
            # plt.suptitle('z='+str(k0*dx[2])+'m')
            # fig.savefig(os.path.join(path_out_figs, fig_name))
            # plt.close(fig)




        print 'plotting rim geometry'
        # plotting rcenter
        fig_name = 'rim_geometry_t' + str(np.int(t0)) + '.png'
        ncol = 3
        fig, axes = plt.subplots(1, ncol, sharey='all', figsize=(5 * ncol, 5))
        count_color = 2 * np.double(it) / len(times)
        ax = axes[0]
        ax.plot(rcenter[1,it,:], krange, '-o', label='r center')
        ax.plot(rmin[it,:], krange, '-o', label='r min')
        ax.plot(rmax[it,:], krange, '-o', label='r max')
        ax.legend(fontsize=10)
        ax.set_xlabel('radius r  [m]')
        ax.set_ylabel('height k  (dz='+str(dx[2])+'m)')
        ax = axes[1]
        ax.plot(wmin[it, :], krange, '-o', label='w min')
        ax.plot(wmax[it, :], krange, '-o', label='w max')
        ax.plot([0,0], [krange[0], krange[-1]], 'k-')
        ax.legend(fontsize=10)
        ax.set_xlabel('w [m/s]')
        ax.set_ylabel('height k  (dz=' + str(dx[2]) + 'm)')
        ax = axes[2]
        ax.plot(omega_minus, krange, '-o', label='om+')
        ax.plot(omega_plus, krange, '-o', label='om-')
        # ax.plot([0,0], [krange[0], krange[-1]], 'k-')
        ax.legend(fontsize=10)
        ax.set_xlabel('vorticity = dw/dr')
        ax.set_ylabel('height k  (dz=' + str(dx[2]) + 'm)')
        plt.suptitle('t=' + str(t0) + 's')
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)


    # output wmin, wmax into data_analysis/stats_radial_averaged.nc
    print('!!! name output file')
    file_name = 'stats_radial_averaged_.nc'
    dump_minmax_profiles('wmax', wmax, kmin, kmax, tmin, tmax, file_name)
    dump_minmax_profiles('wmin', wmin, kmin, kmax, tmin, tmax, file_name)
    dump_minmax_profiles('r_wmin', rmin, kmin, kmax, tmin, tmax, file_name)
    dump_minmax_profiles('r_wmax', rmax, kmin, kmax, tmin, tmax, file_name)
    dump_minmax_profiles('r_wcenter', rcenter[1,:,:], kmin, kmax, tmin, tmax, file_name)



    return
# _______________________________
# _______________________________

def dump_minmax_profiles(var_name, variable, kmin, kmax, tmin, tmax, file_name):
    print('do output', var_name, tmin, tmax, kmin, kmax)
    rootgrp = nc.Dataset(os.path.join(path_data, file_name), 'r+')
    # rootgrp = nc.Dataset(os.path.join(path_data, file_name), 'a', format='NETCDF4')
    # root_grp = nc.Dataset(os.path.join(path_data, file_name), 'r+', format='NETCDF4')
    nz = rootgrp.groups['dimensions'].dimensions['nz'].size
    nt = rootgrp.groups['timeseries'].dimensions['nt'].size
    try:
        prof_grp = rootgrp.groups['profiles']
    except:
        prof_grp = rootgrp.createGroup('profiles')
        prof_grp.createDimension('nt', nt)
        prof_grp.createDimension('nz', nz)
    try:
        var = prof_grp.variables[var_name]
    except:
        var = prof_grp.createVariable(var_name, 'f8', ('nt', 'nz'))
    var[:, kmin:kmax+1] = variable

    # try:
    #     var = prof_grp.variables['wmax']
    # except:
    #     var = prof_grp.createVariable('wmax', 'f8', ('nt', 'nz'))
    # # print var.shape, wmax.shape
    # var[:, kmin:kmax+1] = wmax

    rootgrp.close()
    return

# _______________________________
# _______________________________
def set_input_parameters(args):
    print ''' setting parameters '''
    global path, path_data, path_out_figs
    path = args.path

    path_data = os.path.join(path, 'data_analysis')
    path_out_figs = os.path.join(path, 'figs_radial_average')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print ''
    print 'paths:'
    print path_data
    print path_out_figs
    print ''

    global case_name
    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx
    dx = np.ndarray(3, dtype=np.int)
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']

    global kmax
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = nz

    global tmin, tmax
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = 100

    # ''' time range '''
    # times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
    #          and tmin <= np.int(name[:-3]) <= tmax]
    # times.sort()

    # print('')
    # print('times', times)
    print('')
    print('tmin, tmax', tmin, tmax)
    print('kmax ', kmax, 'nx ', nx)
    print('')

    # return times, nml
    return nml

# _______________________________

def define_geometry(nml):
    # global ic_arr, jc_arr

    '''--- define geometry ---'''
    global rstar
    if case_name == 'ColdPoolDry_double_2D':
        rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        # zstar = nml['init']['h']
        isep = 4 * irstar
        ic1 = np.int(nx / 3)
        ic2 = ic1 + isep
        jc1 = np.int(ny / 2)
        jc2 = jc1
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_single_3D':
        rstar = nml['init']['r']
        # irstar = np.int(np.round(rstar / dx))
        # zstar = nml['init']['h']
        # dTh = nml['init']['dTh']
        try:
            print('(ic,jc) from nml')
            ic = nml['init']['ic']
            jc = nml['init']['jc']
        except:
            print('(ic,jc) NOT from nml')
            ic = np.int(nx / 2)
            jc = np.int(ny / 2)
        ic_arr = [ic]
        jc_arr = [jc]
        # xc = Gr.x_half[ic + Gr.dims.gw]  # center of cold-pool
        # yc = Gr.y_half[jc + Gr.dims.gw]  # center of cold-pool
    elif case_name == 'ColdPoolDry_double_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        # zstar = nml['init']['h']
        isep = 4 * irstar
        jsep = 0
        ic1 = np.int(np.round((nx + 2 * gw) / 3)) - gw
        jc1 = np.int(np.round((ny + 2 * gw) / 2)) - gw
        ic2 = ic1 + isep
        jc2 = jc1 + jsep
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_triple_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
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

        isep = dhalf

    # ''' --- auxiliary arrays (since no Grid.pyx) ---'''
    # global nx_, ny_, nz_
    # # test file:
    # var = read_in_netcdf_fields('u', os.path.join(path_fields, files[0]))
    # [nx_, ny_, nz_] = var.shape
    #
    # x_half = np.empty((nx_), dtype=np.double, order='c')
    # y_half = np.empty((ny_), dtype=np.double, order='c')
    # z_half = np.empty((nz_), dtype=np.double, order='c')
    # count = 0
    # for i in xrange(nx_):
    #     x_half[count] = (i + 0.5) * dx
    #     count += 1
    # count = 0
    # for j in xrange(ny_):
    #     y_half[count] = (j + 0.5) * dy
    #     count += 1
    # count = 0
    # for i in xrange(nz_):
    #     z_half[count] = (i + 0.5) * dz
    #     count += 1
    #
    # return x_half, y_half, z_half

    return ic_arr, jc_arr

# _______________________________


if __name__ == '__main__':
    main()