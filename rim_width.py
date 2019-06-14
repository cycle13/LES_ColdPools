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
    path_data = os.path.join(path, 'data_analysis')
    path_out_figs = os.path.join(path, 'figs_radial_average')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print 'paths:'
    print path_data
    print path_out_figs
    print ''


    # IDEA:
    # (1) read in averaged w
    #       path_data = os.path.join(path, 'data_analysis')
    #       file_name = 'stats_radial_averaged.nc'
    # (2) start at rmax, find irmax_plot
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
    file_name_in = 'stats_radial_averaged.nc'
    data_stats = nc.Dataset(os.path.join(path_data, file_name_in), 'r')
    nz_stats = data_stats.groups['dimensions'].dimensions['nz']
    krange_stats = data_stats.groups['dimensions'].variables['krange'][:]
    stats_grp = data_stats.groups['stats'].variables
    times_stats = data_stats.groups['timeseries'].variables['time'][:]
    r_range = stats_grp['r'][:]
    nr = len(r_range)
    ''' ----- w ----- '''
    var_name = 'w'
    var = stats_grp[var_name][:, :, :]  # var(t, r, k)
    data_stats.close()

    # krange = np.arange(kmin, kmax+1)
    krange = [np.int(k) for k in krange_stats if k>=kmin and k<=kmax]
    print('krange radial stats: ', krange_stats)
    print('krange:              ', krange)
    print ''

    ''' time range '''
    times = [np.int(t) for t in times_stats if t>=tmin and t<=tmax]
    itmin = np.where(times_stats == tmin)
    itmax = np.where(times_stats == tmax)
    print 'times radial stats: ', times_stats
    print 'times:              ', times
    print('tmin, tmax: ', tmin, tmax)
    print ''

    # declare arrays
    nt = len(times)
    nk = kmax+1 - kmin
    wmin = np.zeros((nt, nk))            # minimum(w)
    wmax = np.zeros((nt, nk))            # maximum(w)
    wint = np.zeros((nt, nk))            # w at inner rim edge (defined as inner edge of downdraft zone)
    wout = np.zeros((nt, nk))            # w at outer rim edge (defined as outer edge of updraft zone)
    wcrit = 1e-2
    rcenter = np.zeros((2, nt, nk))    # 'center' of vortex: radius wher w approx zero
    rmin = np.zeros((nt, nk))          # radius of wmin
    rmax = np.zeros((nt, nk))          # radius of wmax
    rint = np.zeros((nt, nk))            # radius of inner rim edge, defined as maximum(w(r<rmin))
    rout = np.zeros((nt, nk))            # radius at outer rim edge, defined as point of w(r>rmax)<wcrit

    omega_plus = np.zeros((nt, nk))    # vorticity from linear fit on (wmax-wcenter)/(rmax-rcenter)
    omega_minus = np.zeros((nt, nk))   # vorticity from linear fit on (wcenter-wmin)/(rcenter-wmin)

    rmax_plot = 7e3
    irmax_plot = np.where(r_range == rmax_plot)[0][0]

    for it, t0 in enumerate(times):
        print('---- t0='+str(t0)+', '+str(it)+' ----')
        for k0 in range(kmin, kmax+1):
            # print('   ---- k0='+str(k0)+' ----')
            # indices of max / min
            delta = np.int(500./dx[0])
            wmin[it, k0] = np.amin(var[it, delta:, k0])     # search min(w) in interval r=[500,rmax]
            imin = np.argmin(var[it, delta:, k0]) + delta
            rmin[it, k0] = r_range[imin]
            wmax[it, k0] = np.amax(var[it, delta:, k0])
            imax = np.argmax(var[it, delta:, k0]) + delta
            rmax[it, k0] = r_range[imax]

            # find center of vortex (zero point of w)
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

            # find inner edge of rim: defined as maximum(w(r<min))
            delta_int = 3
            wint[it, k0] = np.amax(var[it, delta_int:imin, k0])
            iint = np.argmax(var[it, delta_int:imin, k0]) + delta_int
            rint[it, k0] = r_range[iint]
            # find outer edge of rim: defined as point of w(r>rmax)<wcrit
            i = imax
            while (i < nr and var[it, i - 1, k0] > wcrit):
                i += 1
            iout = i
            wout[it, k0] = var[it, iout, k0]
            rout[it, k0] = r_range[iout]

            # print ''
            # print('wmin', wmin[it, k0])
            # print('wcenter', wcenter)
            # print('wmax', wmax[it, k0])
            # # print 'imax', imax
            # # print 'icenter', icenter
            # # print 'imin', imin
            # print ''


            # compute linear fitting functions
            # solid body rotation: u = r*omega
            omega_plus[it, k0] = (wmax[it, k0] - wcenter) / (rmax[it, k0] - rcenter[1, it, k0])
            omega_minus[it, k0] = (wcenter - wmin[it, k0]) / (rcenter[1, it, k0] - rmin[it, k0])
            dyplus = -omega_plus[it, k0]*rcenter[1,it, k0]
            dyminus = -omega_minus[it, k0]*rcenter[1,it, k0]
            lin_plus = omega_plus[it, k0]*r_range + dyplus
            lin_minus = omega_minus[it, k0]*r_range + dyminus
            # print('mplus', omega_plus[it, k0])
            # print('mminus', omega_minus[it, k0])
            # print('dyplus', dyplus)
            # print('dyminus', dyminus)
            # print ''


            # # # plotting test_fig: w(t) for each k and show rmin, rmax, rcenter
            # fig_name = 'test_fig_t' + str(t0) + '_z' + str(k0 * dx[2]) + 'm.png'
            # plot_test_fig(var, wmin, wmax, rmin, rcenter, rmax, imin, imax,
            #               wint, rint, wout, rout, wcrit,
            #             r_range, nr, irmax_plot,
            #             lin_plus, lin_minus, omega_plus, omega_minus,
            #             it, t0, times, k0, fig_name, path_out_figs)



    ''' plotting '''
    fig_name = 'w_rmin_rmax_test.png'
    plot_w_rmin_rmax_rc(var, wmin, wmax, wcenter, r_range, rmin, rmax, rcenter, rint, rout,
                        delta, times, fig_name, path_out_figs)

    fig_name = 'rim_width' + '.png'
    plot_rim_width(rmin, rmax, rcenter, rint, rout, krange, times, times_stats, fig_name, path_out_figs)


    ''' make output '''
    file_name_out = 'stats_radial_averaged_rimwidth.nc'
    create_output_file(times, krange, nk, file_name_out, path_data)

    rootgrp_out = nc.Dataset(os.path.join(path_data, file_name_out), 'r+')
    nz_out = rootgrp_out.groups['dimensions'].dimensions['nz'].size
    nt_out = rootgrp_out.groups['timeseries'].dimensions['nt'].size
    rootgrp_out.close()

    rootgrp_in = nc.Dataset(os.path.join(path_data, file_name_in), 'r+')
    nz_in = rootgrp_in.groups['dimensions'].dimensions['nz'].size
    nt_in = rootgrp_in.groups['timeseries'].dimensions['nt'].size
    rootgrp_in.close()
    print ''

    dump_minmax_profiles('wmax', wmax, kmin, kmax, tmin, tmax, file_name_out, path_data)
    dump_minmax_profiles('wmin', wmin, kmin, kmax, tmin, tmax, file_name_out, path_data)
    dump_minmax_profiles('r_wmin', rmin, kmin, kmax, tmin, tmax, file_name_out, path_data)
    dump_minmax_profiles('r_wmax', rmax, kmin, kmax, tmin, tmax, file_name_out, path_data)
    dump_minmax_profiles('r_wcenter', rcenter[1,:,:], kmin, kmax, tmin, tmax, file_name_out, path_data)
    # dump_minmax_profiles('r_int', rint, kmin, kmax, tmin, tmax, file_name_out, path_data)
    # dump_minmax_profiles('r_out', rout, kmin, kmax, tmin, tmax, file_name_out, path_data)
    dump_minmax_profiles('omega_plus', omega_plus, kmin, kmax, tmin, tmax, file_name_out, path_data)
    dump_minmax_profiles('omega_minus', omega_minus, kmin, kmax, tmin, tmax, file_name_out, path_data)

    return
# _______________________________
# _______________________________
def plot_w_rmin_rmax_rc(w, wmin, wmax, wcenter,
                        r_range, rmin, rmax, rcenter, rint, rout,
                        delta, times, fig_name, path_out_figs):
    ncol = 1
    nrow = 2
    k0 = 0
    irmax = np.int((np.amax(rmax[:,k0])+2000.)/dx[0])
    print(irmax, np.amax(rmax[:,k0]))
    fig, axis = plt.subplots(nrow, ncol, sharey='none', figsize=(10*ncol, 5 * nrow))
    ax1 = axis[0]
    ax1.plot([r_range[delta], r_range[delta]], [-2, 2], 'k--', linewidth=1, label='delta')
    for it, t0 in enumerate(times[1::2]):
        count_color = 2 * np.double(it) / len(times)
        ax1.plot(r_range[:], w[2*it+1, :, k0], color=cm.jet(count_color), label='t='+str(t0))
        ax1.plot(rmax[2*it+1,k0], wmax[2*it+1,k0], 'ko')
        ax1.plot(rmin[2*it+1,k0], wmin[2*it+1,k0], 'kd')
        ax1.set_xlim(0, np.amax(rmax[:,k0])+1000)
    ax1.legend(loc='best', fontsize=8, ncol=4)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    return

def plot_rim_width(rmin, rmax, rcenter, rint, rout, krange, times, times_stats, fig_name, path_out_figs):
    ncol = 2
    nrow = 2
    fig, axes = plt.subplots(nrow, ncol, sharey='none', figsize=(5 * ncol, 5 * nrow))
    ax = axes[0, 0]
    for it, t0 in enumerate(times[2:]):
        count_color = np.double(it) / len(times)
        ax.plot(rmin[it, :] - rcenter[1, it, :], krange, '-', label='t=' + str(t0),
                color=plt.cm.get_cmap('coolwarm')(count_color))
        ax.plot(rmax[it, :] - rcenter[1, it, :], krange, '-', color=plt.cm.get_cmap('coolwarm')(count_color))
    ax.set_xlim(-3e3, 1.2e3)
    ax.set_title('rmin-rcenter; rmax-rcenter')
    ax.set_xlabel('r-rcenter')
    ax.set_ylabel('height k  (dz=' + str(dx[2]) + 'm)')
    ax.legend(loc='upper left', bbox_to_anchor=(0, 1),
              fancybox=True, shadow=True, ncol=2, fontsize=6)
    ax = axes[0, 1]
    ax.set_title('rmax-rmin')
    ax.fill_between(times, 200, 500, color='0.8')
    ax.plot([400, 400], [0, 3e3], 'k', linewidth=1)
    for k0 in krange:
        count_color = np.double(k0) / len(krange)
        ax.plot(times, rmax[:, k0] - rmin[:, k0], '-', label='z=' + str(k0 * dx[2]) + 'm',
                color=plt.cm.get_cmap('coolwarm')(count_color), linewidth=2)
    ax.set_xlim(0, 3500)
    ax.set_ylim(0, 3e3)
    ax.legend(loc='upper left', bbox_to_anchor=(0, 1),
              fancybox=True, shadow=True, ncol=3, fontsize=8)
    ax.set_ylabel('radius r  [m]')
    ax.set_xlabel('times')
    ax = axes[1, 0]
    for it, t0 in enumerate(times[2:]):
        count_color = np.double(it) / len(times)
        ax.plot(rout[it, :] - rcenter[1, it, :], krange, '-', label='t=' + str(t0),
                color=plt.cm.get_cmap('bone')(count_color))
        ax.plot(rint[it, :] - rcenter[1, it, :], krange, '-', color=plt.cm.get_cmap('coolwarm')(count_color))
    # ax.set_xlim(-3e3,1.2e3)
    ax.set_title('rint-rcenter; rout-rcenter')
    ax.set_xlabel('r-rcenter')
    ax.set_ylabel('height k  (dz=' + str(dx[2]) + 'm)')
    ax.legend(loc='upper left', bbox_to_anchor=(-0.3, 1),
              fancybox=True, shadow=True, ncol=3, fontsize=6)
    ax.set_xlabel('vorticity = dw/dr')
    ax = axes[1, 1]
    ax.set_title('rout-rcenter')
    ax.fill_between(times, 1000, 2000, color='0.8')
    ax.plot([400, 400], [0, 3e3], 'k', linewidth=1)
    for k0 in krange:
        count_color = np.double(k0) / len(krange)
        ax.plot(times, rout[:, k0] - rcenter[1, :, k0], '-', label='z=' + str(k0 * dx[2]) + 'm',
                color=plt.cm.get_cmap('coolwarm')(count_color), linewidth=2)
    ax.set_xlim(0, 3500)
    ax.set_ylim(0, 3e3)
    # ax.legend(loc='upper left', bbox_to_anchor=(0, 1),
    #           fancybox=True, shadow=True, ncol=3, fontsize=8)
    ax.set_ylabel('radius r  [m]')
    ax.set_xlabel('times')
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    return
# _______________________________
# _______________________________
def create_output_file(timerange, krange, nk, filename_out, path_out):
    print ''
    print('-------- create statistics file -------- ')
    print(path_out + ', ' + filename_out)
    print('')

    nt = len(timerange)

    rootgrp = nc.Dataset(os.path.join(path_out, filename_out), 'w', format='NETCDF4')

    dims_grp = rootgrp.createGroup('dimensions')
    dims_grp.createDimension('dx', dx[0])
    dims_grp.createDimension('dy', dx[1])
    dims_grp.createDimension('dz', dx[2])
    dims_grp.createDimension('nz', nk)
    var = dims_grp.createVariable('krange', 'f8', ('nz'))
    # var[:] = np.arange(kmin, kmax+1, dtype=np.int)
    var[:] = krange

    ts_grp = rootgrp.createGroup('timeseries')
    ts_grp.createDimension('nt', nt)
    var = ts_grp.createVariable('time', 'f8', ('nt'))
    var.unit = "s"
    var[:] = timerange

    prof_grp = rootgrp.createGroup('rim_width')
    prof_grp.createDimension('nt', nt)
    prof_grp.createDimension('nz', nk)
    # prof_grp.createVariable(var_name, 'f8', ('nt', 'nz'))

    rootgrp.close()
    return


def dump_minmax_profiles(var_name, variable, kmin, kmax, tmin, tmax, file_name, path_data):
    # print('do output', var_name, tmin, tmax, kmin, kmax)
    rootgrp = nc.Dataset(os.path.join(path_data, file_name), 'r+')
    # rootgrp = nc.Dataset(os.path.join(path_data, file_name), 'a', format='NETCDF4')
    # root_grp = nc.Dataset(os.path.join(path_data, file_name), 'r+', format='NETCDF4')
    nz = rootgrp.groups['dimensions'].dimensions['nz'].size
    nt = rootgrp.groups['timeseries'].dimensions['nt'].size

    prof_grp = rootgrp.groups['rim_width']
    try:
        var = prof_grp.variables[var_name]
    except:
        var = prof_grp.createVariable(var_name, 'f8', ('nt', 'nz'))

    var[:, kmin:kmax+1] = variable

    # try:
    #     var = prof_grp.variables['wmax']
    # except:
    #     var = prof_grp.createVariable('wmax', 'f8', ('nt', 'nz'))
    # var[:, kmin:kmax+1] = wmax

    rootgrp.close()
    return

# _______________________________
# _______________________________

def plot_test_fig(var, wmin, wmax, rmin, rcenter, rmax, imin, imax,
                  wint, rint, wout, rout, wcrit,
                  r_range, nr, irmax,
                  lin_plus, lin_minus, omega_plus, omega_minus,
                  it, t0, times, k0, fig_name, path_out_figs):

    ncol = 3
    fig, axes = plt.subplots(1, ncol, sharey='none', figsize=(5 * ncol, 5))
    count_color = 2 * np.double(it) / len(times)
    ax = axes[0]
    ax.plot([0, r_range[irmax]], [0., 0.], 'k')
    ax.plot(rmin[it, k0], wmin[it, k0], 'kx')
    ax.plot(rmax[it, k0], wmax[it, k0], 'kx')
    ax.plot(rint[it, k0], wint[it, k0], 'rx')
    ax.plot(rout[it, k0], wout[it, k0], 'gx')
    ax.plot([rmin[it, k0], rmin[it, k0]], [wmin[it, k0], wmax[it, k0]], '0.5')
    ax.plot([rmax[it, k0], rmax[it, k0]], [wmin[it, k0], wmax[it, k0]], '0.5')
    ax.plot([rint[it, k0], rint[it, k0]], [wmin[it, k0], wmax[it, k0]], 'r')
    ax.plot([rout[it, k0], rout[it, k0]], [wmin[it, k0], wmax[it, k0]], 'g')
    ax.plot([rcenter[1, it, k0], rcenter[1, it, k0]], [wmin[it, k0], wmax[it, k0]], 'k')
    ax.plot(r_range[:irmax], var[it, :irmax, k0], color=cm.jet(count_color), label='t=' + str(np.int(t0)))
    ax = axes[1]
    ax.plot([0, r_range[irmax]], [0., 0.], 'k')
    ax.plot([0, r_range[irmax]], [wmax[it, k0], wmax[it, k0]], '0.5')
    ax.plot([0, r_range[irmax]], [wmin[it, k0], wmin[it, k0]], '0.5')
    ax.plot([rmin[it, k0], rmin[it, k0]], [wmin[it, k0], wmax[it, k0]], '0.5')
    ax.plot([rmax[it, k0], rmax[it, k0]], [wmin[it, k0], wmax[it, k0]], '0.5')
    ax.plot([rcenter[1, it, k0], rcenter[1, it, k0]], [wmin[it, k0], wmax[it, k0]], 'k')
    ax.plot(r_range[:irmax], var[it, :irmax, k0], color=cm.jet(count_color), label='t=' + str(np.int(t0)))
    ax.plot(r_range[:irmax], lin_plus[:irmax], label='om+=' + str(np.round(omega_plus[k0], 3)))
    ax.plot(r_range[:irmax], lin_minus[:irmax], label='om-=' + str(np.round(omega_minus[k0], 3)))
    ax.legend(fontsize=10)
    ax.set_xlim(r_range[np.maximum(imin - 10, 0)], r_range[np.minimum(imax + 10, nr - 1)])
    ax.set_ylim([wmin[it, k0] - 1e-1, wmax[it, k0] + 1e-1])
    ax.legend(loc=4, fontsize=10)
    ax = axes[2]
    min = np.maximum(imin - 30, 0)
    max = np.minimum(imax + 30, nr - 1)
    ax.plot([r_range[min], r_range[irmax]], [wcrit, wcrit], 'r')
    ax.plot([r_range[min], r_range[irmax]], [0., 0.], 'k')
    ax.plot(rmax[it, k0], wmax[it, k0], 'kx')
    ax.plot(rmin[it, k0], wmin[it, k0], 'kx')
    ax.plot(rout[it, k0], wout[it, k0], 'rx')
    ax.plot([rmin[it, k0], rmin[it, k0]], [wmin[it, k0], wmax[it, k0]], '0.5')
    ax.plot([rmax[it, k0], rmax[it, k0]], [wmin[it, k0], wmax[it, k0]], '0.5')
    ax.plot([rcenter[1, it, k0], rcenter[1, it, k0]], [wmin[it, k0], wmax[it, k0]], 'k')
    ax.plot([rout[it, k0], rout[it, k0]], [wout[it, k0], wout[it, k0]], 'r')
    ax.plot(r_range, var[it, :, k0], '-',
            color=cm.jet(count_color), label='t=' + str(np.int(t0)))
    # ax.set_xlim(r_range[np.maximum(imin - 10, 0)], r_range[np.minimum(imax + 10, nr - 1)])
    ax.set_xlim(r_range[imax - 10], r_range[irmax])
    ax.set_ylim(-0.1,0.1)
    plt.suptitle('z=' + str(k0 * dx[2]) + 'm')
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    return

# _______________________________
# _______________________________
def set_input_parameters(args):
    print ''' setting parameters '''
    global path
    path = args.path

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

    global kmin, kmax
    ''' krange '''
    if args.kmin:
        kmin = np.int(args.kmin)
    else:
        kmin = 0
    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = 0

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
    print('kmin: ', kmin)
    print('kmax:', kmax)
    print('nx ', nx)
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
        irstar = np.int(np.round(rstar / dx[0]))
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
        irstar = np.int(np.round(rstar / dx[0]))
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