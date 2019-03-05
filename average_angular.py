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
    parser.add_argument("--kmax")
    args = parser.parse_args()

    times, nml = set_input_parameters(args)
    nt = len(times)
    ic_arr, jc_arr = define_geometry(nml)

    r_field = np.zeros((nx, ny), dtype=np.int)    # radius
    th_field = np.zeros((nx, ny), dtype=np.int)   # angle

    ic = ic_arr[0]
    jc = jc_arr[0]
    irange = np.minimum(nx-ic, ic)
    jrange = np.minimum(ny-jc, jc)
    rmax = np.int(np.ceil(np.sqrt(irange**2 + jrange**2)))

    # compute radius
    for i in range(irange):
        for j in range(jrange):
            r_field[ic+i,jc+j] = np.round(np.sqrt(i**2+j**2))
            r_field[ic-i,jc+j] = r_field[ic+i,jc+j]
            r_field[ic-i,jc-j] = r_field[ic+i,jc+j]
            r_field[ic+i,jc-j] = r_field[ic+i,jc+j]
            if i == 0:
                i = 1e-9
            aux = np.arctan(j/i)
            th_field[ic+i, jc+j] = aux
            th_field[ic-i, jc+j] = aux
            th_field[ic-i, jc+j] = -aux
            th_field[ic+i, jc-j] = -aux


    # compute radial velocity
    print ''
    print('-------- compute radial velocity -------- ')
    uv_list = ['u', 'v']
    t0 = tmin
    it = 0
    k0 = 0
    v_rad = np.zeros((nt, nx, ny, kmax))
    for it, t0 in enumerate(times):
        print('t=' + str(t0))
        fullpath_in = os.path.join(path, 'fields', str(t0) + '.nc')
        rootgrp = nc.Dataset(fullpath_in, 'r')
        for k0 in range(kmax):
            v_hor = np.zeros((2, nx, ny))
            v_hor[0,:,:] = rootgrp.groups['fields'].variables['u'][:,:,k0]
            v_hor[1,:,:] = rootgrp.groups['fields'].variables['v'][:,:,k0]
            v_rad[it, :, :, k0] = compute_radial_vel(v_hor, th_field)
        rootgrp.close()
    create_vrad_field(v_rad, kmax)





    # compute angular average
    var_list = ['w', 's', 'phi', 'temperature']

    file_name = 'stats_radial_averaged.nc'
    create_statistics_file(var_list, file_name, times, rmax)

    data_dict_av = {}
    for var_name in var_list:
        data_dict_av[var_name] = np.zeros((nt, rmax, kmax))
    v_rad_av = np.zeros((nt, rmax, kmax))

    # t0 = tmin
    # it = 0
    # k0 = 0
    for it, t0 in enumerate(times):
        fullpath_in = os.path.join(path, 'fields', str(t0) + '.nc')
        data_dict = read_in_vars(fullpath_in, var_list)
        print('running t='+str(t0))
        for k in range(kmax):
            v_rad_av[it, :, k] = compute_average_var(v_rad[it, :, :, k], rmax, r_field)
            for var_name in var_list:
                data_dict_av[var_name][it, :, k] = compute_average_var(data_dict[var_name][:,:,k], rmax, r_field)

    dump_statistics_file(data_dict_av, v_rad_av, var_list, file_name)






    ''' plotting '''
    data_stats = nc.Dataset(os.path.join(path_out_data, file_name), 'r')
    stats_grp = data_stats.groups['stats'].variables
    times = data_stats.groups['timeseries'].variables['time'][:]
    # var = np.array('nt', 'nr', 'nz')
    r_range = stats_grp['r'][:]

    var_list = ['w', 'v_rad', 's']
    ncol = len(var_list)
    rmax_plot = 10e3
    irmax = np.where(r_range == rmax_plot)[0]

    for k0 in range(kmax):
        fig_name = 'radial_average_k'+str(k0)+'.png'
        fig, axes = plt.subplots(1, ncol, sharey='none', figsize=(5*ncol, 5))
        for i, ax in enumerate(axes):
            var = stats_grp[var_list[i]][:,:,:]
            if len(times) >= 2:
                for it, t0 in enumerate(times[1::2]):
                    count_color = 2 * np.double(it) / len(times)
                    if var_list[i] == 's':
                        ax.plot(r_range[:-1], var[2*it+1, :-1, 0], color=cm.jet(count_color), label='t='+str(t0))
                    else:
                        ax.plot(r_range, var[2*it+1, :, 0], color=cm.jet(count_color), label='t='+str(t0))
            else:
                for it, t0 in enumerate(times):
                    count_color = np.double(it) / len(times)
                    if var_list[i] == 's':
                        ax.plot(r_range[:rmax_plot], var[it, :rmax_plot, 0], color=cm.jet(count_color), label='t='+str(t0))
                    else:
                        ax.plot(r_range[:rmax_plot], var[it, :rmax_plot, 0], color=cm.jet(count_color), label='t='+str(t0))
            ax.set_title(var_list[i])
        axes[2].legend(loc='upper center', bbox_to_anchor=(1.2, 1.),
                   fancybox=True, shadow=True, ncol=1, fontsize=10)
        plt.suptitle('z='+str(k0*dx[2])+'m')
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)

    data_stats.close()

    return


# _______________________________

def create_vrad_field(v_rad, kmax):
    path_out = os.path.join(path, 'fields_v_rad')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    file_name = 'v_rad.nc'
    rootgrp = nc.Dataset(os.path.join(path_out, file_name), 'w', format='NETCDF4')

    rootgrp.createDimension('time', None)
    rootgrp.createDimension('nx', nx)
    rootgrp.createDimension('ny', ny)
    rootgrp.createDimension('nz', kmax)

    var = rootgrp.createVariable('v_rad', 'f8', ('time', 'nx', 'ny', 'nz'))
    var[:,:,:,:] = v_rad
    rootgrp.close()

    return
# _______________________________

def create_statistics_file(var_list, file_name, timerange, rmax):
    print('-------- create statistics file -------- ')
    print(path_out_data + ', ' + file_name)
    print('')

    nt = len(timerange)

    rootgrp = nc.Dataset(os.path.join(path_out_data, file_name), 'w', format='NETCDF4')

    dims_grp = rootgrp.createGroup('dimensions')
    dims_grp.createDimension('dx', dx[0])
    dims_grp.createDimension('dy', dx[1])
    dims_grp.createDimension('dz', dx[2])
    dims_grp.createDimension('nz', kmax)
    var = dims_grp.createVariable('krange', 'f8', ('nz'))
    var[:] = np.arange(0, kmax, dtype=np.int)

    ts_grp = rootgrp.createGroup('timeseries')
    ts_grp.createDimension('nt', nt)
    var = ts_grp.createVariable('time', 'f8', ('nt'))
    var.unit = "s"
    var[:] = timerange

    stats_grp = rootgrp.createGroup('stats')
    stats_grp.createDimension('nt', nt)
    stats_grp.createDimension('nr', rmax)
    stats_grp.createDimension('nz', kmax)
    ri_range = np.arange(0, rmax, dtype=np.int)
    r_range = np.arange(0, rmax, dtype=np.double)*dx[0]
    var = stats_grp.createVariable('r', 'f8', ('nr'))
    var.unit = 'm'
    var[:] = r_range
    var = stats_grp.createVariable('ri', 'f8', ('nr'))
    var.unit = '-'
    var[:] = ri_range

    var = stats_grp.createVariable('v_rad', 'f8', ('nt', 'nr', 'nz'))
    var.long_name = 'radial velocity'
    var.units = "m/s"

    for var_name in var_list:
        stats_grp.createVariable(var_name, 'f8', ('nt', 'nr', 'nz'))
        # var.units = "m"

    rootgrp.close()

    return



def dump_statistics_file(data_dictionary, v_rad_av, var_list, file_name):
    print('-------- dump statistics file -------- ')
    rootgrp = nc.Dataset(os.path.join(path_out_data, file_name), 'r+', format='NETCDF4')

    # ts_grp = rootgrp.groups['timeseries']
    # var = ts_grp.variables['r_av']
    # var[it,ik] = rim_vel_av[0]
    # var = ts_grp.variables['U_av']
    # var[it,ik] = rim_vel_av[1]

    stats_grp = rootgrp.groups['stats']
    for var_name in var_list:
        var = stats_grp.variables[var_name]
        var[:,:,:] = data_dictionary[var_name][:,:,:]
    var = stats_grp.variables['v_rad']
    var[:, :, :] = v_rad_av[:, :, :]

    # prof_grp = rootgrp.groups['profiles']
    # var = prof_grp.variables['k_dumped']
    # var[ik] = 1

    rootgrp.close()
    return


def read_in_vars(fullpath_in, var_list):

    var_dict = {}

    rootgrp = nc.Dataset(fullpath_in, 'r')
    for var_name in var_list:
        var = rootgrp.groups['fields'].variables[var_name]
        data = var[:, :, :]
        var_dict[var_name] = data
    rootgrp.close()

    return var_dict


# _______________________________

def compute_average_var(var1, rmax, r_field):

    var1_av = np.zeros(rmax, dtype=np.double)
    count = np.zeros(rmax, dtype=np.int)
    for i in range(nx):
        for j in range(ny):
            r = r_field[i, j]
            count[r] += 1
            var1_av[r] += var1[i, j]

    for r in range(rmax):
        if count[r] > 0:
            var1_av[r] /= count[r]

    return var1_av
# _______________________________

def compute_radial_vel(uv, th_field):
    ur = np.zeros((nx,ny), dtype=np.double)
    for i in range(nx):
        for j in range(ny):
            th = th_field[i,j]
            ur[i,j] = uv[0,i,j]*np.cos(th) + uv[1,i,j]*np.sin(th)
    return ur

# _______________________________
# _______________________________
def set_input_parameters(args):
    print ''' setting parameters '''
    global path, path_fields, path_out_data, path_out_figs
    path = args.path
    path_fields = os.path.join(path, 'fields')

    path_out_data = os.path.join(path, 'data_analysis')
    if not os.path.exists(path_out_data):
        os.mkdir(path_out_data)
    path_out_figs = os.path.join(path, 'figs_radial_average')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print ''
    print 'paths:'
    print path_out_data
    print path_out_figs
    print path_fields
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

    ''' time range '''
    times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
             and tmin <= np.int(name[:-3]) <= tmax]
    times.sort()

    print('')
    print('times', times)
    print('')
    print('kmax ', kmax, 'nx ', nx)
    print('')

    return times, nml

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
