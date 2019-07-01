import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import netCDF4 as nc
import argparse
import json as simplejson
import os
import sys
import time

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

    times, nml, ID = set_input_parameters(args)
    print('id: '+ ID)
    radius = np.int(np.double(ID[12:]) / dx[0])

    path_out_data = os.path.join(path, 'data_analysis')
    if not os.path.exists(path_out_data):
        os.mkdir(path_out_data)
    path_out_data_2D = os.path.join(path, 'fields_v_rad')
    if not os.path.exists(path_out_data_2D):
        os.mkdir(path_out_data_2D)
    path_out_figs = os.path.join(path, 'figs_radial_average')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    print ''
    print 'paths:'
    print 'data:    ', path_out_data
    print 'data 2D: ', path_out_data_2D
    print 'figs:    ', path_out_figs
    print 'fields:  ', path_fields
    print ''

    ic_arr, jc_arr = define_geometry(nml)
    ic = ic_arr[0]
    jc = jc_arr[0]
    irange = np.minimum(nx-ic, ic)
    jrange = np.minimum(ny-jc, jc)
    rmax = np.int(np.ceil(np.sqrt(irange**2 + jrange**2)))

    # plot configuration test file
    plot_configuration(ic, jc, path_out_data_2D)


    print ''
    print('----- compute radius ----------- ')
    # OUTPUT: th_field[nx, ny], r_field[nx, ny]
    file_name_rthfield = 'r_th_field.nc'    # in 'fields_v_rad'
    th_field, r_field = compute_radius(ic, jc, irange, jrange, file_name_rthfield, path_out_data_2D)


    print ''
    print('----- compute radial velocity ----------- ')
    # creates output-file with v_rad[nt, nx, ny, kmax], v_ran[nt, nx, ny, kmax] and r_field[nx, ny]
    file_name_vradfield = 'v_rad.nc'    # in 'fields_v_rad'
    compute_radial_velocity(th_field, r_field, times, file_name_vradfield, ic, jc, 0.5*np.amax(r_field), path_out_data_2D)


    print ''
    print('----- compute angular average ----------- ')
    # OUTPUT: file with angular averaged statistics, e.g. v_rad[nt, nr, nz]
    # file_name_stats = 'stats_radial_averaged_test.nc'
    file_name_stats = 'stats_radial_averaged.nc'
    compute_angular_average(rmax, times, file_name_stats, path_out_data, path_out_data_2D)


    print ''
    print('----- compute angular average CP height ----------- ')
    sth = 0.5
    file_name_CP_height = 'CP_height_' + ID + '_sth' + str(sth) + '.nc' # in path_out_data
    compute_CP_height_radial_av(rmax, times, file_name_CP_height, path_out_data, path_out_data_2D)



    print ''
    print('----- plotting ----------- ')
    file_name_stats = 'stats_radial_averaged.nc'
    plot_radially_averaged_vars(times, file_name_stats, path_out_data, path_out_figs)
    # ----- plot CP height radially averaged
    file_name_CP_height = 'CP_height_' + id + '_sth' + str(sth) + '.nc' # in path_out_data
    plot_radially_averaged_CP_height(times, file_name_CP_height, path_out_data)

    return


# _______________________________
# _______________________________
def compute_radius(ic, jc, irange, jrange, file_name, path_out_data_2D):
    r_field = np.zeros((nx, ny), dtype=np.int)  # radius
    th_field = np.zeros((nx, ny), dtype=np.double)  # angle
    # irange = 10
    # jrange = 10
    for i in range(irange):
        for j in range(jrange):
            r_field[ic+i, jc+j] = np.round(np.sqrt(i**2+j**2))
            r_field[ic-i, jc+j] = r_field[ic+i,jc+j]
            r_field[ic-i, jc-j] = r_field[ic+i,jc+j]
            r_field[ic+i, jc-j] = r_field[ic+i,jc+j]
            if i == 0:
                i = 1e-9
            aux = np.arctan(np.double(j)/i)
            # print(aux, np.pi-aux, np.pi + aux, 2*np.pi - aux)
            th_field[ic+i, jc+j] = aux
            th_field[ic-i, jc+j] = np.pi - aux
            th_field[ic-i, jc-j] = np.pi + aux
            th_field[ic+i, jc-j] = 2*np.pi - aux

    # th_range = np.zeros((2, 4 * nx), dtype=np.double)
    # th_range_aux = np.zeros((nx, ny), dtype=np.int)
    # r0 = 45
    # count = 0
    # for i in range(nx):
    #     for j in range(ny):
    #         if r_field[i,j] < r0 + 1 and r_field[i,j] > r0 -1:
    #             th_range_aux[i, j] = 1
    #             th_range[0, count] = r_field[i,j]
    #             th_range[1, count] = th_field[i,j]
    #             count += 1
    #

    fullpath_out = os.path.join(path_out_data_2D, file_name)



    fig, axis = plt.subplots(3, 2)
    ax1 = axis[0,0]
    ax2 = axis[0,1]
    ax3 = axis[1,0]
    ax4 = axis[1,1]
    ax5 = axis[2,0]
    ax6 = axis[2,1]
    cf = ax1.imshow(r_field[ic-irange:ic+irange+1,jc-jrange:jc+jrange+1].T, origin='lower')
    plt.colorbar(cf, ax=ax1)
    cf = ax2.imshow(th_field[ic-irange:ic+irange+1,jc-jrange:jc+jrange+1].T, origin='lower')
    plt.colorbar(cf, ax=ax2)
    cf = ax3.imshow(r_field[ic-irange:ic+irange+1,jc-jrange:jc+jrange+1].T, origin='lower')
    # cf = ax3.imshow(th_range_aux[ic-irange:ic+irange+1,jc-jrange:jc+jrange+1].T, origin='lower')
    plt.colorbar(cf, ax=ax3)
    cf = ax5.imshow(np.cos(th_field)[ic-irange:ic+irange+1,jc-jrange:jc+jrange+1].T, origin='lower')
    cf = ax6.imshow(np.sin(th_field)[ic-irange:ic+irange+1,jc-jrange:jc+jrange+1].T, origin='lower')
    # plt.colorbar(cf, ax=ax4)
    # ax4.plot(th_range[1, :])
    plt.savefig(os.path.join(path_out_data_2D, 'test_field.png'))
    plt.close()


    # rootgrp = nc.Dataset(fullpath_out, 'w', format='NETCDF4')
    # nx_ = 200
    # rootgrp.createDimension('nx', nx_)
    # rootgrp.createDimension('ny', nx_)
    # # var = rootgrp.createVariable('r_field', 'f8', ('nx', 'ny'))
    # # var[:, :] = r_field[:,:]
    # # var = rootgrp.createVariable('th_field', 'f8', ('nx', 'ny'))
    # # print('aaa', var.shape, th_field.shape)
    # # var[:,:] = np.ones(shape=var.shape)
    # # var[:, :] = th_field[:,:]
    # rootgrp.close()

    return th_field, r_field

# _______________________________

def compute_radial_velocity(th_field, r_field, times, filename, ic, jc, rmax, path_out_data_2D):
    nt = len(times)

    uv_list = ['u', 'v']
    # v_rad = np.zeros((nt, nx, ny, kmax))
    # v_tan = np.zeros((nt, nx, ny, kmax))
    v_rad_int = np.zeros((nt, nx, ny, kmax))
    v_tan_int = np.zeros((nt, nx, ny, kmax))
    for it, t0 in enumerate(times):
        print('t=' + str(t0))
        fullpath_in = os.path.join(path, 'fields', str(t0) + '.nc')
        rootgrp = nc.Dataset(fullpath_in, 'r')
        for k0 in range(kmax):
            print('   k=' + str(k0))
            v_hor_int = np.zeros((2, nx, ny))
            u = rootgrp.groups['fields'].variables['u'][:,:,k0]
            for i in range(1,nx-1):
                v_hor_int[0,i,:] = 0.5*(u[i,:]+u[i-1,:])
            del u
            v = rootgrp.groups['fields'].variables['v'][:,:,k0]
            for j in range(1,ny-1):
                v_hor_int[1,:,j] = 0.5*(v[:,j]+v[:,j-1])
            del v
            v_rad_int[it, :, :, k0], v_tan_int[it, :, :, k0] = compute_radial_vel(v_hor_int, th_field)
            # v_hor = np.zeros((2, nx, ny))
            # v_hor[0,:,:] = rootgrp.groups['fields'].variables['u'][:,:,k0]
            # v_hor[1,:,:] = rootgrp.groups['fields'].variables['v'][:,:,k0]
            # v_rad[it, :, :, k0], v_tan[it, :, :, k0] = compute_radial_vel(v_hor, th_field)
        rootgrp.close()

        # v_hor_diff = v_hor - v_hor_int

        fig, axis = plt.subplots(2, 2, figsize=(10,10))
        ax11 = axis[0, 0]
        ax12 = axis[0, 1]
        ax21 = axis[1, 0]
        ax22 = axis[1, 1]
        circle1 = plt.Circle((ic, jc), rmax/2, fill=False, color='r', linewidth=2)
        cf = ax11.imshow(v_rad_int[it, :, :, 0].T, origin='lower')
        plt.colorbar(cf, ax=ax11, shrink=0.8)
        ax11.set_title('radial velocity')
        ax11.add_artist(circle1)
        ax11.plot(ic-0.5, jc-0.5, 'ow', markersize=7)
        cf = ax12.imshow(v_tan_int[it, :, :, 0].T, origin='lower')
        ax12.set_title('tangential velocity')
        plt.colorbar(cf, ax=ax12, shrink=0.8)
        ax21.imshow(v_hor_int[0, :, :].T, origin='lower')
        plt.colorbar(cf, ax=ax21, shrink=0.8)
        ax21.set_title('u')
        ax22.imshow(v_hor_int[1, :, :].T, origin='lower')
        plt.colorbar(cf, ax=ax22, shrink=0.8)
        ax22.set_title('v')
        ax11.set_xlim(0, nx)
        ax11.set_ylim(0, ny)
        ax12.set_xlim(0, nx)
        ax12.set_ylim(0, ny)
        plt.tight_layout()
        plt.savefig(os.path.join(path_out_data_2D, 'test_field_vrad_vtan_t'+str(t0)+'.png'))
        plt.close()


    create_vrad_field(v_rad_int, v_tan_int, r_field, kmax, filename, path_out_data_2D)
    # del v_rad, v_tan
    del v_rad_int, v_tan_int
    return

# _______________________________

def compute_angular_average(rmax, times, file_name, path_out_data, path_out_data_2D):
    t0 = time.time()
    var_list = ['w', 's', 'phi', 'temperature']
    # file_name = 'stats_radial_averaged_test.nc'
    create_statistics_file(var_list, file_name, times, rmax, path_out_data)

    print('1', time.time() - t0)
    t0 = time.time()
    data_dict_av = {}
    for var_name in var_list:
        data_dict_av[var_name] = np.zeros((rmax, kmax))
    v_rad_av = np.zeros((rmax, kmax))
    v_tan_av = np.zeros((rmax, kmax))
    print('2', time.time() - t0)
    t0 = time.time()
    # read in v_rad-field
    file_name_vradfield = 'v_rad.nc'
    fullpath_in = os.path.join(path_out_data_2D, file_name_vradfield)
    rootgrp = nc.Dataset(fullpath_in)
    v_rad = rootgrp.variables['v_rad'][:, :, :, :]
    v_tan = rootgrp.variables['v_tan'][:, :, :, :]
    r_field = rootgrp.variables['r_field'][:, :]
    rootgrp.close()
    if v_rad.shape[3] < kmax:
        print('v_rad-fields contain less levels than kmax=' + str(kmax))
        sys.exit()
    print('3', time.time() - t0)
    t0 = time.time()
    for it, t0 in enumerate(times):
        fullpath_in = os.path.join(path, 'fields', str(t0) + '.nc')
        data_dict = read_in_vars(fullpath_in, var_list)
        print('running t=' + str(t0))
        v_rad_av[:, :] = compute_average_var(v_rad[it, :, :, :kmax], rmax, kmax, r_field)
        v_tan_av[:, :] = compute_average_var(v_tan[it, :, :, :kmax], rmax, kmax, r_field)
        for var_name in var_list:
            data_dict_av[var_name][:, :] = compute_average_var(data_dict[var_name][:, :, :], rmax, kmax, r_field)

        dump_statistics_file(data_dict_av, v_rad_av, v_tan_av, var_list, it, file_name, path_out_data)
    return

# _______________________________

def create_vrad_field(v_rad, v_tan, r_field, kmax, file_name, path_out_data_2D):
    # file_name = 'v_rad.nc'
    rootgrp = nc.Dataset(os.path.join(path_out_data_2D, file_name), 'w', format='NETCDF4')

    rootgrp.createDimension('time', None)
    rootgrp.createDimension('nx', nx)
    rootgrp.createDimension('ny', ny)
    rootgrp.createDimension('nz', kmax)

    var = rootgrp.createVariable('v_rad', 'f8', ('time', 'nx', 'ny', 'nz'))
    var[:,:,:,:] = v_rad
    var = rootgrp.createVariable('v_tan', 'f8', ('time', 'nx', 'ny', 'nz'))
    var[:,:,:,:] = v_tan
    var = rootgrp.createVariable('r_field', 'f8', ('nx', 'ny'))
    var[:,:] = r_field
    rootgrp.close()

    return
# _______________________________

def create_statistics_file(var_list, file_name, timerange, rmax, path_out_data):
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
    var = stats_grp.createVariable('v_tan', 'f8', ('nt', 'nr', 'nz'))
    var.long_name = 'tangential velocity'
    var.units = "m/s"

    for var_name in var_list:
        stats_grp.createVariable(var_name, 'f8', ('nt', 'nr', 'nz'))
        # var.units = "m"

    rootgrp.close()

    return



def dump_statistics_file(data_dictionary, v_rad_av, v_tan_av, var_list, it, file_name, path_out_data):
    print('-------- dump statistics file -------- ')
    print(os.path.join(path_out_data, file_name))
    rootgrp = nc.Dataset(os.path.join(path_out_data, file_name), 'r+', format='NETCDF4')

    # ts_grp = rootgrp.groups['timeseries']
    # var = ts_grp.variables['r_av']
    # var[it,ik] = rim_vel_av[0]
    # var = ts_grp.variables['U_av']
    # var[it,ik] = rim_vel_av[1]

    stats_grp = rootgrp.groups['stats']
    for var_name in var_list:
        var = stats_grp.variables[var_name]
        var[it, :, :] = data_dictionary[var_name][:, :]
    var = stats_grp.variables['v_rad']
    var[it, :, :] = v_rad_av[:, :]
    var = stats_grp.variables['v_tan']
    var[it, :, :] = v_tan_av[:, :]

    rootgrp.close()
    return

# _______________________________

def read_in_vars(fullpath_in, var_list):

    var_dict = {}

    rootgrp = nc.Dataset(fullpath_in, 'r')
    for var_name in var_list:
        var = rootgrp.groups['fields'].variables[var_name]
        data = var[:, :, :kmax]     # x, y, z
        var_dict[var_name] = data
    rootgrp.close()

    return var_dict

# _______________________________

def compute_average_var(var, rmax, kmax_, r_field):
    count = np.zeros(rmax, dtype=np.int)

    if kmax_ == 1:
        var_av = np.zeros((rmax), dtype=np.double)
        for i in range(nx):
            for j in range(ny):
                r = r_field[i, j]
                count[r] += 1
                var_av[r] += var[i, j]

    elif kmax_ > 1:
        var_av = np.zeros((rmax, kmax_), dtype=np.double)
        for i in range(nx):
            for j in range(ny):
                r = r_field[i, j]
                count[r] += 1
                var_av[r, :] += var[i, j, :]

    for r in range(rmax):
        if count[r] > 0:
            var_av[r] /= count[r]

    return var_av
# _______________________________

def compute_radial_vel(uv, th_field):
    ur = np.zeros((nx,ny), dtype=np.double)
    uth = np.zeros((nx,ny), dtype=np.double)
    for i in range(nx):
        for j in range(ny):
            th = th_field[i,j]
            # # # clockwise rotation
            # # ur[i,j] = uv[0,i,j]*np.cos(th) + uv[1,i,j]*np.sin(th)
            # counter-clockwise rotation
            ur[i,j] = uv[0,i,j]*np.cos(th) + uv[1,i,j]*np.sin(th)
            uth[i,j] = -uv[0,i,j]*np.sin(th) + uv[1,i,j]*np.cos(th)
    return ur, uth

# _______________________________

def compute_CP_height_radial_av(rmax, times, file_name_in, path_out_data, path_out_data_2D):
    # (1) compute angular average (compute_average_var)
    # (2) append CP_height(t,r) to stats-file data_analysis/CP_height_dTh3_z2000_r600_sth0.5.nc

    # ----- read in CP height 2D field
    print(os.path.join(path_out_data, file_name_in))
    file = nc.Dataset(os.path.join(path_out_data, file_name_in))
    time_CP_height = file.groups['timeseries'].variables['time'][:]

    if len(times) != len(time_CP_height) or (time_CP_height - times).any() != 0:
        if time_CP_height[0] != times[0] or len(times) > len(time_CP_height):
            print('time dimensions not fitting')
            print times[0], time_CP_height[0]
            sys.exit()
        else:
            itmax = np.where(time_CP_height == times[-1])[0][0] + 1
    else:
        itmax = len(time_CP_height)
    CP_height_2d = file.groups['fields_2D'].variables['CP_height_2d'][:itmax, :, :]
    file.close()

    # ----- read in r_field
    file_name_vradfield = 'v_rad.nc'
    fullpath_in = os.path.join(path_out_data_2D, file_name_vradfield)
    rootgrp = nc.Dataset(fullpath_in)
    r_field = rootgrp.variables['r_field'][:, :]
    rootgrp.close()

    # ----- azimuthally average CP_height_field
    CP_height_av = np.zeros((itmax, rmax))
    for it, t0 in enumerate(times):
        CP_height_av[it, :] = compute_average_var(CP_height_2d[it, :, :], rmax, 1, r_field)

    # ----- dump
    file = nc.Dataset(os.path.join(path_out_data, file_name_in), 'r+')
    if not 'stats' in file.groups.keys():
        stats_grp = file.createGroup('stats')
        stats_grp.createDimension('nt', len(time_CP_height))
        stats_grp.createDimension('nr', rmax)
    else:
        stats_grp = file.groups['stats']
    if not 'r' in stats_grp.variables.keys():
        ri_range = np.arange(0, rmax, dtype=np.int)
        r_range = np.arange(0, rmax, dtype=np.double) * dx[0]
        var = stats_grp.createVariable('r', 'f8', ('nr'))
        var.unit = 'm'
        var[:] = r_range
        var = stats_grp.createVariable('ri', 'f8', ('nr'))
        var.unit = '-'
        var[:] = ri_range
    if not 'CP_height_rad' in stats_grp.variables.keys():
        var = stats_grp.createVariable('CP_height_rad', 'f8', ('nt', 'nr'))
        var.long_name = 'CP height, azimuthally averaged'
        var.units = "m/s"
    else:
        var = stats_grp.variables['CP_height_rad'][:, :]
    var[:itmax, :] = CP_height_av[:, :]
    file.close()
    return
# _______________________________
# _______________________________

def plot_radially_averaged_vars(times, file_name, path_out_data, path_out_figs):

    print path_out_figs
    print ''
    # read in file
    # file_name = 'stats_radial_averaged.nc'
    data_stats = nc.Dataset(os.path.join(path_out_data, file_name), 'r')
    stats_grp = data_stats.groups['stats'].variables
    r_range = stats_grp['r'][:]
    times_ = data_stats.groups['timeseries'].variables['time'][:]
    ta = np.where(times_ == tmin)[0][0]
    tb = np.where(times_ == tmax)[0][0]
    print ''
    print 'times'
    print times
    print times_
    print ta, tb

    # var = np.array('nt', 'nr', 'nz')
    var_list = ['w', 'v_rad', 'v_tan', 's']
    ncol = len(var_list)
    rmax_plot = 10e3
    irmax = np.where(r_range == rmax_plot)[0][0]

    for k0 in range(kmax):
        print('-- k='+str(k0))
        fig_name = 'radial_average_k'+str(k0)+'.png'
        fig, axes = plt.subplots(1, ncol, sharey='none', figsize=(5*ncol, 5))
        for i, ax in enumerate(axes):
            var = stats_grp[var_list[i]][:,:,:]
            max = np.amax(var[:, :irmax, :kmax])
            min = np.amin(var[:, :irmax, :kmax])
            if len(times) >= 2:
                for it, t0 in enumerate(times_[1::2]):
                    if it >= ta and it <= tb:
                        count_color = 2 * np.double(it) / len(times)
                        ax.plot(r_range[:irmax], var[2*it+1, :irmax, k0], color=cm.jet(count_color), label='t='+str(np.int(t0)))
            else:
                for it, t0 in enumerate(times_):
                    if it >= ta and it <= tb:
                        count_color = np.double(it) / len(times)
                        ax.plot([0, r_range[irmax]], [0.,0.], 'k')
                        ax.plot(r_range[:irmax], var[it, :irmax, k0], color=cm.jet(count_color), label='t='+str(t0))
            if var_list[i] == 's':
                ax.set_ylim(6850, max+10)
            elif var_list[i] == 'w':
                ax.set_ylim(-5, max)
            else:
                ax.set_ylim(min, max)
            ax.set_title(var_list[i])
            ax.set_xlabel('radius r  [m]')
            ax.set_ylabel(var_list[i])
        # axes[0].legend(loc='lower left', bbox_to_anchor=(0, 0.),
        #            fancybox=True, shadow=True, ncol=4, fontsize=6)
        axes[-1].legend(loc='upper center', bbox_to_anchor=(1.2, 1.),
                   fancybox=True, shadow=True, ncol=1, fontsize=10)
        plt.subplots_adjust(bottom=0.12, right=.9, left=0.04, top=0.9, wspace=0.15)
        plt.suptitle('z='+str(k0*dx[2])+'m')
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)

    data_stats.close()

    return



def plot_radially_averaged_CP_height(times, file_name_in, path_out_data):
    data = nc.Dataset(os.path.join(path_out_data, file_name_in), 'r')
    time_CP_height = data.groups['timeseries'].variables['time'][:]
    CP_height_av = data.groups['stats'].variables['CP_height_rad'][:, :]
    r_range = data.groups['stats'].variables['r'][:]
    data.close()
    ta = np.where(time_CP_height == tmin)[0][0]
    tb = np.where(time_CP_height == tmax)[0][0]
    rmax_plot = 10e3
    irmax = np.where(r_range == rmax_plot)[0][0]
    fig_name = 'CP_height_radial_average.png'
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 5))
    print CP_height_av.shape
    print time_CP_height.shape
    if len(times) >= 2:
        for it, t0 in enumerate(time_CP_height[1::2]):
            if it >= ta and it <= tb:
                count_color = 2 * np.double(it) / len(times)
                ax0.plot(r_range[:irmax], CP_height_av[2 * it + 1, :irmax], color=cm.jet(count_color),
                         label='t=' + str(np.int(t0)))
    else:
        for it, t0 in enumerate(time_CP_height):
            if it >= ta and it <= tb:
                count_color = np.double(it) / len(times)
                ax0.plot([0, r_range[irmax]], [0., 0.], 'k')
                ax0.plot(r_range[:irmax], CP_height_av[it, :irmax], color=cm.jet(count_color), label='t=' + str(t0))
    ax0.set_xlabel('radius r  [m]')
    ax0.set_ylabel('CP height  [m]')
    ax0.legend(loc='upper center', bbox_to_anchor=(1.2, 1.),
               fancybox=True, shadow=True, ncol=1, fontsize=10)
    plt.subplots_adjust(bottom=0.12, right=.9, left=0.04, top=0.9, wspace=0.15)
    fig.savefig(os.path.join(path, 'figs_CP_height', fig_name))
    plt.close(fig)
    return

# _______________________________
def plot_configuration(ic, jc, path_out):
    fig_name = 'test_config.png'
    fullpath_in = os.path.join(path, 'fields', '0.nc')
    rootgrp = nc.Dataset(fullpath_in, 'r')
    s = rootgrp.groups['fields'].variables['s'][:, :, :]
    fig, (ax1,ax2,ax3) = plt.subplots(1, 3, sharey='none', figsize=(16, 5))
    ax1.imshow(s[:,:,0].T, origin='lower')
    ax1.plot(ic, jc, 'ko', markersize=3)
    ax1.plot([ic,ic], [0,ny], 'k-')
    ax1.plot([0,nx], [jc,jc], 'k-')
    di = 20
    # print imin, imax, jc, s.shape
    ax2.imshow(s[:,jc,:100].T, origin='lower')
    ax2.plot([jc,jc], [0,100], 'k-', linewidth=1)
    ax3.imshow(s[ic, :,:100].T, origin='lower')
    ax3.plot([ic,ic], [0,100], 'k-', linewidth=1)

    radius = np.int(np.double(os.path.basename(path[:-1])[12:])/dx[0])
    circle1 = plt.Circle((ic, jc), radius, fill=False, color='r', linewidth=2)
    circle2 = plt.Circle((ic+1, jc+1), radius, fill=False, color='r', linewidth=2)
    circle3 = plt.Circle((ic-1, jc-1), radius, fill=False, color='r', linewidth=2)
    ax1.add_artist(circle1)
    ax1.add_artist(circle2)
    ax1.add_artist(circle3)

    ax1.set_xlim(ic-di,ic+di+1)
    ax1.set_ylim(jc-di,jc+di+1)
    ax2.set_xlim(ic-di,ic+di+1)
    ax2.set_ylim(0,40)
    ax3.set_xlim(jc-di,jc+di+1)
    ax3.set_ylim(0,40)
    fig.savefig(os.path.join(path_out, fig_name))
    plt.close(fig)
    del s
    return

# _______________________________
# _______________________________
def set_input_parameters(args):
    print ''' setting parameters '''
    global path, path_fields
    path = args.path
    path_fields = os.path.join(path, 'fields')

    global case_name
    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    id = os.path.basename(path[:])
    if id == '':
        id = os.path.basename(path[:-1])
    global nx, ny, nz, dx, gw
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
        kmax = np.int(args.kmax)+1
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

    return times, nml, id

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
