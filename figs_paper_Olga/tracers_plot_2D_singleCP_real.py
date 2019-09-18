import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os

import plot_CP_singleCP
from plot_CP_singleCP import get_radius_vel

label_size = 12
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['text.usetex'] = 'true'

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    # # parser.add_argument("--k0")
    # parser.add_argument("--tmin")
    # parser.add_argument("--tmax")
    # parser.add_argument("--timerange", nargs='+', type=int)
    # parser.add_argument("--dx")
    # parser.add_argument("--imin")
    # parser.add_argument("--shift_x")
    # parser.add_argument("--shift_t")
    args = parser.parse_args()

    path_olga = '/nbi/ac/conv1/henneb/results/coldpool/lindp2K_13/output/cp/'
    file_olga = 'coldpool_tracer_out_all.txt'
    CP_id_olga = 3
    col_id = 4  # column in textfile of CP-ID

    return
# ----------------------------------------------------------------------

def plot_each_timestep(coordinates, n_tracers, var_name, var, v_rad, lvls_var, lvls_v_rad,
                       k0, shift_t, shift, times, path_out_figs, colmap):
    '''------- plot each timestep -------'''
    for t0 in times[:-1]:
        it = np.int(t0 / dt_fields)
        print('plot: t='+str(t0) + ', '+str(it))
        fig_name = var_name + '_t' + str(np.int(t0)) + '_k'+str(k0) + '.png'
        fig, (ax0, ax1) = plt.subplots(1,2, figsize=(20,10))
        cf = ax0.contourf(var[it,:,:].T, cmap=colmap, levels=lvls_var)
        plt.colorbar(cf, ax=ax0, shrink=0.5)
        cf = ax1.contourf(v_rad[it,:,:].T, cmap=colmap, levels=lvls_v_rad)
        plt.colorbar(cf, ax=ax1, shrink=0.5)
        for ax in [ax0, ax1]:
            ax.set_aspect('equal')
            for i in range(n_tracers):
                ax.plot(coordinates[it+shift_t, i, 0] + shift, coordinates[it+shift_t, i, 1] + shift, 'ok', markersize=2)
        ax0.set_title('w')
        ax1.set_title('v_{rad}')
        plt.suptitle('t-shift '+str(shift_t)+', x-shift '+str(shift))
        plt.tight_layout()
        plt.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)



def plot_each_timestep_circles(coordinates, cp_id, n_tracers, n_cps,
                               var_name, var, v_rad, lvls_var, lvls_v_rad,
                               k0, shift_t, shift, dx, times, path_data, path_tracer_file, path_out_figs, colmap):
    '''------- plot each timestep incl. circle of average radius -------'''
    ''' a) read in CP radius r_torus(t)=r_av and CP spreading velocity u_torus(t)=U_rad_av from tracer statistics '''
    print('path tracers: ', path_tracer_file)
    nt = len(times)
    nk = 1
    dist_av = np.zeros((nt, nk))
    U_rad_av = np.zeros((nt, nk))
    for it, t0 in enumerate(times):
        # print('---t0: ' + str(t0) + '---', it)
        dist_av[it, k0], U_rad_av[it, k0] = get_radius_vel(path_tracer_file, it, cp_id, n_tracers, n_cps)
    r_tracers_av = dist_av * dx
    del dist_av

    ''' b) get azimuthally averaged fields '''
    rad_stats_file = nc.Dataset(os.path.join(path_data, 'data_analysis', 'stats_radial_averaged.nc'))
    v_rad_av = rad_stats_file.groups['stats'].variables['v_rad'][:, :, :]
    r_av = rad_stats_file.groups['stats'].variables['r'][:]  # nr
    rad_stats_file.close()
    v_rad_gradient = np.zeros(shape=v_rad_av.shape)
    for ir, r in enumerate(r_av[1:-1]):
        v_rad_gradient[:, ir, :] = (v_rad_av[:, ir + 1, :] - v_rad_av[:, ir - 1, :]) / (r_av[ir + 1] - r_av[ir - 1])
    ir_vrad_grad_max = np.argmin(v_rad_gradient, axis=1)  # nt, nk

    print ''
    for t0 in times[:-1]:
        if t0 < 1000:
            imin = 100
        elif t0 < 2000:
            imin = 50
        else:
            imin = 0
        it = np.int(t0 / dt_fields)
        print('plot: t=' + str(t0) + ', ' + str(it))
        fig_name = 'circle_' + var_name + '_t' + str(np.int(t0)) + '_k' + str(k0) + '.png'
        fig, axis = plt.subplots(2, 2, figsize=(12, 10))
        ax11 = axis[0, 0]
        ax12 = axis[0, 1]
        ax21 = axis[1, 0]
        ax22 = axis[1, 1]
        radius_tracers = r_tracers_av[it+shift_t, k0] / dx
        radius_gradient  = ir_vrad_grad_max[it,k0]
        circle1_grad = plt.Circle((ic+0.5, jc+0.5), radius_gradient, fill=False, color='b', linewidth=2)
        circle2_grad = plt.Circle((ic+0.5, jc+0.5), radius_gradient, fill=False, color='b', linewidth=2)
        circle1_tr = plt.Circle((ic+0.5, jc+0.5), radius_tracers, fill=False, color='g', linewidth=3)
        circle2_tr = plt.Circle((ic+0.5, jc+0.5), radius_tracers, fill=False, color='g', linewidth=3)
        circle3 = plt.Circle((ic, jc), radius_tracers, fill=False, color='g', linewidth=3)
        circle4 = plt.Circle((ic, jc), radius_tracers, fill=False, color='g', linewidth=3)
        cf = ax11.contourf(var[it, :, :].T, cmap=colmap, levels=lvls_var)
        plt.colorbar(cf, ax=ax11, shrink=0.5)
        cf = ax12.contourf(v_rad[it, :, :].T, cmap=colmap, levels=lvls_v_rad)
        plt.colorbar(cf, ax=ax12, shrink=0.5)
        ax11.set_title('vertical velocity')
        ax11.add_artist(circle1_grad)
        ax11.add_artist(circle1_tr)
        ax12.add_artist(circle2_grad)
        ax12.add_artist(circle2_tr)
        ax11.plot(ic - 0.5, jc - 0.5, 'ow', markersize=7)
        ax12.set_title('radial velocity')
        cf = ax21.contourf(var[it, :, :].T, cmap=colmap, levels=lvls_var)
        plt.colorbar(cf, ax=ax21, shrink=0.5)
        cf = ax22.contourf(v_rad[it, :, :].T, cmap=colmap, levels=lvls_v_rad)
        plt.colorbar(cf, ax=ax22, shrink=0.5)

        for ax in [ax11, ax12, ax21, ax22]:
            # ax.set_aspect('equal')
            for i in range(n_tracers):
                ax.plot(coordinates[it + shift_t, i, 0] + shift, coordinates[it + shift_t, i, 1] + shift, 'ok',
                        markersize=1)

        ax21.add_artist(circle3)
        ax22.add_artist(circle4)
        ax22.set_title('v')
        for ax in axis.flat:
            ax.set_xlim(imin, nx-imin)
            ax.set_ylim(imin, nx-imin)
        plt.tight_layout()
        plt.suptitle('t-shift ' + str(shift_t) + ', x-shift ' + str(shift))
        plt.savefig(os.path.join(path_out_figs, fig_name))
        plt.close()

# ----------------------------------------------------------------------

def plot_tracer_input_fields(coordinates, n_tracers, shift, path_tracers, path_out_figs):
    print(os.path.join(path_tracers, 'input', 'uv_alltimes.nc'))
    rootgrp = nc.Dataset(os.path.join(path_tracers, 'input', 'uv_alltimes.nc'))
    for it, t0 in enumerate(times):
        if it > 0:
            print('-plot time: ' + str(t0))
            fig_name = 'uv_input_fields' + '_t' + str(t0) + '_tracers.png'
            fig, ax = plt.subplots(1, 1, figsize=(7, 6), sharey='all')
            var1 = rootgrp.variables['u'][it, :, :]
            var2 = rootgrp.variables['v'][it, :, :]
            max = np.amax(var1)
            min = -max
            # min = np.amin(var)
            ax.contourf(var1.T, levels=np.linspace(min, max, 1e2), cmap=cm_bwr)
            # ax.contour(var2.T, levels=np.linspace(min, max, 2e1), colors='k', linewidth=0.5)
            ax.contour(var2.T, colors='k', linewidths=0.5)
            ax.set_title('u and v')
            ax.set_aspect('equal')
            # for i in range(n_tracers):
            #     ax.plot(coordinates[it, i, 0] + shift, coordinates[it, i, 1] + shift, 'ok', markersize=2)
            plt.tight_layout()
            fig.savefig(os.path.join(path_out_figs, fig_name))
            plt.close(fig)
    rootgrp.close()

    return



# ----------------------------------------------------------------------


def get_tracer_coords(cp_id, n_cps, n_tracers, times, dt_fields, fullpath_in):
    print('get tracer coordinates')
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    column = lines[0].split()

    nt = len(times)
    coords = np.zeros((nt, n_tracers, 2))

    for it, t0 in enumerate(times):
        print('----t0='+str(t0), it, '----')
        i = 0
        # count = t0 * n_cps * n_tracers + (cp_id - 1) * n_tracers
        # count = it * n_cps * n_tracers + (cp_id - 1) * n_tracers
        count = np.int(t0/dt_fields) * n_cps * n_tracers + (cp_id - 1) * n_tracers
        # while CP age is 0 and CP ID is cp_id
        timestep = int(lines[count].split()[0])
        cp_ID = int(lines[count].split()[3])
        # while (timestep - 1 == it and cp_ID == cp_id):
        while (timestep - 1 == t0/dt_fields and cp_ID == cp_id):
            columns = lines[count].split()
            coords[it,i,0] = float(columns[4])
            coords[it,i,1] = float(columns[5])
            i += 1
            count += 1
            cp_ID = int(lines[count].split()[3])
            timestep = int(lines[count].split()[0])

    f.close()
    # print ''
    return coords


def get_number_cps(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    cp_number = int(lines[-1].split()[3])
    f.close()

    return cp_number


def get_number_tracers(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    count = 0
    # while CP age is 0 and CP ID is 1
    cp_age = int(lines[count].split()[0])
    cp_ID = int(lines[count].split()[3])
    print('cp_age', cp_age)
    while (cp_age == 1 and cp_ID == 1):
        count += 1
        cp_age = int(lines[count].split()[0])
        cp_ID = int(lines[count].split()[3])
    n_tracers = count
    f.close()

    return n_tracers
# ----------------------------------------------------------------------

def set_input_parameters(case_name, path_data):
    nml = simplejson.loads(open(os.path.join(path_data, case_name + '.in')).read())
    global nx, ny
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    # nz = nml['grid']['nz']
    # dx = np.zeros(3, dtype=np.int)
    # dx[0] = nml['grid']['dx']
    # dx[1] = nml['grid']['dy']
    # dx[2] = nml['grid']['dz']
    # gw = nml['grid']['gw']
    # dV = dx[0] * dx[1] * dx[2]

    global ic, jc
    try:
        ic = nml['init']['ic']
        jc = nml['init']['jc']
    except:
        ic = np.int(nx/2)
        jc = np.int(ny/2)
    print('ic, jc', ic, jc)

    return

# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()
