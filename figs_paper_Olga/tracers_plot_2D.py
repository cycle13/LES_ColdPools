import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os

label_size = 10
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.labelsize'] = 15

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    # parser.add_argument("casename")
    # parser.add_argument("path")
    # parser.add_argument("--k0")
    # parser.add_argument("--tmin")
    # parser.add_argument("--tmax")
    # parser.add_argument("--path_tracers")
    parser.add_argument("--t_shift")
    parser.add_argument("--x_shift")
    args = parser.parse_args()
    # set_input_parameters(args)

    path_olga = '/nbi/ac/conv1/henneb/results/coldpool/lindp2K_13/output/cp/'
    file_olga = 'coldpool_tracer_out_all.txt'
    CP_id_olga = 3
    col_id = 4  # column in textfile of CP-ID
    # tracer file:
    # - col=0:      timestep of simulation
    # - col=1:      age of CP (timestep since beginning of first CP)
    # - col=2:      tracer id
    # - col=3:      CP id
    # - col=4,5:
    # - col=14,15:  CP center (xc,yc)

    path_out_figs = os.path.join('/nbi/home/meyerbe/paper_olga/figs_2D_real/')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)

    path_tracer_file = os.path.join(path_olga, file_olga)

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')

    print('???? dt_fields')
    dt_fields = 100

    if args.t_shift:
        t_shift = np.int(args.t_shift)
    else:
        t_shift = 0
    if args.x_shift:
        x_shift = np.int(args.x_shift)
    else:
        x_shift = 0
    print('t_shift: '+str(t_shift))
    print('x_shift: '+str(x_shift))


    print('')
    print('path figs:    ' + path_out_figs)
    print('path tracers: ' + path_tracer_file)
    print('')




    cp_id = CP_id_olga
    print('CP ID: ' + str(cp_id))
    # n_cps = get_number_cps(path_tracer_file)
    n_cps = 2265
    print('number of CPs: ' + str(n_cps))
    # n_tracers = get_number_tracers(path_tracer_file)
    n_tracers = 999
    print('number of tracers: ' + str(n_tracers))

    tau, t0, t1 = get_cp_lifetime(cp_id, path_tracer_file)
    xc, yc = get_cp_center(cp_id, tau, n_tracers, path_tracer_file)
    # tau = 5
    # t0 = 47
    # t1 = 51
    lifetime = np.arange(t0, t1+1)
    print('CP lifetime: ' + str(tau))
    print('CP center: '+ str(xc)+', '+str(yc))
    print('')

    ''' tracer coordinates '''
    coordinates = get_tracer_coords(cp_id, n_cps, n_tracers, lifetime, path_tracer_file)
    print('')
    print('COORDINATES')
    print(coordinates[-1,:,:])
    print('')

    ''' averaging '''
    aux = np.zeros(shape=coordinates.shape)
    # for it in range(tau):
    #     aux[it,:,:] = coordinates[it,:,:] -


    '''' plot input fields '''
    # plot_tracers(coordinates, cp_id, n_tracers, x_shift, t_shift, lifetime, path_out_figs)
    path_fields = path_out_figs
    plot_input_fields(coordinates, xc, yc, cp_id, n_tracers, x_shift, t_shift, lifetime,
                      'u', path_fields, path_out_figs)





    # var_list = ['s', 'w']
    # for it,t0 in enumerate(times):
    #     print('-plot time: '+str(t0))
    #     fig_name = 's_w' + '_t' + str(t0) + '_tracers.png'
    #     fig, axis = plt.subplots(1, 2, figsize=(11, 6), sharey='all')
    #     rootgrp = nc.Dataset(os.path.join(path_fields, str(t0)+'.nc'))
    #     for j, var_name in enumerate(var_list):
    #         print var_name, j
    #         var = rootgrp.groups['fields'].variables[var_name][:,:,k0]
    #         max = np.amax(var)
    #         if var_name in ['w', 'v_rad', 'v_tan']:
    #             min = -max
    #             cm_ = cm_bwr
    #         else:
    #             min = np.amin(var)
    #             cm_ = cm_hsv
    #         # axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_)
    #         axis[j].imshow(var.T, vmin=min, vmax=max, cmap=cm_, origin='lower')
    #         axis[j].set_title(var_name)
    #         axis[j].set_aspect('equal')
    #     rootgrp.close()
    #     for i in range(n_tracers):
    #         for j in range(len(var_list)):
    #             axis[j].plot(coordinates[it,i,0]+shift, coordinates[it,i,1]+shift, 'ok', markersize=2)
    #     plt.tight_layout()
    #     fig.savefig(os.path.join(path_out_figs, fig_name))
    #     plt.close(fig)
    #
    #
    # var_list = ['v_rad', 'v_tan']
    # rootgrp = nc.Dataset(os.path.join(path, 'fields_v_rad', 'v_rad.nc'))
    # for it, t0 in enumerate(times):
    #     print('-plot time: ' + str(t0))
    #     fig_name = 'v_rad_tan' + '_t' + str(t0) + '_tracers.png'
    #     fig, axis = plt.subplots(1, 2, figsize=(11, 6), sharey='all')
    #     for j, var_name in enumerate(var_list):
    #         print var_name, j
    #         var = rootgrp.variables[var_name][it, :, :, k0]
    #         max = np.amax(var)
    #         if var_name in ['w', 'v_rad', 'v_tan']:
    #             min = -max
    #         else:
    #             min = np.amin(var)
    #         axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_bwr)
    #         axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_bwr)
    #         axis[j].set_title(var_name)
    #         axis[j].set_aspect('equal')
    #     for i in range(n_tracers):
    #         for j in range(len(var_list)):
    #             axis[j].plot(coordinates[it,i,0]+shift, coordinates[it,i,1]+shift, 'ok', markersize=2)
    #     plt.tight_layout()
    #     fig.savefig(os.path.join(path_out_figs, fig_name))
    #     plt.close(fig)
    # rootgrp.close()
    #
    #
    #
    # var_list = ['u', 'v']
    # rootgrp = nc.Dataset(os.path.join(path_tracers, 'input', 'uv_alltimes.nc'))
    # for it, t0 in enumerate(times):
    #     print('-plot time: ' + str(t0))
    #     fig_name = 'uv_alltimes' + '_t' + str(t0) + '_tracers.png'
    #     fig, axis = plt.subplots(1, 2, figsize=(11, 6), sharey='all')
    #     for j, var_name in enumerate(var_list):
    #         print var_name, j
    #         var = rootgrp.variables[var_name][it, :, :]
    #         max = np.amax(var)
    #         if var_name in ['w', 'v_rad', 'v_tan']:
    #             min = -max
    #         else:
    #             min = np.amin(var)
    #         axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_bwr)
    #         axis[j].contourf(var.T, levels=np.linspace(min, max, 1e2), cmap=cm_bwr)
    #         axis[j].set_title(var_name)
    #         axis[j].set_aspect('equal')
    #     for i in range(n_tracers):
    #         for j in range(len(var_list)):
    #             axis[j].plot(coordinates[it,i,0]+shift, coordinates[it,i,1]+shift, 'ok', markersize=2)
    #     plt.tight_layout()
    #     fig.savefig(os.path.join(path_out_figs, fig_name))
    #     plt.close(fig)
    # rootgrp.close()


    return
# ----------------------------------------------------------------------

def plot_tracers(coordinates, cp_id, n_tracers, x_shift, t_shift, lifetime, path_out_figs):
    print('plot tracers: ', n_tracers)
    lvls = np.linspace(0,1,n_tracers)
    for it, t0 in enumerate(lifetime):
        print('-plot time: ' + str(t0))
        fig_name = 'cp'+str(cp_id)+'_tracers_t' + str(t0) + '_tracers.png'
        fig, axis = plt.subplots(2, 2, figsize=(7, 6), sharey='none')
        ax0 = axis[0,0]
        ax1 = axis[0,1]
        ax2 = axis[1,0]
        ax0.set_aspect('equal')
        for i in range(n_tracers):
            ax0.plot(coordinates[it, i, 0] + shift, coordinates[it, i, 1] + shift,
                    'o', color=cm_hsv(lvls[i]), markeredgecolor='w', markersize=5)
            ax1.plot(i, coordinates[it, i, 1],
                     'o', color=cm_hsv(lvls[i]), markeredgecolor='w', markersize=5)
            ax2.plot(coordinates[it,i,0], i,
                     'o', color=cm_hsv(lvls[i]), markeredgecolor='w', markersize=5)
        plt.tight_layout()
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)
    return




def plot_input_fields(coordinates, xc, yc, cp_id, n_tracers, x_shift, t_shift, lifetime,
                      field_name, path_fields, path_out_figs):
    k0 = 0
    rootgrp = nc.Dataset(os.path.join(path_fields, 'input_u.nc'))
    var1 = rootgrp.variables['var1'][:, k0, :, :]
    rootgrp.close()
    rootgrp = nc.Dataset(os.path.join(path_fields, 'input_v.nc'))
    var2 = rootgrp.variables['var1'][:, k0, :, :]
    rootgrp.close()
    print('plot input fields: ')
    print(lifetime)

    for it, t0 in enumerate(lifetime):
        print('-plot time: ' + str(t0), it)
        fig_name = 'uv_input_fields' + '_t' + str(t0) + '_tracers.png'
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 5), sharey='all')
        max = np.amax(var1[t0,:,:])
        min = -max
        lvls = np.linspace(min, max, 10)
        cf = ax0.contourf(var1[t0,:,:], levels=lvls, cmap=cm_bwr)
        ax0.set_title('u')
        plt.colorbar(cf, ax=ax0)
        max = np.amax(var2[t0,:,:])
        min = -max
        lvls = np.linspace(min, max, 10)
        cf = ax1.contourf(var2[t0,:,:], levels=lvls, cmap=cm_bwr)
        ax1.set_title('v')
        plt.colorbar(cf, ax=ax1)
        # # ax.contour(var2.T, levels=np.linspace(min, max, 2e1), colors='k', linewidth=0.5)
        # ax.contour(var2.T, colors='k', linewidths=0.5)
        # ax.set_title('u and v')

        # print(coordinates[it,:,0])
        print(np.amin(coordinates[it,:,0]), np.amax(coordinates[it,:,0]))
        print(np.amin(coordinates[it,:,1]), np.amax(coordinates[it,:,1]))
        for i in range(n_tracers):
            ax0.plot(coordinates[it, i, 0] + x_shift, coordinates[it, i, 1] + x_shift, 'ok', markersize=2)
            ax1.plot(coordinates[it, i, 0] + x_shift, coordinates[it, i, 1] + x_shift, 'ok', markersize=2)
        for ax in [ax0, ax1]:
            ax.set_aspect('equal')
            # ax.set_xlim(np.amin(coordinates[it, :, 0] - 1), np.amax(coordinates[it, :, 0] + 1))
            # ax.set_ylim(np.amin(coordinates[it, :, 1] - 1), np.amax(coordinates[it, :, 1] + 1))
            ax.set_xlim(630, 670)
            ax.set_ylim(580, 620)
            ax.plot(xc[it],yc[it], 'xk', markersize=20)
        plt.suptitle('t='+str(t0)+'s, t-shift '+str(t_shift)+', x-shift '+str(x_shift))
        plt.tight_layout()
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)

    return



# ----------------------------------------------------------------------


def get_tracer_coords(cp_id, n_cps, n_tracers, times, fullpath_in):
    print('get tracer coordinates')
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    column = lines[0].split()

    nt = len(times)
    coords = np.zeros((nt, n_tracers, 2))

    count_start = 0
    ID = int(lines[count_start].split()[3])
    while (ID < cp_id):
        count_start += 1
        ID = int(lines[count_start].split()[3])
    print('count start', count_start)

    for it, t0 in enumerate(times):
        print('----t0='+str(t0), it, '----')
        i = 0
        count = count_start + it * n_tracers
        # count = it * n_cps * n_tracers + (cp_id - 1) * n_tracers
        print('count: ', count)
        print(int(lines[count].split()[3]),
              int(lines[count-1].split()[0]), int(lines[count].split()[0]),
              int(lines[count].split()[2]))
        # count = t0/dt_fields * n_cps * n_tracers + (cp_id - 1) * n_tracers
        # while CP age is 0 and CP ID is cp_id
        timestep = int(lines[count].split()[0])
        print('timestep: ', timestep, t0)
        # cp_ID = int(lines[count].split()[3])
        tracer_ID = int(lines[count].split()[2])
        # while (timestep - 1 == it and cp_ID == cp_id):
        # while (timestep - 1 == t0/dt_fields and cp_ID == cp_id):
        while (timestep == t0 and tracer_ID <= n_tracers):
            columns = lines[count].split()
            coords[it,i,0] = float(columns[4])
            coords[it,i,1] = float(columns[5])
            i += 1
            count += 1
            tracer_ID = int(lines[count].split()[2])
            timestep = int(lines[count].split()[0])
        print('i', i)
    f.close()
    print ''
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
    cp_age = int(lines[count].split()[0])
    cp_ID = int(lines[count].split()[3])
    # print('cp_age', cp_age)
    # print('cp_ID', cp_ID)
    age = cp_age
    ID = cp_ID
    while (age == cp_age and ID == cp_ID):
        count += 1
        age = int(lines[count].split()[0])
        ID = int(lines[count].split()[3])
    n_tracers = count
    f.close()

    return n_tracers



def get_cp_lifetime(cp_ID, fullpath_in):
    f = open(fullpath_in, 'r')
    lines = f.readlines()
    count = 0
    ID = int(lines[count].split()[3])
    while (ID < cp_ID):
        count += 1
        ID = int(lines[count].split()[3])
    t0 = int(lines[count].split()[0])
    ID = int(lines[count].split()[3])

    while (ID == cp_ID):
        count += 1
        ID = int(lines[count].split()[3])
    t1 = int(lines[count-1].split()[0])
    tau = t1-t0+1


    f.close()
    print('lifetime: tau='+str(tau), t0, t1)
    return tau, t0, t1



def get_cp_center(cp_ID, tau, n_tracers, fullpath_in):
    xc = np.zeros(tau)
    yc = np.zeros(tau)

    f = open(fullpath_in, 'r')
    lines = f.readlines()
    count = 0
    ID = int(lines[count].split()[3])
    while (ID < cp_ID):
        count += 1
        ID = int(lines[count].split()[3])
    count_start = count

    for it in range(tau):
        xc[it] = float(lines[count_start+it*n_tracers].split()[14])
        yc[it] = float(lines[count_start+it*n_tracers].split()[15])

    f.close()
    print('xc: ', xc)
    print('yc: ', yc)
    return xc, yc
# ----------------------------------------------------------------------

def set_input_parameters(args):
    print('--- set input parameters ---')
    global path, path_fields, case_name
    path = args.path
    path_fields = os.path.join(path, 'fields')
    case_name = args.casename

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
    global times, files
    times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
             and tmin <= np.int(name[:-3]) <= tmax]
    times.sort()
    files = [str(t) + '.nc' for t in times]
    print('tmin, tmax: ', tmin, tmax)
    print('times: ', times)
    print('')

    return

# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()
