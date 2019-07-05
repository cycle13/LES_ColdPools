import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os
import scipy
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

label_size = 10
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.labelsize'] = 15

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--path_tracers")
    args = parser.parse_args()
    set_input_parameters(args)
    nt = len(times)

    dt_fields = 100


    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')


    path_out_figs = os.path.join(path, 'figs_tracers')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    if args.path_tracers:
        path_tracers = os.path.join(path, args.path_tracers)
        path_tracer_file = os.path.join(path, args.path_tracers, 'output')
    else:
        k_tracers = 0
        path_tracers = os.path.join(path, 'tracer_k' + str(k_tracers))
        path_tracer_file = os.path.join(path, 'tracer_k' + str(k_tracers), 'output')

    print('')
    print('path figs:    ' + path_out_figs)
    print('path tracers: ' + path_tracer_file)
    print('')

    ''' read in coordinates '''
    cp_id = 1
    n_cps = get_number_cps(path_tracer_file)
    n_tracers = get_number_tracers(path_tracer_file)
    print('number of CPs: ' + str(n_cps))
    print('CP ID: ' + str(cp_id))
    coordinates_xy, coordinates_pol = get_tracer_coords(cp_id, n_cps, n_tracers, times, dt_fields, path_tracer_file)
    dist_mean = np.mean(coordinates_pol[:,:,0], axis=1)
    print(coordinates_xy.shape, coordinates_pol.shape)
    print('')


    ''' compute polar coordinates separately '''
    shift = 0
    ic = nx/2
    jc = ny/2
    shifts = np.arange(-1,2)
    coordinates_pol_new = np.zeros((len(shifts), nt, n_tracers, 2))
    coordinates_pol_new_sorted = np.zeros((nt, n_tracers, 2))
    it = 0
    for it,t0 in enumerate(times):
        print('t0: '+ str(t0))
        for si,s in enumerate(shifts):
            for i in range(n_tracers):
                # DELTAx = coordinates_xy[i, 0] - dsize_x / 2. !(COMx(i) + dx(i)) - 1.
                # DELTAy = coordinates_xy[i, 1] - dsize_y / 2. !(COMy(i) + dy(i)) - 1 !dsize_y / 2.
                DELTAx = coordinates_xy[it, i, 0] - ic + s
                DELTAy = coordinates_xy[it, i, 1] - jc + s
                # coordinates_pol_new[]
                if (DELTAx > 0 and DELTAy > 0):
                    phi = np.arctan(DELTAy / DELTAx)
                elif (DELTAx < 0 and DELTAy > 0):
                    phi = np.arctan(abs(DELTAx) / DELTAy) + 1. / 2. * np.pi
                elif (DELTAx < 0 and DELTAy < 0):
                    phi = np.arctan(DELTAy / DELTAx) + np.pi
                elif (DELTAx > 0. and DELTAy < 0.):
                    phi = np.arctan(abs(DELTAx) / abs(DELTAy)) + 3. / 2. * np.pi
                elif (DELTAx > 0. and DELTAy == 0.):
                    phi = 0.
                elif (DELTAx < 0. and DELTAy == 0):
                    phi = np.pi
                elif (DELTAx == 0. and DELTAy > 0):
                    phi = 1. / 2. * np.pi
                elif (DELTAx == 0. and DELTAy < 0):
                    phi = 3. / 2. * np.pi

                coordinates_pol_new[si,it,i,1] = phi
                coordinates_pol_new[si,it,i,0] = np.sqrt(DELTAy ** 2. + DELTAx ** 2.)
        corr = np.amax(coordinates_pol_new[0,it,:20,1])
        print('corr', corr)
        coordinates_pol_new_sorted[it,:,0] = coordinates_pol_new[0,it,:,0]
        coordinates_pol_new_sorted[it,:,1] = np.mod(coordinates_pol_new[0,it,:,1] - corr, 2*np.pi)
            # print('...')
            # print(coordinates_pol_new[si,it,:10,1])
            # print(coordinates_pol_new[si,it,:10,1])
    # r = np.sqrt( (coordinates_xy[:,:,0] - ic)**2*dx[0]**2 + (coordinates_xy[:,:,1] - jc)**2*dx[1]**2 )
    # print('r', r.shape)

    shift = 0
    for it,t0 in enumerate(times):
        if it >= 2:
            shift = 1
        print('-plot time: '+str(t0))
        fig_name = 'coordinates' + '_t' + str(t0) + '_tracers.png'
        fig, axis = plt.subplots(1, 4, figsize=(24, 6), sharey='none')
        axis[0].set_aspect('equal')
        radius = np.int(np.double(os.path.basename(path[:])[12:]) / dx[0]) + 5
        # for s in np.arange(0,1):
        s = 1
        circle1 = plt.Circle((ic, jc), radius + 2*s*(it-1)**2, fill=False, color='k', linewidth=7)
        circle2 = plt.Circle((ic+s, jc+s), radius + 2*s*(it-1)**2, fill=False, color='gray', linewidth=7)
        axis[0].add_artist(circle1)
        axis[0].add_artist(circle2)
        axis[0].plot(ic,jc, 'ko', markersize=10)
        for i in range(n_tracers):
            count_col = np.double(i)/n_tracers
            axis[0].plot(coordinates_xy[it,i,0], coordinates_xy[it,i,1],
                         'o', color=cm_hsv(count_col), markersize=3, markeredgecolor='w')
            axis[1].plot([0,2*np.pi],[dist_mean[it],dist_mean[it]], 'k', linewidth=0.5)
            axis[1].plot(coordinates_pol[it,i,1], coordinates_pol[it,i,0],
                         'o', color=cm_hsv(count_col), markersize=4, markeredgecolor='w')
            axis[1].plot(coordinates_pol_new[1,it,i,1], coordinates_pol_new[1,it,i,0],
                         'ko', markersize=2, markeredgecolor='w')
            axis[2].plot([0, 2 * np.pi], [dist_mean[it], dist_mean[it]], 'k', linewidth=0.5)
            axis[2].plot(coordinates_pol_new[0,it, i, 1], coordinates_pol_new[0,it, i, 0],
                         'o', color=cm_hsv(count_col), markersize=4, markeredgecolor='w')
            axis[3].plot([0, 2 * np.pi], [dist_mean[it], dist_mean[it]], 'k', linewidth=0.5)
            # axis[3].plot(coordinates_pol_new[0, it, i, 1], coordinates_pol_new[0, it, i, 0],
            #              'o', color=cm_hsv(count_col), markersize=4, markeredgecolor='w')
        axis[3].plot(coordinates_pol_new_sorted[it, :, 1], coordinates_pol_new_sorted[it, :, 0],
                         'o-', color='gray', linewidth=1, markersize=2)
        # axis[3].plot(coordinates_pol_new_sorted[it, :, 1],
        #              'o-', color='gray', linewidth=1, markersize=2)

        # angles_ = coordinates_pol[it,:,1]
        # radii_ = coordinates_pol[it,:,0]
        # linear_interp = interp1d(angles_, radii_)
        # angles = np.linspace(1, 6, 50)
        # # print(coordinates_pol_new_sorted[it,:10,1])
        # print(np.amin(angles), np.amax(angles), np.amin(angles_), np.amax(angles_))
        # linear_results = linear_interp(angles)

        axis[0].set_xlabel('x')
        axis[0].set_title('circle centere shifted by '+str(shift))
        axis[0].set_ylabel('y')
        axis[1].set_title('output tracers')
        axis[1].set_xlabel('angle')
        axis[1].set_ylabel('radius')
        axis[1].set_xlim(0, 2*np.pi)
        axis[2].set_title('tracers shifted by '+str(shifts[0]))
        axis[2].set_xlabel('angle')
        axis[2].set_ylabel('radius')
        axis[2].set_xlim(0, 2*np.pi)
        axis[3].set_xlim(0, 2*np.pi)
        plt.tight_layout()
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)

    # scipy
    # # - finding the roots of a scalar function (1.5.5.3: scipy.optimize.root())
    # # - interpolate between tracers that are sorted by increasing angle: scipy.interpolate.interp1d


    return
# ----------------------------------------------------------------------


def get_tracer_coords(cp_id, n_cps, n_tracers, times, dt_fields, fullpath_in):
    print('get tracer coordinates')
    f = open(fullpath_in + '/coldpool_tracer_out.txt', 'r')
    lines = f.readlines()
    column = lines[0].split()

    nt = len(times)
    coords = np.zeros((nt, n_tracers, 2))
    coords_pol = np.zeros((nt, n_tracers, 2))

    for it, t0 in enumerate(times):
        print('----t0='+str(t0), it, '----')
        i = 0
        # count = t0 * n_cps * n_tracers + (cp_id - 1) * n_tracers
        # count = it * n_cps * n_tracers + (cp_id - 1) * n_tracers
        count = t0/dt_fields * n_cps * n_tracers + (cp_id - 1) * n_tracers
        # while CP age is 0 and CP ID is cp_id
        timestep = int(lines[count].split()[0])
        cp_ID = int(lines[count].split()[3])
        while (timestep - 1 == t0/dt_fields and cp_ID == cp_id):
            columns = lines[count].split()
            coords[it,i,0] = float(columns[4])
            coords[it,i,1] = float(columns[5])
            coords_pol[it,i,0] = float(columns[8])
            coords_pol[it,i,1] = float(columns[9])
            i += 1
            count += 1
            cp_ID = int(lines[count].split()[3])
            timestep = int(lines[count].split()[0])

    f.close()
    # print ''
    return coords, coords_pol


def get_number_cps(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in + '/coldpool_tracer_out.txt', 'r')
    lines = f.readlines()
    # count = 0
    # # while CP age is 0 and CP ID is 1
    # while (int(lines[count].split()[0]) == 1):
    #     count += 1
    # cp_number = int(lines[count-1].split()[3])
    cp_number = int(lines[-1].split()[3])

    f.close()

    return cp_number


def get_number_tracers(fullpath_in):
    # get number of tracers in each CP
    f = open(fullpath_in + '/coldpool_tracer_out.txt', 'r')
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