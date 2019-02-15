import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35
# plt.rcParams['savefig.edgecolor'] = 'white'


def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    # parser.add_argument("--tmin")
    # parser.add_argument("--tmax")
    args = parser.parse_args()

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')

    nml = set_input_parameters(args)
    # nt = len(times)

    file_name = 'field_vort_yz.nc'
    rootgrp = nc.Dataset(os.path.join(path_in, file_name), 'r')
    grp_descr =rootgrp.groups['description'].variables
    ic_arr = grp_descr['ic_arr'][:]
    jc_arr = grp_descr['jc_arr'][:]

    grp_fields = rootgrp.groups['fields']
    print grp_fields.dimensions.keys()
    nt = grp_fields.dimensions['time'].size
    nx = grp_fields.dimensions['nx'].size
    ny = grp_fields.dimensions['ny'].size
    nz = grp_fields.dimensions['nz'].size
    print('shape fields: ', nx, ny, nz)
    vort_yz = grp_fields.variables['vort_yz'][:,:,:]
    time_range = grp_fields.variables['time'][:]
    rootgrp.close()


    it = 0
    t0 = 100

    for it, t0 in enumerate(time_range):
        # location of maximum / minimum
        vort_max = np.amax(vort_yz[it,:,:])
        vort_min = np.amin(vort_yz[it,:,:])
        vort_max_i = np.unravel_index(np.argmax(vort_yz[it,:,:]), (ny,nz))
        vort_min_i = np.unravel_index(np.argmin(vort_yz[it,:,:]), (ny,nz))
        # print vort_max, vort_yz[it, vort_max_i]
        # print('location maximum: ', vort_max_i)

        jmin = 50
        jmax = 350
        kmax = 20

        # computing circulation
        # !!!!!!!! need to define limits of search area?! (dynamically as it moves with time, maybe read in radius)
        circ = np.zeros(nt)
        vort_c = 2e-2
        test_field = np.zeros((jmax-jmin, kmax), dtype = np.int)
        for j in range(jmin, jmax):
            for k in range(kmax):
                if vort_yz[it, j,k] >= vort_c:
                    circ[it] += vort_yz[it, j,k]
                    test_field[j-jmin, k] = 1

        # plotting
        lvls = np.linspace(vort_min, vort_max, 1e3)
        fig, axes = plt.subplots(1, 2, figsize=(15, 4))
        ax = axes[0]
        ax.contourf(vort_yz[it, jmin:jmax, :kmax].T, levels=lvls)
        ax.plot(vort_max_i[0] - jmin, vort_max_i[1], 'kx', markersize=20)
        # ax.contourf(vort_yz[it, :, :kmax].T, levels=lvls)
        # ax.plot(vort_max_i[0], vort_max_i[1], 'kx', markersize=20)
        # ax.plot(vort_min_i[0]- jmin, vort_min_i[1], 'kx', markersize=20)
        ax = axes[1]
        ax.imshow(test_field[vort_max_i[0]-50-jmin:vort_max_i[0]+25-jmin,:].T, origin='lower')
        ax.set_title('circulation: C='+str(circ[it]))
        # for ax in axes:
        #     ax.legend(loc='best')
        plt.tight_layout
        plt.suptitle('t='+str(t0)+'s')
        fig_name = 'vort_max_t' + str(np.int(t0)) + '.png'
        plt.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)

    return

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

def set_input_parameters(args):
    print('--- set input parameters ---')
    global case_name, times
    global path, path_in, path_fields
    global path_out_figs, path_out_data

    path = args.path
    path_in = os.path.join(path, 'fields_vorticity')
    path_fields = os.path.join(path, 'fields')
    path_out_figs = os.path.join(path, 'figs_vorticity')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)
    path_out_data = os.path.join(path, 'fields_vorticity')
    if not os.path.exists(path_out_data):
        os.mkdir(path_out_data)

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

    ''' determine time range '''
    # if args.tmin:
    #     tmin = np.int(args.tmin)
    # else:
    #     tmin = 100
    # if args.tmax:
    #     tmax = np.int(args.tmax)
    # else:
    #     tmax = tmin
    # times = np.arange(tmin, tmax + 100, 100)
    # nt = len(times)
    # print('timerange', times)
    #
    # return nml, times
    return nml

# ----------------------------------
# ----------------------------------

if __name__ == '__main__':
    main()