import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import json as simplejson
import os


from define_cp_rim_plottingfct import set_colorbars
from define_cp_rim_plottingfct import plot_yz_crosssection, plot_w_field, plot_s, \
    plot_outlines, plot_rim_mask, plot_angles, plot_cp_outline_alltimes, \
    plot_cp_rim_velocity, plot_cp_rim_averages, plot_rim_thickness

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--casename")
    parser.add_argument("--path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    global path_fields, path_out
    if args.path:
        path = args.path
    else:
        path = '/Users/bettinameyer/polybox/ClimatePhysics/Copenhagen/Projects/LES_ColdPool/' \
               'triple_3D_noise/Out_CPDry_triple_dTh2K/'
        # path = '/nbi/ac/cond1/meyerbe/ColdPools/triple_3D_noise/Out_CPDry_triple_Th3K/'
    if os.path.exists(os.path.join(path, 'fields')):
        path_fields = os.path.join(path, 'fields')
    elif os.path.exists(os.path.join(path, 'fields_k120')):
        path_fields = os.path.join(path, 'fields_k120')
    path_out = os.path.join(path, 'figs_cp_rim')
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    global case_name
    if args.casename:
        case_name = args.casename
    else:
        case_name = 'ColdPoolDry_triple_3D'

    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 400
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = tmin
    timerange = np.arange(tmin, tmax + 100, 100)
    nt = len(timerange)

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dy, dz
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    gw = nml['grid']['gw']

    global cm_bwr, cm_grey, cm_vir
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_vir = plt.cm.get_cmap('viridis')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    set_colorbars(cm_bwr, cm_vir, cm_grey)  # to set colorbars as global functions in define_cp_rim_plottingfct.py



    k0 = 5
    perc = 95
    fullpath_in = os.path.join(path_out, 'rimmask_perc'+str(perc)+'th_t'+str(tmin)+'.nc')

    ''' read in mask file '''
    mask = read_in_netcdf_file('mask', 'fields', fullpath_in)[:,:,:]
    plt.figure()
    plt.imshow(mask[:,:,k0].T, origin='lower')
    plt.savefig(os.path.join(path_out, 'test_read_in_mask.png'))
    plt.close()
    del mask

    rim_int = read_in_netcdf_file('rim_inner', 'fields', fullpath_in)
    rim_out = read_in_netcdf_file('rim_outer', 'fields', fullpath_in)
    plt.figure()
    plt.subplot(1,2,1)
    plt.imshow(rim_int[:,:,k0].T, origin='lower')
    plt.title('inner rim')
    plt.subplot(1,2,2)
    plt.imshow(rim_out[:,:,k0].T, origin='lower')
    plt.title('outer rim')
    plt.savefig(os.path.join(path_out, 'test_read_in_rim.png'))
    plt.close()

    return



def read_in_netcdf_file(variable_name, group_name, fullpath_in):
    print(fullpath_in)
    print('')
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups[group_name].variables[variable_name]
    if group_name == 'fields':
        data = var[:, :, :]
        rootgrp.close()
        return data
    else:
        rootgrp.close()
        return

if __name__ == '__main__':
    main()