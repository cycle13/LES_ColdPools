import netCDF4 as nc
import argparse
import json as simplejson
import os
import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import matplotlib.cbook as cbook
import time

execfile('settings.py')
label_size = 15
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['axes.labelsize'] = 18

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("--level")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    case_name_1CP = 'ColdPoolDry_single_3D'
    case_name_2CP = 'ColdPoolDry_double_3D'
    case_name_3CP = 'ColdPoolDry_triple_3D'

    dTh = 5
    zstar = 1000
    rstar = 2000
    rst = str(rstar)
    sep = d_range[rst][0]
    print('Parameters: ')
    print('d-range: ' + str(d_range[rst]))
    print('')

    case = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
    case_xCP = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar) + '_d' + str(sep) + 'km'
    path_1CP = '/nbi/ac/conv4/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/single_3D_noise/run5_PE_scaling_dx100m/'
    path_2CP = '/nbi/ac/coag/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/double_3D/'
    path_3CP = '/nbi/ac/cond2/meyerbe/ColdPools/3D_sfc_fluxes_off/triple_3D/'
    print('')
    print('Case: ' + case)
    print('path 1CP:   ' + path_1CP)
    print('path 2CP:   ' + path_2CP)
    print('path 3CP:   ' + path_3CP)
    path_data = os.path.join(path_3CP, 'data_analysis')
    # path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single'
    path_out_figs = '/nbi/home/meyerbe/paper_CP'
    print('Path Figures: ' + path_out_figs)
    print('')

    if args.level:
        zlevel = np.int(args.level)
    else:
        zlevel = 1000
    print('')

    # plotting parameters
    xmin_3CP = np.zeros(3, dtype=np.int)
    xmax_3CP = np.zeros(3, dtype=np.int)
    ymin_3CP = np.zeros(3, dtype=np.int)
    ymax_3CP = np.zeros(3, dtype=np.int)
    if rstar == 1100:
        xmin_3CP[0] = 0
        xmax_3CP[0] = 400
        ymin_3CP[0] = 0
        xmin_3CP[1] = 0
        xmax_3CP[1] = 400
        ymin_3CP[1] = 0
        xmin_3CP[2] = 0
        xmax_3CP[2] = 400
        ymin_3CP[2] = 0


    for i, sep in enumerate(d_range[rst][:]):
        print('--- d=' + str(sep) + 'km ---')
        case_xCP = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar) + '_d' + str(sep) + 'km'

        nx_1CP, nx_2CP, nx_3CP, dt_fields_1CP, dt_fields_2CP, dt_fields_3CP = set_input_output_parameters(args, case_name_1CP, case_name_2CP, case_name_3CP,
                                                                     path_1CP, path_2CP, path_3CP,
                                                                     case, case_xCP)

        ic_arr_1CP, jc_arr_1CP, ic_arr_2CP, jc_arr_2CP, \
                ic_arr_3CP, jc_arr_3CP, \
                ic_2CP, jc_2CP, ic_3CP, jc_3CP = define_geometry(nx_1CP, nx_2CP, nx_3CP, sep)


        ''' Mass Flux '''
        filename_data = 'mass_flux_z' + str(zlevel) + '.nc'
        path_in = os.path.join(path_1CP, case, 'data_analysis', filename_data)
        print(path_in)
        root_in = nc.Dataset(path_in, 'r')
        mass_flux_1CP = root_in.groups['fields_2D'].variables['mass_flux_2D'][:, :, :]
        mass_flux_pos_1CP = root_in.groups['fields_2D'].variables['mass_flux_2D_positive'][:, :, :]
        root_in.close()
        # path_in = os.path.join(path_2CP, case_xCP, 'data_analysis', filename_data)
        # print(path_in)
        # root_in = nc.Dataset(path_in, 'r')
        # mass_flux_2CP = root_in.groups['fields_2D'].variables['mass_flux_2D'][:, :, :]
        # mass_flux_pos_2CP = root_in.groups['fields_2D'].variables['mass_flux_2D_positive'][:, :, :]
        # time_data_2CP = root_in.groups['timeseries'].variables['time'][:]
        # root_in.close()
        # path_in = os.path.join(path_3CP, case_xCP, 'data_analysis', filename_data)
        # print(path_in)
        # root_in = nc.Dataset(path_in, 'r')
        # mass_flux_3CP = root_in.groups['fields_2D'].variables['mass_flux_2D'][:, :, :]
        # mass_flux_pos_3CP = root_in.groups['fields_2D'].variables['mass_flux_2D_positive'][:, :, :]
        # time_data_3CP = root_in.groups['timeseries'].variables['time'][:]
        # root_in.close()
    return



# _______________________________________________________

def set_input_output_parameters(args, case_name_1CP, case_name_2CP, case_name_3CP,
                                path_1CP, path_2CP, path_3CP,
                                case, case_xCP):
    print('')
    print('--- set input parameters ---')
    nml_1CP = simplejson.loads(open(os.path.join(path_1CP, case, case_name_1CP + '.in')).read())
    nml_2CP = simplejson.loads(open(os.path.join(path_2CP, case_xCP, case_name_2CP + '.in')).read())
    nml_3CP = simplejson.loads(open(os.path.join(path_3CP, case_xCP, case_name_3CP + '.in')).read())
    nx_1CP = np.zeros(3, dtype=np.int)
    nx_2CP = np.zeros(3, dtype=np.int)
    nx_3CP = np.zeros(3, dtype=np.int)
    nx_1CP[0] = nml_1CP['grid']['nx']
    nx_1CP[1] = nml_1CP['grid']['ny']
    nx_1CP[2] = nml_1CP['grid']['nz']
    nx_2CP[0] = nml_2CP['grid']['nx']
    nx_2CP[1] = nml_2CP['grid']['ny']
    nx_2CP[2] = nml_2CP['grid']['nz']
    nx_3CP[0] = nml_3CP['grid']['nx']
    nx_3CP[1] = nml_3CP['grid']['ny']
    nx_3CP[2] = nml_3CP['grid']['nz']
    global dx
    dx = np.zeros(3, dtype=np.int)
    dx[0] = nml_3CP['grid']['dx']
    dx[1] = nml_3CP['grid']['dy']
    dx[2] = nml_3CP['grid']['dz']
    print('grid 1CP: ', nx_1CP)
    print('grid 2CP: ', nx_2CP)
    print('grid 3CP: ', nx_3CP)

    dt_fields_1CP = np.int(nml_1CP['fields_io']['frequency'])
    dt_fields_2CP = np.int(nml_2CP['fields_io']['frequency'])
    dt_fields_3CP = np.int(nml_3CP['fields_io']['frequency'])

    return nx_1CP, nx_2CP, nx_3CP, dt_fields_1CP, dt_fields_2CP, dt_fields_3CP


def define_geometry(nx_1CP, nx_2CP, nx_3CP, sep):
    # 1-CP simulation
    ic_arr_1CP = [np.int(nx_1CP[0]/2)]
    jc_arr_1CP = [np.int(nx_1CP[1]/2)]

    # 2-CP simulation
    isep = np.int(np.round(sep / dx[0]))
    jsep = 0
    ic_2CP = np.int(np.round(nx_2CP[0] / 2))
    jc_2CP = np.int(np.round(nx_2CP[1] / 2))
    ic1 = ic_2CP - np.int(np.round(isep / 2))
    jc1 = jc_2CP
    ic2 = ic1 + isep
    jc2 = jc1 + jsep
    ic_arr_2CP = [ic1, ic2]
    jc_arr_2CP = [jc1, jc2]

    # 3-CP simulation
    i_d = np.int(np.round(sep / dx[0]))
    idhalf = np.int(np.round(i_d / 2))
    a = np.int(np.round(i_d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
    r_int = np.int(np.sqrt(3.) / 6 * i_d)  # radius of inscribed circle
    # point of 3-CP collision (ic, jc)
    ic_3CP = np.int(np.round(nx_3CP[0] / 2))
    jc_3CP = np.int(np.round(nx_3CP[1] / 2))
    ic1 = ic_3CP - r_int
    ic2 = ic1
    ic3 = ic_3CP + (a - r_int)
    jc1 = jc_3CP - idhalf
    jc2 = jc_3CP + idhalf
    jc3 = jc_3CP
    ic_arr_3CP = [ic1, ic2, ic3]
    jc_arr_3CP = [jc1, jc2, jc3]
    print('')
    print('CP configuration: ')
    print('1CP: ', ic_arr_1CP, jc_arr_1CP)
    print('2CP: ', ic_arr_2CP, jc_arr_2CP)
    print('3CP: ', ic_arr_3CP, jc_arr_3CP)
    print('CP collision points: ')
    print('2CP: ', ic_2CP, jc_2CP)
    print('3CP: ', ic_3CP, jc_3CP)
    print('')

    return ic_arr_1CP, jc_arr_1CP, ic_arr_2CP, jc_arr_2CP, ic_arr_3CP, jc_arr_3CP, ic_2CP, jc_2CP, ic_3CP, jc_3CP

# _______________________________________________________

if __name__ == '__main__':
    main()
