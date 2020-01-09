import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse
import scipy
from scipy import stats


execfile('settings.py')

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    case_name = 'ColdPoolDry_single_3D'
    path_root = '/nbi/ac/conv4/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/single_3D_noise/run5_PE_scaling_dx100m/'
    path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single/figs_run5'
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)

    # loop through all cases




    return


# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()