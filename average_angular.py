import numpy as np
import matplotlib.pyplot as plt
#import netCDF4 as nc
import argparse
import json as simplejson
import os
import sys

def main():
    global nx, ny
    nx = 20
    ny = 20
    r_field = np.zeros((nx, ny), dtype=np.int)

    ic = np.int(nx/2)
    jc = np.int(ny/2)
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

    # t0 = 500
    # fullpath_in = os.path.join(path, fields, str(t0) + '.nc')
    # rootgrp = nc.Dataset(fullpath_in, 'r')
    # for var_name in var_list:
    #     var = rootgrp.groups['fields'].variables[var_name]
    #     var_av = compute_average_var(var, rmax, r_field)
    # rootgrp.close()
    # var_list = ['w', 's', 'phi']
    # data_dict = read_in_vars(fullpath_in, var_list)
    # for var_name in var_list:
    #     var_av = compute_average_var(data_dict[var_name], rmax, r_field)
    
    var1 = np.ones((nx, ny))  # should be from 3D LES fields, read in
    var1_av = compute_average_var(var1, rmax, r_field)

    return


# _______________________________

def read_in_vars(fullpath_in, var_list):

    var_dict = {}

    rootgrp = nc.Dataset(fullpath_in, 'r')
    for var_name in var_list:
        var = rootgrp.groups['fields'].variables[var_name]
        data = var[:, :, :]
        # var_dict[var_name] = data

    rootgrp.close()

    return data
    # return var_dict


# _______________________________

def compute_average_var(var1, rmax, r_field):

    var1_av = np.zeros(rmax, dtype=np.double)
    count = np.zeros(rmax, dtype=np.int)
    for i in range(nx):
        for j in range(ny):
            r = r_field[i, j]
            count[r] += 1
            var1_av[r] += var1[i, j]
    print count

    for r in range(rmax):
        if count[r] > 0:
            var1_av[r] /= count[r]

    return var1_av
# _______________________________


if __name__ == '__main__':
    main()