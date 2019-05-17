import numpy as np
import netCDF4 as nc
import argparse
import os
import matplotlib.pyplot as plt
import json as simplejson




def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()
    set_input_output_parameters(args)

    t0 = 200
    k0 = 1
    rootgrp_in = nc.Dataset(os.path.join(path_fields, str(t0)+'.nc'))
    s = rootgrp_in.groups['fields'].variables['s'][:,:,k0]

    s_fft = np.fft.fft2(s)
    s_fft_real = s_fft.real
    print('s_fft')
    print type(s_fft), type(s_fft[0,0]), s_fft[0,0].real, s_fft[0,0].imag
    print 'fft(s): ', np.amax(s_fft.real), np.amin(s_fft.real)

    max = np.amax(s_fft_real)
    max_list = []
    max_loc = []
    max_loc_2d = []
    max_loc_x = []
    max_loc_y = []
    s_fft_real_cp = np.array(s_fft_real, copy=True)
    while (max > 1e3):
        print 'max: ', max
        max_list.append(max)
        max_loc.append(np.argmax(s_fft_real_cp))
        max_loc_2d.append(np.unravel_index(np.argmax(s_fft_real_cp, axis=None), s_fft_real_cp.shape))
        max_loc_x.append(np.unravel_index(np.argmax(s_fft_real_cp, axis=None), s_fft_real_cp.shape)[0])
        max_loc_y.append(np.unravel_index(np.argmax(s_fft_real_cp, axis=None), s_fft_real_cp.shape)[1])
        s_fft_real_cp = np.ma.masked_greater_equal(s_fft_real_cp, max)
        max = np.amax(s_fft_real_cp)
    print ' '
    print max_loc_2d
    print max_loc_x
    max_x = np.amax([max for max in max_loc_x if max < nx/2])
    max_y = np.amax([max for max in max_loc_y if max < ny/2])
    print 'max x, y: ', max_x, max_y
    # times1 = [np.int(name[:-3]) for name in os.listdir(path_fields1) if name[-2:] == 'nc'
    #           and np.int(name[:-3]) >= tmin and np.int(name[:-3]) <= tmax]

    max_loc_mask = np.zeros(shape=s_fft_real.shape)
    for [i,j] in max_loc_2d:
        max_loc_mask[i,j] = 1



    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    ax1 = axes[0]
    ax1.imshow(s.T, origin='lower')
    ax2 = axes[1]
    im = ax2.imshow(s_fft[:10, :10].real, origin='lower')
    plt.colorbar(im, ax=ax2, shrink=0.5)
    ax2.set_title('real')
    ax3 = axes[2]
    im = ax3.imshow(s_fft.imag.T, origin='lower')
    plt.colorbar(im, ax=ax3, shrink=0.5)

    ax4 = axes[3]
    im = ax4.imshow(max_loc_mask.T, origin='lower', cmap=plt.cm.get_cmap('Reds'))
    plt.colorbar(im, ax=ax4, shrink=0.5)
    plt.savefig('./fft_s_t'+str(t0)+'s.png')
    plt.close()

    return


# ----------------------------------------------------------------------

def set_input_output_parameters(args):
    path = args.path

    global path_fields
    path_fields = os.path.join(path, 'fields')

    global case_name
    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dy, dz, dV, gw
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = np.zeros(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx[0] * dx[1] * dx[2]

    return


# ----------------------------------

if __name__ == '__main__':
    main()