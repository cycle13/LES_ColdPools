import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import netCDF4 as nc
import argparse
import json as simplejson
import os

# (a) vertical velocity is zero from divergence to convergence zone (behind front)
#       >> velocity should be zero when averaged over a larger area
# (b) radially average


def main():
    path = '/Users/bettinameyer/polybox/ClimatePhysics/Copenhagen/Projects/LES_ColdPool/' \
           'triple_3D_noise/Out_CPDry_triple_dTh2K/'
    path_fields = os.path.join(path, 'fields')
    path_out = os.path.join(path, 'figs_cp_rim')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    print('path out: ', path_out)

    global case_name
    case_name = 'ColdPoolDry_triple_3D'
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dy, dz
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    gw = nml['grid']['gw']

    # define subdomain to scan
    # --- for triple coldpool ---
    d = np.int(np.round(ny / 2)) + gw
    a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
    rstar = 5000.0  # half of the width of initial cold-pools [m]
    irstar = np.int(np.round(rstar / dx))
    ic = np.int(np.round(a / 2))
    jc = np.int(np.round(d / 2))
    ishift = 20
    jshift = ishift
    id = irstar + ishift
    jd = id
    nx_ = 2 * id
    ny_ = 2 * jd
    print(ic,jc,id,jd)

    # t0 = 200
    cm = plt.get_cmap('bwr')
    for t0 in np.arange(100,500,100):
        print('time: '+ str(t0))

        k0 = 2

        w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0) + '.nc'))
        w_roll = np.roll(w[:,:,k0],[ishift,jshift],[0,1])
        w_ = w_roll[ic-id+ishift:ic+id+ishift, jc-jd+jshift:jc+jd+jshift]

        # ic + ishift + id  = ic + id - ic + id = 2*id
        # ic + ishift - id = ic + id - ic - id = 0

        cm_bwr = plt.cm.get_cmap('bwr')
        max = np.amax(w)
        plt.figure(figsize=(20,6))
        ax1 = plt.subplot(1,5,1)
        ax1.set_title('w')
        # plt.contourf(w[:, :, k0].T, cmap=cm_bwr, aspect=1.)
        # plt.colorbar()
        ax2 = plt.subplot(1,5,2)
        ax2.set_title('np.roll(w)')
        ax3 = plt.subplot(1,5,3)
        ax3.set_title('np.roll(w)[:ic,:jc,k0]')
        ax4 = plt.subplot(1,5,4)
        ax4.set_title('masked array')
        ax5 = plt.subplot(1,5,5)
        ax5.set_title('mask')
        ax1.imshow(w[:,:,k0].T, cmap=cm_bwr, origin='lower', vmin=-max, vmax=max)
        ax1.plot([ic,ic],[0,ny-1])
        ax1.plot([0,nx-1],[jc,jc])
        ax2.imshow(w_roll.T, cmap=cm_bwr, origin='lower', vmin=-max,vmax=max)
        ax2.plot([ic+ishift,ic+ishift],[0,ny-1])
        ax2.plot([0,nx-1],[jc+jshift,jc+jshift])
        ax3.imshow(w_.T, cmap=cm_bwr, origin='lower', vmin=-max,vmax=max)


        w_c = 5e-1
        # w_mask = np.ma.masked_less(w_, w_c)
        w_mask = np.ma.masked_less(w_, w_c)
        print('...', w.shape, w_.shape, w_roll.shape, w_mask.mask.shape)
        if w_mask.mask.any():
            w_bin = np.asarray(
                [np.int(w_mask.mask.reshape(nx_ * ny_)[i]) for i in range(nx_ * ny_)]).reshape(nx_, ny_)
        w_mask_r = np.ma.masked_where(w_ > w_c, w_)
        if not w_mask_r.mask.any():
            print('STOP (t='+str(t0)+')' )
            continue
        else:
            w_bin_r = np.asarray(
                [np.int(w_mask_r.mask.reshape(nx_ * ny_)[i]) for i in range(nx_ * ny_)]).reshape(nx_, ny_)

        ax4.imshow(w_mask.T, cmap=cm_bwr, origin='lower', vmin=-max,vmax=max)
        ax5.imshow(w_mask.mask.T, origin='lower', vmin=-max,vmax=max)

        plt.suptitle('w masked on w='+str(w_c)+'m/s (z=' + str(k0 * dz) + ')')
        plt.savefig(os.path.join(path_out, 'w_masked_k'+str(k0)+'_thresh'+str(w_c)+'_t'+str(t0)+'.png'))
        plt.close()

    #     s = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))
    #     # s_roll = np.roll(w[:, :, k0], [ishift, jshift], [0, 1])
    #     s_mask = np.ma.masked_where(w[:,:,k0] <= w_c, s[:,:,k0])
    #     plt.figure()
    #     ax1 = plt.subplot(1,2,1)
    #     ax2 = plt.subplot(1,2,2)
    #     ax1.imshow(s[:,:,k0].T, origin='lower')
    #     ax2.imshow(s_mask.T, origin='lower')
    #     # plt.colorbar(ax2)
    #     plt.suptitle('s masked on w='+str(w_c)+'m/s (z='+str(k0*dz)+')')
    #     plt.savefig(os.path.join(path_out, 's_masked_k'+str(k0)+'_thresh'+str(w_c)+'_t'+str(t0)+'.png'))
    #     plt.close()


        # i = 0
        # j = 0
        # j = jc
        # while (w_mask_r[i,j] and i<2*id-1):
        #     i += 1

        del w_roll

        # w_mask = True, if w<w_c
        # w_mask_r = True, if w>w_c
        rim = np.zeros((nx_,ny_),dtype=np.int)
        rim_r = np.zeros((nx_,ny_), dtype=np.int)
        # print('limits: ', ic-id, ic+id, jc-jd, jc+jd, ishift, jshift)
        # print(w_mask_r.mask.shape, nx_, ny_)
        for i in range(ic - id - ishift, ic + id - ishift):
            for j in range(jc - jd - jshift, jc + jd - jshift):
                # print('i,j', i,j, w_mask_r.mask.shape)
                if w_mask_r.mask[i,j]:
                    a = np.count_nonzero(w_bin_r[i-1:i+2,j-1:j+2])
                    print('first step', i, j, ic, jc)
                    # print(w_bin_r[i-1:i+2,j-1:j+2], a)

                    if a > 4 and a < 9:
                        rim[i,j] = 1
                        print('second step')

        # i0 = 44
        # j0 = 27
        # i1 = i0
        # j1 = 10
        plt.figure()
        plt.subplot(1,3,1)
        plt.imshow(w_mask.T, origin='lower')
        plt.subplot(1,3,2)
        plt.imshow(w_mask.mask.T, origin='lower')
        plt.subplot(1,3,3)
        plt.imshow(rim.T, origin='lower')
        plt.colorbar()
        plt.savefig(os.path.join(path_out, 'rim_searching_test_t0' + str(t0) + '.png'))


    return


# ----------------------------------
def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:,:,:]
    rootgrp.close()
    return data

if __name__ == '__main__':
    main()