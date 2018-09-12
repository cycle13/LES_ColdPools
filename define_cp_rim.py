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

    cm_bwr = plt.cm.get_cmap('bwr')

    # define subdomain to scan
    # --- for triple coldpool ---
    d = np.int(np.round(ny / 2)) + gw
    a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
    rstar = 5000.0  # half of the width of initial cold-pools [m]
    irstar = np.int(np.round(rstar / dx))
    ic = np.int(np.round(a / 2))
    jc = np.int(np.round(d / 2))
    shift = 20
    ishift = shift
    jshift = shift
    id = irstar + ishift
    jd = irstar + jshift
    nx_ = 2 * id
    ny_ = 2 * jd

    print(ic,jc,id,jd)

    # t0 = 200
    for t0 in np.arange(100,1000,100):
        print('time: '+ str(t0))

        k0 = 5

        w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0) + '.nc'))
        w_roll = np.roll(w[:, :, k0], [ishift, jshift], [0, 1])
        w_ = w_roll[ic - id + ishift:ic + id + ishift, jc - jd + jshift:jc + jd + jshift]
        icshift = id
        jcshift = jd

        max = np.amax(w)
        plt.figure()
        cax = plt.imshow(w[ic,:,:100].T, cmap=cm_bwr, origin='lower', vmin=-max, vmax=max)
        plt.colorbar(cax, ticks=np.arange(-np.floor(max),np.floor(max)+1))
        plt.savefig(os.path.join(path_out,'w_crosssection_t'+str(t0)+'.png'))
        plt.close()

        plt.figure(figsize=(20,6))
        ax1 = plt.subplot(1,5,1)
        ax1.set_title('w')
        # plt.contourf(w[:, :, k0].T, cmap=cm_bwr, aspect=1.)
        # plt.colorbar()
        ax2 = plt.subplot(1,5,2)
        ax2.set_title('np.roll(w)')
        ax3 = plt.subplot(1,5,3)
        ax3.set_title('w_')
        ax4 = plt.subplot(1,5,4)
        ax4.set_title('masked array')
        ax5 = plt.subplot(1,5,5)
        ax5.set_title('mask')
        cax = ax1.imshow(w[:,:,k0].T, cmap=cm_bwr, origin='lower', vmin=-max, vmax=max)
        cbar = plt.colorbar(cax, ticks=np.arange(-np.floor(max), np.floor(max) + 1), shrink=0.5)
        ax1.plot([ic,ic],[0,ny-1])
        ax1.plot([0,nx-1],[jc,jc])
        ax2.imshow(w_roll.T, cmap=cm_bwr, origin='lower', vmin=-max,vmax=max)
        ax2.plot([ic+ishift,ic+ishift],[0,ny-1])
        ax2.plot([ic-id+ishift,ic-id+ishift],[0,ny-1], '--')
        ax2.plot([ic+id+ishift,ic+id+ishift],[0,ny-1], '--')
        ax2.plot([0,nx-1],[jc+jshift,jc+jshift])
        ax3.imshow(w_.T, cmap=cm_bwr, origin='lower', vmin=-max,vmax=max)
        ax3.plot([icshift, icshift], [0, ny_ - 1])
        ax3.plot([0, nx_-1], [jcshift, jcshift])

        perc = 90
        w_c = np.percentile(w_, perc)
        w_mask = np.ma.masked_less(w_, w_c)
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

        plt.suptitle('w masked on '+str(perc)+'th percentile: w='+str(np.round(w_c,2))+'m/s (z=' + str(k0 * dz) + ')')
        plt.savefig(os.path.join(path_out, 'w_masked_k'+str(k0)+'_perc'+str(perc)+'_t'+str(t0)+'.png'))
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
        del w_roll

        ''' Rim based on outline '''
        rim_test = np.zeros((nx_, ny_), dtype=np.int)
        rim_test1 = np.zeros((nx_, ny_), dtype=np.int)
        rim_test2 = np.zeros((nx_, ny_), dtype=np.int)
        j = 0
        while(j<ny_-1):
            i = 0
            while (w_mask.mask[i,j] and i<nx_-1):
                i += 1
            if i <= ic:
                rim_test1[i,j] = 1
                rim_test[i, j] = 1
            i = nx_-1
            while (w_mask.mask[i,j] and i>=0):
                i -= 1
            if i >= ic:
                rim_test1[i, j] = 1
                rim_test[i, j] = 1
            j += 1

        i = 0
        while (i<nx_-1):
            j = 0
            while (w_mask.mask[i,j] and j<ny_-1):
                j+=1
            if j <= jc:
                rim_test[i,j] = 1
                rim_test2[i,j] = 1
            j = ny_-1
            while (w_mask.mask[i,j] and j>=0):
                j-=1
            if j >= jc:
                rim_test[i, j] = 1
                rim_test2[i,j] = 1
            i+=1



        ''' Rim based on number of neighbours '''
        # w_mask = True, if w<w_c
        # w_mask_r = True, if w>w_c
        rim2 = np.zeros((nx_,ny_),dtype=np.int)
        for i in range(nx_):
            for j in range(ny_):
                if w_mask_r.mask[i,j]:
                    a = np.count_nonzero(w_bin_r[i-1:i+2,j-1:j+2])
                    if a > 5 and a < 9:
                        rim2[i,j] = 1
        rim = np.zeros((nx_,ny_),dtype=np.int)
        for i in range(nx_):
            for j in range(jcshift+1):
                if w_mask_r.mask[i,j]:
                    a = np.count_nonzero(w_bin_r[i-1:i+2,j-1:j+2])
                    if a > 5 and a < 9:
                        rim[i,j] = 1
                        a = np.count_nonzero(w_bin_r[i-1:i+2, j:j+3])
                        if a<=5 or a>=9:
                            break
            for j in range(ny_-1,jcshift,-1):
                if w_mask_r.mask[i,j]:
                    a = np.count_nonzero(w_bin_r[i-1:i+2,j-1:j+2])
                    if a > 5 and a < 9:
                        rim[i,j] = 1
                        a = np.count_nonzero(w_bin_r[i-1:i+2, j-2:j+1])
                        if a <= 5 or a >= 9:
                            break


        # i0 = 44
        # j0 = 27
        # i1 = i0
        # j1 = 10
        nx_plots = 4
        ny_plots = 2
        plt.figure(figsize=(24,8))
        plt.subplot(ny_plots,nx_plots,1)
        plt.imshow(w_mask.T, cmap=cm_bwr, origin='lower', vmin=-max, vmax=max)
        plt.plot([nx_-1,nx_-1],[0,ny_-1],'r')
        plt.plot([0, nx_ - 1], [ny_-1, ny_ - 1], 'r')
        plt.title('w masked')
        plt.subplot(ny_plots,nx_plots,2)
        plt.imshow(w_mask.mask.T, origin='lower')
        plt.title('mask')
        # plt.colorbar(shrink=0.5)
        plt.subplot(ny_plots,nx_plots,3)
        plt.title('rim - #neighbours')
        plt.imshow(rim2.T, origin='lower')
        plt.plot([0, nx_ - 1], [jcshift, jcshift], 'w')

        plt.subplot(ny_plots,nx_plots,4)
        plt.title('rim2 - #neighbours')
        plt.imshow(rim.T, origin='lower')
        plt.plot([0,nx_-1],[jcshift, jcshift],'w')

        plt.subplot(ny_plots,nx_plots,5)
        plt.title('rim')
        imin = 15
        imax = 75
        plt.imshow(rim_test1[imin:imax,imin:imax].T, origin='lower')
        ax = plt.gca()
        ax.set_xticks(np.arange(-0.5, imax - imin - 0.5, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, imax - imin - 0.5, 1), minor=True)
        ax.grid(which='minor', color='w', linewidth=0.2)#, linestyle='-', linewidth=2)
        plt.subplot(ny_plots,nx_plots,6)
        plt.title('rim y')
        plt.imshow(rim_test2.T, origin='lower')
        plt.subplot(ny_plots,nx_plots,7)
        plt.title('rim - outermost')
        plt.imshow(rim_test.T, origin='lower')

        plt.subplot(ny_plots,nx_plots,8)
        plt.title('rim - outermost vs #neighbours')
        plt.imshow(rim_test.T, origin='lower')
        plt.imshow(rim.T, origin='lower',alpha =0.5)
        plt.plot([0, nx_ - 1], [jcshift, jcshift], 'w')

        plt.savefig(os.path.join(path_out, 'rim_searching_perc'+str(perc)+'_t0' + str(t0) + '.png'))
        plt.close()

        plt.figure()
        imin = 20
        imax = 60
        plt.imshow(rim_test1[imin:imax, :].T, origin='lower')
        # plt.grid()
        ax = plt.gca()
        # Minor ticks
        ax.set_xticks(np.arange(-0.5, imax - imin - 0.5, 1), minor=True)
        ax.grid(which='minor', color='w', linewidth=0.2)  # , linestyle='-', linewidth=2)
        # plt.show()


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