import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import netCDF4 as nc
import argparse
import json as simplejson
import os

def main():
    define_geometry()


    th_g = 300.0  # temperature for neutrally stratified background (value from Soares Surface)

    dTh_range = [2, 3, 4]
    # dTh_range = [3]
    for dTh in dTh_range:
        if dTh == 2:
            zstar_range = [2450, 1225, 815]
            rstar_range = [815, 1225, 2450]
        elif dTh == 3:
            zstar_range = [4000, 2000, 1500, 1000, 670, 500, 250]
            rstar_range = [250, 500, 670, 1000, 1500, 2000, 4000]
            # zstar_range = [1500, 1000]
            # rstar_range = [670, 1000]
        elif dTh == 4:
            zstar_range = [1730, 870, 430]
            rstar_range = [430, 870, 1730]
        n_thermo = len(dTh_range)
        n_geom = len(zstar_range)

        for ng in range(n_geom):
            zstar = zstar_range[ng]
            rstar = rstar_range[ng]
            irstar = np.int(np.round(rstar / dx))
            kstar = np.int(np.round(zstar / dz))
            print('r*, ir*, z*, k*: ', rstar, irstar, zstar, kstar)

            k_max_arr = (-1)*np.ones((2, nlg[0], nlg[1]), dtype=np.double)
            z_max_arr = np.zeros((2, nlg[0], nlg[1]), dtype=np.double)
            theta = th_g * np.ones(shape=(nlg[0], nlg[1], nlg[2]))
            theta_pert = np.random.random_sample(npg)
            entropy = np.empty((npl), dtype=np.double, order='c')
            for i in xrange(nlg[0]):
                ishift = i * nlg[1] * nlg[2]
                for j in xrange(nlg[1]):
                    jshift = j * nlg[2]
                    r = np.sqrt((x_half[i + indx_lo[0]] - xc) ** 2 +
                                (y_half[j + indx_lo[1]] - yc) ** 2)
                    if r <= rstar:
                        k_max = kstar * (np.cos(r / rstar * np.pi / 2)) ** 2
                        k_max_arr[0, i, j] = np.int(np.round(k_max))
                        z_max = zstar * ( np.cos(r/rstar * np.pi/2) **2 )
                        z_max_arr[0, i, j] = z_max
                    if r <= (rstar + marg):
                        # count_1 += 1
                        k_max = (kstar + marg_i) * (np.cos(r / (rstar + marg) * np.pi / 2)) ** 2
                        k_max_arr[1, i, j] = np.int(np.round(k_max))
                        z_max = (zstar + marg) * (np.cos(r / (rstar + marg) * np.pi / 2) ** 2)
                        z_max_arr[1, i, j] = z_max
            # count_0 = 0
            # for i in xrange(nlg[0]):
            #     for j in xrange(nlg[1]):
            #         for k in xrange(gw,nlg[2]-gw):
                    for k in xrange(nlg[2]-gw):
                        ijk = ishift + jshift + k
                        if k <= k_max_arr[0, i, j] + gw:
                            if (k-gw) > kstar:
                                print('k', k)
                            # count_0 += 1
                            theta[i, j, k] = th_g - dTh
                        elif (k-gw) <= k_max_arr[1, i, j]:
                            th = th_g - dTh * np.sin(( (k-gw) - k_max_arr[1, i, j]) / (k_max_arr[1, i, j] - k_max_arr[0, i, j]) * np.pi/2) ** 2
                            # th = th_g - dTh * np.sin((k - k_max_arr[1, i, j]) / (k_max_arr[1, i, j] - k_max_arr[0, i, j])) ** 2
                            # th = th_g - dTh * np.cos((k - k_max_arr[1, i, j]) / (k_max_arr[1, i, j] - k_max_arr[0, i, j])) ** 2
                            theta[i, j, k] = th
                        if k <= kstar + 2:
                            theta_pert_ = (theta_pert[ijk] - 0.5) * 0.1
                        else:
                            theta_pert_ = 0.0
                        # PV.values[s_varshift + ijk] = entropy_from_thetas_c(theta[i, j, k] + theta_pert_, 0.0)
                        # entropy[ijk] = entropy_from_thetas_c(theta[i, j, k] + theta_pert_, 0.0)

            theta_z = th_g * np.ones(shape=(nlg[0], nlg[1], nlg[2]))
            for i in xrange(nlg[0]):
                for j in xrange(nlg[1]):
                    for k in xrange(nlg[2]):
                        ijk = ishift + jshift + k
                        if z_half[k] <= z_max_arr[0, i, j]:
                            theta_z[i, j, k] = th_g - dTh
                        elif z_half[k] <= z_max_arr[1, i, j]:
                            th = th_g - dTh * np.sin((z_half[k] - z_max_arr[1, i, j]) / (z_max_arr[1, i, j] - z_max_arr[0, i, j]) * np.pi/2) ** 2
                            theta_z[i, j, k] = th
                        if k <= kstar + 2:
                            theta_pert_ = (theta_pert[ijk] - 0.5) * 0.1
                        else:
                            theta_pert_ = 0.0
                        # PV.values[s_varshift + ijk] = entropy_from_thetas_c(theta_z[i, j, k] + theta_pert_, 0.0)
                        # entropy[ijk] = entropy_from_thetas_c(theta_z[i, j, k] + theta_pert_, 0.0)

            icg = ic + gw
            jcg = jc + gw
            print('max(k_max[0,:,:])', np.amax(k_max_arr[0,:,:]))
            print('max(k_max[0,ic+gw-1,:])', np.amax(k_max_arr[0,icg-1,:]))
            print('max(k_max[0,ic+gw,:])', np.amax(k_max_arr[0,icg,:]))
            print('max(k_max[0,ic+gw+1,:])', np.amax(k_max_arr[0,icg+1,:]))
            print('')


            ''' plot theta[k=0]'''
            theta_ = theta[gw:-gw,gw:-gw,gw:-gw]
            theta_z_ = theta_z[gw:-gw, gw:-gw, gw:-gw]
            fig, axes = plt.subplots(2, 4, figsize=(20, 10))
            # ax = axes[0, 0]
            ax1 = axes[0,0]
            ax2 = axes[1,0]
            im = ax1.imshow(theta_[:,:,0].T, origin='lower', cmap = cm.bwr)
            plt.colorbar(im, ax=ax1, shrink=0.5)
            ax1.plot(ic, jc, 'or')
            ax1.plot([ic, ic], [0, ny], 'k')
            ax1.plot([ic+irstar, ic+irstar], [0, ny], 'w:')
            ax1.plot([ic-irstar, ic-irstar], [0, ny], 'w:')
            ax1.plot([jc - irstar - marg_i, jc - irstar - marg_i], [0, ny], '--', color='w', linewidth=1,
                     label='jc-irstar-marg_i')
            ax1.plot([jc + irstar + marg_i, jc + irstar + marg_i], [0, ny], '--', color='w', linewidth=1)
            ax1.plot([xc/dx, xc/dx], [0, ny], '--k')
            ax1.plot([0, nx], [jc, jc], 'k')
            circle1 = plt.Circle((ic, jc), irstar, fill=False, color='lime', linewidth=2)
            circle2 = plt.Circle((ic, jc), irstar + marg_i, fill=False, color='lime', linewidth=1)
            ax1.add_artist(circle1)
            ax1.add_artist(circle2)
            ax1.set_xlim([0, nx])
            ax1.set_ylim([0, ny])
            # ax1.set_title(str(count_0))
            ax2.contourf(x_half[gw:-gw], y_half[gw:-gw], theta_[:,:,0])
            ax2.plot([xc,xc], [y_half[gw],y_half[-gw]], 'k')
            ax2.plot([x_half[gw], x_half[-gw]], [yc, yc], 'k')


            ax1 = axes[0, 1]
            ax2 = axes[1, 1]
            im = ax1.imshow(theta_[ic,:,:].T, origin='lower', cmap = cm.bwr)
            ax1.plot(k_max_arr[0,icg,gw:-gw], 'gold', label='k_max[0]', linewidth=3)
            ax1.plot(k_max_arr[1,icg,gw:-gw], 'lime', label='k_max[1]', linewidth=3)
            ax1.plot([jc,jc],[0,nz],'k')
            ax1.plot([jc - irstar, jc - irstar], [0, nz], ':', color='lightgray', linewidth=2, label='jc+irstar')
            ax1.plot([jc + irstar, jc + irstar], [0, nz], ':', color='lightgray', linewidth=2)
            ax1.plot([jc - irstar - marg_i, jc - irstar - marg_i], [0, nz], '--', color='lightgray', linewidth=2,
                     label='jc-irstar-marg_i')
            ax1.plot([jc + irstar + marg_i, jc + irstar + marg_i], [0, nz], '--', color='lightgray', linewidth=2)
            ax1.plot([0, ny], [kstar, kstar], color='aqua', label='kstar', linewidth=2)
            ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
                       fancybox=True, shadow=True, ncol=2, fontsize=10)
            ax1.set_title('max(kmax0)='+str(np.amax(k_max_arr[0,:,:])) +
                            ', max(kmax0[ic])=' + str(np.round(np.amax(k_max_arr[0, icg, :]))), fontsize = 10)
            ax1.grid()
            ax1.set_xlim([0,ny])
            ax1.set_ylim([0,nz])
            ax2.contourf(theta_[ic,:,:].T, cmap=cm.bwr)
            ax2.plot(k_max_arr[0,icg,gw:-gw], 'gold', label='k_max[0]', linewidth=3)
            ax2.plot(k_max_arr[1,icg,gw:-gw], 'lime', label='k_max[1]', linewidth=3)
            ax2.plot([jc,jc],[0,nz],'k')
            ax2.plot([jc-irstar,jc-irstar],[0,nz], ':', color='lightgray', linewidth=2, label='jc+irstar')
            ax2.plot([jc+irstar,jc+irstar],[0,nz],':', color='lightgray', linewidth=2)
            ax2.plot([jc-irstar-marg_i,jc-irstar-marg_i],[0,nz],'--', color='lightgray', linewidth=2, label='jc-irstar-marg_i')
            ax2.plot([jc+irstar+marg_i,jc+irstar+marg_i],[0,nz],'--', color='lightgray', linewidth=2)
            ax2.grid()
            ax2.set_xlim([0,nx])

            ax1 = axes[0, 2]
            ax2 = axes[1, 2]
            ax1.imshow(theta_z_[ic,:,:].T, origin='lower', cmap = cm.bwr)
            ax1.plot(z_max_arr[0,icg,gw:-gw]/dz, 'gold', label='z_max[0]/dz', linewidth=3)
            ax1.plot(z_max_arr[1,icg,gw:-gw]/dz, 'lime', label='z_max[1]/dz', linewidth=3)
            ax1.plot([jc,jc],[0,nz],'k')
            ax1.plot([jc - irstar, jc - irstar], [0, nz], ':', color='lightgray', linewidth=2, label='jc+irstar')
            ax1.plot([jc + irstar, jc + irstar], [0, nz], ':', color='lightgray', linewidth=2)
            ax1.plot([jc - irstar - marg_i, jc - irstar - marg_i], [0, nz], '--', color='lightgray', linewidth=2,
                     label='jc-irstar-marg_i')
            ax1.plot([jc + irstar + marg_i, jc + irstar + marg_i], [0, nz], '--', color='lightgray', linewidth=2)
            ax1.plot([0, ny], [kstar, kstar], color='aqua', label='kstar', linewidth=2)
            ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
               fancybox=True, shadow=True, ncol=2, fontsize=10)
            ax1.set_title('max(zmax0)=' + str(np.amax(z_max_arr[0, :, :])) +
                          ', max(zmax0[ic])=' + str(np.round(np.amax(z_max_arr[0, icg, :]))), fontsize=10)
            ax1.grid()
            ax1.set_xlim([0,ny])
            ax1.set_ylim([0,nz])
            ax2.contourf(y_half[gw:-gw], z_half[gw:-gw], theta_z_[ic,:,:].T, cmap = cm.bwr)
            ax2.plot(y_half[gw:-gw], z_max_arr[0,icg,gw:-gw], 'gold', label='z_max[0]', linewidth=3)
            ax2.plot(y_half[gw:-gw], (np.round(z_max_arr[0,icg,gw:-gw]/dz))*dz, '--', color='k', label='z_max[0]', linewidth=3)
            ax2.plot(y_half[gw:-gw], z_max_arr[1,icg,gw:-gw], 'lime', label='z_max[1]', linewidth=3)
            ax2.plot(y_half[gw:-gw], (np.round(z_max_arr[1, icg, gw:-gw] / dz)) * dz, '--', color='k', label='z_max[0]', linewidth=3)
            ax2.plot([jc,jc],[0,nz],'k')
            ax2.plot([jc-irstar,jc-irstar],[0,nz],'w:')
            ax2.plot([jc+irstar,jc+irstar],[0,nz],'w:')
            ax2.plot([jc-irstar-marg_i,jc-irstar-marg_i],[0,nz],'w--')
            ax2.plot([jc+irstar+marg_i,jc+irstar+marg_i],[0,nz],'w--')
            # ax4.legend()
            # legend = plt.legend(frameon=1)
            # frame.set_facecolor('green')
            # frame.set_edgecolor('red')

            ax1 = axes[0, 3]
            ax2 = axes[1, 3]
            im = ax1.imshow(theta_z_[:, :, 0].T, origin='lower', cmap=cm.bwr)
            plt.colorbar(im, ax=ax1, shrink=0.5)
            ax1.plot(ic, jc, 'or')
            ax1.plot([ic, ic], [0, ny], 'k')
            ax1.plot([ic + irstar, ic + irstar], [0, ny], 'w:')
            ax1.plot([ic - irstar, ic - irstar], [0, ny], 'w:')
            ax1.plot([jc - irstar - marg_i, jc - irstar - marg_i], [0, ny], '--', color='w', linewidth=1,
                     label='jc-irstar-marg_i')
            ax1.plot([jc + irstar + marg_i, jc + irstar + marg_i], [0, ny], '--', color='w', linewidth=1)
            ax1.plot([xc / dx, xc / dx], [0, ny], '--k')
            ax1.plot([0, nx], [jc, jc], 'k')
            circle1 = plt.Circle((ic, jc), irstar, fill=False, color='lime', linewidth=2)
            circle2 = plt.Circle((ic, jc), irstar+marg_i, fill=False, color='lime', linewidth=1)
            ax1.add_artist(circle1)
            ax1.add_artist(circle2)
            ax1.set_xlim([0, nx])
            ax1.set_ylim([0, ny])
            ax2.contourf(x_half[gw:-gw], y_half[gw:-gw], theta_z_[:, :, 0])
            ax2.plot([xc, xc], [y_half[gw], y_half[-gw]], 'k')
            ax2.plot([x_half[gw], x_half[-gw]], [yc, yc], 'k')

            plt.tight_layout
            plt.suptitle('r*='+str(rstar)+' (irstar='+str(irstar)+'), z*='+str(zstar) +' (kstar='+str(kstar)+')')
            fig.savefig('./initialization_k_dTh'+str(dTh)+'_z'+str(zstar)+'_r'+str(rstar)+'.png')
            plt.close(fig)





    return




def define_geometry():

    global nx, ny, nz, dx, dy, dz
    nx = 80
    ny = 80
    nz = 30
    dx = 100
    dy = 100
    dz = 100
    global gw, nl, nlg, npl, npg
    gw = 5
    # gw = 1
    # Gr.dims.nlg[i], i=0,1,2
    mpi_dims = [1, 1, 1]
    nl = np.ndarray(3, dtype=np.int)
    nl[0] = nx / mpi_dims[0]
    nl[1] = ny / mpi_dims[1]
    nl[2] = nz / mpi_dims[2]
    nlg = np.ndarray(3, dtype=np.int)
    nlg[0] = nx + 2*gw
    nlg[1] = ny + 2*gw
    nlg[2] = nz + 2*gw
    npl = nl[0] * nl[1] * nl[2]
    npg = nlg[0] * nlg[1] * nlg[2]
    global indx_lo
    # Gr.dims.indx_lo[i]
    indx_lo = np.zeros(3, dtype=np.int)

    global x_half, y_half, z_half
    x_half = np.empty((nx+2*gw), dtype=np.double, order='c')
    y_half = np.empty((ny+2*gw), dtype=np.double, order='c')
    z_half = np.empty((nz+2*gw), dtype=np.double, order='c')
    count = 0
    for i in xrange(-gw, nx+gw):
        x_half[count] = (i + 0.5) * dx
        count += 1
    count = 0
    for j in xrange(-gw, ny+gw):
        y_half[count] = (j + 0.5) * dy
        count += 1
    count = 0
    for i in xrange(-gw, nz+gw):
        z_half[count] = (i + 0.5) * dz
        count += 1

    global ic, jc, xc, yc, marg_i, marg
    # zstar_range = [4000, 2000, 1500, 1000, 670, 500, 250]
    # rstar_range = [250, 500, 670, 1000, 1500, 2000, 4000]
    # rstar_range = [2000]
    # zstar_range = [2000]
    ic = np.int(nx / 2)
    jc = np.int(ny / 2)
    xc = x_half[ic + gw]  # center of cold-pool !!!
    yc = y_half[jc + gw]  # center of cold-pool !!!
    marg_i_ = np.int(500. / np.round(dx))  # width of margin
    marg_i = np.int(np.round(500. / dx))  # width of margin
    marg_ = marg_i_ * dx  # width of margin
    marg = marg_i * dx  # width of margin


    return

if __name__ == '__main__':
    main()