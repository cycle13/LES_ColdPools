import numpy as np
import matplotlib.pyplot as plt
import os

global cm_bwr, cm_grey, cm_vir
cm_bwr = plt.cm.get_cmap('bwr')
cm_vir = plt.cm.get_cmap('jet')
cm_grey = plt.cm.get_cmap('gist_gray_r')
# ----------------------------------
def plot_cp_rim_averages(rim_vel, rim_vel_av, timerange):
    nt = rim_vel.shape[1]
    plt.figure(figsize=(12, 5))
    plt.subplot(121)
    plt.plot(timerange[:nt], rim_vel_av[0,:nt],'-o')
    plt.xlabel('t  [s]')
    plt.ylabel('r  [m]')
    plt.title('average radius')


    plt.subplot(122)
    plt.plot(timerange[:nt], rim_vel_av[1,:],'-o')
    plt.xlabel('t  [s]')
    plt.ylabel('U  [m/s]')
    plt.title('average rim velocity')

    plt.savefig(os.path.join(path_out, 'rim_velocity_av.png'))
    plt.close()
    return


def plot_cp_rim_velocity(rim_vel, rim_vel_av, timerange):
    nt = rim_vel.shape[1]
    plt.figure(figsize=(12,5))
    ax = plt.subplot(122)
    for it, t0 in enumerate(timerange[0:nt]):
        ax.plot(rim_vel[0, it, :], rim_vel[2, it, :],
                label='t=' + str(t0) + 's', color=cm_vir(np.double(it) / nt))
    # plt.legend()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.0), ncol=3, fontsize=12)
                   # fancybox=True, shadow=True, ncol=5)
    plt.xlabel('phi  [deg]')
    plt.ylabel('U(phi)  [m/s]')
    plt.ylim([20,np.amax(rim_vel[2,:nt,:])+2])

    ax = plt.subplot(121, projection='polar')
    for it, t0 in enumerate(timerange[0:nt]):
        ax.plot(rim_vel[1, it, :], rim_vel[2, it, :],
                label='t=' + str(t0) + 's', color=cm_vir(np.double(it) / nt))
    ax.set_rmax(45.0)
    plt.suptitle('radial velocity of CP expansion (dr/dt)')
    plt.savefig(os.path.join(path_out, 'rim_velocity.png'))
    plt.close()
    return


def plot_cp_outline_alltimes(rim_intp_all, timerange):
    nt = rim_intp_all.shape[1]
    plt.figure()
    ax = plt.subplot(111, projection='polar')
    for it, t0 in enumerate(timerange[0:nt]):
        ax.plot(rim_intp_all[1,it,:], rim_intp_all[2,it,:],'-o', label='t='+str(t0)+'s',
                color = cm_vir(np.double(it)/nt))
    plt.legend()
    plt.suptitle('outline CP (all times)')
    plt.savefig(os.path.join(path_out, 'rim_cp1_alltimes.png'))
    plt.close()

    return




def plot_angles(rim_list, rim_intp, t0):
    plt.figure(figsize=(20,5))
    nx_plots = 4
    plt.subplot(1, nx_plots, 1)
    for i in range(len(rim_list)):
        plt.plot(rim_list[i][1][1], rim_list[i][1][0], 'x', color=cm_vir(float(i) / len(rim_list)))
    plt.xlabel('th')
    plt.ylabel('r')
    plt.title('rim list')
    plt.subplot(1, nx_plots, 2)
    for i in range(len(rim_list)):
        plt.plot(rim_list[i][1][1], rim_list[i][0][0], 'x', color=cm_vir(float(i) / len(rim_list)))
    plt.xlabel('th')
    plt.ylabel('x')
    plt.title('rim list')
    plt.subplot(1, nx_plots, 3)
    plt.plot(rim_intp[0,:], rim_intp[2,:], '-o')
    plt.title('rim interpolated')
    plt.xlabel('th')
    plt.ylabel('r')

    plt.subplot(144, projection='polar')
    plt.plot(rim_intp[1, :], rim_intp[2, :], '-o')

    plt.suptitle('t='+str(t0)+'s', fontsize=28)
    plt.savefig(os.path.join(path_out, 'angles_t'+str(t0)+'s.png'))
    plt.close()
    return



def plot_rim_mask(w_mask, rim, rim_list, rim_list_backup,
                  icshift, jcshift, nx_, ny_, t0):
    max = np.amax(w_mask)
    nx_plots = 4
    ny_plots = 2
    plt.figure(figsize=(5*nx_plots, 6*ny_plots))
    plt.subplot(ny_plots, nx_plots, 1)
    plt.imshow(w_mask.T, cmap=cm_bwr, origin='lower', vmin=-max, vmax=max)
    plt.title('w')
    plt.subplot(ny_plots, nx_plots, 2)
    plt.imshow(rim.T, origin='lower', cmap=cm_vir)
    plt.title('rim')
    plt.subplot(ny_plots, nx_plots, 3)
    plt.imshow(rim.T, origin='lower')
    for i in range(len(rim_list)):
        plt.plot(rim_list_backup[i][0], rim_list_backup[i][1], 'yx', markersize=2)
    plt.title('rim + rim_list')
    plt.xlim([0, nx_ - 1])
    plt.ylim([0, ny_ - 1])
    plt.subplot(ny_plots, nx_plots, 4)
    plt.plot(rim_list_backup)
    plt.title('orange=y, blue=x')
    plt.subplot(ny_plots, nx_plots, 5)
    for i in range(len(rim_list_backup)):
        plt.plot(rim_list_backup[i][0],rim_list_backup[i][1], 'x', color=cm_vir(float(i)/len(rim_list)))
    plt.plot(rim_list_backup[0][0], rim_list_backup[0][1], 'ko')
    plt.title('before sort')
    plt.subplot(ny_plots, nx_plots, 6)
    for i in range(len(rim_list_backup)):
        plt.plot(rim_list_backup[i][0]-icshift,rim_list_backup[i][1]-jcshift,
                 'x', color=cm_vir(float(i)/len(rim_list)))
    plt.plot(rim_list_backup[0][0]-icshift, rim_list_backup[0][1]-jcshift, 'ko')
    plt.title('shifted, before sort')
    plt.subplot(ny_plots, nx_plots, 7)
    for i in range(len(rim_list)):
        plt.plot(rim_list[i][0][0], rim_list[i][0][1], 'x', color=cm_vir(float(i)/len(rim_list)))
    plt.title('after sort (c=order)')
    plt.subplot(ny_plots, nx_plots, 8)
    for i in range(len(rim_list)):
        plt.plot(rim_list[i][0][0], rim_list[i][0][1],
                 'x', color=cm_vir(rim_list[i][1][1]/360))
    plt.title('after sort (c=angle)')
    plt.suptitle('cold pool outline - t='+str(t0)+'s',fontsize=28)
    plt.savefig(os.path.join(path_out,'rim_mask_t'+str(t0)+'s.png'))
    plt.close()

    return


def plot_outlines(perc, w_mask, rim, rim_list, rim_aux, icshift, jcshift, nx_, ny_, t0, path_out):
    max = np.amax(w_mask)
    nx_plots = 5
    ny_plots = 2
    plt.figure(figsize=(24,8))
    plt.subplot(ny_plots,nx_plots,1)
    plt.imshow(w_mask.T, cmap=cm_bwr, origin='lower', vmin=-max, vmax=max)
    plt.plot([nx_-1,nx_-1],[0,ny_-1],'r')
    plt.plot([0, nx_ - 1], [ny_-1, ny_ - 1], 'r')
    plt.title('w masked')
    plt.xlim([0, nx_-1])
    plt.ylim([0, ny_-1])

    plt.subplot(ny_plots,nx_plots,2)
    plt.imshow(w_mask.mask.T, origin='lower', cmap=cm_grey, interpolation='nearest')
    plt.title('mask')
    # plt.colorbar(shrink=0.5)

    plt.subplot(ny_plots, nx_plots, 3)
    plt.imshow(rim_aux.T, origin='lower', cmap=cm_grey)
    plt.title('covered ij (rim_aux)')

    plt.subplot(ny_plots,nx_plots,4)
    plt.title('rim int')
    plt.plot([0,nx_-1],[jcshift, jcshift],'b')
    plt.imshow(rim.T, origin='lower', cmap=cm_grey)
    for tup in rim_list[-2:]:
        plt.plot(tup[0], tup[1],'xr')
    plt.xlim([0, nx_ - 1])
    plt.ylim([0, ny_ - 1])

    plt.subplot(ny_plots,nx_plots,5)
    plt.title('mask + rim int')
    plt.imshow(w_mask.mask.T, origin='lower', cmap=cm_vir)
    plt.imshow(rim.T, origin='lower', cmap=cm_grey, alpha=0.5)
    plt.plot([icshift,icshift],[0,ny_-1], 'k')
    plt.plot([0,nx_-1],[jcshift,jcshift],'k')
    plt.xlim([0, nx_ - 1])
    plt.ylim([0, ny_ - 1])

    # plt.subplot(ny_plots,nx_plots,6)

    # plt.subplot(ny_plots,nx_plots,7)

    #



    plt.savefig(os.path.join(path_out, 'rim_searching_perc'+str(perc)+'_t' + str(t0) + '_v2.png'))
    plt.close()
    return


def plot_s(w,w_c,t0,k0):
    s = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))
    # s_roll = np.roll(w[:, :, k0], [ishift, jshift], [0, 1])
    s_mask = np.ma.masked_where(w[:,:,k0] <= w_c, s[:,:,k0])
    plt.figure()
    ax1 = plt.subplot(1,2,1)
    ax2 = plt.subplot(1,2,2)
    ax1.imshow(s[:,:,k0].T, origin='lower')
    ax2.imshow(s_mask.T, origin='lower')
    # plt.colorbar(ax2)
    plt.suptitle('s masked on w='+str(w_c)+'m/s (z='+str(k0*dz)+')')
    plt.savefig(os.path.join(path_out, 's_masked_k'+str(k0)+'_thresh'+str(w_c)+'_t'+str(t0)+'.png'))
    plt.close()
    del s
    return



def plot_w_field(w_c, perc, w, w_roll, w_, w_mask,
                 ishift, jshift, id, jd, ic, jc, icshift, jcshift,
                 k0, t0, dz, gw, nx_ , ny_, nx, ny, path_out):
    max = np.amax(w)
    plt.figure(figsize=(20, 6))
    ax1 = plt.subplot(1, 5, 1)
    ax1.set_title('w')
    # plt.contourf(w[:, :, k0].T, cmap=cm_bwr, aspect=1.)
    # plt.colorbar()
    ax2 = plt.subplot(1, 5, 2)
    ax2.set_title('np.roll(w)')
    ax3 = plt.subplot(1, 5, 3)
    ax3.set_title('w_')
    ax4 = plt.subplot(1, 5, 4)
    ax4.set_title('masked array')
    ax5 = plt.subplot(1, 5, 5)
    ax5.set_title('mask')
    cax = ax1.imshow(w[:, :, k0].T, cmap=cm_bwr, origin='lower', vmin=-max, vmax=max)
    cbar = plt.colorbar(cax, ticks=np.arange(-np.floor(max), np.floor(max) + 1), shrink=0.5)
    ax1.plot([ic, ic], [0, ny - 1], 'b')
    ax1.plot([ic + gw, ic + gw], [0, ny - 1], 'b--')
    ax1.plot([0, nx - 1], [jc, jc], 'r')
    ax1.plot([0, nx - 1], [jc + gw, jc + gw], 'r--')
    ax2.imshow(w_roll.T, cmap=cm_bwr, origin='lower', vmin=-max, vmax=max)
    ax2.plot([ic + ishift, ic + ishift], [0, ny - 1])
    ax2.plot([ic - id + ishift, ic - id + ishift], [0, ny - 1], '--')
    ax2.plot([ic + id + ishift, ic + id + ishift], [0, ny - 1], '--')
    ax2.plot([0, nx - 1], [jc + jshift, jc + jshift])
    ax2.plot([0, nx - 1], [jc + jshift - jd, jc + jshift - jd], 'k--')
    ax2.plot([0, nx - 1], [jc + jshift + jd, jc + jshift + jd], 'k--')
    ax3.imshow(w_.T, cmap=cm_bwr, origin='lower', vmin=-max, vmax=max)
    ax3.plot([icshift, icshift], [0, ny_ - 1])
    ax3.plot([0, nx_ - 1], [jcshift, jcshift])
    ax4.imshow(w_mask.T, cmap=cm_bwr, origin='lower', vmin=-max, vmax=max)
    ax5.imshow(w_mask.mask.T, origin='lower', vmin=-max, vmax=max)

    plt.suptitle(
        'w masked on ' + str(perc) + 'th percentile: w=' + str(np.round(w_c, 2)) + 'm/s (z=' + str(k0 * dz) + ')')
    plt.savefig(os.path.join(path_out, 'w_masked_k' + str(k0) + '_perc' + str(perc) + '_t' + str(t0) + '.png'))
    plt.close()
    return



def plot_yz_crosssection(w,ic,path_out,t0):
    max = np.amax(w)
    plt.figure()
    cax = plt.imshow(w[ic, :, :100].T, cmap=cm_bwr, origin='lower', vmin=-max, vmax=max)
    plt.colorbar(cax, ticks=np.arange(-np.floor(max), np.floor(max) + 1))
    plt.savefig(os.path.join(path_out, 'w_crosssection_t' + str(t0) + '.png'))
    plt.close()

    return
