import numpy as np
import matplotlib.pyplot as plt
import os

# global cm_bwr, cm_grey, cm_vir
# cm_bwr = plt.cm.get_cmap('bwr')
# cm_vir = plt.cm.get_cmap('viridis')
# cm_grey = plt.cm.get_cmap('gist_gray_r')

def set_colorbars(cm_bwr_, cm_vir_, cm_grey_):

    global cm_bwr, cm_grey, cm_vir
    cm_bwr = cm_bwr_
    cm_grey = cm_grey_
    cm_vir = cm_vir_
    return



# ----------------------------------
def plot_cp_rim_averages(rim_vel, rim_vel_av, timerange, path_out):
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

    plt.savefig(os.path.join(path_out, 'rim_velocity_av_v2.png'))
    plt.close()
    return


def plot_cp_rim_velocity(rim_vel, rim_vel_av, timerange, path_out):
    nt = rim_vel.shape[1]
    plt.figure(figsize=(12,5))
    ax = plt.subplot(121, projection='polar')
    for it, t0 in enumerate(timerange[0:nt]):
        ax.plot(rim_vel[1, it, :], rim_vel[3, it, :],
                label='t=' + str(t0) + 's', color=cm_vir(np.double(it) / nt))
    # ax.set_rmax(np.minimum(np.amax(rim_vel[3,:,:]),0.1))

    ax = plt.subplot(122)
    for it, t0 in enumerate(timerange[0:nt]):
        ax.plot(rim_vel[0, it, :], rim_vel[3, it, :], '-', linewidth=2,
                label='t=' + str(t0) + 's', color=cm_vir(np.double(it) / nt))
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    # fancybox=True, shadow=True, ncol=5)
    plt.xlabel('phi  [deg]')
    plt.ylabel('U(phi)  [m/s]')
    # plt.ylim([0,np.minimum(np.amax(rim_vel[3,:,:]),0.1)])

    plt.suptitle('radial velocity of CP expansion (dr/dt)')
    plt.savefig(os.path.join(path_out, 'rim_velocity_v2.png'))
    plt.close()
    return


def plot_rim_thickness(rim_intp_all, timerange, dx, path_out):
    nt = rim_intp_all.shape[1]
    plt.figure(figsize=(12,6))
    # thickness = rim_intp_all[2, :, :] - rim_intp_all[3, :, :]
    thickness = rim_intp_all[4, :, :]
    th_av = np.average(thickness[:,:],axis=1)
    ax = plt.subplot(122)
    for it, t0 in enumerate(timerange[0:nt]):
        ax.plot(rim_intp_all[1,it,:],thickness[it,:], label='t='+str(t0)+'s',
                color = cm_vir(np.double(it)/nt))
    plt.ylim([100,2000])
    plt.xlabel('angle phi  [deg]')
    plt.ylabel('thickness D  [m] (dx=' + str(dx) + ')')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    plt.title('rim thickness')
    ax = plt.subplot(121)
    ax.plot(timerange, th_av, '-o')
    plt.xlabel('time')
    plt.ylabel('average thickness D  [m] (dx=' + str(dx) + ')')
    plt.title('average rim thickness D')

    plt.savefig(os.path.join(path_out, 'rim_cp1_thickness.png'))
    plt.close()
    return


def plot_cp_outline_alltimes(rim_intp_all, timerange, dx, path_out):
    nt = rim_intp_all.shape[1]
    plt.figure(figsize=(12,10))

    ax = plt.subplot(121, projection='polar')
    for it, t0 in enumerate(timerange[0:nt]):
        ax.plot(rim_intp_all[1,it,:], rim_intp_all[2,it,:],'-o', label='t='+str(t0)+'s',
                color = cm_vir(np.double(it)/nt))
    # plt.legend()
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)

    ax = plt.subplot(122)
    for it, t0 in enumerate(timerange[0:nt]):
        ax.plot(rim_intp_all[1,it,:], rim_intp_all[2,it,:], '-o', label='t='+str(t0)+'s',
                color = cm_vir(np.double(it)/nt))
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    plt.xlabel('angle phi  [deg]')
    plt.ylabel('radius r  [m] (dx='+str(dx)+')')
    plt.suptitle('outline CP (outer)', fontsize=28)
    plt.savefig(os.path.join(path_out, 'rim_cp1_alltimes_v2.png'))
    plt.close()

    plt.figure(figsize=(12,10))
    ax = plt.subplot(121, projection='polar')
    for it, t0 in enumerate(timerange[0:nt]):
        ax.plot(rim_intp_all[1,it,:], rim_intp_all[3,it,:],'-o', label='t='+str(t0)+'s',
                color = cm_vir(np.double(it)/nt))

    ax = plt.subplot(122)
    for it, t0 in enumerate(timerange[0:nt]):
        ax.plot(rim_intp_all[1,it,:], rim_intp_all[3,it,:], '-o', label='t='+str(t0)+'s',
                color = cm_vir(np.double(it)/nt))
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    plt.xlabel('angle phi  [deg]')
    plt.ylabel('radius r  [m] (dx=' + str(dx) + ')')
    plt.suptitle('outline CP (inner)', fontsize=28)
    plt.savefig(os.path.join(path_out, 'rim_cp1_alltimes_int_v2.png'))
    plt.close()







    return


def plot_angles(rim_list, rim_list_int, rim_intp, t0, path_out):
    plt.figure(figsize=(20,10))
    nx_plots = 4
    ny_plots = 2
    plt.subplot(ny_plots, nx_plots, 1)
    for i in range(len(rim_list)):
        plt.plot(rim_list[i][1][1], rim_list[i][1][0], 'x', color=cm_vir(float(i) / len(rim_list)))
    plt.xlabel('th')
    plt.ylabel('r')
    plt.title('rim list')

    plt.subplot(ny_plots, nx_plots, 2)
    for i in range(len(rim_list)):
        plt.plot(rim_list[i][1][1], rim_list[i][0][0], 'x', color=cm_vir(float(i) / len(rim_list)))
    plt.xlabel('th')
    plt.ylabel('x')
    plt.title('rim list')

    plt.subplot(ny_plots, nx_plots, 3)
    plt.plot(rim_intp[0,:], rim_intp[2,:], '-o')
    plt.title('rim interpolated')
    plt.xlabel('th')
    plt.ylabel('r')

    plt.subplot(244, projection='polar')
    plt.plot(rim_intp[1, :], rim_intp[2, :], '-o')

    plt.subplot(ny_plots, nx_plots, 5)
    for i in range(len(rim_list_int)):
        plt.plot(rim_list_int[i][1][1], rim_list_int[i][1][0], 'x', color=cm_vir(float(i) / len(rim_list_int)))
    plt.xlabel('th')
    plt.ylabel('r')
    plt.title('rim list')

    plt.subplot(ny_plots, nx_plots, 6)
    for i in range(len(rim_list_int)):
        plt.plot(rim_list_int[i][1][1], rim_list_int[i][0][0], 'x', color=cm_vir(float(i) / len(rim_list_int)))
    plt.xlabel('th')
    plt.ylabel('x')
    plt.title('rim list')

    plt.subplot(ny_plots, nx_plots, 7)
    plt.plot(rim_intp[0,:], rim_intp[3 ,:], '-o')
    plt.title('rim interpolated')
    plt.xlabel('th')
    plt.ylabel('r')

    plt.subplot(248, projection='polar')
    plt.plot(rim_intp[1, :], rim_intp[3, :], '-o')

    plt.suptitle('t='+str(t0)+'s', fontsize=28)
    plt.savefig(os.path.join(path_out, 'angles_t'+str(t0)+'s_v2.png'))
    plt.close()
    return



def plot_rim_mask(w, w_mask, rim_out, rim_int, rim_list, rim_list_int,
                  icshift, jcshift, nx_, ny_, t0, path_out):
    max = np.amax(w)
    nx_plots = 3
    ny_plots = 2
    plt.figure(figsize=(5*nx_plots, 6*ny_plots))
    plt.subplot(ny_plots, nx_plots, 1)
    plt.imshow(w.T, cmap=cm_bwr, origin='lower', vmin=-max, vmax=max)
    plt.title('w')
    plt.subplot(ny_plots, nx_plots, 2)
    plt.imshow(rim_out.T, origin='lower', cmap=cm_vir)
    plt.title('rim out')
    plt.subplot(ny_plots, nx_plots, 3)
    for i in range(len(rim_list)):
        plt.plot(rim_list[i][0][0], rim_list[i][0][1],
                 'x', color=cm_vir(rim_list[i][1][1]/360))
    plt.title('outer: after sort (c=angle)')

    plt.subplot(ny_plots, nx_plots, 4)
    plt.imshow(w_mask.T, cmap=cm_bwr, origin='lower', vmin=-max, vmax=max)
    plt.title('w mask')
    plt.subplot(ny_plots, nx_plots, 5)
    plt.imshow(rim_int.T, origin='lower', cmap=cm_vir)
    plt.title('rim int')
    plt.subplot(ny_plots, nx_plots, 6)
    for i in range(len(rim_list_int)):
        plt.plot(rim_list_int[i][0][0], rim_list_int[i][0][1],
                 'x', color=cm_vir(rim_list_int[i][1][1]/360))
    plt.title('inner: after sort (c=angle)')
    plt.suptitle('cold pool outline - t='+str(t0)+'s',fontsize=28)
    plt.savefig(os.path.join(path_out,'rim_mask_t'+str(t0)+'s_v2.png'))
    plt.close()

    return


def plot_outlines(perc, w_mask, rim_int, rim_out, rim_list, rim_aux, rmax2, icshift, jcshift, imin, imax, jmin, jmax,
                  nx_, ny_, t0, path_out):
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

    ax = plt.subplot(ny_plots,nx_plots,2)
    plt.imshow(w_mask.mask.T, origin='lower', cmap=cm_grey, interpolation='nearest')
    circle1 = plt.Circle((icshift, jcshift), np.sqrt(rmax2), fill=False, color='y', linewidth=2)
    ax.add_artist(circle1)
    plt.plot([nx_ - 1, nx_ - 1], [jmin, jmin], 'r', linewidth=2)
    plt.plot([nx_ - 1, nx_ - 1], [jmax, jmax], 'r', linewidth=2)
    plt.plot([imin, imin], [ny_ - 1, ny_ - 1], 'r', linewidth=2)
    plt.plot([imax, imax], [ny_ - 1, ny_ - 1], 'r', linewidth=2)
    plt.xlim([-20,150])
    plt.title('mask')

    plt.subplot(ny_plots, nx_plots, 3)
    plt.imshow(rim_aux.T, origin='lower', cmap=cm_grey)
    plt.title('covered ij (rim_aux)')

    plt.subplot(ny_plots,nx_plots,4)
    plt.title('rim int')
    plt.plot([0,nx_-1],[jcshift, jcshift],'b')
    plt.imshow(rim_int.T, origin='lower', cmap=cm_grey)
    plt.xlim([0, nx_ - 1])
    plt.ylim([0, ny_ - 1])

    plt.subplot(ny_plots,nx_plots,5)
    plt.title('mask + rim int')
    plt.imshow(w_mask.mask.T, origin='lower', cmap=cm_vir)
    plt.imshow(rim_int.T, origin='lower', cmap=cm_grey, alpha=0.5)
    plt.plot([icshift,icshift],[0,ny_-1], 'k')
    plt.plot([0,nx_-1],[jcshift,jcshift],'k')
    plt.xlim([0, nx_ - 1])
    plt.ylim([0, ny_ - 1])


    # plt.subplot(ny_plots,nx_plots,6)

    plt.subplot(ny_plots,nx_plots,9)
    plt.title('rim out')
    plt.plot([0, nx_ - 1], [jcshift, jcshift], 'b')
    plt.imshow(rim_out.T, origin='lower', cmap=cm_grey)
    plt.xlim([0, nx_ - 1])
    plt.ylim([0, ny_ - 1])
    # for tup in rim_list[-2:]:
    #     plt.plot(tup[0], tup[1],'xr')

    plt.subplot(ny_plots,nx_plots,10)
    plt.title('mask + rim int')
    plt.imshow(w_mask.mask.T, origin='lower', cmap=cm_vir)
    plt.imshow(rim_out.T, origin='lower', cmap=cm_grey, alpha=0.5)
    plt.plot([icshift,icshift],[0,ny_-1], 'k')
    plt.plot([0,nx_-1],[jcshift,jcshift],'k')
    plt.xlim([0, nx_ - 1])
    plt.ylim([0, ny_ - 1])

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
    plt.savefig(os.path.join(path_out, 'w_masked_k' + str(k0) + '_perc' + str(perc) + '_t' + str(t0) + '_v2.png'))
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
