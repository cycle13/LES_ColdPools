import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# from matplotlib.widgets import TextBox
import netCDF4 as nc
import argparse
import json as simplejson
import os

execfile('settings.py')

# label_size = 12
# plt.rcParams['xtick.labelsize'] = label_size
# plt.rcParams['ytick.labelsize'] = label_size
# plt.rcParams['lines.linewidth'] = 1
# plt.rcParams['lines.markersize'] = 6
# plt.rcParams['legend.fontsize'] = 12
# plt.rcParams['axes.labelsize'] = 15
plt.rcParams['font.sans-serif'] = 'Helvetica'
# plt.rcParams['text.usetex'] = 'true'
plt.rcParams['legend.numpoints'] = 1

from thermodynamic_functions import thetas_c

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    args = parser.parse_args()

    res = 25
    run = 'run4'
    dTh = 3
    rstar = 1000
    zstar = 1000
    case = 'dTh'+str(dTh)+'_z'+str(zstar)+'_r'+str(rstar)
    # case_name = args.casename
    case_name = 'ColdPoolDry_single_3D'

    path_root = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/'
    path = os.path.join(path_root, run + '_dx'+str(res)+'m', case)
    path_fields = os.path.join(path, 'fields')
    # path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single'
    path_out_figs = '/nbi/home/meyerbe/paper_CP'

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    ic = nml['init']['ic']
    jc = nml['init']['jc']
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    dx = np.ndarray(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    kmax = np.int(1.0e3 / dx[2])
    print('res: '+str(dx))
    print('CP centre: '+ str(ic)+', '+str(jc))
    print('')

    print('read in vorticity')
    # root_vort = nc.Dataset(os.path.join(path, 'fields_vorticity/field_vort_xz.nc'))
    # time_vort = root_vort.groups['fields'].variables['time'][:]
    # vorticity_ = root_vort.groups['fields'].variables['vort_xz'][:, :, :kmax]
    # root_vort.close()
    root_vort = nc.Dataset(os.path.join(path, 'fields_vorticity/field_vort_yz.nc'))
    time_vort = root_vort.groups['fields'].variables['time'][:]
    vorticity_ = root_vort.groups['fields'].variables['vort_yz'][:, :, :kmax]
    root_vort.close()

    print('read in azimuthally averaged fields')
    root_vrad_2D = nc.Dataset(os.path.join(path, 'fields_v_rad/v_rad.nc'))
    # time_rad = root_vrad_2D.variables['time'][:]
    vrad_2D = root_vrad_2D.variables['v_rad'][:, :, jc, :kmax]
    root_vrad_2D.close()
    root_vrad = nc.Dataset(os.path.join(path, 'data_analysis/stats_radial_averaged.nc'))
    time_rad = root_vrad.groups['timeseries'].variables['time'][:]
    v_rad = root_vrad.groups['stats'].variables['v_rad'][:, :, :kmax]
    w_rad = root_vrad.groups['stats'].variables['w'][:, :, :kmax]
    s_rad = root_vrad.groups['stats'].variables['s'][:, :, :kmax]
    root_vrad.close()
    theta_rad = thetas_c(s_rad, 0.0)
    del s_rad
    root_vort_rad = nc.Dataset(os.path.join(path, 'fields_vorticity', 'field_vort_phi.nc'))
    vort_rad = root_vort_rad.groups['fields'].variables['vort_phi'][:,:,:]
    root_vort_rad.close()
    print('')
    print('dimensions: ', w_rad.shape, vort_rad.shape)
    print('')


    jmin_range = [150, 100, 100]
    for i,t0 in enumerate([900, 1200, 1500]):
    # for it,t0 in enumerate([900]):
        it = t0 / 100
        # print('')
        print('time: ', t0)
        # print('')
        jmin = jmin_range[i]
        jmax = ny - jmin
        rmax = np.int(.5*(jmax-jmin))
        fullpath_in = os.path.join(path_fields, str(t0)+'.nc')
        root_field = nc.Dataset(fullpath_in)
        grp = root_field.groups['fields']
        s = grp.variables['s'][ic,:,:kmax]
        w = grp.variables['w'][ic,:,:kmax]
        # v = grp.variables['v'][ic,:,:kmax]
        u = grp.variables['u'][ic,:,:kmax]
        root_field.close()
        theta = thetas_c(s, 0.0)[ic-nx/2:ic+nx/2, :]
        vorticity = vorticity_[np.int(t0/100),:,:]


        # fig_name = 'CP_structure_rad_dx' + str(res) + '_' + case + '_t'+str(t0)+'_test.png'
        # plot_test(theta, theta_rad, u, v_rad, w, w_rad, vorticity, vort_rad,
        #       it, jmin, jmax, rmax, kmax, dx, nx, ny,
        #       fig_name, path_out_figs)
        fig_name = 'CP_structure_rad_dx' + str(res) + '_' + case + '_t' + str(t0) + '.png'
        print(fig_name)
        plot_figure_rad(theta, theta_rad, u, v_rad, w, w_rad, vorticity, vort_rad,
              it, jmin, jmax, rmax, kmax, dx, nx, ny,
              fig_name, path_out_figs)



    return





def plot_figure_rad(theta, theta_rad, u, v_rad, w, w_rad, vorticity, vort_rad,
              it, jmin, jmax, rmax, kmax, dx, nx, ny,
              fig_name, path_out_figs):
    nlev = 2e2
    ncol = 2
    nrow = 4
    # title_pos_x = 850
    title_pos_x = -53
    title_pos_y = 43
    textprops = dict(facecolor='white', alpha=0.9, linewidth=0.)

    fig, axes = plt.subplots(nrow, ncol, sharex='col', figsize=(5 * ncol, 2 * nrow))
    for ax in axes[:, 1].flatten():
        ax.plot([0, jmax], [0, 0], '.5', linewidth=0.5)

    ax = axes[0, :]
    # min = 298
    min = np.round(np.amin(theta[:, :kmax]), 1)
    max = 300
    cf = ax[0].contourf(theta_rad[it, :, :kmax].T, levels=np.linspace(min, max, nlev), cmap=cm_bw_r, extend='max')
    cbar = plt.colorbar(cf, ax=ax[0], shrink=.75, aspect=12, ticks=np.arange(np.ceil(min), max + 0.5, 1.))
    # cbar = plt.colorbar(cf, ax=ax[0], shrink=.75, aspect=12, ticks=np.arange(min, max+0.5, 1.))
    for k in range(0, 500/dx[2]+1, 4):
        ax[1].plot(theta_rad[it, :, k], str((1.7*np.double(k)/kmax)), label='z='+str(k*dx[2])+'m')
        print('k', k, k*dx[2], str(1.8*np.double(k) / kmax))
    # ax[1].plot(theta_rad[it,:,1], '.5')
    ax[1].legend(loc='best', fontsize=8)
    ax[1].set_ylabel(r'$\theta$ [K]')
    ax[1].set_ylim(min, max)
    y_ticks = [np.int(ti) for ti in ax[1].get_yticks()]
    ax[1].set_yticklabels(y_ticks)
    for label in ax[1].yaxis.get_ticklabels():
        label.set_visible(False)
    if it == 9:
        for label in ax[1].yaxis.get_ticklabels()[4::5]:
            label.set_visible(True)
    elif it == 12:
        for label in ax[1].yaxis.get_ticklabels()[2::5]:
            label.set_visible(True)
    elif it == 15:
        for label in ax[1].yaxis.get_ticklabels()[1::5]:
            label.set_visible(True)
    ax[0].text(title_pos_x, title_pos_y, 'a) potential temperature', fontsize=15, horizontalalignment='left',
               bbox=textprops)

    ax = axes[1, :]
    min = -2.
    max = -min
    cf = ax[0].contourf(w_rad[it, :, :kmax].T, levels=np.linspace(min, max, nlev), cmap=cm_bwr, extend='both')
    cbar = plt.colorbar(cf, ax=ax[0], shrink=0.75, aspect=12, ticks=np.arange(min, max + 0.5, 1.))
    for k in range(0, 500/dx[2]+1, 4):
        ax[1].plot(w_rad[it, :, k], str((1.7*np.double(k)/kmax)), label='z='+str(k*dx[2])+'m')
    ax[1].set_ylabel(r'$w$ [m/s]')
    y_ticks = [np.int(ti) for ti in ax[1].get_yticks()]
    ax[1].set_yticklabels(y_ticks)
    for label in ax[1].yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    ax[0].text(title_pos_x, title_pos_y, 'c) vertical velocity', fontsize=15, horizontalalignment='left',
               bbox=textprops)

    ax = axes[2, :]
    min = -6.
    max = -min
    cf = ax[0].contourf(v_rad[it, :, :kmax].T, levels=np.linspace(min, max, nlev), cmap=cm_bwr, extend='both')
    cbar = plt.colorbar(cf, ax=ax[0], shrink=0.75, aspect=12, ticks=np.arange(min, max + 0.5, 2))
    # cbar.set_label(cont_var_name, rotation=90)
    for k in range(0, 500 / dx[2] + 1, 4):
        ax[1].plot(v_rad[it, :, k], str((1.7 * np.double(k) / kmax)), label='z=' + str(k * dx[2]) + 'm')
    ax[1].set_ylabel(r'$v_r$ [m/s]')
    # ax[1].set_ylim(-7,max)
    y_ticks = [ti for ti in ax[1].get_yticks()]
    ax[1].set_yticklabels(y_ticks)
    # for label in ax[1].yaxis.get_ticklabels()[1::2]:
    #     label.set_visible(False)
    ax[0].text(title_pos_x, title_pos_y, 'b) radial velocity', fontsize=15, horizontalalignment='left', bbox=textprops)

    ax = axes[3, :]
    min = -0.15
    max = np.amax(-vort_rad[it,:,:kmax])
    cf = ax[0].contourf(-vort_rad[it, :, :kmax].T, levels=np.linspace(0, max, nlev), cmap=cm_bw, extend='max')
    cbar = plt.colorbar(cf, ax=ax[0], shrink=0.75, aspect=12, ticks=np.arange(min, max + 0.02, 0.05))
    for k in range(0, 500/dx[2]+1, 4):
        ax[1].plot(-vort_rad[it, :, k], str((1.7*np.double(k)/kmax)), label='z='+str(k*dx[2])+'m')
    ax[1].set_ylabel(r'$\omega$ [m/s$^2$]')
    y_ticks = [ti for ti in ax[1].get_yticks()]
    ax[1].set_yticklabels(y_ticks)
    for label in ax[1].yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    ax[1].set_ylim(0,max)
    ax[0].text(title_pos_x, title_pos_y, 'd) horizontal vorticity', fontsize=15, horizontalalignment='left', bbox=textprops)

    speed_rad = np.sqrt(v_rad[it, :rmax, :] ** 2 + w_rad[it, :rmax, :] ** 2)
    if speed_rad.max() > 0.:
        lw_rad = 5 * (speed_rad+.5) / speed_rad.max()
    else:
        lw_rad = 5 * speed_rad / 1.
    ax[0].streamplot(np.arange(rmax), np.arange(kmax), v_rad[it, :rmax, :kmax].T, w_rad[it, :rmax, :kmax].T,
                     linewidth=lw_rad.T, color='k')

    for ax in axes.flat:
        ax.set_xlim(0, rmax)
        x_ticks = [ti*dx[0] * 1e-3 for ti in ax.get_xticks()]
        ax.set_xticklabels(x_ticks)
    for ax in axes[:3, 0].flat:
        y_ticks = [ti*dx[2] * 1e-3 for ti in ax.get_yticks()]
        ax.set_yticklabels(y_ticks)
        for label in ax.yaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
        ax.set_ylabel('Height  [km]')
    ax = axes[3,0]
    y_ticks = [ti * dx[2] * 1e-3 for ti in ax.get_yticks()]
    ax.set_yticklabels(y_ticks)
    for label in ax.yaxis.get_ticklabels()[0::2]:
        label.set_visible(False)
    ax.set_ylabel('Height  [km]')
    for ax in axes[-1, :]:
        ax.set_xlabel('Radius [km]')

    # axes[-1].legend(loc='upper center', bbox_to_anchor=(1.2, 1.),
    #                 fancybox=True, shadow=True, ncol=1, fontsize=10)
    plt.subplots_adjust(bottom=0.06, right=.95, left=0.1, top=0.95, wspace=0.2, hspace=0.4)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    return








def plot_test(theta, theta_rad, u, v_rad, w, w_rad, vorticity, vort_rad,
              it, jmin, jmax, rmax, kmax, dx, nx, ny,
              fig_name, path_out_figs):
    nlev = 2e2
    ncol = 3
    nrow = 4
    # title_pos_x = 850
    title_pos_x = jmin
    title_pos_y = 43
    textprops = dict(facecolor='white', alpha=0.9, linewidth=0.)

    fig, axes = plt.subplots(nrow, ncol, sharex='col', figsize=(5 * ncol, 2 * nrow))
    for ax in axes[:, 2].flatten():
        ax.plot([0, jmax], [0, 0], '.5', linewidth=0.5)

    ax = axes[0, :]
    # min = 298
    min = np.round(np.amin(theta[:, :kmax]), 1)
    max = 300
    cf = ax[0].contourf(theta[:, :kmax].T, levels=np.linspace(min, max, nlev), cmap=cm_bw_r, extend='max')
    # cbar = plt.colorbar(cf, ax=ax[0], shrink=.75, aspect=12, ticks=np.arange(min, max+0.5, 1.))
    cbar = plt.colorbar(cf, ax=ax[0], shrink=.75, aspect=12, ticks=np.arange(np.ceil(min), max + 0.5, 1.))
    cf = ax[1].contourf(theta_rad[it, :, :kmax].T, levels=np.linspace(min, max, nlev), cmap=cm_bw_r, extend='max')
    cbar = plt.colorbar(cf, ax=ax[1], shrink=.75, aspect=12, ticks=np.arange(np.ceil(min), max + 0.5, 1.))
    for k in range(0,kmax,4):
        ax[2].plot(theta_rad[it, :, k], str(np.double(k) / kmax))
        print('k', k, str(np.double(k) / kmax))
    # ax[2].plot(theta_rad[it,:,1], '.5')
    ax[2].set_ylabel(r'$\theta$ [K]')
    ax[2].set_ylim(min, max)
    y_ticks = [np.int(ti) for ti in ax[1].get_yticks()]
    ax[2].set_yticklabels(y_ticks)
    for label in ax[2].yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    ax[0].text(title_pos_x, title_pos_y, 'a) potential temperature', fontsize=15, horizontalalignment='left',
               bbox=textprops)

    ax = axes[1, :]
    min = -6.
    max = -min
    cf = ax[0].contourf(u[:, :kmax].T, levels=np.linspace(min, max, nlev), cmap=cm_bwr, extend='both')
    cbar = plt.colorbar(cf, ax=ax[0], shrink=0.75, aspect=12, ticks=np.arange(min, max + 0.5, 2))
    # cbar.set_label(cont_var_name, rotation=90)
    cf = ax[1].contourf(v_rad[it, :, :kmax].T, levels=np.linspace(min, max, nlev), cmap=cm_bwr, extend='both')
    cbar = plt.colorbar(cf, ax=ax[1], shrink=0.75, aspect=12, ticks=np.arange(min, max + 0.5, 2))
    ax[2].plot(v_rad[it, :, 0], '0.')
    ax[2].plot(v_rad[it, :, 1], '0.5')
    ax[2].set_ylabel(r'$v_r$ [m/s]')
    # ax[1].set_ylim(-7,max)
    y_ticks = [ti for ti in ax[2].get_yticks()]
    ax[2].set_yticklabels(y_ticks)
    # for label in ax[1].yaxis.get_ticklabels()[1::2]:
    #     label.set_visible(False)
    ax[0].text(title_pos_x, title_pos_y, 'b) radial velocity', fontsize=15, horizontalalignment='left', bbox=textprops)

    ax = axes[2, :]
    min = -2.
    max = -min
    cf = ax[0].contourf(w[:, :kmax].T, levels=np.linspace(min, max, nlev), cmap=cm_bwr, extend='both')
    cbar = plt.colorbar(cf, ax=ax[0], shrink=0.75, aspect=12, ticks=np.arange(min, max + 0.5, 1.))
    cf = ax[1].contourf(w_rad[it, :, :kmax].T, levels=np.linspace(min, max, nlev), cmap=cm_bwr, extend='both')
    cbar = plt.colorbar(cf, ax=ax[1], shrink=0.75, aspect=12, ticks=np.arange(min, max + 0.5, 1.))
    ax[2].plot(w_rad[it, :, 0], '0.')
    ax[2].plot(w_rad[it, :, 1], '0.5')
    ax[2].set_ylabel(r'$w$ [m/s]')
    y_ticks = [np.int(ti) for ti in ax[2].get_yticks()]
    ax[2].set_yticklabels(y_ticks)
    for label in ax[2].yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    ax[0].text(title_pos_x, title_pos_y, 'c) vertical velocity', fontsize=15, horizontalalignment='left',
               bbox=textprops)

    ax = axes[3, :]
    min = -0.15
    max = -min
    cf = ax[0].contourf(vorticity[:, :kmax].T, levels=np.linspace(min, max, nlev), cmap=cm_bwr, extend='both')
    cbar = plt.colorbar(cf, ax=ax[0], shrink=0.75, aspect=12, ticks=np.arange(min, max + 0.02, 0.05))
    cf = ax[1].contourf(vort_rad[it, :, :kmax].T, levels=np.linspace(min, max, nlev), cmap=cm_bwr, extend='both')
    cbar = plt.colorbar(cf, ax=ax[1], shrink=0.75, aspect=12, ticks=np.arange(min, max + 0.02, 0.05))
    ax[2].plot(vort_rad[it, :, 0], '0.')
    ax[2].plot(vort_rad[it, :, 1], '0.5')
    ax[2].set_ylabel(r'$\omega$ [m/s$^2$]')
    y_ticks = [ti for ti in ax[2].get_yticks()]
    ax[2].set_yticklabels(y_ticks)
    for label in ax[2].yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    ax[0].text(title_pos_x, title_pos_y, 'd) vorticity', fontsize=15, horizontalalignment='left', bbox=textprops)

    # speed = np.sqrt(v**2 + w**2)
    speed = np.sqrt(u ** 2 + w ** 2)
    if speed[:, :kmax].max() > 0.:
        lw = 5 * (speed[:, :kmax]+.5) / speed[:, :kmax].max()
    else:
        lw = 5 * speed[:, :kmax] / 1.
    speed_rad = np.sqrt(v_rad[it, :rmax, :kmax] ** 2 + w_rad[it, :rmax, :kmax] ** 2)
    if speed_rad.max() > 0.:
        lw_rad = 5 * speed_rad / speed_rad.max()
    else:
        lw_rad = 5 * speed_rad / 1.
    ax[0].streamplot(np.arange(ny), np.arange(kmax), u[:, :kmax].T, w[:, :kmax].T, density=1., linewidth=lw[:,:].T, color='k')
    ax[1].streamplot(np.arange(rmax), np.arange(kmax), v_rad[it, :rmax, :kmax].T, w_rad[it, :rmax, :kmax].T, linewidth=lw_rad.T, color='k')


    for ax in axes[:, 0].flat:
        ax.set_xlim(jmin, jmax)
        x_ticks = [ti * dx[0] * 1e-3 for ti in ax.get_xticks()]
        ax.set_xticklabels(x_ticks)
        for label in ax.yaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
        ax.set_ylabel('Height  [km]')
    # for label in axes[-1,0].yaxis.get_ticklabels()[1::2]:
    #     label.set_visible(False)
    for n in range(nrow):
        for c in range(2):
            ax = axes[n, c]
            y_ticks = [ti * dx[2] * 1e-3 for ti in ax.get_yticks()]
            ax.set_yticklabels(y_ticks)
            for label in ax.yaxis.get_ticklabels()[1::2]:
                label.set_visible(False)

    for ax in axes[:, 1:3].flat:
        # ax.set_xlim(jmin,jmax)
        ax.set_xlim(0, rmax)
        x_ticks = [ti * dx[0] * 1e-3 for ti in ax.get_xticks()]
        ax.set_xticklabels(x_ticks)
    axes[-1, 0].set_xlabel('distance [km]')
    for ax in axes[-1, 1:]:
        ax.set_xlabel('Radius [km]')

    # axes[-1].legend(loc='upper center', bbox_to_anchor=(1.2, 1.),
    #                 fancybox=True, shadow=True, ncol=1, fontsize=10)
    plt.subplots_adjust(bottom=0.06, right=.95, left=0.1, top=0.95, wspace=0.25, hspace=0.45)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    return



# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()
