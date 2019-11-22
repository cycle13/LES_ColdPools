import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# from matplotlib.widgets import TextBox
import netCDF4 as nc
import argparse
import json as simplejson
import os

label_size = 12
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['lines.markersize'] = 6
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['text.usetex'] = 'true'
plt.rcParams['legend.numpoints'] = 1

execfile('settings.py')

from thermodynamic_functions import thetas_c

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    # parser.add_argument("--tmin")
    # parser.add_argument("--tmax")
    # parser.add_argument("--dx")
    # parser.add_argument("--shift")
    # parser.add_argument("--irmax")
    # parser.add_argument("--nsub")
    args = parser.parse_args()

    res = 50
    run = 'run3'
    dTh = 3
    rstar = 1000
    zstar = 1000
    case = 'dTh'+str(dTh)+'_z'+str(zstar)+'_r'+str(rstar)
    # case_name = args.casename
    case_name = 'ColdPoolDry_single_3D'

    path_root = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/'
    path = os.path.join(path_root, run + '_dx'+str(res)+'m', case)
    path_fields = os.path.join(path, 'fields')
    path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single'

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    ic = nml['init']['ic']
    jc = nml['init']['jc']
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    dx = np.ndarray(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    kmax = np.int(1.2e3 / dx[2])
    print('res: '+str(dx))
    print('CP centre: '+ str(ic)+', '+str(jc))
    print('')

    root_vort = nc.Dataset(os.path.join(path, 'fields_vorticity/field_vort_xz.nc'))
    time_vort = root_vort.groups['fields'].variables['time'][:]
    vorticity_ = root_vort.groups['fields'].variables['vort_xz'][:, :, :kmax]
    root_vort.close()

    root_vrad_2D = nc.Dataset(os.path.join(path, 'fields_v_rad/v_rad.nc'))
    # time_rad = root_vrad_2D.variables['time'][:]
    vrad_2D = root_vrad_2D.variables['v_rad'][:, :, jc, :kmax]
    root_vrad_2D.close()
    root_vrad = nc.Dataset(os.path.join(path, 'data_analysis/stats_radial_averaged.nc'))
    time_rad = root_vrad.groups['timeseries'].variables['time'][:]
    vrad = root_vrad.groups['stats'].variables['v_rad'][:, :, :kmax]
    root_vrad.close()


    for t0 in [900, 1200, 1500]:
    # for t0 in [1200]:
        # print('')
        # print('time: ', time_vort[np.int(t0 / 100)])
        # print('')
        imin = 50
        imax = nx-imin
        fullpath_in = os.path.join(path_fields, str(t0)+'.nc')
        root_field = nc.Dataset(fullpath_in)
        grp = root_field.groups['fields']
        s = grp.variables['s'][:,jc,:kmax]
        w = grp.variables['w'][:,jc,:kmax]
        v = grp.variables['v'][:,jc,:kmax]
        root_field.close()

        theta = thetas_c(s, 0.0)[ic-nx/2:ic+nx/2, :]

        vorticity = vorticity_[np.int(t0/100),:,:]



        fig_name = 'CP_structure_dx' + str(res) + '_' + case + '_t'+str(t0)+'.png'
        cm = plt.cm.get_cmap('rainbow')
        cm_bwr = plt.cm.get_cmap('bwr')
        nlev = 2e2
        ncol = 2
        nrow = 5
        title_pos_x = 460
        title_pos_y = 26
        textprops = dict(facecolor='white', alpha=0.9, linewidth=0.)

        fig, axes = plt.subplots(nrow, ncol, sharex='col', figsize=(5*ncol, 2*nrow))
        ax = axes[0,:]
        min = 298.5
        max = 300
        cf = ax[0].contourf(theta[:,:kmax].T, levels=np.linspace(min, max,nlev), extend='max')
        # plt.colorbar(cf, ax = ax[0])
        cbar = plt.colorbar(cf, ax=ax[0], shrink=.75, aspect=12, ticks=np.arange(min, max+0.5, 0.5))
        ax[1].plot(theta[:,0], 'k')
        ax[1].set_ylabel(r'$\theta$ [K]')
        ax[1].set_ylim(min,max)
        ax[0].text(title_pos_x, title_pos_y, 'potential temperature', fontsize=15, horizontalalignment='center', bbox=textprops)
        # ax[0].text(0.73, 0.9, textstr, transform=ax.transAxes, fontsize=24,
        #         verticalalignment='top', bbox=props)

        ax = axes[1, :]
        min = -0.25
        max = -min
        cf = ax[0].contourf(v[:, :kmax].T, levels=np.linspace(min, max, nlev), cmap=cm_bwr)
        # plt.colorbar(cf, ax=ax[0])
        cbar = plt.colorbar(cf, ax=ax[0], shrink=0.75, aspect=12, ticks=np.arange(-0.2, 0.2+0.1, 0.1))
        # cbar.set_label(cont_var_name, rotation=90)
        ax[1].plot(v[:, 1], 'k')
        ax[1].set_ylabel(r'$v_r$ [m/s]')
        ax[0].text(title_pos_x, title_pos_y, 'radial velocity', fontsize=15, horizontalalignment='center', bbox=textprops)

        ax = axes[2, :]
        min = -2.
        max = -min
        cf = ax[0].contourf(w[:, :kmax].T, levels=np.linspace(-2., 2., nlev), cmap=cm_bwr)
        cbar = plt.colorbar(cf, ax=ax[0], shrink=0.75, aspect=12, ticks=np.arange(min, max+0.5, 0.5))
        ax[1].plot(w[:, 0], 'k')
        ax[1].set_ylabel(r'$w$ [m/s]')
        ax[0].text(title_pos_x, title_pos_y, 'vertical velocity', fontsize=15, horizontalalignment='center', bbox=textprops)

        ax = axes[3, :]
        min = -0.1
        max = -min
        cf = ax[0].contourf(vorticity[:, :kmax].T, levels=np.linspace(-0.1, 0.1, nlev))
        cbar = plt.colorbar(cf, ax=ax[0], shrink=0.75, aspect=12, ticks=np.arange(min, max+0.02, 0.02))
        ax[1].plot(vorticity[:, 0], 'k')
        ax[1].set_ylabel(r'$\omega$ [m/s$^2$]')
        ax[0].text(title_pos_x, title_pos_y, 'vorticity', fontsize=15, horizontalalignment='center', bbox=textprops)

        ax = axes[4, :]
        speed = np.sqrt(v**2 + w**2)
        if speed[:, :kmax].max() > 0.:
            lw = 5 * speed[:, :kmax] / speed[:, :kmax].max()
        else:
            lw = 5 * speed[:, :kmax] / 1.
        cf = ax[0].contourf(speed[:, :kmax].T, levels=np.linspace(0., 2., nlev))
        plt.colorbar(cf, ax=ax[0])
        # ax[1].streamplot(np.arange(nx), np.arange(kmax), v[:, :kmax].T, w[:, :kmax].T,
        #                  color='k', density=1.5, linewidth=lw[:, :].T)
        # ax[1].streamplot(np.arange(nx), np.arange(kmax), v[:, :kmax].T, w[:, :kmax].T, color='w')
        ax[1].streamplot(np.arange(nx), np.arange(kmax), w[:, :kmax].T, v[:, :kmax].T, density=1.5, linewidth=lw[:, :].T)
        # ax[1].plot(vorticity[:, 0], 'k')
        ax[0].text(title_pos_x, title_pos_y, 'streamlines', fontsize=15, horizontalalignment='center', bbox=textprops)

        for ax in axes[:,0].flat:
            ax.set_xlim(imin,imax)
            # print('x ticks:', ax.get_xticks())
            x_ticks = [ti*dx[0]*1e-3 for ti in ax.get_xticks()]
            # print('x ticks:', x_ticks)
            ax.set_xticklabels(x_ticks)
            # y_ticks = [np.round(ti*dx[2]*1e-3,1) for ti in ax.get_yticks()]
            y_ticks = [ti*dx[2]*1e-3 for ti in ax.get_yticks()]
            ax.set_yticklabels(y_ticks)
            # for label in ax.yaxis.get_ticklabels()[1::2]:
            #     label.set_visible(dFalse)
            ax.set_ylabel('height  [km]')
        for label in axes[-1,0].yaxis.get_ticklabels()[1::2]:
            # print('lable', label)
            label.set_visible(False)
        for ax in axes[:,1].flat:
            ax.set_xlim(imin,imax)
            x_ticks = [ti*dx[0]*1e-3 for ti in ax.get_xticks()]
            ax.set_xticklabels(x_ticks)
            y_ticks = [ti for ti in ax.get_yticks()]
            ax.set_yticklabels(y_ticks)
            # for label in ax.yaxis.get_ticklabels()[1::2]:
            #     label.set_visible(False)
        for ax in axes[-1,:]:
            ax.set_xlabel('distance [km]')

        # axes[-1].legend(loc='upper center', bbox_to_anchor=(1.2, 1.),
        #                 fancybox=True, shadow=True, ncol=1, fontsize=10)
        plt.subplots_adjust(bottom=0.075, right=.95, left=0.1, top=0.9, wspace=0.3, hspace=0.5)
        # plt.suptitle('z=' + str(k0 * dx[2]) + 'm')
        fig.savefig(os.path.join(path_out_figs, fig_name))
        plt.close(fig)


    return



# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()
