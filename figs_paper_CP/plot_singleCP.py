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
label_size = 13
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 6
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['text.usetex'] = 'true'
plt.rcParams['legend.numpoints'] = 1

# tracer file:
# - col=0:      timestep of simulation
# - col=1:      age of CP (timestep since beginning of first CP)
# - col=2:      tracer id
# - col=3:      CP id
# - col=4,5:    tracer position (coordinates)
# - col=6,7:    tracer position (coordinates), rounded
# - col=8,9:    tracer radius & angle
# - col=10,11:  u, v wind components
# - col=12,13:  radial & tangential wind components
# - col=14,15:  x-, y-distance of tracer to center
# - col=16,17:  CP center (xc,yc) (COM)

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

    res = 25
    run = 'run4'
    dTh = 3
    rstar = 1000
    zstar = 1000
    case = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
    # case_name = args.casename
    case_name = 'ColdPoolDry_single_3D'

    time_range = [900, 1200, 1500]
    imin = 100
    k0 = 0

    path_root = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/'
    path = os.path.join(path_root, run + '_dx' + str(res) + 'm', case)
    print(path)
    path_fields = os.path.join(path, 'fields')
    # path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single'
    path_out_figs = '/nbi/home/meyerbe/paper_CP_single/figs_run4'
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    ic = nml['init']['ic']
    jc = nml['init']['jc']
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    dx = np.ndarray(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    # kmax = np.int(1.0e3 / dx[2])
    dt_fields = nml['fields_io']['frequency']
    print('res: ' + str(dx))
    print('CP centre: ' + str(ic) + ', ' + str(jc))
    print('')

    ''' radial velocity '''
    root_vrad_2D = nc.Dataset(os.path.join(path, 'fields_v_rad/v_rad.nc'))
    # time_rad = root_vrad_2D.variables['time'][:]
    vrad_2D_ = root_vrad_2D.variables['v_rad'][:, :, :, k0]  # v_rad[nt,nx,ny,nz]
    root_vrad_2D.close()
    # root_vrad = nc.Dataset(os.path.join(path, 'data_analysis/stats_radial_averaged.nc'))
    # time_rad = root_vrad.groups['timeseries'].variables['time'][:]
    # vrad = root_vrad.groups['stats'].variables['v_rad'][:, :, k0]         # v_rad[nt,nr,nz]
    # root_vrad.close()

    # ''' vorticity '''
    # root_vort = nc.Dataset(os.path.join(path, 'fields_vorticity/field_vort_yz.nc'))
    # time_vort = root_vort.groups['fields'].variables['time'][:]
    # vorticity_ = root_vort.groups['fields'].variables['vort_yz'][:, :, :kmax]     # vort_yz[t,y,z]
    # root_vort.close()


    # Figure with potential temperature, vertical velocity, radial velocity ?
    # zoom in to see clefts

    fig_name = 'CP_crosssection_dx' + str(res) + '_' + case + '.png'
    nlev = 2e2
    ncol = 3
    nrow = 4
    # # title_pos_x = 850





    fig, axes = plt.subplots(nrow, ncol, sharex='col', figsize=(5 * ncol, 2 * nrow))

    # jmin_range = [150, 100, 100]
    for i, t0 in enumerate(time_range):
        it = np.int(t0/dt_fields)
        print('time: ', t0)

        fullpath_in = os.path.join(path_fields, str(t0) + '.nc')
        root_field = nc.Dataset(fullpath_in)
        grp = root_field.groups['fields']
        s = grp.variables['s'][:, :, k0]
        w = grp.variables['w'][:, :, k0]
        # v = grp.variables['v'][:,:,k0]
        # u = grp.variables['u'][:, :, k0]
        root_field.close()
        theta = thetas_c(s, 0.0)#[ic - nx / 2:ic + nx / 2, :]
        # vorticity = vorticity_[it, :, :]
        vrad_2D = vrad_2D_[it, :, :]

        ax = axes[0, :]
        # min = 298
        min = np.round(np.amin(theta[:, :]))
        max = 300
        cf = ax[i].contourf(theta[:, :].T, levels=np.linspace(min, max, nlev), cmap=cm_bw_r, extend='max')

        ax = axes[1, :]
        min = np.round(np.amin(w[:, :]))
        max = np.round(np.amax(w[:, :]))
        cf = ax[i].contourf(w[:, :].T, levels=np.linspace(min, max, nlev), cmap=cm_bw_r, extend='max')

        ax = axes[2, :]
        min = np.round(np.amin(vrad_2D[:, :]))
        max = np.round(np.amax(vrad_2D[:, :]))
        max = 300
        cf = ax[i].contourf(vrad_2D[:, :].T, levels=np.linspace(min, max, nlev), cmap=cm_bw_r, extend='max')

    textprops = dict(facecolor='white', alpha=0.9, linewidth=0.)
    title_pos_x = 10
    title_pos_y = 0
    txt = 'a) potential temperature'
    axes[0,0].text(title_pos_x, title_pos_y, txt, fontsize=15, horizontalalignment='left', bbox=textprops)
    txt = 'b) vertical velocity'
    # title_pos_y = 0
    axes[1,0].text(title_pos_x, title_pos_y, txt, fontsize=15, horizontalalignment='left', bbox=textprops)
    txt = 'c) radial velocity'
    # title_pos_y = 0
    axes[2,0].text(title_pos_x, title_pos_y, txt, fontsize=15, horizontalalignment='left', bbox=textprops)

    # axes[-1].legend(loc='upper center', bbox_to_anchor=(1.2, 1.),
    #                 fancybox=True, shadow=True, ncol=1, fontsize=10)
    plt.subplots_adjust(bottom=0.06, right=.95, left=0.1, top=0.95, wspace=0.25, hspace=0.45)
    fig.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)

    return


# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()