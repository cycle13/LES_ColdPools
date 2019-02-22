import netCDF4 as nc
import argparse
import os, sys
import numpy as np
import json as  simplejson
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import matplotlib.pylab as plt


label_size = 15
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
# plt.rcParams['xtick.top']='False'
plt.rcParams['legend.fontsize'] = label_size
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['grid.linewidth'] = 20
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['xtick.minor.width'] = 1.
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['ytick.minor.width'] = 1.
plt.rcParams['pdf.fonttype'] = 42         # Output Type 3 (Type3) or Type 42 (TrueType)

def main():
    path_r1km = '/nbi/ac/cond2/meyerbe/ColdPools_dry/3D_sfc_fluxes_off/double_3D_noise/dTh3K_z2000_r1000_double_sep5r/'
    path_r2km = '/nbi/ac/cond2/meyerbe/ColdPools_dry/3D_sfc_fluxes_off/double_3D_noise/dTh3K_z2000_r2000_double/'
    imin = 70
    imax = 177
    times = [0, 200, 300, 400, 500, 700, 900, 1200, 1400]
    plot_double_collision(path_r1km, times, imin, imax)
    imin = 78
    imax = 200
    times = [0, 200, 400, 500, 600, 800,1000,1400]
    plot_double_collision(path_r2km, times, imin, imax)


    path = '/nbi/ac/cond2/meyerbe/ColdPools_dry/3D_sfc_fluxes_off/triple_3D_noise/dTh3K_triple/'
    plot_streamlines_xy_collision()

    return









def plot_double_collision(path_r1km, times, imin, imax):
    case_name = 'ColdPoolDry_double_3D'

    nml_r1km = simplejson.loads(open(os.path.join(path_r1km, case_name + '.in')).read())
    dx = nml_r1km['grid']['dx']
    dz = nml_r1km['grid']['dz']

    file = nc.Dataset(os.path.join(path_r1km, 'fields_merged', 'fields_allt_xz_j150.nc'))
    theta = file.variables['theta'][:,:,:]
    time_array = file.variables['time'][:]
    file.close()

    fig_name = 'crosssections.png'
    fullpath_out = os.path.join(path_r1km, fig_name)
    print fullpath_out



    nt = len(times)

    kmax = 30
    min = 296.9
    max = 300.1
    lvls = np.linspace(min, max, 1e2)
    fig, axes = plt.subplots(nrows=1, ncols=nt, figsize=(nt*5,4), sharey='all')
    for i, t0 in enumerate(times):
        it = np.where(time_array == t0)[0][0]
        ax = axes[i]
        # cmaps: cm.bone, cm.bwr, cm.gnuplot
        cf = ax.contourf(theta[it, imin:imax,:kmax].T, levels=lvls, cmap = cm.bone)
        ax.set_title('t='+str(t0)+'s')
    cbar_ax = fig.add_axes([0.92, 0.18, 0.012, 0.7])
    fig.colorbar(cf, cax=cbar_ax, ticks=np.arange(np.floor(min), np.floor(max)+1, 1))
    ax = axes[0]
    ax.set_ylabel('height z  [km]')
    y_ticks = ax.get_yticks()*dz * 1e-3
    ax.set_yticklabels(y_ticks)

    x_ticks = axes[0].get_xticks()
    for i,ti in enumerate(x_ticks):
        x_ticks[i] = np.round(ti * dx * 1e-3, 0)
    for i in range(nt):
        axes[i].set_xlabel('x  [km]')
        axes[i].tick_params(top=False, right=False)
        axes[i].set_xticklabels(x_ticks)

    plt.subplots_adjust(bottom=0.15, right=.9, left=0.05, top=0.9, wspace=0.1)
    fig.savefig(fullpath_out)
    plt.close(fig)
    return



def plot_double_collision(path_r1km, times, imin, imax):
    case_name = 'ColdPoolDry_double_3D'

    nml_r1km = simplejson.loads(open(os.path.join(path_r1km, case_name + '.in')).read())
    dx = nml_r1km['grid']['dx']
    dz = nml_r1km['grid']['dz']

    file = nc.Dataset(os.path.join(path_r1km, 'fields_merged', 'fields_allt_xz_j150.nc'))
    theta = file.variables['theta'][:,:,:]
    time_array = file.variables['time'][:]
    file.close()

    fig_name = 'crosssections.png'
    fullpath_out = os.path.join(path_r1km, fig_name)
    print fullpath_out



    nt = len(times)

    kmax = 30
    min = 296.9
    max = 300.1
    lvls = np.linspace(min, max, 1e2)
    fig, axes = plt.subplots(nrows=1, ncols=nt, figsize=(nt*5,4), sharey='all')
    for i, t0 in enumerate(times):
        it = np.where(time_array == t0)[0][0]
        ax = axes[i]
        # cmaps: cm.bone, cm.bwr, cm.gnuplot
        cf = ax.contourf(theta[it, imin:imax,:kmax].T, levels=lvls, cmap = cm.bone, norm=LogNorm(min, maxs))
        ax.set_title('t='+str(t0)+'s')
    cbar_ax = fig.add_axes([0.92, 0.18, 0.012, 0.7])
    # fig.colorbar(cf, cax=cbar_ax, ticks=np.arange(np.log(np.floor(min)), np.log(np.floor(max))+1, 1))
    fig.colorbar(cf, cax=cbar_ax)
    ax = axes[0]
    ax.set_ylabel('height z  [km]')
    y_ticks = ax.get_yticks()*dz * 1e-3
    ax.set_yticklabels(y_ticks)

    x_ticks = axes[0].get_xticks()
    for i,ti in enumerate(x_ticks):
        x_ticks[i] = np.round(ti * dx * 1e-3, 0)
    for i in range(nt):
        axes[i].set_xlabel('x  [km]')
        axes[i].tick_params(top=False, right=False)
        axes[i].set_xticklabels(x_ticks)

    plt.subplots_adjust(bottom=0.15, right=.9, left=0.05, top=0.9, wspace=0.1)
    fig.savefig(fullpath_out)
    plt.close(fig)
    return


if __name__ == '__main__':
    main()