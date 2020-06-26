import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import netCDF4 as nc
import argparse
import json as simplejson
import os

execfile('settings.py')


def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("path_root")
    parser.add_argument("--casename")
    parser.add_argument("--dTh", type=int)
    parser.add_argument("--zparams", nargs='+', type=int)
    parser.add_argument('--rparams', nargs='+', type=int)
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    parser.add_argument("--dz")
    args = parser.parse_args()

    if args.dz:
        dz=np.int(args.dz)
    else:
        dz = 100

    zstar_ref = 1000
    rstar_ref = 1000
    if dz == 100:
        path_ref = '/nbi/ac/conv4/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/dTh3_z1000_r1000'
    elif dz == 50:
        path_ref = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run3_dx50m/dTh3_z1000_r1000'
    print('')
    print('reference case: '+path_ref)
    print('')
    kmax = 50

    # root = nc.Dataset(os.path.join(path_ref, 'fields', '0.nc'), 'r')
    # s = root.groups['fields'].variables['s'][:,:,:kmax]
    # [nx,ny,nz] = s.shape
    # root.close()
    #
    root = nc.Dataset(os.path.join(path_ref, 'stats', 'Stats.ColdPoolDry_single_3D.nc'))
    rho0_ref = root.groups['reference'].variables['rho0'][:kmax]
    z_ref = root.groups['reference'].variables['z'][:kmax]
    root.close()

    g = 9.80665

    # theta_ref = theta_s(s)
    # th_env_ref = np.mean(np.mean(theta_ref[:50,:,:], axis=0), axis=0)
    # th_g_ref = th_env[0]
    # th_g_ = 300.
    #
    #
    # ''' testing '''
    # # PE1 = np.zeros(kmax, dtype=np.double)
    # # PE2 = np.zeros(kmax, dtype=np.double)
    # # PE3 = np.zeros(kmax, dtype=np.double)
    # # PE1_ = np.zeros(kmax, dtype=np.double)
    # # PE2_ = np.zeros(kmax, dtype=np.double)
    # # PE3_ = np.zeros(kmax, dtype=np.double)
    # # PE1_a = np.zeros(kmax, dtype=np.double)
    # # PE2_a = np.zeros(kmax, dtype=np.double)
    # # PE3_a = np.zeros(kmax, dtype=np.double)
    # # PE1_b = np.zeros(kmax, dtype=np.double)
    # # PE2_b = np.zeros(kmax, dtype=np.double)
    # # PE3_b = np.zeros(kmax, dtype=np.double)
    # # # PE1[k] = np.sum(theta[:, :, k].reshape(nx ** 2, kmax), axis=0)
    # # for k in range(kmax):
    # #     PE1_[k] = np.sum((th_env[k] - theta[:, :, k]) / th_env[k])
    # #     PE2_[k] = np.sum((th_g - theta[:, :, k]) / th_g)
    # #     PE3_[k] = np.sum((th_g_ - theta[:, :, k]) / th_g_)
    # #     PE1_a[k] = np.sum(z[k]*(th_env[k] - theta[:, :, k])/th_env[k])
    # #     PE2_a[k] = np.sum(z[k]*(th_g - theta[:, :, k])/th_g)
    # #     PE3_a[k] = np.sum(z[k]*(th_g_ - theta[:, :, k])/th_g_)
    # #     PE1_b[k] = np.sum(rho0[k]*(th_env[k] - theta[:, :, k]) / th_env[k])
    # #     PE2_b[k] = np.sum(rho0[k]*(th_g - theta[:, :, k]) / th_g)
    # #     PE3_b[k] = np.sum(rho0[k]*(th_g_ - theta[:, :, k]) / th_g_)
    # #     PE1[k] = np.sum(g*rho0[k]*z[k] * (th_env[k] - theta[:, :, k]) / th_env[k])
    # #     PE2[k] = np.sum(g*rho0[k]*z[k] * (th_g - theta[:, :, k]) / th_g)
    # #     PE3[k] = np.sum(g*rho0[k]*z[k] * (th_g_ - theta[:, :, k]) / th_g_)
    # #
    # #
    # # fig,(ax0, ax1, ax2, ax3, ax4, ax5) = plt.subplots(1,6, figsize=(14,5))
    # # ax0.imshow(theta[:,:,0])
    # # ax1.plot(th_env, np.arange(kmax))
    # # ax1.plot(th_g*np.ones(kmax), np.arange(kmax))
    # # ax1.plot(th_g_*np.ones(kmax), np.arange(kmax))
    # # ax2.plot(PE1_, np.arange(kmax), '-d', label='PE1 (th0[z])')
    # # ax2.plot(PE2_, np.arange(kmax), '-x', label='PE2 (th0[0])')
    # # ax2.plot(PE3_, np.arange(kmax), '-', label='PE3 (300K)')
    # # ax2.set_title('sum[ (th-th0)/th0 ]')
    # # ax2.legend()
    # # # ax.set_title('sum[ (th-th0)/th0 ]')
    # # ax3.plot(PE1_a, np.arange(kmax), '-d', label='PE1')
    # # ax3.plot(PE2_a, np.arange(kmax), '-x', label='PE2')
    # # ax3.plot(PE3_a, np.arange(kmax), '-', label='PE3')
    # # ax3.set_title('sum[ z*(th-th0)/th0 ]')
    # # ax4.plot(PE1_b, np.arange(kmax), '-d', label='PE1')
    # # ax4.plot(PE2_b, np.arange(kmax), '-x', label='PE2')
    # # ax4.plot(PE3_b, np.arange(kmax), '-', label='PE3')
    # # ax4.set_title('sum[ rho0*(th-th0)/th0 ]')
    # # ax5.plot(PE1, np.arange(kmax), '-d', label='PE1')
    # # ax5.plot(PE2, np.arange(kmax), '-x', label='PE2')
    # # ax5.plot(PE3, np.arange(kmax), '-', label='PE3')
    # # ax5.set_title('PE')
    # # for ax in [ax2, ax3, ax4, ax5]:
    # #     ax.plot([0,0], [0,kmax], 'k', linewidth=1)
    # # plt.savefig('/nbi/home/meyerbe/to_cp_to_laptop/test_energy.png')
    #
    #
    #
    #
    ''' looping '''
    path_root = args.path_root
    if args.dTh:
        dTh = np.int(args.dTh)
    else:
        dTh = 3
    if args.zparams:
        z_params = args.zparams
    else:
        z_params = [1000]
    if args.rparams:
        r_params = args.rparams
    else:
        r_params = [1000]
    print('z*: ', z_params)
    print('r*: ', r_params)
    n_params = len(z_params)

    if args.casename:
        case_name = args.casename
    else:
        case_name = 'ColdPoolDry_single_3D'
    id0 = 'dTh' + str(dTh) + '_z' + str(z_params[0]) + '_r' + str(r_params[0])
    nml = simplejson.loads(open(os.path.join(path_root, id0, case_name + '.in')).read())
    global dx, dV, gw
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = np.zeros(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx[0] * dx[1] * dx[2]

    ''' determine file range '''
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 0
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = tmin
    timerange = np.arange(tmin, tmax + 100, 100)
    nt = len(timerange)
    th_g = 300.

    PE_computed_1 = np.zeros((n_params, nt, kmax))
    PE_computed_2 = np.zeros((n_params, nt, kmax))
    PE_ref_1 = np.zeros((nt, kmax))
    PE_ref_2 = np.zeros((nt, kmax))
    for it, t0 in enumerate(timerange):
        print('-- t=' + str(t0) + '(it=' + str(it) + ') --')
        root = nc.Dataset(os.path.join(path_ref, 'fields', '0.nc'), 'r')
        s = root.groups['fields'].variables['s'][:, :, :kmax]
        root.close()
        theta_ref = theta_s(s)
        th_env_ref = np.mean(np.mean(theta_ref[:50, :, :], axis=0), axis=0)
        for k in range(kmax):
            PE_ref_1[it, k] = np.sum(g * rho0_ref[k] * z_ref[k] * (th_env_ref[k] - theta_ref[:, :, k]) / th_env_ref[k])
            PE_ref_2[it, k] = np.sum(g * rho0_ref[k] * z_ref[k] * (th_g - theta_ref[:, :, k]) / th_g)

    for istar in range(n_params):
        zstar = z_params[istar]
        rstar = r_params[istar]
        ID = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        print(ID)
        path_in = os.path.join(path_root, ID)
        root = nc.Dataset(os.path.join(path_in, 'stats', 'Stats.ColdPoolDry_single_3D.nc'))
        rho0 = root.groups['reference'].variables['rho0'][:]
        z = root.groups['reference'].variables['z'][:]
        root.close()
        for it, t0 in enumerate(timerange):
            print('-- t=' + str(t0) + '(it=' + str(it) + ') --')
            root = nc.Dataset(os.path.join(path_in, 'fields', str(t0)+'.nc'), 'r')
            s = root.groups['fields'].variables['s'][:, :, :kmax]
            root.close()
            th_s = theta_s(s)
            del s
            th_env = np.mean(np.mean(th_s[:50,:,:], axis=0), axis=0)
            for k in range(kmax):
                PE_computed_1[istar, it, k] = np.sum(g*rho0[k]*z[k] * (th_env[k] - th_s[:, :, k]) / th_env[k])
                PE_computed_2[istar, it, k] = np.sum(g*rho0[k]*z[k] * (th_g - th_s[:, :, k]) / th_g)


    fig,(ax0, ax1, ax2, ax3) = plt.subplots(1,4, figsize=(14,5))
    for istar in range(n_params):
        lbl = 'z*'+str(z_params[istar])+', r*'+str(r_params[istar])
        ax0.plot(PE_computed_1[istar,0,:], np.arange(kmax), '-d', color=colorlist[istar], label='PE1 (th0[z]), '+lbl)
        ax0.plot(PE_computed_2[istar,0,:], np.arange(kmax), '-', color=colorlist[istar], label='PE1 (300K), '+lbl)
        ax2.plot(timerange, np.sum(PE_computed_1[istar,:,:], axis=1), '-d', color=colorlist[istar], label='PE(th0[z]), '+lbl)
        ax2.plot(timerange, np.sum(PE_computed_2[istar,:,:], axis=1), '-', color=colorlist[istar], label='PE(300K), '+lbl)
        ax3.plot(timerange, np.sum(PE_computed_1[istar, :, :], axis=1)/np.sum(PE_ref_1[:,:], axis=1), '-d', color=colorlist[istar], label='PE(th0[z]), ' + lbl)
        ax3.plot(timerange, np.sum(PE_computed_2[istar, :, :], axis=1)/np.sum(PE_ref_2[:,:], axis=1), '-', color=colorlist[istar], label='PE(300K), ' + lbl)
    ax2.plot(timerange, np.sum(PE_ref_1[:,:], axis=1), '-', color='r', label='PE0(300K)')
    ax2.plot(timerange, np.sum(PE_ref_2[:,:], axis=1), '-', color='r', label='PE0(300K)')
    ax0.set_title('PE(t=0,z)')
    ax0.set_xlabel('PE [J]')
    ax0.set_ylabel('Height z [m]')
    lbl =  np.round(np.sum(PE_computed_1[:, 0, :], axis=-1)/np.sum(PE_ref_1[0,:]), 3)
    ax1.plot(z_params, np.sum(PE_computed_1[:, 0, :], axis=-1), '-d', color='k', label='PE(th0[z]), ' + str(lbl))
    ax1.plot(zstar_ref, np.sum(PE_ref_1[0, :], axis=-1), '-', color='r', label='PE0(th0[z])')
    lbl = np.round(np.sum(PE_computed_2[:, 0, :], axis=-1)/np.sum(PE_ref_2[0,:]), 3)
    ax1.plot(z_params, np.sum(PE_computed_2[:, 0, :], axis=-1), '-', color='k', label='PE(300K), '+str(lbl))
    ax1.plot(zstar_ref, np.sum(PE_ref_2[0, :], axis=-1), '-', color='r', label='PE0(300K)')
    ax1.set_title('PE(t=0)')
    ax1.set_xlabel('z* [m]')
    ax1.set_ylabel('PE [J]')
    ax1.legend(loc='best', fontsize=8)
    ax2.set_title('PE TS')
    ax2.set_xlabel('time [s]')
    ax2.set_ylabel('PE [J]')
    ax3.set_title('PE TS')
    ax3.set_xlabel('time [s]')
    ax3.set_ylabel('PE [J]')
    ax3.set_ylim(0.,1.5)
    # # ax1.legend()
    ax3.legend(loc='upper center', bbox_to_anchor=(1.3, 1.), fontsize = 9)
    for ax in [ax0]:
        ax.plot([0,0], [0,kmax], 'k', linewidth=1)
    ax2.plot([tmin, tmax],[0, 0], 'k', linewidth=1)
    plt.subplots_adjust(left=0.05, right=.85, bottom=0.1, top=0.9, hspace=0.2, wspace=0.3)
    print('saving: '+'/nbi/home/meyerbe/to_cp_to_laptop/test_energy_PE.png')
    plt.savefig('/nbi/home/meyerbe/to_cp_to_laptop/test_energy_PE.png')

    return



# ----------------------------------
def theta_s(s):
    # T_tilde = 298.15
    # thetas_c(s, qt){
    #     return T_tilde * exp((s - (1.0 - qt) * sd_tilde - qt * sv_tilde) / cpm_c(qt));
    # }
    # cpm_c(qt){
    #     return (1.0 - qt) * cpd + qt * cpv;
    # }
    # parameters from pycles
    T_tilde = 298.15
    sd_tilde = 6864.8
    cpd = 1004.0
    th_s = T_tilde * np.exp( (s - sd_tilde)/cpd )
    return th_s


if __name__ == '__main__':
    main()

