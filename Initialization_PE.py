import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import netCDF4 as nc
import argparse
import json as simplejson
import os
from scipy.integrate import odeint
import sys

from thermodynamic_functions import thetas_c
# import thermodynamic_profiles
# from thermodynamic_profiles import alpha_c

def main():
    cm_rain = plt.cm.get_cmap('rainbow')
    cm_rep = plt.cm.get_cmap('prism')
    cm = cm_rain

    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--dx")
    parser.add_argument("--dTh_min")
    parser.add_argument("--dTh_max")
    args = parser.parse_args()

    path_out = './figs_Initialization/'

    'surface values'
    Th_g = 300.0  # temperature for neutrally stratified background (value from Soares Surface)
    if args.dx:
        dx_ = np.int(args.dx)
    else:
        dx_ = 50
    marg = 200.
    z_half = define_geometry(dx_, marg)
    set_parameters()
    print('dx: '+str(dx_))


    ''' Parameter range '''
    if args.dTh_min:
        dTh_min = np.int(args.dTh_min)
    else:
        dTh_min = 2
    if args.dTh_max:
        dTh_max = np.int(args.dTh_max)
    else:
        dTh_max = 4
    dTh_range = np.arange(dTh_min, dTh_max+1)
    dTh_range = [2, 3, 4]
    # for dx=100m, no matches for r*>4200m;
    # for dx=50m, matches for r*=..., 3600, 4900, 5000; no matches for r<400m
    rstar_min = 5e2
    rstar_max = 25e2
    zstar_min = 5e2
    zstar_max = 25e2
    rstar_range = np.arange(rstar_min, rstar_max+100, 100)
    zstar_range = np.arange(zstar_min, zstar_max+100, 1e2)
    n_thermo = len(dTh_range)
    n_geom_z = len(zstar_range)
    n_geom_r = len(rstar_range)
    print('dTh: ' + str(dTh_range))
    print('zstar: ' + str(zstar_range))
    print('rstar: ' + str(rstar_range))
    print('n_geom_z: ' + str(n_geom_z))
    print('n_geom_r: ' + str(n_geom_r))

    ''' PE scaling-range '''
    # scaling = [2**(-1), 2**0, 2**1, 2**2, 2**3]
    scaling = [2**0]
    print('scaling: ' + str(scaling))
    print('')

    ''' Reference case '''
    dTh_ref = 3
    rstar_ref = 1000
    zstar_ref = 1000
    if dx_== 100:
        path_ref = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/dTh3_z1000_r1000/'
    elif dx_==50:
        path_ref = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run3_dx50m/dTh3_z1000_r1000_nz240/'
    case_name = 'ColdPoolDry_single_3D'
    print('Reference case: ' + path_ref)
    rootgrp = nc.Dataset(os.path.join(path_ref, 'stats', 'Stats.' + case_name + '.nc'))
    rho0_stats = rootgrp.groups['reference'].variables['rho0'][:]
    # alpha0_stats = rootgrp.groups['reference'].variables['alpha0'][:]
    zhalf_stats = rootgrp.groups['reference'].variables['z'][:]
    # rho_unit = rootgrp.groups['reference'].variables['rho0'].units
    rootgrp.close()
    print('')
    print('z-arrays:')
    print('z_half:', z_half[:20])
    print('z_half_stats:', zhalf_stats[:20])
    print('')


    ''' reference PE '''
    print('Reference case: ')
    z_max_arr, theta_z = compute_envelope(dTh_ref, Th_g, rstar_ref, zstar_ref, marg, z_half)
    PE_ref = compute_PE(theta_z, Th_g, zstar_ref, rho0_stats, zhalf_stats, z_half)
    PE_ref_approx = compute_PE_density_approx(dTh_ref, zstar_ref, rstar_ref, z_half)
    print('PE_ref: ', PE_ref)
    # test reference numerically
    rootgrp_field = nc.Dataset(os.path.join(path_ref, 'fields', '0.nc'))
    s0 = rootgrp_field.groups['fields'].variables['s'][:,:,:]
    rootgrp_field.close()
    ic_ = s0.shape[0] / 2
    jc_ = s0.shape[1] / 2
    theta = thetas_c(s0, 0.0)[ic_-nx/2:ic_+nx/2, jc_-ny/2:jc_+ny/2,:nz]
    PE_ref_num = compute_PE(theta, Th_g, zstar_ref, rho0_stats, zhalf_stats, z_half)
    del s0
    print('  (from envel.) PE_ref =        ' + str(PE_ref))
    print('  (from field)  PE_ref_num =    ' + str(PE_ref_num))
    print('   difference: '+str((PE_ref-PE_ref_num)/PE_ref))
    print('')
    # print('   (upper limit)  PE_up = ' + str(PE_ref_up))
    # print('   (lower limit)  PE_low = ' + str(PE_ref_low))
    # print('   difference: '+str((PE_ref-PE_ref_up)/PE_ref))
    # print('   difference: '+str((PE_ref-PE_ref_low)/PE_ref))
    print('')


    for iTh, dTh in enumerate(dTh_range):
        rstar_guess = np.sqrt(PE_ref_approx / (dTh * zstar_max ** 2))
        rstar_min = rstar_guess - dx
        rstar_max = np.minimum(np.sqrt(PE_ref_approx / (dTh * zstar_min ** 2)) + 2*dx, 2500)
        print('rstar', rstar_guess, rstar_min, rstar_max)
        print('')
        rstar_range = dx * np.arange(np.floor(rstar_min/dx), np.ceil(rstar_max/dx)+1)
        n_geom_r = len(rstar_range)
        PE = np.zeros((n_geom_z, n_geom_r))
        diff = np.zeros((3, n_thermo, n_geom_z))
        diff[2, :, :] = 2*PE_ref*np.ones((n_thermo, n_geom_z))

        print('dTh: ', dTh, 'r*-range', rstar_range)
        fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(11,5))
        ax2.fill_between([rstar_range[0], rstar_range[-1]], 0.9, 1.1, color='gray', alpha=0.2)
        ax2.fill_between([rstar_range[0], rstar_range[-1]], 0.95, 1.05, color='gray', alpha=0.4)

        for iz, zstar in enumerate(zstar_range):
            PE_guess = 0.0
            for ir, rstar in enumerate(rstar_range):
                print('z*: ', zstar, 'r*:', rstar)

                z_max_arr, theta = compute_envelope(dTh, rstar, zstar, Th_g)
                PE[iz, ir] = compute_PE(theta, Th_g, z_max_arr, zstar, rho0_stats, zhalf_stats)

                diff_aux = np.abs(PE[iz, ir] - PE_ref)
                if diff_aux < np.abs(diff[2, iTh, iz] - PE_ref):
                    diff[0, iTh, iz] = ir
                    diff[1, iTh, iz] = rstar
                    diff[2, iTh, iz] = PE[iz, ir]

                rstar_guess = np.sqrt(PE_ref_approx / (dTh * zstar ** 2))
                if rstar == rstar_guess:
                    PE_guess = PE[iz, ir]

            ax1.plot(rstar_range, PE[iz, :], '-o', color=cm(np.double(iz)/n_geom_z))
            ax2.plot(rstar_range, PE[iz, :]/PE_ref, '-o', color=cm(np.double(iz)/n_geom_z),
                     # label='dTh='+str(dTh)+', z*='+str(zstar)+', r*='+str(diff[1, iTh, iz])
                     #       + ', (PE/PE_ref='+str(diff[2, iTh, iz]) + ')')
                     label='z*='+str(zstar)+', r*='+str(diff[1, iTh, iz])
                           + ', (PE/PE_ref='+str(np.round(diff[2, iTh, iz]/PE_ref, 2)) + ')')
            ax1.plot(diff[1, iTh, iz], PE[iz, diff[0, iTh, iz]], 'ko')
            ax2.plot(diff[1, iTh, iz], PE[iz, diff[0, iTh, iz]]/PE_ref, 'ko')
            ax1.plot(rstar_guess, PE_guess, 'kd')
            ax2.plot(rstar_guess, PE_guess/PE_ref, 'dk')
        ax1.plot([rstar_range[0], rstar_range[-1]], [PE_ref, PE_ref], '-k', linewidth=2,
                 label='PE ref (dTh=' + str(dTh_ref) + ', z*=' + str(zstar_ref) + ', ' + 'r*=' + str(rstar_ref) + ')')
        ax2.plot([rstar_range[0], rstar_range[-1]], [1., 1.], '-k', linewidth=2)

        ax1.set_xlabel('r*   [m]')
        ax2.set_xlabel('r*   [m]')
        ax1.set_ylabel('PE(r*)   [J]')
        ax2.set_ylabel('PE(r*)/PE_ref')
        ax1.set_xlim(rstar_range[0], diff[1, iTh, 0]+2*dx)
        ax1.set_ylim(np.amin(PE), PE_ref*1.3)
        ax2.set_xlim(rstar_range[0], diff[1, iTh, 0]+2*dx)
        ax2.set_ylim(0.5, 1.5)
        ax2.legend(loc='center left', bbox_to_anchor=(1., 0.5),fontsize=8)
        plt.suptitle('dTh='+str(dTh) + ', marg='+str(marg)+'m', fontsize=15)
        plt.subplots_adjust(bottom=0.12, right=.75, left=0.07, top=0.9, wspace=0.25)
        fig1.savefig(os.path.join(path_out, 'initialization_dTh'+str(dTh)+'K_marg'+str(marg)+'m.png'))
        plt.close(fig1)


        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11,5))
        cf = ax1.contourf(PE[:,:])
        plt.colorbar(cf, ax=ax1)
        ax1.set_title('PE(z*,r*)')
        ax1.set_xlabel('z*')
        ax1.set_ylabel('r*')
        lvls = np.linspace(0,2, 1e3)
        cf = ax2.contourf(PE[:,:]/PE_ref, levels=lvls, extend='max')
        plt.colorbar(cf, ax=ax2)
        ax2.set_title('PE(z*,r*)/PE(1000m, 1000m)')
        ax2.set_xlabel('z*')
        ax2.set_ylabel('r*')
        plt.savefig(os.path.join(path_out, 'initialization_dTh' + str(dTh) + 'K_marg'+str(marg)+'m_contourf.png'))
        plt.close(fig)

        file_name = 'PE_init_dTh'+str(dTh)+'K_marg'+str(marg)+'m.nc'
        dump_file(file_name, n_geom_z, n_geom_r, zstar_range, rstar_range, PE, PE_ref, path_out)



    print('')
    for iTh, dTh in enumerate(dTh_range):
        print('dTh='+str(dTh) + ': z*='+str(zstar_range))
        print('r*=' +str(diff[1, iTh, :]))
        print('PE_ref - PE=' +str(diff[2, iTh, :]))
        print('PE/PE_ref = '+str(diff[2, iTh, :]/PE_ref))

    return




#_______________________________

def dump_file(fname, n_zstar, n_rstar, zstar_range, rstar_range, PE, PE_ref, path_out):
    rootgrp = nc.Dataset(os.path.join(path_out, fname), 'w', format='NETCDF4')
    rootgrp.createDimension('n_zstar', n_zstar)
    rootgrp.createDimension('n_rstar', n_rstar)
    var = rootgrp.createVariable('zstar', 'f8', ('n_zstar'))
    var[:] = zstar_range[:]
    var = rootgrp.createVariable('rstar', 'f8', ('n_rstar'))
    var[:] = rstar_range[:]
    var = rootgrp.createVariable('PE_ref', 'f8', )
    var[:] = PE_ref
    var = rootgrp.createVariable('PE', 'f8', ('n_zstar', 'n_rstar'))
    var[:,:] = PE[:,:]
    rootgrp.close()
    return



# -----------------------------------------
def compute_envelope(dTh, th_g, rstar, zstar, marg, z_half):
    # k_max_arr = (-1) * np.ones((2, nlg[0], nlg[1]), dtype=np.double)
    z_max_arr = np.zeros((2, nlg[0], nlg[1]), dtype=np.double)
    # theta = th_g * np.ones(shape=(nlg[0], nlg[1], nlg[2]))
    theta_z = th_g * np.ones(shape=(nlg[0], nlg[1], nlg[2]))
    theta_pert = np.random.random_sample(npg)
    # entropy = np.empty((npl), dtype=np.double, order='c')

    for i in xrange(nlg[0]):
        ishift = i * nlg[1] * nlg[2]
        for j in xrange(nlg[1]):
            jshift = j * nlg[2]
            r = np.sqrt((x_half[i + indx_lo[0]] - xc) ** 2 +
                        (y_half[j + indx_lo[1]] - yc) ** 2)
            if r <= rstar:
                # k_max = kstar * (np.cos(r / rstar * np.pi / 2)) ** 2
                # k_max_arr[0, i, j] = np.int(np.round(k_max))
                z_max = zstar * (np.cos(r / rstar * np.pi / 2) ** 2)
                z_max_arr[0, i, j] = z_max
            if r <= (rstar + marg):
                # k_max = (kstar + marg_k) * (np.cos(r / (rstar + marg) * np.pi / 2)) ** 2
                # k_max_arr[1, i, j] = np.int(np.round(k_max))
                z_max = (zstar + marg) * (np.cos(r / (rstar + marg) * np.pi / 2) ** 2)
                z_max_arr[1, i, j] = z_max

            kstar = np.int(np.round(zstar / dz))
            for k in xrange(gw, nlg[2]-gw):
                ijk = ishift + jshift + k
                if z_half[k] <= z_max_arr[0, i, j]:
                    theta_z[i, j, k] = th_g - dTh
                elif z_half[k] <= z_max_arr[1, i, j]:
                    th = th_g - dTh * np.sin(
                        (z_half[k] - z_max_arr[1, i, j]) / (z_max_arr[1, i, j] - z_max_arr[0, i, j]) * np.pi / 2) ** 2
                    theta_z[i, j, k] = th
                # if k <= kstar + 2:
                #     theta_pert_ = (theta_pert[ijk] - 0.5) * 0.1
                # else:
                #     theta_pert_ = 0.0
                # PV.values[s_varshift + ijk] = entropy_from_thetas_c(theta_z[i, j, k] + theta_pert_, 0.0)
                # entropy[ijk] = entropy_from_thetas_c(theta_z[i, j, k] + theta_pert_, 0.0)
    return z_max_arr, theta_z



def compute_PE_density_approx(dTh, zstar, rstar):

    return dTh * zstar**2 * rstar**2



def compute_PE(theta_z, th_g, z_max, rho0_stats, zhalf_stats):
    g = 9.80665
    dV = dx * dy * dz

    # p0 = np.log(Pg)
    # # Perform the integration at integration points z_half
    # # p = odeint(rhs, p0, z_, hmax=1.0)[:, 0]  # type: # object
    # # p = np.exp(p)
    # p_half = odeint(rhs, p0, z_half[gw:nz+gw+1], hmax=1.0)[1:, 0]
    # p_half = np.exp(p_half)
    # # temperature[k], ql[k], qi[k] = Thermodynamics.eos(p_[k], self.sg, self.qtg)
    # T, ql, qi = eos(p_half, sg, qtg)
    # # qv[k] = self.qtg - (ql[k] + qi[k])
    # # qv = np.zeros(p.shape)
    # qv = np.zeros(p_half.shape)
    # alpha0 = alpha_c(p_half, T, qtg, qv)
    # rho0 = 1./alpha0
    # # print(rho0.shape, rho0_stats.shape, z_half.shape)
    # # print('diff rho: ', np.amax(np.abs(rho0[:nz] - rho0_stats[:nz])), np.amax(rho0[:nz]))

    # plt.figure()
    # plt.subplot(1,2,1)
    # plt.plot(zhalf_stats[:nz], label='stats', linewidth=3)
    # plt.plot(z_half[gw:nz+gw], '-x', label='calc (gw)')
    # plt.subplot(1,2,2)
    # plt.plot(rho0_stats[:nz], label='stats', linewidth=3)
    # plt.plot(rho0[gw:nz+gw], '-x', label='calc (gw)')
    # plt.legend()
    # plt.show()

    kmax = np.int(z_max/dz) + 20
    if not kmax <= nz:
        print('in PE computation: looping outwards of vertical domain extent')
        sys.exit()

    PE = 0.0
    PE_av = 0.0
    theta_av = np.average(np.average(theta_z, axis=0), axis=0)
    for i in range(nx):
        for j in range(ny):
            for k in range(kmax):
                delta_th = th_g - theta_z[i,j,k]
                PE += zhalf_stats[k] * delta_th * dV * rho0_stats[k]
                delta_th_av = theta_z[i, j, k] - theta_av[k]
                PE_av += z_half[k + gw] * delta_th_av * dV * rho0_stats[k]
    PE = g / th_g * PE
    PE_av = g / th_g * PE_av

    # print('test: ', PE_test / PE)

    return PE

#_______________________________
# compute Temperature from pressure and entropy
def eos(pd, s, qt):
    ql = np.zeros(pd.shape)
    qi = np.zeros(pd.shape)
    eos_c = T_tilde * (np.exp( (s - sd_tilde + Rd*np.log(pd / p_tilde) ) / cpd ) )
    return eos_c, ql, qi

def alpha_c(p0, T, qt, qv):
    return (Rd * T)/p0 * (1.0 - qt + eps_vi * qv)

def rhs(p, z):
    ql = 0.0
    qi = 0.0
    # given sfc values for pressure, temperature and moisture
    # >> compute sfc entropy (= constant value throught atmosphere for reference profile being defines as constant-entropy profile)
    # compute temperature from pressure and entropy (at any given height)
    T, ql, qi = eos(np.exp(p), sg, qtg)
    T, ql, qi = eos(np.exp(p), sg, qtg)
    rhs_ = -g / (Rd * T * (1.0 - qtg + eps_vi * (qtg - ql - qi)))
    return rhs_
#_______________________________



def set_parameters():
    global g, Rd, cpd, T_tilde, sd_tilde, p_tilde, eps_vi
    g = 9.80665
    Rd = 287.1
    Rv = 461.5
    eps_v = 0.62210184182
    eps_vi = 1.60745384883
    cpd = 1004.0
    cpv = 1859.0
    cvd = 718.
    gamma = cpd / cvd
    T_tilde = 298.15
    p_tilde = 100000.0
    pv_star_t = 611.7
    sd_tilde = 6864.8
    sv_tilde = 10513.6

    'surface values'
    global Pg, Tg, qtg, sg
    Pg = 1.0e5
    Tg = 300.0
    qtg = 0.0
    # sg = entropy_dry(Pg, Tg, qtg, 0.0, 0.0)
    sg = 6e3
    return




def define_geometry():

    global nx, ny, nz, dx, dy, dz
    nx = 80
    ny = 80
    nz = 50
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

    global ic, jc, xc, yc, marg_i, marg_k, marg
    # zstar_range = [4000, 2000, 1500, 1000, 670, 500, 250]
    # rstar_range = [250, 500, 670, 1000, 1500, 2000, 4000]
    # rstar_range = [2000]
    # zstar_range = [2000]
    ic = np.int(nx / 2)
    jc = np.int(ny / 2)
    xc = x_half[ic + gw]  # center of cold-pool !!!
    yc = y_half[jc + gw]  # center of cold-pool !!!
    marg = 200.
    marg_i = np.int(np.round( marg / dx ))
    marg_k = np.int(np.round( marg / dz))  # width of margin


    return

if __name__ == '__main__':
    main()