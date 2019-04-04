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
from scipy.integrate import odeint

# import thermodynamic_profiles
# from thermodynamic_profiles import alpha_c

def main():
    define_geometry()

    th_g = 300.0  # temperature for neutrally stratified background (value from Soares Surface)

    dTh_range = [2, 3, 4]
    dTh_range = [3]
    for dTh in dTh_range:
        if dTh == 2:
            # zstar_range = [2450, 1225, 815]
            # rstar_range = [815, 1225, 2450]
            zstar_range = [500, 1000, 1600, 2100, 2500]
            rstar_range = [2400, 1200, 800, 600, 500]
        elif dTh == 3:
            # zstar_range = [4000, 2000, 1500, 1000, 670, 500, 250]
            # rstar_range = [250, 500, 670, 1000, 1500, 2000, 4000]
            # zstar_range = [500, 1000, 1700, 2100, 2600]
            # rstar_range = [2000, 1000, 600, 500, 400]
            zstar_range = [1000]
            rstar_range = [1000]
        elif dTh == 4:
            # zstar_range = [1730, 870, 430]
            # rstar_range = [430, 870, 1730]
            zstar_range = [500, 1100, 1500, 2200]
            rstar_range = [1700, 800, 600, 400]
        n_thermo = len(dTh_range)
        n_geom = len(zstar_range)

        for ng in range(n_geom):
            zstar = zstar_range[ng]
            rstar = rstar_range[ng]
            irstar = np.int(np.round(rstar / dx))
            kstar = np.int(np.round(zstar / dz))
            print('r*, ir*, z*, k*: ', rstar, irstar, zstar, kstar)
            rstar2 = rstar ** 2
            rstar_marg2 = (rstar + marg) ** 2

            k_max_arr = (-1)*np.ones((2, nlg[0], nlg[1]), dtype=np.double)
            z_max_arr = np.zeros((2, nlg[0], nlg[1]), dtype=np.double)
            theta = th_g * np.ones(shape=(nlg[0], nlg[1], nlg[2]))
            theta_z = th_g * np.ones(shape=(nlg[0], nlg[1], nlg[2]))
            theta_pert = np.random.random_sample(npg)
            entropy = np.empty((npl), dtype=np.double, order='c')
            for i in xrange(nlg[0]):
                ishift = i * nlg[1] * nlg[2]
                for j in xrange(nlg[1]):
                    jshift = j * nlg[2]
                    r = np.sqrt((x_half[i + indx_lo[0]] - x_half[ic1 + gw]) ** 2 +
                                (y_half[j + indx_lo[1]] - y_half[jc1 + gw]) ** 2)
                    r2 = ((x_half[i + indx_lo[0]] - x_half[ic1 + gw]) ** 2 +
                          (y_half[j + indx_lo[1]] - y_half[jc1 + gw]) ** 2)
                    if r2 <= rstar_marg2:
                        k_max = (kstar + marg_k) * (np.cos(r / (rstar + marg) * np.pi / 2)) ** 2
                        k_max_arr[1, i, j] = np.int(np.round(k_max))
                        z_max = (zstar + marg) * (np.cos(r / (rstar + marg) * np.pi / 2) ** 2)
                        z_max_arr[1, i, j] = z_max
                        z_max_arr[1, i + (ic2 - ic1), j + (jc2 - jc1)] = z_max
                        z_max_arr[1, i + (ic3 - ic1), j + (jc3 - jc1)] = z_max
                        if r2 <= rstar2:
                            k_max = kstar * (np.cos(r / rstar * np.pi / 2)) ** 2
                            k_max_arr[0, i, j] = np.int(np.round(k_max))
                            k_max_arr[0, i + (ic2 - ic1), j + (jc2 - jc1)] = k_max
                            k_max_arr[0, i + (ic3 - ic1), j + (jc3 - jc1)] = k_max
                            z_max = zstar * ( np.cos(r / rstar * np.pi/2) **2 )
                            z_max_arr[0, i, j] = z_max
                            z_max_arr[0, i + (ic2 - ic1), j + (jc2 - jc1)] = z_max
                            z_max_arr[0, i + (ic3 - ic1), j + (jc3 - jc1)] = z_max


                    # for k in xrange(gw,nlg[2]-gw):
                    for k in xrange(nlg[2]):
                        ijk = ishift + jshift + k
                        if k <= k_max_arr[0, i, j] + gw:
                            if (k-gw) > kstar:
                                print('k', k)
                            # count_0 += 1
                            theta[i, j, k] = th_g - dTh
                        elif (k-gw) <= k_max_arr[1, i, j]:
                            th = th_g - dTh * np.sin(( (k-gw) - k_max_arr[1, i, j]) / (k_max_arr[1, i, j] - k_max_arr[0, i, j]) * np.pi/2) ** 2
                            # th = th_g - dTh * np.sin((k - k_max_arr[1, i, j]) / (k_max_arr[1, i, j] - k_max_arr[0, i, j])) ** 2
                            theta[i, j, k] = th
                        if k <= kstar + 2:
                            theta_pert_ = (theta_pert[ijk] - 0.5) * 0.1
                        else:
                            theta_pert_ = 0.0
                        # PV.values[s_varshift + ijk] = entropy_from_thetas_c(theta[i, j, k] + theta_pert_, 0.0)
                        # entropy[ijk] = entropy_from_thetas_c(theta[i, j, k] + theta_pert_, 0.0)

                    for k in xrange(nlg[2]):
                        ijk = ishift + jshift + k
                        if z_half[k] <= z_max_arr[0, i, j]:
                            theta_z[i, j, k] = th_g - dTh
                        elif z_half[k] <= z_max_arr[1, i, j]:
                            th = th_g - dTh * np.sin((z_half[k] - z_max_arr[1, i, j]) / (z_max_arr[0, i, j] - z_max_arr[1, i, j]) * np.pi/2) ** 2
                            theta_z[i, j, k] = th

                        # --- adding noise ---
                        if k <= kstar + 2:
                            theta_pert_ = (theta_pert[ijk] - 0.5) * 0.1
                        else:
                            theta_pert_ = 0.0
                        # PV.values[s_varshift + ijk] = entropy_from_thetas_c(theta_z[i, j, k] + theta_pert_, 0.0)
                        # entropy[ijk] = entropy_from_thetas_c(theta_z[i, j, k] + theta_pert_, 0.0)

            icg = ic1 + gw
            icg = ic1
            jcg = jc1 + gw
            jcg = jc1
        #     print('max(k_max[0,:,:])', np.amax(k_max_arr[0,:,:]))
        #     print('max(k_max[0,ic+gw-1,:])', np.amax(k_max_arr[0,icg-1,:]))
        #     print('max(k_max[0,ic+gw,:])', np.amax(k_max_arr[0,icg,:]))
        #     print('max(k_max[0,ic+gw+1,:])', np.amax(k_max_arr[0,icg+1,:]))
        #
        #
            plotting(dTh, rstar, irstar, zstar, kstar, theta, theta_z, k_max_arr, z_max_arr, icg)
        #
        #     'surface values'
        #     set_parameters()
        #     compute_PE(theta_z, th_g, z_max_arr)
        #
        #     print('')

    return



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
def compute_PE(theta_z, th_g, z_max_arr):
    g = 9.80665
    dV = dx * dy * dz

    path = '/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/dTh3_z1000_r1000/stats/'
    case_name = 'ColdPoolDry_single_3D'
    rootgrp = nc.Dataset(os.path.join(path, 'Stats.' + case_name + '.nc'))
    rho0_stats = rootgrp.groups['reference'].variables['rho0'][:]
    alpha0_stats = rootgrp.groups['reference'].variables['alpha0'][:]
    zhalf_stats = rootgrp.groups['reference'].variables['z'][:]
    rho_unit = rootgrp.groups['reference'].variables['rho0'].units
    rootgrp.close()

    p0 = np.log(Pg)
    # Perform the integration at integration points z_half
    # p = odeint(rhs, p0, z_, hmax=1.0)[:, 0]  # type: # object
    # p = np.exp(p)
    p_half = odeint(rhs, p0, z_half[gw:nz+gw+1], hmax=1.0)[1:, 0]
    p_half = np.exp(p_half)
    # temperature[k], ql[k], qi[k] = Thermodynamics.eos(p_[k], self.sg, self.qtg)
    T, ql, qi = eos(p_half, sg, qtg)
    # qv[k] = self.qtg - (ql[k] + qi[k])
    # qv = np.zeros(p.shape)
    qv = np.zeros(p_half.shape)
    alpha0 = alpha_c(p_half, T, qtg, qv)
    rho0 = 1./alpha0
    print(rho0.shape, rho0_stats.shape, z_half.shape)
    print('diff: ', np.amax(np.abs(rho0[:nz] - rho0_stats[:nz])), np.amax(rho0[:nz]))

    # plt.figure()
    # plt.subplot(1,2,1)
    # plt.plot(zhalf_stats[:nz], label='stats', linewidth=3)
    # plt.plot(z_half[gw:nz+gw], '-x', label='calc (gw)')
    # plt.subplot(1,2,2)
    # plt.plot(rho0_stats[:nz], label='stats', linewidth=3)
    # plt.plot(rho0[gw:nz+gw], '-x', label='calc (gw)')
    # plt.legend()
    # plt.show()

    theta_av = np.average(np.average(theta_z, axis=0), axis=0)
    PE = 0.0
    PE_av = 0.0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                # delta_th = theta_z[i,j,k] - th_g
                # PE += z_half[k+gw] * delta_th * dV * 1./alpha0[k]
                if z_half[k] <= z_max_arr[0, i, j]:
                    delta_th = theta_z[i,j,k] - th_g
                    # PE += z_half[k+gw] * delta_th * dV * rho0_stats[k]
                    PE += zhalf_stats[k] * delta_th * dV * rho0_stats[k]
                    delta_th_av = theta_z[i,j,k] - theta_av[k]
                    PE_av += z_half[k+gw] * delta_th_av * dV * rho0_stats[k]
    print(np.amax(theta_z), np.amin(theta_z), th_g)
    print(delta_th, delta_th_av)
    PE_av = g / th_g * PE_av
    PE = g / th_g * PE
    print('PE: '+str(PE), PE_av)

    return


def plotting(dTh, rstar, irstar, zstar, kstar, theta, theta_z, k_max_arr, z_max_arr, icg):
    ''' plot theta[k=0]'''
    theta_ = theta[gw:-gw, gw:-gw, gw:-gw]
    theta_z_ = theta_z[gw:-gw, gw:-gw, gw:-gw]

    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    # ax = axes[0, 0]
    ax1 = axes[0, 0]
    ax2 = axes[1, 0]
    # im = ax1.imshow(theta_[:, :, 0].T, origin='lower', cmap=cm.bwr)
    # plt.colorbar(im, ax=ax1, shrink=0.5)
    ax1.plot(ic1, jc1, 'or')
    ax1.plot([ic1, ic1], [0, ny], 'k')
    ax1.plot([ic1 + irstar, ic1 + irstar], [0, ny], 'w:')
    ax1.plot([ic1 - irstar, ic1 - irstar], [0, ny], 'w:')
    ax1.plot([jc1 - irstar - marg_i, jc1 - irstar - marg_i], [0, ny], '--', color='w', linewidth=1,
             label='jc1-irstar-marg_i')
    ax1.plot([jc1 + irstar + marg_i, jc1 + irstar + marg_i], [0, ny], '--', color='w', linewidth=1)
    ax1.plot([x_half[ic1] / dx, x_half[ic1] / dx], [0, ny], '--k')
    ax1.plot([0, nx], [jc1, jc1], 'k')
    circle1 = plt.Circle((ic1, jc1), irstar, fill=False, color='lime', linewidth=2)
    circle2 = plt.Circle((ic1, jc1), irstar + marg_i, fill=False, color='lime', linewidth=1)
    ax1.add_artist(circle1)
    ax1.add_artist(circle2)
    ax1.set_xlim([0, nx])
    ax1.set_ylim([0, ny])
    # ax1.set_title(str(count_0))
    ax2.contourf(x_half[gw:-gw], y_half[gw:-gw], theta_[:, :, 0])
    ax2.plot([x_half[ic1], x_half[ic1]], [y_half[gw], y_half[-gw]], 'k')
    ax2.plot([x_half[gw], x_half[-gw]], [y_half[jc1], y_half[jc1]], 'k')

    ax1 = axes[0, 1]
    ax2 = axes[1, 1]
    im = ax1.imshow(theta_[ic1, :, :].T, origin='lower', cmap=cm.bwr)
    ax1.plot(k_max_arr[0, icg, gw:-gw], 'gold', label='k_max[0]', linewidth=3)
    ax1.plot(k_max_arr[1, icg, gw:-gw], 'lime', label='k_max[1]', linewidth=3)
    ax1.plot([jc1, jc1], [0, nz], 'k')
    ax1.plot([jc1 - irstar, jc1 - irstar], [0, nz], ':', color='lightgray', linewidth=2, label='jc+irstar')
    ax1.plot([jc1 + irstar, jc1 + irstar], [0, nz], ':', color='lightgray', linewidth=2)
    ax1.plot([jc1 - irstar - marg_i, jc1 - irstar - marg_i], [0, nz], '--', color='lightgray', linewidth=2,
             label='jc-irstar-marg_i')
    ax1.plot([jc1 + irstar + marg_i, jc1 + irstar + marg_i], [0, nz], '--', color='lightgray', linewidth=2)
    ax1.plot([0, ny], [kstar, kstar], color='aqua', label='kstar', linewidth=2)
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
               fancybox=True, shadow=True, ncol=2, fontsize=10)
    ax1.set_title('max(kmax0)=' + str(np.amax(k_max_arr[0, :, :])) +
                  ', max(kmax0[ic])=' + str(np.round(np.amax(k_max_arr[0, icg, :]))), fontsize=10)
    ax1.grid()
    ax1.set_xlim([0, ny])
    ax1.set_ylim([0, nz])
    ax2.contourf(theta_[ic1, :, :].T, cmap=cm.bwr)
    ax2.plot(k_max_arr[0, icg, gw:-gw], 'gold', label='k_max[0]', linewidth=3)
    ax2.plot(k_max_arr[1, icg, gw:-gw], 'lime', label='k_max[1]', linewidth=3)
    ax2.plot([jc1, jc1], [0, nz], 'k')
    ax2.plot([jc1 - irstar, jc1 - irstar], [0, nz], ':', color='lightgray', linewidth=2, label='jc+irstar')
    ax2.plot([jc1 + irstar, jc1 + irstar], [0, nz], ':', color='lightgray', linewidth=2)
    ax2.plot([jc1 - irstar - marg_i, jc1 - irstar - marg_i], [0, nz], '--', color='lightgray', linewidth=2,
             label='jc-irstar-marg_i')
    ax2.plot([jc1 + irstar + marg_i, jc1 + irstar + marg_i], [0, nz], '--', color='lightgray', linewidth=2)
    ax2.grid()
    ax2.set_xlim([0, nx])

    ax1 = axes[0, 2]
    ax2 = axes[1, 2]
    ax1.imshow(theta_z_[ic1, :, :].T, origin='lower', cmap=cm.bwr)
    ax1.plot(z_max_arr[0, icg, gw:-gw] / dz, 'gold', label='z_max[0]/dz', linewidth=3)
    ax1.plot(z_max_arr[1, icg, gw:-gw] / dz, 'lime', label='z_max[1]/dz', linewidth=3)
    ax1.plot([jc1, jc1], [0, nz], 'k')
    ax1.plot([jc1 - irstar, jc1 - irstar], [0, nz], ':', color='lightgray', linewidth=2, label='jc+irstar')
    ax1.plot([jc1 + irstar, jc1 + irstar], [0, nz], ':', color='lightgray', linewidth=2)
    ax1.plot([jc1 - irstar - marg_i, jc1 - irstar - marg_i], [0, nz], '--', color='lightgray', linewidth=2,
             label='jc-irstar-marg_i')
    ax1.plot([jc1 + irstar + marg_i, jc1 + irstar + marg_i], [0, nz], '--', color='lightgray', linewidth=2)
    ax1.plot([0, ny], [kstar, kstar], color='aqua', label='kstar', linewidth=2)
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
               fancybox=True, shadow=True, ncol=2, fontsize=10)
    ax1.set_title('max(zmax0)=' + str(np.amax(z_max_arr[0, :, :])) +
                  ', max(zmax0[ic])=' + str(np.round(np.amax(z_max_arr[0, icg, :]))), fontsize=10)
    ax1.grid()
    ax1.set_xlim([0, ny])
    # ax1.set_ylim([0, nz])
    ax2.contourf(y_half[gw:-gw], z_half[gw:-gw], theta_z_[ic1, :, :].T, cmap=cm.bwr)
    ax2.plot(y_half[gw:-gw], z_max_arr[0, icg, gw:-gw], 'gold', label='z_max[0]', linewidth=3)
    ax2.plot(y_half[gw:-gw], (np.round(z_max_arr[0, icg, gw:-gw] / dz)) * dz, '--', color='k', label='z_max[0]',
             linewidth=3)
    ax2.plot(y_half[gw:-gw], z_max_arr[1, icg, gw:-gw], 'lime', label='z_max[1]', linewidth=3)
    ax2.plot(y_half[gw:-gw], (np.round(z_max_arr[1, icg, gw:-gw] / dz)) * dz, '--', color='k', label='z_max[0]',
             linewidth=3)
    ax2.plot([jc1, jc1], [0, nz], 'k')
    ax2.plot([jc1 - irstar, jc1 - irstar], [0, nz], 'w:')
    ax2.plot([jc1 + irstar, jc1 + irstar], [0, nz], 'w:')
    ax2.plot([jc1 - irstar - marg_i, jc1 - irstar - marg_i], [0, nz], 'w--')
    ax2.plot([jc1 + irstar + marg_i, jc1 + irstar + marg_i], [0, nz], 'w--')
    # ax4.legend()
    # legend = plt.legend(frameon=1)
    # frame.set_facecolor('green')
    # frame.set_edgecolor('red')

    ax1 = axes[0, 3]
    ax2 = axes[1, 3]
    im = ax1.imshow(theta_z_[:, :, 0].T, origin='lower', cmap=cm.bwr)
    plt.colorbar(im, ax=ax1, shrink=0.5)
    ax1.plot(ic1, jc1, 'or')
    ax1.plot([ic1, ic1], [0, ny], 'k')
    ax1.plot([ic1 + irstar, ic1 + irstar], [0, ny], 'w:')
    ax1.plot([ic1 - irstar, ic1 - irstar], [0, ny], 'w:')
    ax1.plot([jc1 - irstar - marg_i, jc1 - irstar - marg_i], [0, ny], '--', color='w', linewidth=1,
             label='jc-irstar-marg_i')
    ax1.plot([jc1 + irstar + marg_i, jc1 + irstar + marg_i], [0, ny], '--', color='w', linewidth=1)
    ax1.plot([xc1 / dx, xc1 / dx], [0, ny], '--k')
    ax1.plot([0, nx], [jc1, jc1], 'k')
    circle1 = plt.Circle((ic1, jc1), irstar, fill=False, color='lime', linewidth=2)
    circle2 = plt.Circle((ic1, jc1), irstar + marg_i, fill=False, color='lime', linewidth=1)
    ax1.add_artist(circle1)
    ax1.add_artist(circle2)
    ax1.set_xlim([0, nx])
    ax1.set_ylim([0, ny])
    ax2.contourf(x_half[gw:-gw], y_half[gw:-gw], theta_z_[:, :, 0])
    ax2.plot([xc1, xc1], [y_half[gw], y_half[-gw]], 'k')
    ax2.plot([x_half[gw], x_half[-gw]], [y_half[jc1], y_half[jc1]], 'k')

    plt.tight_layout
    plt.suptitle('r*=' + str(rstar) + ' (irstar=' + str(irstar) + '), z*=' + str(zstar) + ' (kstar=' + str(kstar) + ')')
    fig.savefig('./figs_Initialization/triple_initialization_k_dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar) + '.png')
    plt.close(fig)

    return



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
    nx = 150
    ny = 150
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

    global ic_arr, jc_arr, xc1, yc1, marg_i, marg_k, marg
    global ic1, ic2, ic3, jc1, jc2, jc3
    # d = namelist['init']['d']
    d = 6000.
    id = np.int(np.round(d / dx))
    idhalf = np.int(np.round(id / 2))
    a = np.int(np.round(id * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
    r_int = np.int(np.sqrt(3.) / 6 * id)  # radius of inscribed circle
    # point of 3-CP collision (ic, jc)
    ic = np.int(np.round(nx / 2))
    jc = np.int(np.round(ny / 2))
    ic1 = ic - r_int
    ic2 = ic1
    ic3 = ic + r_int
    jc1 = jc - idhalf
    jc2 = jc + idhalf
    jc3 = jc
    print ic1, ic2, jc1, jc2, r_int
    ic_arr = np.asarray([ic1, ic2, ic3])
    jc_arr = np.asarray([jc1, jc2, jc3])
    xc1 = x_half[ic1 + gw]  # center of cold-pool 1
    yc1 = y_half[jc1 + gw]  # center of cold-pool 1
    # xc1 = x_half[ic + gw]  # center of cold-pool !!!
    # yc1 = y_half[jc + gw]  # center of cold-pool !!!
    marg_k = np.int(np.round(0. / dz))  # width of margin
    marg_i_ = np.int(0. / np.round(dx))  # width of margin
    marg_i = np.int(np.round(0. / dx))  # width of margin
    marg_ = marg_i_ * dx  # width of margin
    marg = marg_i * dx  # width of margin

    print('')
    # print('nx: ' + str(Gr.dims.n[0]), str(Gr.dims.n[1]))
    # print('nyg: ' + str(Gr.dims.ng[0]), str(Gr.dims.ng[1]))
    # print('gw: ' + str(Gr.dims.gw))
    print('d: ' + str(d) + ', id: ' + str(id))
    print('Cold Pools:')
    print('collision point: ', ic, jc)
    print('cp1: [' + str(ic1) + ', ' + str(jc1) + ']')
    print('cp2: [' + str(ic2) + ', ' + str(jc2) + ']')
    print('cp3: [' + str(ic3) + ', ' + str(jc3) + ']')
    print('')


    return

if __name__ == '__main__':
    main()