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
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--dz")
    parser.add_argument("--nz")
    args = parser.parse_args()
    define_geometry(args)


    th_g = 300.0  # temperature for neutrally stratified background (value from Soares Surface)

    dTh_range = [2, 3, 4]
    dTh_range = [3]
    for dTh in dTh_range:
        if dTh == 2:
            zstar_range = [500, 1000, 1600, 2100, 2500]
            rstar_range = [2400, 1200, 800, 600, 500]
        elif dTh == 3:
            zstar_range = [500, 1000, 1700, 2100, 2600]
            rstar_range = [2000, 1000, 600, 500, 400]
            zstar_range = [1000]
            rstar_range = [1000]
        elif dTh == 4:
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

            # k_max_arr = (-1)*np.ones((2, nlg[0], nlg[1]), dtype=np.double)
            z_max_arr = np.zeros((2, nlg[0], nlg[1]), dtype=np.double)
            theta_bg = np.empty(shape=(nlg[2]), dtype=np.double)
            theta_neutral = th_g * np.ones(shape=(nlg[0], nlg[1], nlg[2]))
            theta_strat = th_g * np.ones(shape=(nlg[0], nlg[1], nlg[2]))
            theta_pert = np.random.random_sample(npg)
            entropy = np.empty((npl), dtype=np.double, order='c')

            # crate background profile
            Nv = 5e-5
            g = 9.81
            for k in xrange(nlg[2]):
                if z_half[k] <= 1000.:
                    theta_bg[k] = th_g
                else:
                    print(z_half[k], k)
                    theta_bg[k] = th_g * np.exp(Nv/g*(z_half[k]-1000.))

            # Cold Pool anomaly
            for i in xrange(nlg[0]):
                ishift = i * nlg[1] * nlg[2]
                for j in xrange(nlg[1]):
                    jshift = j * nlg[2]
                    r = np.sqrt((x_half[i + indx_lo[0]] - xc) ** 2 +
                                (y_half[j + indx_lo[1]] - yc) ** 2)
                    if r <= rstar:
                        z_max = zstar * ( np.cos(r/rstar * np.pi/2) ** 2)
                        z_max_arr[0, i, j] = z_max
                    if r <= (rstar + marg):
                        z_max = (zstar + marg) * (np.cos(r / (rstar + marg) * np.pi / 2) ** 2)
                        z_max_arr[1, i, j] = z_max

                    for k in xrange(nlg[2]):
                        ijk = ishift + jshift + k
                        theta_strat[i, j, k] = theta_bg[k]
                        if z_half[k] <= z_max_arr[0, i, j]:
                            theta_neutral[i, j, k] = th_g - dTh
                            theta_strat[i, j, k] = theta_bg[k] - dTh
                        elif z_half[k] <= z_max_arr[1, i, j]:
                            th = dTh * np.sin((z_half[k] - z_max_arr[1, i, j]) / (z_max_arr[1, i, j] - z_max_arr[0, i, j]) * np.pi/2) ** 2
                            theta_neutral[i, j, k] = th_g - th
                            theta_strat[i, j, k] = theta_bg[k] - th
                        if k <= kstar + 2:
                            theta_pert_ = (theta_pert[ijk] - 0.5) * 0.1
                        else:
                            theta_pert_ = 0.0
                        # PV.values[s_varshift + ijk] = entropy_from_thetas_c(theta_neutral[i, j, k] + theta_pert_, 0.0)
                        # entropy[ijk] = entropy_from_thetas_c(theta_neutral[i, j, k] + theta_pert_, 0.0)

            icg = ic + gw
            jcg = jc + gw
            print('max(z_max[0,:,:])      ', np.amax(z_max_arr[0,:,:]))
            print('max(z_max[0,ic+gw-1,:])', np.amax(z_max_arr[0,icg-1,:]))
            print('max(z_max[0,ic+gw,:])  ', np.amax(z_max_arr[0,icg,:]))
            print('max(z_max[0,ic+gw+1,:])', np.amax(z_max_arr[0,icg+1,:]))
            print('')
            print(theta_bg[5+gw:15+gw])
            print(theta_strat[ic,jc,5+gw:15+gw])


            plotting(dTh, rstar, irstar, zstar, kstar, theta_neutral, theta_strat, z_max_arr, theta_bg, icg)
        #
        #     'surface values'
        #     set_parameters()
        #     compute_PE(theta_neutral, th_g, z_max_arr)

            # print('')

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



def plotting(dTh, rstar, irstar, zstar, kstar, theta_n, theta_s, z_max_arr, theta_bg, icg):
    ''' plot theta[k=0]'''
    theta_n_ = theta_n[gw:-gw, gw:-gw, gw:-gw]
    theta_s_ = theta_s[gw:-gw, gw:-gw, gw:-gw]
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))

    ax1 = axes[0, 0]
    ax2 = axes[1, 0]
    im = ax1.imshow(theta_n_[:, :, 0].T, origin='lower', cmap=cm.bwr)
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
    circle2 = plt.Circle((ic, jc), irstar + marg_i, fill=False, color='lime', linewidth=1)
    ax1.add_artist(circle1)
    ax1.add_artist(circle2)
    ax1.set_xlim([0, nx])
    ax1.set_ylim([0, ny])
    ax1.set_title(str(np.amin(theta_n_)) + ', ' + str(np.amax(theta_n_)))


    ax1 = axes[0, 1]
    ax2 = axes[1, 1]
    im = ax1.imshow(theta_n_[ic, :, :].T, origin='lower', cmap=cm.bwr)
    plt.colorbar(im, ax=ax1, shrink=0.5)
    ax1.plot(z_max_arr[0, icg, gw:-gw] / dz, 'gold', label='z_max[0]/dz', linewidth=3)
    ax1.plot(z_max_arr[1, icg, gw:-gw] / dz, 'lime', label='z_max[1]/dz', linewidth=3)
    ax1.plot([jc, jc], [0, nz], 'k')
    ax1.plot([jc - irstar, jc - irstar], [0, nz], ':', color='lightgray', linewidth=2, label='jc+irstar')
    ax1.plot([jc + irstar, jc + irstar], [0, nz], ':', color='lightgray', linewidth=2)
    ax1.plot([jc - irstar - marg_i, jc - irstar - marg_i], [0, nz], '--', color='lightgray', linewidth=2,
             label='jc-irstar-marg_i')
    ax1.plot([jc + irstar + marg_i, jc + irstar + marg_i], [0, nz], '--', color='lightgray', linewidth=2)
    ax1.plot([0, ny], [kstar, kstar], color='aqua', label='kstar', linewidth=2)
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
               fancybox=True, shadow=True, ncol=2, fontsize=10)
    ax1.grid()
    ax1.set_xlim([0, ny])
    ax1.set_ylim([0, nz])


    ax1 = axes[0, 2]
    ax2 = axes[1, 2]
    ax1.imshow(theta_s_[ic, :, :].T, origin='lower', cmap=cm.bwr)
    plt.colorbar(im, ax=ax1, shrink=0.5)
    ax1.plot(z_max_arr[0, icg, gw:-gw] / dz, 'gold', label='z_max[0]/dz', linewidth=3)
    ax1.plot(z_max_arr[1, icg, gw:-gw] / dz, 'lime', label='z_max[1]/dz', linewidth=3)
    ax1.plot([jc, jc], [0, nz], 'k')
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
    ax1.set_xlim([0, ny])
    ax1.set_ylim([0, nz])



    ax1 = axes[0, 3]
    im = ax1.imshow(theta_s_[:, :, 0].T, origin='lower', cmap=cm.bwr)
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
    circle2 = plt.Circle((ic, jc), irstar + marg_i, fill=False, color='lime', linewidth=1)
    ax1.add_artist(circle1)
    ax1.add_artist(circle2)
    ax1.set_xlim([0, nx])
    ax1.set_ylim([0, ny])
    ax1.set_title(str(np.amin(theta_s_)) + ', ' + str(np.amax(theta_s_)))


    ax1 = axes[1, 2]
    ax2 = axes[1, 3]
    ax1.plot(theta_bg[gw:-gw], np.arange(0,nz), linewidth=2)
    ax2.plot(theta_bg[gw:-gw], np.arange(0,nz)*dz, linewidth=2)

    # plt.tight_layout
    # plt.suptitle('r*=' + str(rstar) + ' (irstar=' + str(irstar) + '), z*=' + str(zstar) + ' (kstar=' + str(kstar) + ')')
    fig.savefig('./figs_Initialization/initialization_strat_dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar) + '.png')
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

def define_geometry(args):

    global nx, ny, nz, dx, dy, dz
    if args.dz:
        dz = np.int(args.dz)
    else:
        dz = 100
    if args.nz:
        nz = np.int(args.nz)
    else:
        nz = 30
    nx = 80
    ny = 80
    dx = dz
    dy = dz
    print('nz', nz, 'dz', dz)
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
    marg_k = np.int(np.round(0. / dz))  # width of margin
    marg_i_ = np.int(0. / np.round(dx))  # width of margin
    marg_i = np.int(np.round(0. / dx))  # width of margin
    marg_ = marg_i_ * dx  # width of margin
    marg = marg_i * dx  # width of margin


    return

if __name__ == '__main__':
    main()