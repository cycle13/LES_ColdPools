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

    path_out = './figs_Initialization/'

    Th_g = 300.0  # temperature for neutrally stratified background (value from Soares Surface)
    dx_ = 50
    marg = 200.
    define_geometry(dx_, marg)
    set_parameters()

    # reference PE
    dTh_ref = 3
    rstar_ref = 1000
    zstar_ref = 1000

    # parameter range
    dTh_min = 3
    dTh_max = 3#10
    # dTh_max = 2
    dTh_range = np.arange(dTh_min, dTh_max+1)
    rstar_min = 200
    # for dx=100m, no matches for r*>4200m;
    # for dx=50m, matches for r*=..., 3600, 4900, 5000; no matches for r<400m
    rstar_max = 5e3
    zstar_min = 4e2
    zstar_max = 3e3
    rstar_range = np.arange(rstar_min, rstar_max+100, 100)
    zstar_range = np.arange(zstar_min, zstar_max+100, 2e2)
    # zstar_range = [zstar_ref]
    n_thermo = len(dTh_range)
    n_geom_z = len(zstar_range)
    n_geom_r = len(rstar_range)
    print('zstar', zstar_range)
    print('rstar', rstar_range)
    print('n_geom_r', n_geom_r)
    print('n_geom_z', n_geom_z)

    # scaling-range
    # scaling = [2**(-1), 2**0, 2**1, 2**2, 2**3]
    scaling = [2**(-1), 2**0, 2**1]
    print('scaling: ' + str(scaling))
    print('')

    # path to reference PE
    if dx_== 100:
        path = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/dTh3_z1000_r1000/'
    elif dx_==50:
        path = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run3_dx50m/dTh3_z1000_r1000/'
    case_name = 'ColdPoolDry_single_3D'
    rootgrp = nc.Dataset(os.path.join(path, 'stats', 'Stats.' + case_name + '.nc'))
    rho0_stats = rootgrp.groups['reference'].variables['rho0'][:]
    # alpha0_stats = rootgrp.groups['reference'].variables['alpha0'][:]
    zhalf_stats = rootgrp.groups['reference'].variables['z'][:]
    rootgrp.close()


    ''' reference PE '''
    # reference PE
    z_max_arr, theta_z = compute_envelope(dTh_ref, rstar_ref, zstar_ref, Th_g, marg)
    PE_ref = compute_PE(theta_z, Th_g, zstar_ref, rho0_stats, zhalf_stats)
    PE_ref_approx = compute_PE_density_approx(dTh_ref, zstar_ref, rstar_ref)
    #
    # # test reference numerically
    # rootgrp_field = nc.Dataset(os.path.join(path, 'fields', '0.nc'))
    # s0 = rootgrp_field.groups['fields'].variables['s'][:,:,:]
    # rootgrp_field.close()
    # ic_ = s0.shape[0] / 2
    # jc_ = s0.shape[1] / 2
    # theta = thetas_c(s0, 0.0)[ic_-nx/2:ic_+nx/2, jc_-ny/2:jc_+ny/2,:nz]
    # PE_ref_num = compute_PE(theta, Th_g, zstar_ref, rho0_stats, zhalf_stats)
    # del s0
    #
    # theta_diff = theta_z - theta
    # PE_ref_diff = compute_PE(theta_diff+Th_g, Th_g, zstar_ref, rho0_stats, zhalf_stats)
    # print('theta noise: ', np.mean(theta))
    # print('theta diff: ', np.amax(theta_diff), np.amin(theta_diff), np.mean(theta_diff))
    # print('PE_ref:      ', PE_ref)
    # print('PE_ref_num:  ', PE_ref_num)
    # print('PE_ref_diff: ', PE_ref_diff)
    # print('diff PE:     ', PE_ref - PE_ref_num)

    # fig_name = 'PE_ref.png'
    # fig, axes = plt.subplots(2,3, figsize=(12,10))
    # kmax = np.int(zstar_ref/dz) + 10
    # y0 = np.int(theta_z.shape[1]/2)
    # lvls = np.linspace(296.5, 300.5, 1e2)
    # ax = axes[0,0]
    # cf = ax.contourf(theta_z[:,y0,:kmax].T, levels=lvls, vmin=296.5, vmax=300.1)
    # plt.colorbar(cf, ax=ax)
    # ax.set_title(str(theta_z[1,1,50])+', '+str(y0))
    # ax = axes[0,1]
    # cf = ax.contourf(theta[:,y0,:kmax].T, levels=lvls, vmin=296.5, vmax=300.1)
    # ax.set_title(str(theta[1,1,50])+', '+str(y0))
    # plt.colorbar(cf, ax=ax)
    # ax = axes[0,2]
    # cf = ax.contourf(theta_diff[:,y0,:kmax].T)
    # plt.colorbar(cf, ax=ax)
    #
    # ax = axes[1,0]
    # ax.contourf(theta_z[:,:,0].T)
    # ax.plot([0,nx], [y0, y0], 'k-')
    # ax.plot([ic,ic],[0,ny], 'k-')
    # ax = axes[1,1]
    # ax.contourf(theta[:,:,0].T)
    # ax.plot([0,nx], [y0, y0], 'k-')
    # ax.plot([ic,ic],[0,ny], 'k-')
    # plt.suptitle('y='+str(y0))
    # plt.savefig(os.path.join(path_out, fig_name))
    # plt.close()
    # del theta
    # del z_max_arr, theta_z

    for dTh in dTh_range:
        print('--- dTh='+str(dTh))
        PE_range = np.zeros(len(scaling), dtype=np.single)
        fig1, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 10))
        fig2, (ax1_, ax2_, ax3_) = plt.subplots(1, 3, figsize=(20, 10))
        ax1.plot(zstar_range, PE_ref * np.ones(n_geom_z), 'k-', linewidth=3)
        ax2.plot(zstar_range, np.ones(n_geom_z), 'k-', linewidth=3)
        ax1_.plot(rstar_range, PE_ref * np.ones(n_geom_r), 'k-', linewidth=3)
        ax2_.plot(rstar_range, PE_ref * np.ones(n_geom_r), 'k-', linewidth=3)

        PE = np.zeros((n_geom_z, n_geom_r))
        # params = np.zeros(4)  # PE[dTh,r*], dTh, r*, res
        params_dict = {}
        count = np.ones(len(scaling))
        for i, s in enumerate(scaling):
            print('s=' + str(s), i)
            params_dict[str(s)] = np.zeros((1,4))
            PE_range[i] = PE_ref * s
            ax1.plot(zstar_range, PE_range[i] * np.ones(n_geom_z), 'k-', linewidth=1.5)
            ax2.plot(zstar_range, PE_range[i] * np.ones(n_geom_z) / PE_ref, 'k-', linewidth=1.5)
            ax1.fill_between([zstar_range[0], zstar_range[-1]],
                             scaling[i]*PE_ref*0.9, scaling[i]*PE_ref*1.1,
                             color='gray', alpha=0.2, label=r'$\pm 10$%$')
            ax1.fill_between([zstar_range[0], zstar_range[-1]],
                             scaling[i]*PE_ref*0.95, scaling[i]*PE_ref*1.05,
                             color='gray', alpha=0.4, label=r'$\pm 5$%$')
            ax2.fill_between([zstar_range[0], zstar_range[-1]],
                             scaling[i]-0.1, scaling[i]+0.1,
                             color='gray', alpha=0.4, label=r'$\pm$ 0.1 PE_ref')
            ax3.fill_between([zstar_range[0], zstar_range[-1]],
                             scaling[i] - 0.1, scaling[i] + 0.1,
                             color='gray', alpha=0.4, label=r'$\pm$ 0.1 PE_ref')
            ax1_.plot(rstar_range, PE_range[i] * np.ones(n_geom_r), 'k-', linewidth=1.5)
            ax2_.plot(rstar_range, PE_range[i] * np.ones(n_geom_r) / PE_ref, 'k-', linewidth=1.5)
            ax1_.fill_between([rstar_range[0], rstar_range[-1]],
                             scaling[i] * PE_ref * 0.9, scaling[i] * PE_ref * 1.1,
                             color='gray', alpha=0.2, label=r'$\pm 10$%$')
            ax1_.fill_between([rstar_range[0], rstar_range[-1]],
                             scaling[i] * PE_ref * 0.95, scaling[i] * PE_ref * 1.05,
                             color='gray', alpha=0.4, label=r'$\pm 5$%$')
            ax2_.fill_between([rstar_range[0], rstar_range[-1]],
                             scaling[i] - 0.1, scaling[i] + 0.1,
                             color='gray', alpha=0.4, label=r'$\pm$ 0.1 PE_ref')
            ax3_.fill_between([rstar_range[0], rstar_range[-1]],
                             scaling[i] - 0.1, scaling[i] + 0.1,
                             color='gray', alpha=0.4, label=r'$\pm$ 0.1 PE_ref')
            if s == scaling[-1]:
                ax2.fill_between([zstar_range[0], zstar_range[-1]],
                                 scaling[i] - 0.25, scaling[i] + 0.25,
                                 color='gray', alpha=0.2, label=r'$\pm$ 0.25 PE_ref')
                ax2_.fill_between([rstar_range[0], rstar_range[-1]],
                                 scaling[i] - 0.25, scaling[i] + 0.25,
                                 color='gray', alpha=0.2, label=r'$\pm$ 0.25 PE_ref')
        for ir, rstar in enumerate(rstar_range):
            icol = np.double(np.mod(ir, np.double(n_geom_r) / 3)) / (np.double(n_geom_r) / 3)
            print('icol', ir, n_geom_r, icol)
            for iz, zstar in enumerate(zstar_range):
                print('r*='+str(rstar)+', z*='+str(zstar))
                z_max_arr, theta = compute_envelope(dTh, rstar, zstar, Th_g, marg)
                print('shapes: ', z_max_arr.shape, theta.shape, zstar, np.double(zstar)/dx_)
                PE[iz, ir] = compute_PE(theta, Th_g, zstar, rho0_stats, zhalf_stats)
                for i, s in enumerate(scaling):
                    res = np.abs(PE[iz, ir] - PE_range[i])/PE_range[i]
                    if res < 0.05:
                        count[i] += 1
                        params_dict[str(s)] = np.append(params_dict[str(s)], [PE[iz, ir], zstar, rstar, res]).reshape(count[i], 4)
        #         if PE[iz, ir] > np.amax(scaling)*1.1*PE_ref:
        #             PE[iz+1:, ir] = 20*PE_ref
        #             break
            ax1.plot(zstar_range, PE[:, ir], '-o', color=cm(icol))
            ax2.plot(zstar_range, PE[:, ir] / PE_ref, '-o', color=cm(icol),
                     label='r*=' + str(np.int(rstar)))
            ax3.plot(zstar_range, PE[:, ir] / PE_ref, '-', color='gray')
        for iz, zstar in enumerate(zstar_range):
            icol2 = np.double(iz) / n_geom_z
            ax1_.plot(rstar_range, PE[iz, :], '-o', color=cm(icol2))
            ax2_.plot(rstar_range, PE[iz, :] / PE_ref, '-o', color=cm(icol2),
                     label='z*=' + str(np.int(zstar)))
            ax3_.plot(rstar_range, PE[iz, :] / PE_ref, '-', color='gray')

        for i, s in enumerate(scaling):
            for j in range(1,np.int(count[i])):
                pe  = params_dict[str(s)][j, 0]
                zs = params_dict[str(s)][j, 1]
                rs = params_dict[str(s)][j, 2]
                res = params_dict[str(s)][j, 3]
                ax3.plot(zs, pe/PE_ref, 'o', color=cm(res * 10),
                     label='s=' + str(np.round(s, 1)) + ', r*=' + str(np.int(rs))
                           + ', z*=' + str(np.int(zs)) + ', res=' + str(np.round(res, 2)))
                ax3_.plot(rs, pe/PE_ref, 'o', color=cm(res * 10),
                      label='s=' + str(np.round(s, 1)) + ', r*=' + str(np.int(rs))
                            + ', z*=' + str(np.int(zs)) + ', res=' + str(np.round(res, 2)))

        fsize = 15
        ax1.set_xlabel('z*   [m]', fontsize=fsize)
        ax2.set_xlabel('z*   [m]', fontsize=fsize)
        ax3.set_xlabel('z*   [m]', fontsize=fsize)
        ax1.set_ylabel('PE(z*)   [J]', fontsize=fsize)
        ax2.set_ylabel('PE(z*)/PE_ref', fontsize=fsize)
        ax3.set_ylabel('PE(z*)/PE_ref', fontsize=fsize)
        ax1.set_ylim(0, np.amax(scaling)*1.15*PE_ref)
        ax2.set_ylim(0, np.amax(scaling)+1)
        ax3.set_ylim(0, np.amax(scaling)+1)
        ax1.set_xlim(0, np.amax(zstar_range))
        ax2.set_xlim(0, np.amax(zstar_range))
        ax3.set_xlim(0, np.amax(zstar_range))
        ax1.set_xticks(zstar_range, minor=True)
        ax2.set_xticks(zstar_range, minor=True)
        ax2.set_yticks(np.arange(0, np.amax(scaling)+1), minor=True)
        ax3.set_yticks(np.arange(0, np.amax(scaling)+1), minor=True)
        ax1.grid(which='major', axis='x', linestyle='-', alpha=0.6)
        ax1.grid(which='minor', axis='x', linestyle='-', alpha=0.2)
        ax1.grid(which='major', axis='y', linestyle='-', alpha=0.6)
        ax2.grid(which='major', axis='x', linestyle='-', alpha=0.6)
        ax2.grid(which='minor', axis='x', linestyle='-', alpha=0.2)
        ax2.grid(which='both', axis='y', linestyle='-', alpha=0.6)
        ax3.grid(which='major', axis='x', linestyle='-', alpha=0.6)
        ax3.grid(which='minor', axis='x', linestyle='-', alpha=0.2)
        ax3.grid(which='both', axis='y', linestyle='-', alpha=0.6)
        # ax2.legend(loc='center left', bbox_to_anchor=(1., 0.5), fontsize=9) # for two plots
        ax2.legend(loc='center left', bbox_to_anchor=(-1.65, 0.5), fontsize=9)
        ax3.legend(loc='center left', bbox_to_anchor=(1., .5), fontsize=9)#, ncol=2)
        ax1.set_title('error bars: 5%, 10%')
        ax2.set_title('error bars: 0.1*PE_ref')
        fig1.suptitle('z*=' + str(zstar_ref) + 'm, marg=' + str(marg) + 'm (dx='+str(dx)+'m)', fontsize=15)
        fig1.subplots_adjust(bottom=0.12, right=.85, left=0.1, top=0.9, wspace=0.2)
        fig1.savefig(os.path.join(path_out, 'PE_const_initialization_dTh'+str(dTh)+'_marg'
                                  + str(np.int(marg)) + 'm_dx' + str(dx) + 'm.png'))
        plt.close(fig1)

        ax1_.set_xlabel('r*   [m]', fontsize=fsize)
        ax2_.set_xlabel('r*   [m]', fontsize=fsize)
        ax3_.set_xlabel('r*   [m]', fontsize=fsize)
        ax1_.set_ylabel('PE(dTh)   [J]', fontsize=fsize)
        ax2_.set_ylabel('PE(dTh)/PE_ref', fontsize=fsize)
        ax3_.set_ylabel('PE(dTh)/PE_ref', fontsize=fsize)
        ax1_.set_ylim(0, np.amax(scaling) * 1.15 * PE_ref)
        ax2_.set_ylim(0, np.amax(scaling) + 1)
        ax3_.set_ylim(0, np.amax(scaling) + 1)
        ax1_.set_xlim(0, rstar_max)
        ax2_.set_xlim(0, rstar_max)
        ax3_.set_xlim(0, rstar_max)
        ax1_.set_xticks(rstar_range, minor=True)
        ax2_.set_xticks(rstar_range, minor=True)
        ax2_.set_yticks(np.arange(0, np.amax(scaling) + 1), minor=True)
        ax3_.set_yticks(np.arange(0, np.amax(scaling) + 1), minor=True)
        ax1_.grid(which='major', axis='x', linestyle='-', alpha=0.6)
        ax1_.grid(which='minor', axis='x', linestyle='-', alpha=0.2)
        ax1_.grid(which='major', axis='y', linestyle='-', alpha=0.6)
        ax2_.grid(which='major', axis='x', linestyle='-', alpha=0.6)
        ax2_.grid(which='minor', axis='x', linestyle='-', alpha=0.2)
        ax2_.grid(which='both', axis='y', linestyle='-', alpha=0.6)
        ax3_.grid(which='major', axis='x', linestyle='-', alpha=0.6)
        ax3_.grid(which='minor', axis='x', linestyle='-', alpha=0.2)
        ax3_.grid(which='both', axis='y', linestyle='-', alpha=0.6)
        # ax2.legend(loc='center left', bbox_to_anchor=(1., 0.5), fontsize=9) # for two plots
        ax2_.legend(loc='center left', bbox_to_anchor=(-1.65, 0.5), fontsize=9)
        ax3_.legend(loc='center left', bbox_to_anchor=(1., .5), fontsize=9)  # , ncol=2)
        ax1_.set_title('error bars: 5%, 10%')
        ax2_.set_title('error bars: 0.1*PE_ref')
        fig2.suptitle('z*=' + str(zstar_ref) + 'm, marg=' + str(marg) + 'm (dx='+str(dx)+'m)', fontsize=15)
        fig2.subplots_adjust(bottom=0.12, right=.8, left=0.1, top=0.9, wspace=0.2)
        fig2.savefig(os.path.join(path_out, 'PE_const_initialization_dTh'+str(dTh)+'_marg'
                                  + str(np.int(marg)) + 'm_dx' + str(dx) +'m_rstar.png'))
        plt.close(fig2)

#     file_name = 'PE_scaling_marg'+str(np.int(marg))+'m_dx'+str(dx)+'m.nc'
#     dump_file(file_name, dTh_range, zstar_range, rstar_range, PE_ref, scaling, PE_range, PE,
#               params_dict, count, path_out)
#
#
#
#     # print('')
#     # for iTh, dTh in enumerate(dTh_range):
#     #     print('dTh='+str(dTh) + ': z*='+str(zstar_range))
#     #     print('r*=' +str(diff[1, iTh, :]))
#     #     print('PE_ref - PE=' +str(diff[2, iTh, :]))
#     #     print('PE/PE_ref = '+str(diff[2, iTh, :]/PE_ref))
#
#     return
#
#
#
#
# #_______________________________
#
# def dump_file(fname, dTh_range, zstar_range, rstar_range,
#               PE_ref, scaling, PE_ref_range, PE_numerical,
#               params_dict, count, path_out):
#     rootgrp = nc.Dataset(os.path.join(path_out, fname), 'w', format='NETCDF4')
#     rootgrp.createDimension('n_zstar', len(zstar_range))
#     rootgrp.createDimension('n_rstar', len(rstar_range))
#     rootgrp.createDimension('n_dTh', len(dTh_range))
#     rootgrp.createDimension('n_scal', len(scaling))
#     rootgrp.createDimension('n', np.amax(count))
#     rootgrp.createDimension('params', 3)
#     var = rootgrp.createVariable('zstar', 'f8', ('n_zstar'))
#     var[:] = zstar_range[:]
#     var = rootgrp.createVariable('rstar', 'f8', ('n_rstar'))
#     var[:] = rstar_range[:]
#     var = rootgrp.createVariable('dTh', 'f8', ('n_dTh'))
#     var[:] = dTh_range[:]
#     var = rootgrp.createVariable('scaling', 'f8', ('n_scal'))
#     var[:] = scaling[:]
#     var = rootgrp.createVariable('PE_ref', 'f8', )
#     var[:] = PE_ref
#     var = rootgrp.createVariable('PE_ref_scaled', 'f8', ('n_scal'))
#     var[:] = PE_ref_range[:]
#     var = rootgrp.createVariable('PE_numerical', 'f8', ('n_dTh', 'n_rstar'))
#     var[:,:] = PE_numerical[:,:]
#
#     for i, s in enumerate(scaling):
#         var = rootgrp.createVariable('parameters_s'+str(s), 'f8', ('n', 'params'))
#         n = params_dict[str(s)].shape[0]-1
#         print('s', s, params_dict[str(s)].shape)
#         if n > 0:
#             var[:n,:] = params_dict[str(s)][1:,1:]
#     rootgrp.close()
#     return
#
#
#
#_______________________________
def compute_envelope(dTh, rstar, zstar, th_g, marg):
    z_max_arr = np.zeros((2, nx, ny), dtype=np.double)
    theta_z = th_g * np.ones(shape=(nx, ny, nz))
    theta_pert = np.random.random_sample(npg)
    # entropy = np.empty((npl), dtype=np.double, order='c')

    for i in xrange(nx):
        ishift = i * ny * nz
        for j in xrange(ny):
            jshift = j * nz
            r = np.sqrt((x_half[i] - xc) ** 2 +
                        (y_half[j] - yc) ** 2)
            if r <= rstar:
                z_max = zstar * (np.cos(r / rstar * np.pi / 2) ** 2)
                z_max_arr[0, i, j] = z_max
            if r <= (rstar + marg):
                z_max = (zstar + marg) * (np.cos(r / (rstar + marg) * np.pi / 2) ** 2)
                z_max_arr[1, i, j] = z_max

            for k in xrange(nz):
                # ijk = ishift + jshift + k
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
    [nx,  ny, nz] = theta_z.shape

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
    print('shape 2: ', theta_z.shape, kmax, z_max/dz)
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
                delta_th_av = theta_av[k] - theta_z[i, j, k]
                PE_av += z_half[k + gw] * delta_th_av * dV * rho0_stats[k]
    PE = g / th_g * PE
    PE_av = g / th_g * PE_av
    # print('test: ', PE_test / PE)
    return PE

#_______________________________
# # compute Temperature from pressure and entropy
# def eos(pd, s, qt):
#     ql = np.zeros(pd.shape)
#     qi = np.zeros(pd.shape)
#     eos_c = T_tilde * (np.exp( (s - sd_tilde + Rd*np.log(pd / p_tilde) ) / cpd ) )
#     return eos_c, ql, qi
#
# def alpha_c(p0, T, qt, qv):
#     return (Rd * T)/p0 * (1.0 - qt + eps_vi * qv)
#
# def rhs(p, z):
#     ql = 0.0
#     qi = 0.0
#     # given sfc values for pressure, temperature and moisture
#     # >> compute sfc entropy (= constant value throught atmosphere for reference profile being defines as constant-entropy profile)
#     # compute temperature from pressure and entropy (at any given height)
#     T, ql, qi = eos(np.exp(p), sg, qtg)
#     T, ql, qi = eos(np.exp(p), sg, qtg)
#     rhs_ = -g / (Rd * T * (1.0 - qtg + eps_vi * (qtg - ql - qi)))
#     return rhs_
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




def define_geometry(dx_, marg):

    global nx, ny, nz, dx, dy, dz
    nx = 80
    ny = 80
    nz = 100
    dx = dx_
    dy = dx_
    dz = dx_
    print''
    print('resolution: dx=dz='+str(dx))
    print('margin: marg='+str(marg))
    print''
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
    x_half = np.empty((nx), dtype=np.double, order='c')
    y_half = np.empty((ny), dtype=np.double, order='c')
    z_half = np.empty((nz), dtype=np.double, order='c')
    count = 0
    for i in xrange(nx):
        x_half[count] = (i + 0.5) * dx
        count += 1
    count = 0
    for j in xrange(ny):
        y_half[count] = (j + 0.5) * dy
        count += 1
    count = 0
    for i in xrange(nz):
        z_half[count] = (i + 0.5) * dz
        count += 1

    global ic, jc, xc, yc, marg_i, marg_k
    # zstar_range = [4000, 2000, 1500, 1000, 670, 500, 250]
    # rstar_range = [250, 500, 670, 1000, 1500, 2000, 4000]
    # rstar_range = [2000]
    # zstar_range = [2000]
    ic = np.int(nx / 2)
    jc = np.int(ny / 2)
    xc = x_half[ic]  # center of cold-pool !!!
    yc = y_half[jc]  # center of cold-pool !!!
    marg_i = np.int(np.round( marg / dx ))
    marg_k = np.int(np.round( marg / dz))  # width of margin

    return

if __name__ == '__main__':
    main()