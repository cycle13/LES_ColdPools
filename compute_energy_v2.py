import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import netCDF4 as nc
import argparse
import json as simplejson
import os

# todo: PE computed for total domain


# compute potential temperature by integrating over anomaly
#   PE = \int dz g * (th_anomaly(z) - th_env(z)) * z


def main():

    # ??? reference density???

    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("--casename")
    parser.add_argument("--path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    # parser.add_argument("--k0", nargs = '+', type = int)
    parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    global cm_bwr, cm_grey, cm_vir
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_vir = plt.cm.get_cmap('viridis')
    cm_grey = plt.cm.get_cmap('gist_gray_r')

    global path_in, path_out, path_stats, path_fields
    if args.path:
        path = args.path
    else:
        path = '/Users/bettinameyer/polybox/ClimatePhysics/Copenhagen/Projects/LES_ColdPool/' \
                'triple_3D_noise/Out_CPDry_triple_dTh2K/'
        # path = '/nbi/ac/cond1/meyerbe/ColdPools/triple_3D_noise/Out_CPDry_triple_Th3K/'
    path_in = os.path.join(path, 'fields_cp_rim')
    path_fields = os.path.join(path, 'fields')
    path_out = os.path.join(path, 'figs_energy')
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    global case_name
    if args.casename:
        case_name = args.casename
    else:
        case_name = 'ColdPoolDry_triple_3D'
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dy, dz, dV
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx * dy * dz

    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = tmin
    timerange = np.arange(tmin, tmax + 100, 100)
    nt = len(timerange)

    ''' (A) Potential Energy (PE) '''
    ''' (A1) for LES gravity current '''
    # 1. read in initial s-field
    # 2. convert entropy to potential temperature
    # 3. ??? convert potential temperature to density
    # 4. define background profile (here: take profile at any point outside the anomaly region)
    # 5. integrate
    # compute_PE(ic,jc,id,jd,nx_,ny_,case_name,path,path_fields)

    ''' (B) Kinetic Energy (KE) '''
    # 1. read in velocity fields
    # 2. read in reference rho
    # 3. read in mask cold pool rim (from: define_cp_rim.py)
    # 4. integrate: KE = 0.5*sum_i(rho_i*dV*v_i**2) from center (ic,jc) to rim

    perc = 95
    global nx_, ny_, nk
    global ic, jc, ishift, jshift, id, jd

    # READ IN density profile
    rootgrp = nc.Dataset(os.path.join(path, 'Stats.' + case_name + '.nc'))
    rho0 = rootgrp.groups['reference'].variables['rho0'][:]
    rho_unit = rootgrp.groups['reference'].variables['rho0'].units
    z_half = rootgrp.groups['reference'].variables['z'][:]
    rootgrp.close()

    ''' compute initial PE '''
    mask_file_name = 'rimmask_perc' + str(perc) + 'th' + '_t' + str(np.int(timerange[0])) + '.nc'
    rootgrp = nc.Dataset(os.path.join(path_in, mask_file_name), 'r')
    mask = rootgrp.groups['fields'].variables['mask'][:, :, :]
    # make sure that only non-zero values in krange if nc-files contain this level
    contr = rootgrp.groups['profiles'].variables['k_dumped'][:]
    zrange_ = rootgrp.groups['profiles'].variables['zrange'][:]
    krange_ = rootgrp.groups['profiles'].variables['krange'][:]
    ic = np.int(rootgrp.groups['description'].variables['ic'][:])
    jc = np.int(rootgrp.groups['description'].variables['jc'][:])
    ishift = np.int(rootgrp.groups['description'].variables['ishift'][:])
    jshift = np.int(rootgrp.groups['description'].variables['jshift'][:])
    [nx_, ny_, nk] = mask.shape
    id = np.int(nx_ / 2)
    jd = np.int(ny_ / 2)
    rootgrp.close()
    nk_ = np.int(np.sum(contr))
    nk = len(krange_)
    krange = np.zeros(nk_)
    zrange = np.zeros(nk_)
    krange[:] = krange_[0:nk_]
    zrange[:] = zrange_[0:nk_]
    del contr, krange_, zrange_
    del mask
    print('ic,jc', ic, jc, 'shifts', ishift, jshift)
    PE_init, PEd_init, PE_tot_init, PEd_tot_init = compute_PE_subdomain(rho0, rho_unit, z_half, krange, 0)



    kmax = 21
    # KEd_all = np.zeros((nt, kmax), dtype=np.double)  # kinetic energy density per level
    KE_all = np.zeros((nt, kmax), dtype=np.double)  # kinetic energy per level
    # PEd_all = np.zeros((nt, kmax), dtype=np.double)  # kinetic energy density per level
    PE_all = np.zeros((nt, kmax), dtype=np.double)  # kinetic energy per level
    PE_all_SD = np.zeros((nt, kmax), dtype=np.double)  # kinetic energy per level

    for it,t0 in enumerate(timerange):
        print('------ time: ' + str(t0) + '------')
        # t0 = tmin

        # READ IN mask file from cold pool rim definition
        mask_file_name = 'rimmask_perc' + str(perc) + 'th' + '_t' + str(t0) + '.nc'
        rootgrp = nc.Dataset(os.path.join(path_in, mask_file_name), 'r')
        mask = rootgrp.groups['fields'].variables['mask'][:,:,:]
        # make sure that only non-zero values in krange if nc-files contain this level
        contr = rootgrp.groups['profiles'].variables['k_dumped'][:]
        zrange_ = rootgrp.groups['profiles'].variables['zrange'][:]
        krange_ = rootgrp.groups['profiles'].variables['krange'][:]
        nk_ = np.int(np.sum(contr))
        krange = np.zeros(nk_)
        zrange = np.zeros(nk_)
        krange[:] = krange_[0:nk_]
        zrange[:] = zrange_[0:nk_]
        del contr
        [nx_, ny_, nk] = mask.shape
        ic = np.int(rootgrp.groups['description'].variables['ic'][:])
        jc = np.int(rootgrp.groups['description'].variables['jc'][:])
        ishift = np.int(rootgrp.groups['description'].variables['ishift'][:])
        jshift = np.int(rootgrp.groups['description'].variables['jshift'][:])
        rootgrp.close()

        id = np.int(nx_ / 2)
        jd = np.int(ny_ / 2)
        print('')
        print('read in: ' + os.path.join(path_in, mask_file_name))
        print(nx_, ny_, nk, nk_, len(krange))
        print('')



        ''' COMPUTE PE '''
        PE, PEd, PE_tot, PEd_tot = compute_PE(mask, rho0, rho_unit, z_half, krange, t0)
        PE_all[it, :] = PE[:]
        plt.figure()
        plt.plot(PE, zrange_, '-o')
        plt.xlabel('PE')
        plt.ylabel('z  [m]')
        plt.title('Potential energy in cold pool rim (PE tot: '+str(PE_tot)+'J)')
        plt.savefig(os.path.join(path_out, 'PE_levels_t'+str(t0)+'.png'))
        plt.close()

        plt.figure()
        for it_,t0_ in enumerate(timerange[0:it+1]):
            plt.plot(PE_all[it_,:], z_half[:kmax], '-o', label='t='+str(t0_)+'s')
        plt.xlabel('PE')
        plt.ylabel('z  [m]')
        plt.legend(loc='best')
        plt.title('Potential energy in cold pool rim (PE tot: '+str(PE_tot)+'J)')
        plt.savefig(os.path.join(path_out, 'PE_levels_alltimes.png'))
        plt.close()


        ''' COMPUTE PE subdomain '''
        PE, PEd, PE_tot, PEd_tot = compute_PE_subdomain(rho0, rho_unit, z_half, krange, t0)
        PE_all_SD[it, :] = PE[:]
        plt.figure()
        plt.plot(PE, zrange_, '-o')
        plt.xlabel('PE')
        plt.ylabel('z  [m]')
        plt.title('Potential energy in cold pool rim (PE tot: ' + str(PE_tot) + 'J)')
        plt.savefig(os.path.join(path_out, 'PE_levels_t' + str(t0) + '_subdomain.png'))
        plt.close()

        plt.figure()
        for it_, t0_ in enumerate(timerange[0:it + 1]):
            plt.plot(PE_all_SD[it_, :], z_half[:kmax], '-o', label='t=' + str(t0_) + 's')
        plt.xlabel('PE')
        plt.ylabel('z  [m]')
        plt.legend(loc='best')
        plt.title('Potential energy in subdomain (PE tot: ' + str(PE_tot) + 'J)')
        plt.savefig(os.path.join(path_out, 'PE_levels_alltimes_subdomain.png'))
        plt.close()


        ''' COMPUTE KE '''
        KEd = np.zeros(nk_, dtype=np.double)     # kinetic energy density per level
        KE = np.zeros(nk, dtype=np.double)      # kinetic energy per level
        for ik, k0 in enumerate(krange):
            print('k='+str(k0))

            # plt.figure()
            # plt.contourf(mask[:,:,ik].T)
            # plt.colorbar()
            # plt.savefig(os.path.join(path_out, 'mask_k'+str(k0)+'.png'))

            KEd[ik] = compute_KE(mask[:,:,ik], np.int(k0), t0, path, path_fields)
            KE[ik] = 0.5 * rho0[np.int(k0)] * dV * KEd[ik]
            KE_all[it,ik] = KE[ik]
        KE_tot = np.sum(KE)

        plt.figure()
        plt.plot(KE, zrange_, '-o')
        plt.xlabel('KE')
        plt.ylabel('z  [m]')
        plt.title('Kinetic energy in cold pool rim (KE tot: '+str(KE_tot)+'J)')
        plt.savefig(os.path.join(path_out, 'KE_levels_t'+str(t0)+'.png'))
        plt.close()

        plt.figure()
        for it_,t0_ in enumerate(timerange[0:it+1]):
            plt.plot(KE_all[it_,:], z_half[:kmax], '-o', label='t='+str(t0_)+'s')
        plt.xlabel('KE')
        plt.ylabel('z  [m]')
        plt.legend(loc='best')
        plt.title('Kinetic energy in cold pool rim (KE tot: '+str(KE_tot)+'J)')
        plt.savefig(os.path.join(path_out, 'KE_levels_alltimes.png'))
        plt.close()

        print('')


    return





def compute_KE(mask_, k0, t0, path, path_in):

    u = read_in_netcdf_fields('u', os.path.join(path_in,str(t0)+'.nc'))
    v = read_in_netcdf_fields('v', os.path.join(path_in,str(t0)+'.nc'))
    w = read_in_netcdf_fields('w', os.path.join(path_in,str(t0)+'.nc'))

    # define mask
    u_ = np.roll(np.roll(u[:, :, k0], ishift, axis=0), jshift, axis=1)[ic + ishift - id:ic + ishift + id, jc + jshift - jd:jc + jshift + jd]
    v_ = np.roll(np.roll(v[:, :, k0], ishift, axis=0), jshift, axis=1)[ic + ishift - id:ic + ishift + id, jc + jshift - jd:jc + jshift + jd]
    w_roll = np.roll(np.roll(w[:, :,k0], ishift, axis=0), jshift, axis=1)    # nbi
    # w_roll = np.roll(w[:, :, k0], [ishift, jshift], [0, 1])     # on Laptop gives the same as above
    w_ = w_roll[ic + ishift - id:ic + ishift + id, jc + jshift - jd:jc + jshift + jd]

    w_mask = np.ma.masked_where(mask_==0, w_)
    u_mask = np.ma.masked_where(mask_==0, u_)
    v_mask = np.ma.masked_where(mask_==0, v_)

    plt.figure(figsize=(16,6))
    plt.subplot(2, 4, 2)
    max = np.maximum(np.abs(np.amin(u_)), np.amax(u_))
    plt.contourf(u_[:, :].T, cmap=cm_bwr, vmin=-max, vmax=max)
    plt.colorbar()
    plt.subplot(2, 4, 3)
    max = np.maximum(np.abs(np.amin(v_)), np.amax(v_))
    plt.contourf(v_[:, :].T, cmap=cm_bwr, vmin=-max, vmax=max)
    plt.colorbar()
    plt.subplot(2, 4, 4)
    max = np.maximum(np.abs(np.amin(w_)), np.amax(w_))
    plt.contourf(w_[:, :].T, cmap=cm_bwr, vmin=-max, vmax=max)
    plt.colorbar()

    plt.subplot(2,4,5)
    plt.contourf(mask_[:, :].T)
    plt.colorbar()
    plt.subplot(2,4,6)
    max = np.maximum(np.abs(np.amin(u_)), np.amax(u_))
    plt.contourf(u_mask[:, :].T, cmap=cm_bwr, vmin=-max, vmax=max)
    plt.colorbar()
    plt.subplot(2,4,7)
    max = np.maximum(np.abs(np.amin(v_)), np.amax(v_))
    plt.contourf(v_mask[:, :].T, cmap=cm_bwr, vmin=-max, vmax=max)
    plt.colorbar()
    plt.subplot(2,4,8)
    max = np.maximum(np.abs(np.amin(w_)), np.amax(w_))
    plt.contourf(w_mask[:, :].T, cmap=cm_bwr, vmin=-max, vmax=max)
    plt.colorbar()
    plt.savefig(os.path.join(path_out, 'mask2_k' + str(k0) + '.png'))
    plt.close()

    u2 = u_mask*u_mask
    v2 = v_mask*v_mask
    w2 = w_mask*w_mask
    del u, v, w

    # KE ~ v**2 = (u**2 + v**2 + w**2)
    KEd = 0.5*np.sum(u2+v2+w2)

    #     # for k0 in krange:
    #     for k0 in [0]:
    #         KE[it] += rho0[k0] * dV * np.sum(u2[:,:,k0] + v2[:,:,k0] + w2[:,:,k0])
    #     KE[it] = 0.5 * KE[it]
    #
    # plt.figure(figsize=(12,6))
    # ax1 = plt.subplot(1,2,1)
    # ax2 = plt.subplot(1,2,2)
    # ax1.set_title('KE density')
    # ax2.set_title('KE')
    # ax1.plot(times[1:],KEd[1:])
    # ax2.plot(times[1:],KE[1:])
    # ax1.grid()
    # ax2.grid()
    # ax1.set_xlabel('time [s]')
    # ax2.set_xlabel('time [s]')
    # ax1.set_ylabel('KE density [J/kg]')
    # ax2.set_ylabel('KE [J]')
    # plt.suptitle('kinetic energy in rim (w>0.5m/s)')
    # plt.savefig(os.path.join(path,'KE_density.png'))
    # plt.close()
    return KEd


def compute_PE(mask, rho0, rho_unit, z_half, krange, t0):
    # 1. read in initial s-field
    s0 = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))[:, :, :nk]

    # read in s-field
    s = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))

    # define mask
    s_ = np.roll(np.roll(s[:, :, :], ishift, axis=0), jshift, axis=1)[ic + ishift - id:ic + ishift + id,
         jc + jshift - jd:jc + jshift + jd]
    s_mask = np.ma.masked_where(mask == 0, s_[:, :, :nk])

    # 2. convert entropy to potential temperature
    th_s = theta_s(s_)

    # 4. define background profile (here: take profile of initial distribution outside anomaly)
    # i0 = 0
    # j0 = 0
    # theta_env = th_s[i0,j0,:]
    i0 = 1
    j0 = 1
    theta_env = theta_s(s0[i0, j0, :])
    th_g = theta_env[0]

    # 5. integrate
    g = 9.80665
    # PEd = PE/kg = sum(g*dz*dTh_i) = g*dz*sum(dTh_i)
    # [PE/kg] = m/s^2*m = (m/s)^2
    # PE = m*a*s        >>  [PE] = kg*m/s^2*m = kg*(m/s)^2
    # KE = 0.5*m*v^2    >>  [KE] = kg*(m/s)^2
    # int dz a(z) = sum_i a_i dz_i
    PE = np.zeros(nk, dtype=np.double)
    PEd = np.zeros(nk, dtype=np.double)
    PE_tot = 0.0
    PEd_tot = 0.0
    dV = dx * dy * dz
    for i in range(nx_):
        for j in range(ny_):
            for ik, k in enumerate(krange):
                k = np.int(k)
                if mask[i,j,ik]:
                    PEd[k] += z_half[k] * (theta_env[k] - th_s[i, j, k])
                    PE[k] += z_half[k] * (theta_env[k] - th_s[i, j, k]) * dV * rho0[k]
                    PEd_tot += z_half[k] * (theta_env[k] - th_s[i, j, k])
                    PE_tot += z_half[k] * (theta_env[k] - th_s[i, j, k]) * dV * rho0[k]
    PEd = g / th_g * PEd
    PE = g / th_g * PE
    PEd_tot = g / th_g * PEd_tot
    PE_tot = g / th_g * PE_tot
    # PE_ = g*dz*PE
    print('PE', PE_tot, 'PEd', PEd_tot)
    print('density at 500m: ' + str(rho0[5]) + ' ' + rho_unit)
    print('mass per grid cell at 500m: ' + str(dV * z_half[5]) + ' kg')
    i = ic
    j = jc
    k = 10
    print(z_half[k] * (theta_env[k] - th_s[i, j, k]) * dV * rho0[k])

    # plt.figure(figsize=(12,4))
    # ax1 = plt.subplot(1, 3, 1)
    # plt.imshow(s0[:, jc, :].T, origin='lower')
    # plt.colorbar(shrink=0.5)
    # ax1.set_title('s')
    # ax2 = plt.subplot(1, 3, 2)
    # plt.contourf(th_s[:, np.int(ny_/2), :].T)
    # plt.colorbar(shrink=0.5)
    # ax2.set_title('theta')
    # plt.subplot(1,3,3)
    # plt.plot(rho0,z_half)
    # plt.xlabel('rho0 [' + rho_unit+']')
    # plt.suptitle(case_name + ': PE='+str(np.round(PEd,2)))
    # plt.savefig(os.path.join(path_out,'PE_check_profiles.png'))
    # # plt.show()
    # plt.close()

    del s0, s_

    return PE, PEd, PE_tot, PEd_tot




def compute_PE_subdomain(rho0, rho_unit, z_half, krange, t0):
    # 1. read in initial s-field
    s0 = read_in_netcdf_fields('s', os.path.join(path_fields,'0.nc'))[:,:,:nk]

    # read in s-field
    s = read_in_netcdf_fields('s', os.path.join(path_fields,str(t0)+'.nc'))[:,:,:nk]

    # define mask
    s_ = np.roll(np.roll(s[:, :, :], ishift, axis=0), jshift, axis=1)[ic + ishift - id:ic + ishift + id,
         jc + jshift - jd:jc + jshift + jd]

    # 2. convert entropy to potential temperature
    th_s = theta_s(s_)

    # 4. define background profile (here: take profile of initial distribution outside anomaly)
    # i0 = 0
    # j0 = 0
    # theta_env = th_s[i0,j0,:]
    i0 = 1
    j0 = 1
    theta_env = theta_s(s0[i0, j0, :])
    th_g = theta_env[0]

    # 5. integrate
    g = 9.80665
    # PEd = PE/kg = sum(g*dz*dTh_i) = g*dz*sum(dTh_i)
    # [PE/kg] = m/s^2*m = (m/s)^2
    # PE = m*a*s        >>  [PE] = kg*m/s^2*m = kg*(m/s)^2
    # KE = 0.5*m*v^2    >>  [KE] = kg*(m/s)^2
    # int dz a(z) = sum_i a_i dz_i
    PE = np.zeros(nk, dtype=np.double)
    PEd = np.zeros(nk, dtype=np.double)
    PE_tot = 0.0
    PEd_tot = 0.0
    dV = dx*dy*dz
    for i in range(nx_):
        for j in range(ny_):
            for ik,k in enumerate(krange):
                k = np.int(k)
                PEd[ik] += z_half[k]*(theta_env[k] - th_s[i,j,k])
                PE[ik] += z_half[k]*(theta_env[k] - th_s[i,j,k]) * dV*rho0[k]
                PEd_tot += z_half[k]*(theta_env[k] - th_s[i,j,k])
                PE_tot += z_half[k]*(theta_env[k] - th_s[i,j,k]) * dV*rho0[k]
    PEd = g/th_g * PEd
    PE = g/th_g * PE
    PEd_tot = g/th_g * PEd_tot
    PE_tot = g/th_g * PE_tot
    # PE_ = g*dz*PE
    print('PE', PE_tot, 'PEd', PEd_tot)
    print('density at 500m: ' + str(rho0[5]) + ' ' + rho_unit)
    print('mass per grid cell at 500m: ' + str(dV * z_half[5]) + ' kg')
    i = ic
    j = jc
    k = 10
    print(z_half[k]*(theta_env[k] - th_s[i,j,k]) * dV*rho0[k])


    # plt.figure(figsize=(12,4))
    # ax1 = plt.subplot(1, 3, 1)
    # plt.imshow(s0[:, jc, :].T, origin='lower')
    # plt.colorbar(shrink=0.5)
    # ax1.set_title('s')
    # ax2 = plt.subplot(1, 3, 2)
    # plt.contourf(th_s[:, np.int(ny_/2), :].T)
    # plt.colorbar(shrink=0.5)
    # ax2.set_title('theta')
    # plt.subplot(1,3,3)
    # plt.plot(rho0,z_half)
    # plt.xlabel('rho0 [' + rho_unit+']')
    # plt.suptitle(case_name + ': PE='+str(np.round(PEd,2)))
    # plt.savefig(os.path.join(path_out,'PE_check_profiles.png'))
    # # plt.show()
    # plt.close()

    del s0, s_

    return PE, PEd, PE_tot, PEd_tot


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
    # th_s = T_tilde * np.exp( s*sd_tilde )

    # th_s = np.ndarray(shape=s.shape, dtype=np.double)
    # [nx_,ny_,nz_] = s.shape[:]
    # for i in range(nx_):
    #     for j in range(ny_):
    #         for k in range(nz_):
    #             a = (s[i,j,k]-sd_tilde) / cpd
    #             th_s[i,j,k] = T_tilde * np.exp( a )

    th_s = np.exp( (s - sd_tilde)/cpd )

    return th_s


def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:,:,:]
    rootgrp.close()
    return data

if __name__ == '__main__':
    main()

