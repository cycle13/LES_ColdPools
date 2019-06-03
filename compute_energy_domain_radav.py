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

# compute potential temperature by integrating over anomaly
#   PE = \int dz g * (th_anomaly(z) - th_env(z)) * z
# KE ~ v**2 = (v_rad**2 + v_tan**2 + w**2)

# label_size = 8
# plt.rcParams['xtick.labelsize'] = label_size
# plt.rcParams['ytick.labelsize'] = label_size
# plt.rcParams['lines.linewidth'] = 2
# plt.rcParams['legend.fontsize'] = 8
# plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35

def main():
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    # global cm_bwr, cm_grey, cm_vir, cm_hsv
    # cm_bwr = plt.cm.get_cmap('bwr')
    # cm_grey = plt.cm.get_cmap('gist_gray_r')
    # cm_hsv = plt.cm.get_cmap('hsv')

    filename_in = 'stats_radial_averaged.nc'

    nml, times = set_input_output_parameters(args, filename_in)
    define_geometry(case_name, nml)
    path_fields = os.path.join(path, 'fields')
    path_stats = os.path.join(path, 'data_analysis')
    path_figs = os.path.join(path, 'figs_CP_energy')
    if not os.path.exists(path_figs):
        os.mkdir(path_figs)
    if case_name == 'ColdPoolDry_single_3D':
        ic = ic_arr[0]
        jc = jc_arr[0]

    id = os.path.basename(path)
    print('id: ', id)
    filename_out = 'CP_energy_' + id + '_domain_radav.nc'

    # ''' create output file '''
    # create_output_file(times, filename_in, path_stats, filename_out, path_stats)
    #
    #
    #
    # ''' (A) Domain '''
    # ''' (A1) Potential Energy (PE) '''
    # PE = compute_PE(times, filename_in, filename_out,
    #                 case_name, path, path_fields, path_stats)
    #
    # ''' (A2) Kinetic Energy (KE) '''
    # print times
    # KE, KEd, KE_r = compute_KE(times, id, filename_in, filename_out, path_stats, path_figs)

    ''' plotting '''
    plot_KE_PE(id, filename_out, path_stats, path_figs)
    return

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------



def compute_KE(times, id, filename_in, filename_out, path_stats, path_figs):
    print ''
    print('--- compute KE ---')
    # 1. read in velocity fields
    # 2. read in reference rho
    # 3. define rim of cold pool
    # define_cp_rim()
    # 4. integrate: KE = 0.5*sum_i(rho_i*dV*v_i**2) from center (ic,jc) to rim


    # # th_w = 5e-1


    # 1. read in reference density
    rootgrp = nc.Dataset(os.path.join(path, 'stats', 'Stats.' + case_name + '.nc'))
    rho0 = rootgrp.groups['reference'].variables['rho0'][:]
    # rho_unit = rootgrp.groups['reference'].variables['rho0'].units
    # z_half = rootgrp.groups['reference'].variables['z'][:]
    rootgrp.close()


    # 1. (A) read in azimuthally averaged s-field
    file_radav = nc.Dataset(os.path.join(path, path_stats, filename_in))
    time_in = file_radav.groups['timeseries'].variables['time'][:]
    radius = file_radav.groups['stats'].variables['r'][:]
    krange = file_radav.groups['dimensions'].variables['krange'][:]
    # u_rad = file_radav.groups['stats'].variables['u'][:, :, :]  # u(nt, nr, nz)
    v_rad = file_radav.groups['stats'].variables['v_rad'][:, :, :]  # v(nt, nr, nz)
    v_tan = file_radav.groups['stats'].variables['v_tan'][:, :, :]  # v(nt, nr, nz)
    w_rad = file_radav.groups['stats'].variables['w'][:, :, :]  # w(nt, nr, nz)
    file_radav.close()
    nk = len(krange)
    nr = len(radius)

    v2_rad = v_rad * v_rad
    v2_tan = v_tan * v_tan
    w2_av = w_rad * w_rad




    nt = len(times)
    KE = np.zeros((nt))
    KEd = np.zeros((nt))
    KE_r = np.zeros((nt, nr))           # compute KE[t, r, :] (columnwise integration over z)

    # KE ~ v**2 = (v_rad**2 + v_tan**2 + w**2)
    for it, t0 in enumerate(times):
        print('--t=' + str(t0) + '--')
        aux = np.sum(v2_rad[it,:,:nk]+v2_tan[it,:,:nk]+w2_av[it,:,:nk], 0)
        print 'shapes', v2_rad.shape, aux.shape, rho0.shape
        KEd[it] = 0.5*np.sum(aux)
        KE[it] = 0.5 * dV * np.sum(rho0[:nk] * aux)

        for i in range(nr):
            KE_r[it, i] = 0.5 * dV * np.sum(rho0[:nk] * (v_rad[it,i,:nk] + v_tan[it,i,:nk] + w2_av[it,i,:nk]) )

    fig, axis = plt.subplots(2, 3, figsize=(15, 10), sharey='none')
    rmax = 60
    ax1 = axis[0, 0]
    ax2 = axis[0, 1]
    ax3 = axis[0, 2]
    cf = ax1.contourf(times, radius[:rmax], v2_rad[:, :rmax, 0].T)
    cbar = plt.colorbar(cf, ax=ax1, shrink=0.6, aspect=12)
    cf = ax2.contourf(times, radius[:rmax], v2_tan[:, :rmax, 0].T)
    cbar = plt.colorbar(cf, ax=ax2, shrink=0.6, aspect=12)
    cf = ax3.contourf(times, radius[:rmax], w2_av[:, :rmax, 0].T)
    cbar = plt.colorbar(cf, ax=ax3, shrink=0.6, aspect=12)
    ax1.set_xlabel('time [s]')
    ax2.set_xlabel('time [s]')
    ax3.set_xlabel('time [s]')
    ax1.set_ylabel('radius [m]')

    ax1 = axis[1, 0]
    ax2 = axis[1, 1]
    ax3 = axis[1, 2]
    ax1.plot(times, KEd)
    ax2.plot(times, KE)
    cf = ax3.contourf(times, radius[:rmax], KE_r[:, :rmax].T)
    cbar = plt.colorbar(cf, ax=ax3, shrink=0.6, aspect=12)
    # # ax2.set_ylabel('radius [m]')
    # # ax3.set_ylabel('radius [m]')
    plt.savefig(os.path.join(path_figs, 'KE_test_' + id + '.png'))
    plt.close(fig)


    # # 1. (B) read in 3D fields
    # for it,t0 in enumerate(times):
    #     print('--t='+str(t0)+'--')
    #     u = read_in_netcdf_fields('u', os.path.join(path_fields, str(t0)+'.nc'))[:,:,:kmax]
    #     v = read_in_netcdf_fields('v', os.path.join(path_fields, str(t0)+'.nc'))[:,:,:kmax]
    #     w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0)+'.nc'))[:,:,:kmax]
    #     u2 = u * u
    #     v2 = v * v
    #     w2 = w * w
    #     del u, v, w
    #
    # #     # # define mask
    # #     # u_ = np.roll(u[:, :, :k_max], [ishift, jshift],
    # #     #              [0, 1])[ic + ishift - id:ic + ishift + id, jc + jshift - jd:jc + jshift + jd]
    # #     # v_ = np.roll(v[:, :, :k_max], [ishift, jshift],
    # #     #              [0, 1])[ic + ishift - id:ic + ishift + id, jc + jshift - jd:jc + jshift + jd]
    # #     # w_roll = np.roll(w[:, :, :k_max], [ishift, jshift], [0, 1])
    # #     # w_ = w_roll[ic + ishift - id:ic + ishift + id, jc + jshift - jd:jc + jshift + jd]
    # #     # u2 = u_*u_
    # #     # v2 = v_*v_
    # #     # w2 = w_*w_
    # #     #
    # #     # w_mask = np.ma.masked_where(w_ <= th_w, w_)
    # #     # u_mask = np.ma.masked_where(w_ <= th_w, u_)
    # #     # v_mask = np.ma.masked_where(w_ <= th_w, v_)
    # #     #
    # #     # u2 = u_mask*u_mask
    # #     # v2 = v_mask*v_mask
    # #     # w2 = w_mask*w_mask
    # #     # del u, v, w
    # #
    # #     # KE ~ v**2 = (u**2 + v**2 + w**2)
    # #     aux = np.sum(np.sum(u2+v2+w2, 0), 0)
    # #     KEd[it] = 0.5*np.sum(aux)
    # #     KE[it] = 0.5 * dV * np.sum(rho0[:kmax] * aux)
    # #
    # #     for i in range(nx):
    # #         KE_x[it, i] = 0.5 * dV * np.sum(rho0[:kmax] * (u2[i,jc,:] + v2[i,jc,:] + w2[i,jc,:]) )



    ''' output '''
    rootgrp = nc.Dataset(os.path.join(path_stats, filename_out), 'r+', format='NETCDF4')
    ts_grp = rootgrp.groups['timeseries']
    prof_grp = rootgrp.groups['profiles']
    var = ts_grp.variables['KE']
    var[:] = KE[:]
    var = ts_grp.variables['KEd']
    var[:] = KEd[:]
    var = prof_grp.variables['KE_r']
    var[:, :] = KE_r[:, :]
    rootgrp.close()

    return KE, KEd, KE_r
    # return




def compute_PE(times, filename_in, filename_out,
               case_name, path, path_fields, path_stats):
    print ''
    print('--- compute PE ---')
    # 1. read in initial s-field
    # 2. convert entropy to potential temperature
    # 3. ??? convert potential temperature to density
    # 4. define background profile (here: take profile at any point outside the anomaly region)
    # 5. integrate

    # 1. read in azimuthally averaged s-field
    file_radav = nc.Dataset(os.path.join(path, 'data_analysis', filename_in))
    time_in = file_radav.groups['timeseries'].variables['time'][:]
    radius = file_radav.groups['stats'].variables['r'][:]
    krange = file_radav.groups['dimensions'].variables['krange'][:]
    s_in = file_radav.groups['stats'].variables['s'][:, :, :]       # s(nt, nr, nz)
    file_radav.close()
    nk = len(krange)
    nt = len(times)

    # 2. convert entropy to potential temperature
    th_s = theta_s(s_in)

    # 3. define background profile (here: take profile at any point outside the anomaly region)
    s_ = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
    i0 = 0
    j0 = 0
    theta_env = theta_s(s_[i0,j0,:])
    th_g = theta_env[0]
    del s_

    rootgrp = nc.Dataset(os.path.join(path, 'stats', 'Stats.'+case_name+'.nc'))
    rho0 = rootgrp.groups['reference'].variables['rho0'][:]
    # rho_unit = rootgrp.groups['reference'].variables['rho0'].units
    # z_half = rootgrp.groups['reference'].variables['z'][:]
    rootgrp.close()


    # 4. integrate
    g = 9.80665
    # PEd = PE/kg = sum(g*dz*dTh_i) = g*dz*sum(dTh_i)
    # [PE/kg] = m/s^2*m = (m/s)^2
    # PE = m*a*s        >>  [PE] = kg*m/s^2*m = kg*(m/s)^2
    # KE = 0.5*m*v^2    >>  [KE] = kg*(m/s)^2

    # int dz a(z) = sum_i a_i dz_i
    PE = np.zeros(nt)
    PEd = np.zeros(nt)


    for i,r in enumerate(radius):
        for k in range(nk):
            PEd += z_half[k]*(theta_env[k] - th_s[:nt,i,k])
            PE +=  z_half[k]*(theta_env[k] - th_s[:nt,i,k]) * dV*rho0[k]
    PEd = g/th_g * PEd
    PE = g/th_g * PE
    # # PE_ = g*dz*PE
    # print('PE', PE, 'PEd', PEd)
    # print('density at 500m: ' + str(rho0[5]) + ' ' + rho_unit)
    # print('mass per grid cell at 500m: ' + str(dV * z_half[5]) + ' kg')
    # i = ic
    # j = jc
    # k = 10
    # print(z_half[k]*(theta_env[k] - th_s[i,j,k]) * dV*rho0[k])


    # plt.figure()
    # ax1 = plt.subplot(1, 3, 1)
    # plt.imshow(s[:, jc, :].T, origin='lower')
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
    # plt.savefig(os.path.join(path,'pot_energy_check.png'))
    # plt.savefig('./pot_energy_check.png')
    # # plt.show()
    # plt.close()
    # del s, s_

    rootgrp = nc.Dataset(os.path.join(path_stats, filename_out), 'r+', format='NETCDF4')
    ts_grp = rootgrp.groups['timeseries']
    var = ts_grp.variables['PE']
    # print(time_in)
    # print(PE.shape)
    # print(var.shape)
    var[:] = PE[:]
    var = ts_grp.variables['PEd']
    var[:] = PEd[:]
    rootgrp.close()

    return PE






# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

def plot_KE_PE(id, file_data, path_data, path_figs):

    root = nc.Dataset(os.path.join(path_data, file_data), 'r')
    PE = root.groups['timeseries'].variables['PE'][:]
    KE = root.groups['timeseries'].variables['KE'][:]
    KE_r = root.groups['profiles'].variables['KE_r'][:,:]
    times = root.groups['dimensions'].variables['time'][:]
    radius = root.groups['dimensions'].variables['radius'][:]
    (nt, nr) = KE_r.shape
    root.close()
    print path_figs


    fig, axis = plt.subplots(1, 4, figsize=(15,5))
    ax1 = axis[0]
    ax2 = axis[1]
    ax3 = axis[2]
    ax4 = axis[3]
    ax1.set_title('PE')
    ax2.set_title('KE')
    ax3.set_title('KE, PE')
    ax4.set_title('KE + PE')

    ax2.plot(times, KE)
    ax1.plot(times, PE)
    ax3.plot(times, PE, label='PE')
    ax3.plot(times, KE, label='KE')
    ax3.legend(loc='best')
    ax4.plot(times, PE+KE)

    ax1.set_xlabel('time [s]')
    ax2.set_xlabel('time [s]')
    ax3.set_xlabel('time [s]')
    ax4.set_xlabel('time [s]')
    ax1.set_ylabel('PE [J]')
    ax2.set_ylabel('KE [J]')
    ax4.set_ylabel('KE + PE [J]')

    plt.suptitle('Kinetic and Potential Energy')
    plt.tight_layout()
    plt.savefig(os.path.join(path_figs, 'KE_PE_'+id + '.png'))
    plt.close(fig)




    fig, axis = plt.subplots(1, 2, figsize=(15, 5))
    ax1 = axis[0]
    ax2 = axis[1]
    # ax3 = axis[2]
    # ax4 = axis[3]
    ax1.set_title('KE density')
    ax2.set_title('KE')
    # ax3.set_title('PE')
    # ax4.set_title('KE + PE')

    cf = ax1.contourf(times, radius, KE_r.T)
    # cbar = plt.colorbar(cf, shrink=0.6, ticks=np.arange(min, max + 1, 1), aspect=12)
    cbar = plt.colorbar(cf, ax=ax1, shrink=0.6, aspect=12)
    ax2.plot(times, KE)

    ax1.set_xlabel('time [s]')
    ax2.set_xlabel('time [s]')
    ax1.set_ylabel('radius [m]')
    ax2.set_ylabel('KE [J]')

    plt.suptitle('Kinetic Energy')
    plt.tight_layout()
    plt.savefig(os.path.join(path_figs, 'KE_' + id + '.png'))
    plt.close(fig)

    # # ax1.grid()
    # # ax2.grid()

    # # plt.suptitle('kinetic energy in rim (w>0.5m/s)')
    # # #
    # # # plt.savefig(os.path.join(path,'KE_density.png'))
    # # # plt.close()
    # #
    # # # print '.....', KE_x[:,ic-irstar-100:ic+irstar+100].shape
    # # from matplotlib.colors import LogNorm
    # # fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10,5), sharey='all')
    # # cf = ax0.contourf(KE_x[:,ic-irstar-100:ic+irstar+100])
    # # plt.colorbar(cf, ax=ax0)
    # # # lvls = np.logspace(np.amin(np.log(plot_levels)), np.amax(np.log(plot_levels)), 10)
    # # lvls = np.logspace(4, 10, 100)
    # # cf = ax1.contourf(KE_x[:,ic-irstar-100:ic+irstar+100],norm = LogNorm(), levels=lvls)
    # # plt.colorbar(cf, ax=ax1)
    # # ax0.set_xlabel('x')
    # # ax1.set_xlabel('x')
    # # plt.ylabel('time')
    # # plt.suptitle('Kinetic Energy KE[x,jc,:]')
    # # plt.savefig(os.path.join(path_out,'KE_x_'+id + '.png'))
    # # plt.close()
    return
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

def set_input_output_parameters(args, filename_in):
    print('')
    print('--- set input parameters ---')
    global case_name
    global path
    global files

    if args.casename:
        case_name = args.casename
    else:
        case_name = 'ColdPoolDry_triple_3D'
    path = args.path

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dV, gw
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = np.zeros(3, dtype=np.int)
    dx[0] = nml['grid']['dx']
    dx[1] = nml['grid']['dy']
    dx[2] = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx[0] * dx[1] * dx[2]
    print('nx, ny, nz', nx, ny, nz)

    ''' determine file range '''
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = np.int(100)
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = np.int(10000)
    path_fields = os.path.join(path, 'fields')
    # times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
    #          and np.int(name[:-3]) >= tmin and np.int(name[:-3]) <= tmax]
    # times.sort()
    print('input radav file: ' + filename_in)
    file_radav = nc.Dataset(os.path.join(path, 'data_analysis', filename_in))
    times_ = file_radav.groups['timeseries'].variables['time'][:]
    times = [np.int(t) for t in times_ if np.int(t) >= tmin and np.int(t) <= tmax]
    file_radav.close()
    files = [str(np.int(t)) + '.nc' for t in times]
    # print('times radav', times_)
    print('times', times)
    print ('')
    del times_

    return nml, times




def define_geometry(case_name, nml):
    print('--- define geometry ---')
    global x_half, y_half, z_half
    global ic_arr, jc_arr
    global rstar, irstar, zstar, kstar

    # test file:
    path_fields = os.path.join(path, 'fields')
    var = read_in_netcdf_fields('u', os.path.join(path_fields, files[0]))
    [nx_, ny_, nz_] = var.shape
    del var
    x_half = np.empty((nx_), dtype=np.double, order='c')
    y_half = np.empty((ny_), dtype=np.double, order='c')
    z_half = np.empty((nz_), dtype=np.double, order='c')
    count = 0
    for i in xrange(nx_):
        x_half[count] = (i + 0.5) * dx[0]
        count += 1
    count = 0
    for j in xrange(ny_):
        y_half[count] = (j + 0.5) * dx[1]
        count += 1
    count = 0
    for i in xrange(nz_):
        z_half[count] = (i + 0.5) * dx[2]
        count += 1

    # set coordinates for plots
    if case_name == 'ColdPoolDry_single_3D':
        rstar = nml['init']['r']
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        try:
            ic = nml['init']['ic']
            jc = nml['init']['jc']
            # print('(ic,jc) from nml')
        except:
            ic = np.int(nx/2)
            jc = np.int(ny/2)
            # print('(ic,jc) NOT from nml')
        ic_arr = np.zeros(1)
        jc_arr = np.zeros(1)
        ic_arr[0] = ic
        jc_arr[0] = jc
    elif case_name == 'ColdPoolDry_double_2D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        isep = 4 * irstar
        ic1 = np.int(nx / 3)  # np.int(Gr.dims.ng[0] / 3)
        ic2 = ic1 + isep
        jc1 = np.int(ny / 2)
        jc2 = jc1
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_double_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        kstar = np.int(np.round(zstar / dx[2]))
        isep = 4 * irstar
        jsep = 0
        ic1 = np.int(np.round((nx + 2 * gw) / 3)) - gw
        jc1 = np.int(np.round((ny + 2 * gw) / 2)) - gw
        ic2 = ic1 + isep
        jc2 = jc1 + jsep
        ic_arr = [ic1, ic2]
        jc_arr = [jc1, jc2]
    elif case_name == 'ColdPoolDry_triple_3D':
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx[0]))
        zstar = nml['init']['h']
        kstar = np.int(np.round(zstar / dx[2]))
        d = np.int(np.round(ny / 2))
        dhalf = np.int(np.round(ny / 4))
        a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        ic1 = np.int(np.round(a / 2))  # + gw
        ic2 = ic1
        ic3 = ic1 + np.int(np.round(a))
        jc1 = np.int(np.round(d / 2))  # + gw
        jc2 = jc1 + d
        jc3 = jc1 + np.int(np.round(d / 2))
        ic_arr = [ic1, ic2, ic3]
        jc_arr = [jc1, jc2, jc3]

        isep = dhalf

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

# ----------------------------------
def create_output_file(times, filename_in, path_in, filename_out, path_out):
    # output for each CP:
    # - min, max (timeseries)
    # - CP height (field; max=timeseries)
    # - (ok) CP rim (field)
    nt = len(times)
    print ''
    print('--- create output file ---')
    print(path_out)
    print(times, nt)

    root_in = nc.Dataset(os.path.join(path_in, filename_in), 'r')
    radius = root_in.groups['stats'].variables['r'][:]
    root_in.close()
    nr = len(radius)

    rootgrp = nc.Dataset(os.path.join(path_out, filename_out), 'w', format='NETCDF4')

    dim_grp = rootgrp.createGroup('dimensions')
    dim_grp.createDimension('nt', nt)
    dim_grp.createDimension('nr', nr)
    var = dim_grp.createVariable('time', 'f8', ('nt'))
    var.units = "s"
    var[:] = times
    var = dim_grp.createVariable('radius', 'f8', ('nr'))
    var.units = "m"
    var[:] = radius

    ts_grp = rootgrp.createGroup('timeseries')
    ts_grp.createDimension('nt', nt)
    var = ts_grp.createVariable('KE', 'f8', ('nt'))
    var.units = "J"
    var = ts_grp.createVariable('KEd', 'f8', ('nt'))
    var.units = "(m/s)^2"

    var = ts_grp.createVariable('PE', 'f8', ('nt'))
    var.units = "J"
    var = ts_grp.createVariable('PEd', 'f8', ('nt'))
    var.units = "(m/s)^2"

    prof_grp = rootgrp.createGroup('profiles')
    prof_grp.createDimension('nt', nt)
    prof_grp.createDimension('nr', nr)
    var = prof_grp.createVariable('KE_r', 'f8', ('nt', 'nr'))
    var.units = "J"

    # field_grp = rootgrp.createGroup('fields_2D')
    # field_grp.createDimension('nt', nt)
    # field_grp.createDimension('nx', nx)
    # field_grp.createDimension('ny', ny)
    # var = field_grp.createVariable('w_max', 'f8', ('nt', 'nx', 'ny'))
    # var.units = "m/s"
    # var = field_grp.createVariable('w_max_height', 'f8', ('nt', 'nx', 'ny'))
    # var.units = "m"

    rootgrp.close()

    del radius
    del times
    return


# def dump_output_file(filename, KE, KEd, it):
#     rootgrp = nc.Dataset(os.path.join(path_out, filename), 'r+', format='NETCDF4')
#
#     ts_grp = rootgrp.groups['timeseries']
#     var = ts_grp.variables['KE']
#     var[it] = KE[it]
#     var = ts_grp.variables['KEd']
#     var[it] = KEd[it]
#
#     # field_grp = rootgrp.groups['fields_2D']
#     # var = field_grp.variables['w_max']
#     # var[it,:,:] = w_max[0,:,:]
#     # var = field_grp.variables['w_max_height']
#     # var[it,:,:] = w_max[1,:,:]
#
#     rootgrp.close()
#     return



def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:,:,:]
    rootgrp.close()
    return data


if __name__ == '__main__':
    main()

