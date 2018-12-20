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

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.labelsize'] = 12
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
# plt.rcParams['figure.titlesize'] = 35

def main():
    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path_root")
    parser.add_argument("dTh", type=int)
    parser.add_argument("--zparams", nargs='+', type=int)
    parser.add_argument('--rparams', nargs='+', type=int)
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    args = parser.parse_args()

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')

    nml, dTh, z_params, r_params, tmin, tmax = set_input_parameters(args)
    i0_center, j0_center = define_geometry(case_name, nml)

    ng = len(z_params)
    kmax = np.amax(z_params) + 2000./dx[2]

    # KE = np.ndarray((ng, nt))
    #     , KEd, KE_x

    print ' '
    for istar in range(ng):
        zstar = z_params[istar]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)

        path_in = os.path.join(path_root, id)
        path_fields = os.path.join(path_root, id, 'fields')
        path_out = os.path.join(path_root, id, 'figs_CP_energy')
        if not os.path.exists(path_out):
            os.mkdir(path_out)

        print('id', id)
        print path_in
        print path_out

        ''' create output file '''
        filename = 'CP_energy_' + id + '.nc'
        create_output_file(filename, path_out)

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
        # 3. define rim of cold pool
        # define_cp_rim()
        # 4. integrate: KE = 0.5*sum_i(rho_i*dV*v_i**2) from center (ic,jc) to rim

        ic = ic_arr[0]
        jc = jc_arr[0]
        KE, KEd, KE_x = compute_KE(ic, jc, irstar, tmin, tmax, path_in, path_fields, path_out)
        rootgrp = nc.Dataset(os.path.join(path_out, filename), 'r+', format='NETCDF4')
        ts_grp = rootgrp.groups['timeseries']
        var = ts_grp.variables['KE']
        var[:] = KE[:]
        var = ts_grp.variables['KEd']
        var[:] = KEd[:]
        var = ts_grp.variables['KE_x']
        var[:,:] = KE_x[:,:]
        rootgrp.close()
        del KE, KEd, KE_x


    print ' '
    path_root_out = os.path.join(path_root, 'figs_energy')
    if not os.path.exists(path_root_out):
        os.mkdir(path_root_out)
    print('path_root_out', path_root_out)
    fig = plt.figure()
    for istar in range(ng):
        zstar = z_params[istar]
        rstar = r_params[istar]
        id = 'dTh' + str(dTh) + '_z' + str(zstar) + '_r' + str(rstar)
        path_in = os.path.join(path_root, id, 'figs_CP_energy')
        print('plot ' + id)
        print path_in
        filename = 'CP_energy_' + id + '.nc'
        rootgrp = nc.Dataset(os.path.join(path_in, filename), 'r+', format='NETCDF4')
        ts_grp = rootgrp.groups['timeseries']
        var = ts_grp.variables['KE']
        var[:] = KE[:]
        rootgrp.close()

        plt.plot(times, KE, '-o', label=id)

    plt.legend(loc='best')
    plt.xlabel('time [s]')
    plt.ylabel('KE  [J]')
    plt.savefig(os.path.join(path_root_out, 'dTh_'+str(dTh)+'_KE_domain.png'))
    plt.close(fig)


    return





def compute_KE(ic, jc, irstar, tmin, tmax, path_in, path_fields, path_out):
    times = [np.int(name[:-3]) for name in os.listdir(path_fields) if name[-2:] == 'nc'
             and np.int(name[:-3]) >= tmin and np.int(name[:-3]) <= tmax]
    times.sort()
    # print('times', times)
    files = [str(t) + '.nc' for t in times]
    print 'files', files
    nt = len(times)
    kmax = 100
    krange = np.arange(0,20)

    KE = np.zeros((nt))
    KEd = np.zeros((nt))
    KE_x = np.zeros((nt, nx))       # compute KE[x, jc, :] (columnwise integration over z)

    # ishift = id - irstar
    # jshift = jd - irstar
    # th_w = 5e-1

    print 'path_in', path_in
    print path_fields
    try:
        rootgrp = nc.Dataset(os.path.join(path_in, 'Stats.' + case_name + '.nc'))
    except:
        rootgrp = nc.Dataset(os.path.join(path_in, 'stats', 'Stats.' + case_name + '.nc'))
    rho0 = rootgrp.groups['reference'].variables['rho0'][:]
    rho_unit = rootgrp.groups['reference'].variables['rho0'].units
    # z_half = rootgrp.groups['reference'].variables['z'][:]
    rootgrp.close()

    for it,t0 in enumerate(times):
        print('--t='+str(t0)+'--')
        u = read_in_netcdf_fields('u', os.path.join(path_fields, str(t0)+'.nc'))[:,:,:kmax]
        v = read_in_netcdf_fields('v', os.path.join(path_fields, str(t0)+'.nc'))[:,:,:kmax]
        w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0)+'.nc'))[:,:,:kmax]
        u2 = u * u
        v2 = v * v
        w2 = w * w
        del u, v, w

        # # define mask
        # u_ = np.roll(u[:, :, :k_max], [ishift, jshift],
        #              [0, 1])[ic + ishift - id:ic + ishift + id, jc + jshift - jd:jc + jshift + jd]
        # v_ = np.roll(v[:, :, :k_max], [ishift, jshift],
        #              [0, 1])[ic + ishift - id:ic + ishift + id, jc + jshift - jd:jc + jshift + jd]
        # w_roll = np.roll(w[:, :, :k_max], [ishift, jshift], [0, 1])
        # w_ = w_roll[ic + ishift - id:ic + ishift + id, jc + jshift - jd:jc + jshift + jd]
        # u2 = u_*u_
        # v2 = v_*v_
        # w2 = w_*w_
        #
        # w_mask = np.ma.masked_where(w_ <= th_w, w_)
        # u_mask = np.ma.masked_where(w_ <= th_w, u_)
        # v_mask = np.ma.masked_where(w_ <= th_w, v_)
        #
        # u2 = u_mask*u_mask
        # v2 = v_mask*v_mask
        # w2 = w_mask*w_mask
        # del u, v, w

        # KE ~ v**2 = (u**2 + v**2 + w**2)
        aux = np.sum(np.sum(u2+v2+w2, 0), 0)
        KEd[it] = 0.5*np.sum(aux)
        KE[it] = 0.5 * dV * np.sum(rho0[:kmax] * aux)

        for i in range(nx):
            KE_x[it, i] = 0.5 * dV * np.sum(rho0[:kmax] * (u2[i,jc,:] + v2[i,jc,:] + w2[i,jc,:]) )


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
    #
    # plt.savefig(os.path.join(path,'KE_density.png'))
    # plt.close()

    from matplotlib.colors import LogNorm
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10,5), sharey='all')
    cf = ax0.contourf(KE_x[:,ic-irstar-100:ic+irstar+100])
    plt.colorbar(cf, ax=ax0)
    cf = ax1.contourf(KE_x[:,ic-irstar-100:ic+irstar+100],norm = LogNorm())
    plt.colorbar(cf, ax=ax1)
    ax0.set_xlabel('x')
    ax1.set_xlabel('x')
    plt.ylabel('time')
    plt.suptitle('Kinetic Energy KE[x,jc,:]')
    plt.savefig(os.path.join(path_out,'KE_x.png'))
    plt.close()


    return KE, KEd, KE_x


def compute_PE(ic,jc,id,jd,nx_,ny_,case_name,path,path_fields):
    # 1. read in initial s-field
    s = read_in_netcdf_fields('s', os.path.join(path_fields,'0.nc'))
    s_ = s[ic-id:ic+id,jc-jd:jc+jd,:]
    print('shape s', s_.shape, id, 2*id, jd, 2*jd)
    # 2. convert entropy to potential temperature
    th_s = theta_s(s_)

    # 4. define background profile (here: take profile at any point outside the anomaly region)
    i0 = 0
    j0 = 0
    theta_env = th_s[i0,j0,:]
    th_g = theta_env[0]
    rootgrp = nc.Dataset(os.path.join(path, 'Stats.'+case_name+'.nc'))
    rho0 = rootgrp.groups['reference'].variables['rho0'][:]
    rho_unit = rootgrp.groups['reference'].variables['rho0'].units
    z_half = rootgrp.groups['reference'].variables['z'][:]
    rootgrp.close()


    # 5. integrate
    g = 9.80665
    # PEd = PE/kg = sum(g*dz*dTh_i) = g*dz*sum(dTh_i)
    # [PE/kg] = m/s^2*m = (m/s)^2
    # PE = m*a*s        >>  [PE] = kg*m/s^2*m = kg*(m/s)^2
    # KE = 0.5*m*v^2    >>  [KE] = kg*(m/s)^2
    # int dz a(z) = sum_i a_i dz_i
    PE = 0.0
    PEd = 0.0
    dV = dx*dy*dz
    for i in range(nx_):
        for j in range(ny_):
            for k in range(nz):
                PEd += z_half[k]*(theta_env[k] - th_s[i,j,k])
                PE +=  z_half[k]*(theta_env[k] - th_s[i,j,k]) * dV*rho0[k]
    PEd = g/th_g * PEd
    PE = g/th_g * PE
    # PE_ = g*dz*PE
    print('PE', PE, 'PEd', PEd)
    print('density at 500m: ' + str(rho0[5]) + ' ' + rho_unit)
    print('mass per grid cell at 500m: ' + str(dV * z_half[5]) + ' kg')
    i = ic
    j = jc
    k = 10
    print(z_half[k]*(theta_env[k] - th_s[i,j,k]) * dV*rho0[k])


    plt.figure()
    ax1 = plt.subplot(1, 3, 1)
    plt.imshow(s[:, jc, :].T, origin='lower')
    plt.colorbar(shrink=0.5)
    ax1.set_title('s')
    ax2 = plt.subplot(1, 3, 2)
    plt.contourf(th_s[:, np.int(ny_/2), :].T)
    plt.colorbar(shrink=0.5)
    ax2.set_title('theta')
    plt.subplot(1,3,3)
    plt.plot(rho0,z_half)
    plt.xlabel('rho0 [' + rho_unit+']')
    plt.suptitle(case_name + ': PE='+str(np.round(PEd,2)))
    plt.savefig(os.path.join(path,'pot_energy_check.png'))
    plt.savefig('./pot_energy_check.png')
    # plt.show()
    plt.close()
    del s, s_

    return

# ----------------------------------------------------------------------

def set_input_parameters(args):
    print('--- set input parameters ---')
    global case_name
    global path_root
    global times

    path_root = args.path_root
    path_out_figs = os.path.join(path_root, 'figs_CP_height')
    if not os.path.exists(path_out_figs):
        os.mkdir(path_out_figs)

    dTh = args.dTh
    z_params = args.zparams
    r_params = args.rparams
    print('z*: ', z_params)
    print('r*: ', r_params)

    case_name = args.casename
    id0 = 'dTh' + str(dTh) + '_z' + str(z_params[0]) + '_r' + str(r_params[0])
    nml = simplejson.loads(open(os.path.join(path_root, id0, case_name + '.in')).read())
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


    ''' determine file range '''
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = np.int(100)
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = np.int(100)
    times = np.arange(tmin, tmax + 100, 100)
    # times = [np.int(name[:-3]) for name in files]
    times.sort()
    print('times', times)

    return nml, dTh, z_params, r_params, tmin, tmax




def define_geometry(case_name, nml):
    print('--- define geometry ---')
    global x_half, y_half, z_half
    global ic_arr, jc_arr
    global rstar, irstar, zstar, kstar

    x_half = np.empty((nx), dtype=np.double, order='c')
    y_half = np.empty((ny), dtype=np.double, order='c')
    z_half = np.empty((nz), dtype=np.double, order='c')
    count = 0
    for i in xrange(nx):
        x_half[count] = (i + 0.5) * dx[0]
        count += 1
    count = 0
    for j in xrange(ny):
        y_half[count] = (j + 0.5) * dx[1]
        count += 1
    count = 0
    for i in xrange(nz):
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

    ''' plotting parameters '''
    if case_name == 'ColdPoolDry_single_3D':
        i0_coll = 10
        j0_coll = 10
        i0_center = ic_arr[0]
        j0_center = jc_arr[0]
    elif case_name == 'ColdPoolDry_double_3D':
        i0_coll = 0.5 * (ic_arr[0] + ic_arr[1])
        i0_center = ic_arr[0]
        j0_coll = jc_arr[0]
        j0_center = jc_arr[0]
        # domain boundaries for plotting
    elif case_name == 'ColdPoolDry_triple_3D':
        i0_coll = ic_arr[2]
        i0_center = ic_arr[0]
        j0_coll = jc_arr[2]
        j0_center = jc_arr[0]
        # domain boundaries for plotting

    return i0_center, j0_center



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
def create_output_file(filename, path_out):
    # output for each CP:
    # - min, max (timeseries)
    # - CP height (field; max=timeseries)
    # - (ok) CP rim (field)
    nt = len(times)
    print(path_out)

    rootgrp = nc.Dataset(os.path.join(path_out, filename), 'w', format='NETCDF4')

    ts_grp = rootgrp.createGroup('timeseries')
    ts_grp.createDimension('nt', nt)
    ts_grp.createDimension('nx', nx)
    var = ts_grp.createVariable('time', 'f8', ('nt'))
    var.units = "s"
    var[:] = times

    var = ts_grp.createVariable('KE', 'f8', ('nt'))
    var.units = "J"
    var = ts_grp.createVariable('KEd', 'f8', ('nt'))
    var.units = "(m/s)^2"

    var = ts_grp.createVariable('KE_x', 'f8', ('nt', 'nx'))
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

