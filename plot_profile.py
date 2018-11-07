import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import netCDF4 as nc
import argparse
import json as simplejson
import os

def main():
    parser = argparse.ArgumentParser(prog='LES_CP')
    parser.add_argument("casename")
    parser.add_argument("path")
    parser.add_argument("--tmin")
    parser.add_argument("--tmax")
    # parser.add_argument("--kmin")
    parser.add_argument("--kmax")
    args = parser.parse_args()

    global cm_bwr, cm_grey, cm_vir, cm_hsv
    cm_bwr = plt.cm.get_cmap('bwr')
    cm_grey = plt.cm.get_cmap('gist_gray_r')
    cm_hsv = plt.cm.get_cmap('hsv')
    cm_fall = plt.cm.get_cmap('winter')
    cm_summer = plt.cm.get_cmap('spring')

    timerange, kmax, nml = set_input_output_parameters(args)
    nx_half, ny_half = define_geometry(case_name, nml)
    icshift = nx_half - 1
    jcshift = ny_half

    print('out: ', path_out)

    stats_grp = nc.Dataset(os.path.join(path, 'Stats.'+case_name+'.nc'), 'r')
    zrange = stats_grp.groups['profiles'].variables['z'][:]
    stats_grp.close()


    flag_collison = True
    flag_cp_center = False


    # LOCATIONS
    ienv = np.int((ic3 + ic1) / 2)
    jenv = jc3
    ir = np.int(ic1 + irstar + 2)
    jr = np.int(jc1 + irstar + 2)

    s = read_in_netcdf_fields('s', os.path.join(path_fields, '0.nc'))
    plt.figure()
    plt.imshow(s[:, :, 0].T, origin='lower')
    plt.plot(ic1, jc1, 'o', label='center CP')
    plt.plot(ienv, jenv, 'o', label='env')
    plt.plot(icoll2, jcoll2, 'o', label='2-collision point')
    plt.plot(icoll3, jcoll3, 'o', label='3-collision point')
    plt.plot(ir, jr, 'o', label='rim')
    plt.legend(fontsize=10)
    plt.xlim([0,nx])
    plt.ylim([0,ny])
    plt.savefig(os.path.join(path_out, 'location.png'))
    plt.close()




    if flag_collison:
        # figure for 32-collision pt all times
        fig1a, (ax1a, ax2a) = plt.subplots(1, 2, figsize=(10, 6))  # , sharey=True)
        fig1a.suptitle('Profiles at Collision center')
        ax1a.set_ylabel('z  [m]')
        # ax1b.set_xlabel('s  [J/K]')
        # ax1b.set_title('entropy')
        ax1a.set_xlabel('w  [m/s]')
        ax1a.set_title('vertical velocity')
        ax2a.set_xlabel('theta  [K]')
        ax2a.set_title('potential temperature')

        # figure for 3-collision pt all times
        fig1b, (ax1b, ax2b) = plt.subplots(1, 2, figsize=(10, 6))  # , sharey=True)
        fig1b.suptitle('Profiles at Collision center')
        ax1b.set_ylabel('z  [m]')
        # ax1b.set_xlabel('s  [J/K]')
        # ax1b.set_title('entropy')
        ax1b.set_xlabel('w  [m/s]')
        ax1b.set_title('vertical velocity')
        ax2b.set_xlabel('theta  [K]')
        ax2b.set_title('potential temperature')

        for it, t0 in enumerate(timerange):
            print '--- t=' + str(t0) + ' ---'
            # 1. read in fields
            s = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))
            w = read_in_netcdf_fields('w', os.path.join(path_fields, str(t0) + '.nc'))
            # profile at cold pool center
            s_coll2 = s[icoll2, jcoll2, :]
            # 2. convert entropy to potential temperature
            th_s_coll2 = theta_s(s_coll2)
            w_coll2 = w[icoll2, jcoll2, :]

            s_coll3 = s[icoll3, jcoll3, :]
            # 2. convert entropy to potential temperature
            th_s_coll3 = theta_s(s_coll3)
            w_coll3 = w[icoll3, jcoll3, :]

            alpha_ = np.double(it) / nt
            ax1a.plot(w_coll2[:kmax], zrange[:kmax], color=cm_hsv(alpha_), label='3-coll: t=' + str(t0), linewidth=3)
            # ax1b.plot(s_coll2[:kmax], zrange[:kmax], color=cm_hsv(alpha_), label='3-coll: t=' + str(t0), linewidth=3)
            # ax1b.plot(s_env[:kmax], zrange[:kmax], color=cm_summer(alpha_), linewidth=3)
            ax2a.plot(th_s_coll2[:kmax], zrange[:kmax], color=cm_hsv(alpha_), label='t=' + str(t0), linewidth=3)
            # ax2b.plot(th_s_env[:kmax], zrange[:kmax], color=cm_summer(alpha_), linewidth=3)

            alpha_ = np.double(it) / nt
            ax1b.plot(w_coll3[:kmax], zrange[:kmax], color=cm_hsv(alpha_), label='3-coll: t=' + str(t0), linewidth=3)
            # ax1b.plot(s_coll3[:kmax], zrange[:kmax], color=cm_hsv(alpha_), label='3-coll: t=' + str(t0), linewidth=3)
            # ax1b.plot(s_env[:kmax], zrange[:kmax], color=cm_summer(alpha_), linewidth=3)
            ax2b.plot(th_s_coll3[:kmax], zrange[:kmax], color=cm_hsv(alpha_), label='t=' + str(t0), linewidth=3)
            # ax2b.plot(th_s_env[:kmax], zrange[:kmax], color=cm_summer(alpha_), linewidth=3)

        # ax1a.legend(loc='best', fontsize=8)
        ax2a.legend(loc='right', bbox_to_anchor=(1.2, 0.5), fontsize=8)
        fig1a.savefig(os.path.join(path_out, 'temp_profiles_2-collision.png'))
        plt.close(fig1a)

        # ax1b.legend(loc='best', fontsize=8)
        ax2b.legend(loc='right', bbox_to_anchor=(1.2, 0.5), fontsize=8)
        fig1b.savefig(os.path.join(path_out, 'temp_profiles_3-collision.png'))
        plt.close(fig1b)





    if flag_cp_center:
        # figure for center / env all times
        fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,6))#, sharey=True)
        ax1.set_ylabel('z  [m]')
        ax1.set_xlabel('s  [J/K]')
        ax1.set_title('entropy')
        ax2.set_xlabel('theta  [K]')
        ax2.set_title('potential temperature')

        # figure for rim pt. all times
        fig1a, (ax1a, ax2a) = plt.subplots(1, 2, figsize=(10,6))#, sharey=True)
        ax1a.set_ylabel('z  [m]')
        ax1a.set_xlabel('s  [J/K]')
        ax1a.set_title('entropy')
        ax2a.set_xlabel('theta  [K]')
        ax2a.set_title('potential temperature')

        for it, t0 in enumerate(timerange):
            print '--- t='+str(t0) + ' ---'
            # 1. read in fields
            s = read_in_netcdf_fields('s', os.path.join(path_fields, str(t0) + '.nc'))
            # profile at cold pool center
            s_cp_center = s[ic1, jc1, :]
            # 2. convert entropy to potential temperature
            th_s_cp_center = theta_s(s_cp_center)
            # profile outside cold pool
            s_env = s[ienv,jenv,:]
            th_s_env = theta_s(s_env)
            # profile at rim of cold pool
            s_rim = s[ir, jr, :]
            th_s_rim = theta_s(s_rim)


            alpha_ = np.double(it) / nt
            ax1.plot(s_cp_center[:kmax], zrange[:kmax], color=cm_summer(alpha_), label='CP center: t='+str(t0), linewidth=3)
            ax1.plot(s_env[:kmax], zrange[:kmax], color=cm_fall(alpha_), linewidth=3)
            ax2.plot(th_s_cp_center[:kmax], zrange[:kmax], color=cm_summer(alpha_), linewidth=3)
            ax2.plot(th_s_env[:kmax], zrange[:kmax], color=cm_fall(alpha_), label='env: t='+str(t0), linewidth=3)

            ax1a.plot(s_rim[:kmax], zrange[:kmax], color=cm_summer(alpha_), label='rim: t='+str(t0), linewidth=3)
            # ax1a.plot(s_env[:kmax], zrange[:kmax], color=cm_fall(alpha_), linewidth=3)
            ax2a.plot(th_s_rim[:kmax], zrange[:kmax], color=cm_summer(alpha_), label='t='+str(t0), linewidth=3)
            # ax2a.plot(th_s_env[:kmax], zrange[:kmax], color=cm_fall(alpha_), linewidth=3)


            fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10,6), sharey='all')
            # need to make the axes invisible
            ax1_ = plt.subplot(121)
            ax1_.set_ylabel('z  [m]')
            ax1_.set_xlabel('s  [J/K]')
            ax2_ = plt.subplot(122)
            ax2_.set_xlabel('theta  [K]')
            plt.setp(ax2_.get_yticklabels(), visible=False)
            ax1_.plot(s_cp_center[:kmax], zrange[:kmax], '-bd', alpha = 1., label='CP center')
            ax1_.plot(s_env[:kmax], zrange[:kmax], '-kd',label='env')
            ax1_.plot(s_rim[:kmax], zrange[:kmax], '-dr', label='i=ic+r+5')
            ax1_.legend()
            ax1_.set_title('entropy')
            ax2_.plot(th_s_cp_center[:kmax], zrange[:kmax], '-d',label='CP center')
            ax2_.plot(th_s_env[:kmax], zrange[:kmax], '-d',label='env')
            ax2_.set_title('potential temperature')
            ax2_.legend()
            ax2_.set_xlim([295., 300.2])
            plt.suptitle('temperature Profiles (t='+str(t0)+'s)')
            fig.savefig(os.path.join(path_out, 'temp_profiles_t'+str(t0)+'s.png'))
            plt.close()

        ax1.legend(loc='best', fontsize=10)
        ax2.legend(loc='best', fontsize=10)
        fig1.savefig(os.path.join(path_out, 'temp_profiles.png'))
        plt.close(fig1)

        ax1a.legend(loc='best', fontsize=10)
        ax2a.legend(loc='best', fontsize=10)
        fig1a.savefig(os.path.join(path_out, 'temp_profiles_rim.png'))
        plt.close(fig1a)






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

# ----------------------------------------------------------------------
def set_input_output_parameters(args):
    print('--- set input parameters ---')
    global path, path_in, path_out, path_stats, path_fields
    if args.path:
        path = args.path
    else:
        path = '/Users/bettinameyer/polybox/ClimatePhysics/Copenhagen/Projects/LES_ColdPool/' \
               'triple_3D_noise/Out_CPDry_triple_dTh2K/'
        # path = '/nbi/ac/cond1/meyerbe/ColdPools/triple_3D_noise/Out_CPDry_triple_Th10K/'
    path_in = os.path.join(path, 'fields_cp_rim')
    path_fields = os.path.join(path, 'fields')
    path_out = os.path.join(path, 'figs_profiles')
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    global case_name
    if args.casename:
        case_name = args.casename
    else:
        case_name = 'ColdPoolDry_triple_3D'
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dx, dy, dz, dV, gw
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx * dy * dz

    global nt
    if args.tmin:
        tmin = np.int(args.tmin)
    else:
        tmin = 100
    if args.tmax:
        tmax = np.int(args.tmax)
    else:
        tmax = tmin
    timerange = np.arange(tmin, tmax + 100, 100)
    # timerange = np.append(np.append(np.arange(100,1000,100), np.arange(1000,2900,500)), 2900)
    # timerange = np.append(np.append(np.arange(100,600,500),np.arange(600,1000,200)),  np.arange(1000,2900,200))
    # dTh = 3km, collison:
    timerange = np.append(np.arange(1000, 1500, 200), np.append(np.arange(1500, 2500, 100), np.arange(2700, 3000, 200)))#[1000, 1200, 1400, 1500]
    nt = len(timerange)

    if args.kmax:
        kmax = np.int(args.kmax)
    else:
        kmax = 50

    return timerange, kmax, nml





def define_geometry(case_name, nml):
    global nx_, ny_
    global ic1, ic2, ic3, jc1, jc2, jc3
    global icoll2, jcoll2, icoll3, jcoll3
    global shift, ishift, jshift
    global irstar
    if case_name == 'ColdPoolDry_double_3D':
        flag = 'double'
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]        irstar = np.int(np.round(rstar / dx))
        isep = 4 * irstar
        jsep = 0
        ic1 = np.int(nx / 3)
        jc1 = np.int(ny / 2)
        # ic2 = ic1 + isep
        # jc2 = jc1 + jsep
        shift = 40
        nx_half = irstar + shift
        ny_half = irstar + shift
        ishift = np.max(nx_half - ic, 0)
        jshift = np.max(ny_half - jc, 0)
        nx_ = 2 * nx_half
        ny_ = 2 * ny_half
    elif case_name == 'ColdPoolDry_triple_3D':
        flag = 'triple'
        # d = np.int(np.round(ny / 2))
        d = np.int(np.round((ny + gw) / 2))
        a = np.int(np.round(d * np.sin(60.0 / 360.0 * 2 * np.pi)))  # sin(60 degree) = np.sqrt(3)/2
        try:
            rstar = nml['init']['r']
        except:
            rstar = 5000.0  # half of the width of initial cold-pools [m]
        irstar = np.int(np.round(rstar / dx))
        ic1 = np.int(np.round(a / 2))
        ic2 = ic1
        ic3 = ic1 + np.int(np.round(a))
        jc1 = np.int(np.round(d / 2))   # np.int(np.round(d/2) + gw)
        jc2 = jc1 + d
        jc3 = jc1 + np.int(np.round(d / 2))
        shift = 60
        nx_half = irstar + shift
        ny_half = irstar + shift
        ishift = np.max(nx_half - ic1, 0)
        jshift = np.max(ny_half - jc1, 0)
        nx_ = 2 * nx_half
        ny_ = 2 * ny_half
        #collision points
        icoll2 = ic1
        jcoll2 = np.int(jc1 + np.round(d / 2) -1)
        icoll3 = ic3 - np.int(np.sqrt(3) / 3 * d)
        jcoll3 = np.int(jc1 + np.round(d / 2))

    print('rstar: ' + str(rstar), irstar)
    print('ic1,jc1,id,jd', ic1, jc1, nx_half, ny_half)
    print('nx_,ny_', nx_, ny_)
    print('shift, ishift, jshift', shift, ishift, jshift)
    return nx_half, ny_half

# ----------------------------------------------------------------------

def read_in_netcdf_fields(variable_name, fullpath_in):
    # print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name]
    data = var[:, :, :]
    rootgrp.close()
    return data

if __name__ == '__main__':
    main()