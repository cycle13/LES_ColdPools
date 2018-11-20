import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import netCDF4 as nc
import os
import json as simplejson


'''computing reference profiles in anelastic approximation: hydrostatic, adiabatic(=const entropy) profile'''

def main():
    ''' parameters'''
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

    path = '/cond1/meyerbe/ColdPools/triple_3D_noise/Out_CPDry_triple_dTh2K'
    case_name = 'ColdPoolDry_triple_3D'
    print('')
    print('path: ' + path)
    print('')
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    nz_ = 500
    dz_ = 10
    # nz_ = 10
    # dz = 500
    # global nx, ny, nz, dx, dy, dz, dV
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    gw = nml['grid']['gw']
    dV = dx * dy * dz
    try:
        rstar = nml['init']['r']
    except:
        rstar = 5000.0  # half of the width of initial cold-pools [m]
    zstar = nml['init']['h']




    # compute entropy from pressure and temperature
    # sg = Thermodynamics.entropy(self.Pg, self.Tg, self.qtg, 0.0, 0.0)
    def entropy_dry(pd, T, qt, ql, qi):
        print 'entropy dry'
        qt = 0.0
        ql = 0.0
        qi = 0.0

        sd_c = sd_tilde + cpd * np.log(T / T_tilde) - Rd * np.log(pd / p_tilde)
        return sd_c

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
        rhs_ = -g / (Rd * T * (1.0 - qtg + eps_vi * (qtg - ql - qi)))
        return rhs_



    p0 = np.log(Pg)
    # Construct arrays for integration points
    z_ = np.array(np.arange(nz_)*dz_)
    # Perform the integration
    p = odeint(rhs, p0, z_, hmax=1.0)[:, 0]
    # p_half = odeint(rhs, p0, z_half, hmax=1.0)[1:, 0]
    p = np.exp(p)
    # temperature[k], ql[k], qi[k] = Thermodynamics.eos(p_[k], self.sg, self.qtg)
    T, ql, qi = eos(p, sg, qtg)
    # qv[k] = self.qtg - (ql[k] + qi[k])
    qv = np.zeros(p.shape)
    # alpha[k] = Thermodynamics.alpha(p_[k], temperature[k], self.qtg, qv[k])
    alpha = alpha_c(p, T, qtg, qv)
    rho = 1./alpha



    fig, axes = plt.subplots(1,3, figsize=(12,5))
    ax = axes[0]
    ax.plot(p, z_, '-')
    ax.set_title('pressure')
    ax = axes[1]
    ax.plot(T, z_, '-')
    ax.set_title('temperature')
    ax = axes[2]
    ax.plot(rho, z_, '-')
    ax.plot([0.5, 3.0],[2000, 2000], 'k')
    # ax.plot([3.0, 3.0],[0, 5000], 'k')
    ax.plot([rho[np.int(2000./dz_)],rho[np.int(2000./dz_)]],[0, 5000], 'k')
    ax.set_title('rho')
    plt.savefig('./thermodynamic_profiles_pycles.png')
    plt.close()



    def temp_adiab_hydrostat(Tg,mu,z):
        z0 = Rd * Tg / (mu * g)
        return Tg*(1-(gamma-1)/gamma * z/z0)

    def press_adiab_hydrostat(pg,mu,z):
        z0 = Rd * Tg / (mu * g)
        coeff = gamma / (gamma-1)
        return pg*(1-(gamma-1)/gamma * z/z0)**coeff

    def rho_adiab_hydrostat(pg,Tg,mu,z):
        z0 = Rd * Tg / (mu * g)
        rhog = pg * mu / (Rd * Tg)
        coeff = 1. / (gamma-1)
        return rhog * (1 - (gamma-1)/gamma * z / z0)**coeff

    print('gamma', gamma, gamma / (gamma-1))
    mu = 29.e-3        # molecular weight air; [mu]=kg/m^3
    mu = 1.
    T_ad = temp_adiab_hydrostat(Tg, mu, z_)
    p_ad = press_adiab_hydrostat(Pg, mu, z_)
    rho_ad = rho_adiab_hydrostat(Pg, Tg, mu, z_)


    file_name = 'Stats.ColdPoolDry_triple_3D.nc'
    file = nc.Dataset(os.path.join(path, file_name), 'r')
    rho0_stats = file.groups['reference'].variables['rho0'][:]
    p0_stats = file.groups['reference'].variables['p0'][:]
    T0_stats = file.groups['reference'].variables['temperature0'][:]
    z_stats = file.groups['reference'].variables['z'][:]
    file.close()

    fig, axes = plt.subplots(2, 3, figsize=(12, 10))
    ax = axes[0,0]
    ax.plot(p*1e-5, z_, '-', label='calc')
    ax.plot(p0_stats*1e-5, z_stats, '-', label='stats')
    ax.legend()
    ax.set_ylim([0, 5000])
    ax.set_xlabel('[10^5 Pa]')
    ax.set_title('pressure')
    ax = axes[0,1]
    ax.plot(T, z_, '-', label='calc')
    ax.plot(T0_stats, z_stats, '-', label='stats')
    ax.legend()
    ax.set_ylim([0, 5000])
    ax.set_title('temperature')
    ax = axes[0,2]
    ax.plot(rho, z_, '-', label='calc')
    ax.plot(rho0_stats, z_stats, '-', label='stats')
    ax.plot([0.5, 3.0], [2000, 2000], 'k')
    # ax.plot([3.0, 3.0],[0, 5000], 'k')
    ax.plot([rho[np.int(2000. / dz_)], rho[np.int(2000. / dz_)]], [0, 5000], 'k')
    ax.legend()
    ax.set_xlim([rho[0],rho[-1]])
    ax.set_ylim([0,5000])
    ax.set_title('rho')

    ax = axes[1,0]
    ax.plot(p_ad*1e-5, z_, '-', label='analytic')
    ax.plot(p0_stats*1e-5, z_stats, '-', label='stats')
    ax.legend()
    ax.set_xlim([p_ad[0]*1e-5,p_ad[-1]*1e-5])
    ax.set_ylim([0, 5000])
    ax.set_title('pressure')
    ax.set_xlabel('p  [10^5 Pa]')
    ax.set_ylabel('z [m]')
    ax = axes[1,1]
    ax.plot(T_ad, z_, '-', label='analytic')
    ax.plot(T0_stats, z_stats, '-', label='stats')
    ax.legend()
    ax.set_xlim([T_ad[0],T_ad[-1]])
    ax.set_ylim([0, 5000])
    ax.set_title('temperature')
    ax.set_xlabel('T  [K]')
    ax = axes[1,2]
    ax.plot(rho_ad, z_, '-', label='analytic')
    # ax.plot([0.5, 3.0], [2000, 2000], 'k')
    # ax.plot([3.0, 3.0],[0, 5000], 'k')
    # ax.plot([rho[np.int(2000. / dz_)], rho[np.int(2000. / dz_)]], [0, 5000], 'k')
    ax.plot(rho0_stats, z_stats, '-', label='stats')
    ax.legend()
    ax.set_xlim([rho_ad[0],rho_ad[-1]])
    ax.set_ylim([0, 5000])
    ax.set_title('rho')
    ax.set_xlabel('rho  [kg/m3]')
    plt.savefig('./thermodynamic_profiles.png')
    plt.close()






    def kmax(i_arr,j_arr,ic,jc,rstar,zstar):
        # ir = np.int(np.round(np.sqrt((i - ic) ** 2 + (j - jc) ** 2)))
        irstar = np.int(np.round(rstar / dx))
        kstar = np.int(np.round(zstar / dz))
        if type(i_arr) == np.ndarray:
            k_max = np.zeros((i_arr.shape[0], j_arr.shape[0]))
            for i in i_arr:
                for j in j_arr:
                    r = np.sqrt( (i*dx-ic*dx)**2 + (j*dy-jc*dy)**2 )
                    ir = np.round(np.sqrt((i - ic) ** 2 + (j - jc) ** 2))
                    if r <= rstar:
                        k_max[i,j] = np.int(np.round(kstar * (np.cos(np.double(r) / rstar * np.pi / 2)) ** 2))
                    else:
                        k_max[i,j] = 0
        else:
            ir = np.round(np.sqrt((i_arr - ic) ** 2 + (j_arr - jc) ** 2))
            r = np.round(np.sqrt((i_arr*dx - ic*dx) ** 2 + (j_arr*dy - jc*dy) ** 2))
            if ir <= irstar:
                k_max = np.int(np.round(kstar * (np.cos(np.double(ir) / irstar * np.pi / 2)) ** 2))
            else:
                k_max = 0
        return k_max


    def zmax(x,y,xc,yc,rstar,zstar):
        r = np.sqrt((x - xc) ** 2 + (y - yc) ** 2)
        if r <= rstar:
            return zstar * (np.cos(np.double(r) / rstar * np.pi / 2)) ** 2
        else:
            return 0.0

    nx_half = np.int(nx/2)
    ny_half = np.int(ny/2)
    nz_half = np.int(zstar/dz)+10
    i_arr = np.arange(0,nx_half)
    j_arr = np.arange(0,ny_half)
    x_arr = dx*np.arange(0,nx_half)
    y_arr = dy*np.arange(0,ny_half)
    z_arr = dz*np.arange(0,nz_half)
    ic = np.int(nx/4)
    jc = np.int(ny/4)
    xc = dx*ic
    yc = dy*jc

    # kmax(i_arr,j_arr,ic,jc,1000,1000)

    # varying height; varying radius
    zstar_range = [4000, 2000, 1000, 500, 250]
    rstar_range = [250, 500, 1000, 2000, 4000]
    nzs = len(zstar_range)
    nrs = len(rstar_range)
    V_const = np.zeros(nzs)    # constant rho
    V_rho = np.zeros(nzs)      # reference profile rho
    V_const2 = np.zeros(nzs)    # constant rho
    k_max = np.zeros((nzs, nx_half, ny_half))
    z_max = np.zeros((nzs, nx_half, ny_half))
    rho0_mat = np.tile(rho0_stats, [nx_half, 1])
    min = rho0_stats[0]
    max = rho0_stats[np.int(np.amax(zstar_range)/dz)]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    lvls = np.linspace(rho0_stats[0], rho0_stats[nz_half], 10)
    cf = ax1.contourf(x_arr[10:120], z_arr, rho0_mat[10:120, :np.int(zstar / dz) + 10].T, alpha=0.7)  # , levels=lvls)
    plt.colorbar(cf, ax=ax1, shrink=0.75)
    for count in range(nzs):
        zstar_ = zstar_range[count]
        rstar_ = rstar_range[count]
        # kmax(i_arr,j_arr,ic,jc,1000,1000)
        for i in range(nx_half):
            for j in range(ny_half):
                for k in range(nz_half):
                    k_max[count, i,j] = kmax(i, j, ic, jc, rstar_, zstar_)
                    if k <= k_max[count, i, j]:
                        z_max[count, i,j] = zmax(i*dx, j*dy, xc, yc, rstar_, zstar_)
                        # z = k*dz
                        V_const[count] += (k*dz)*dV
                        V_rho[count] += rho0_stats[k]*(k*dz)*dV
                        # V_const[count] += (k*dz)*dV
                        # V_rho[count] += rho0_stats[k]*(k*dz)*dV
                    if k*dz <= z_max[count, i, j]:
                        V_const2[count] += (k*dz)*dV

        rho_test = np.linspace(min, max, 10)
        cl = str(np.double(count)/(nzs))
        ax1.plot(x_arr[10:120], z_max[count, 10:120, jc], linewidth=3, label='zstar='+str(zstar_) + ', rstar=' + str(rstar_), color=cl)
        ax1.plot(x_arr[10:120], dz*k_max[count, 10:120, jc], '--', linewidth=2, label='zstar='+str(zstar_) + ', rstar=' + str(rstar_), color=cl)
        # if count == 0:
        #     ax2.plot(V_const[count], rho_test, '-x', color=cl, linewidth=2, label= 'V = rho0 * int(z dz)')
        if count == 0:
            ax2.plot(rho_test*V_const[count], rho_test, '-x', color=cl, linewidth=2, label= 'V = rho0 * int(z dz)')
            ax2.plot([V_rho[count], V_rho[count]], [min, max], color=cl, linewidth=2, label='V = int( rho0(z)*z dz )')
            ax3.plot(rho_test*V_const2[count], rho_test, '-x', color=cl, linewidth=2, label= 'V = rho0 * int(z dz)')
        else:
            ax2.plot(rho_test*V_const[count], rho_test, '-x', color=cl, linewidth=2)
            ax2.plot([V_rho[count], V_rho[count]], [min, max], color=cl, linewidth=2)
            ax3.plot(rho_test*V_const2[count], rho_test, '-x', color=cl, linewidth=2)
        rho0_half = rho0_stats[np.int(np.round(np.amax(k_max[count, :, :]) / 2))]
        ax2.plot(rho0_half*V_const[count], rho0_half, 'or')
        ax3.plot(rho0_stats[np.int(zstar_/(2*dz))]*V_const2[count], rho0_stats[np.int(zstar_/dz)], 'ro')
    ax1.legend(loc='best', fontsize=8)
    ax1.set_ylabel('z  [m]')
    ax1.set_title('rstar = '+str(rstar)+'m')
    ax2.legend(loc='best', fontsize=8)
    ax3.legend(loc='best', fontsize=8)
    ax2.set_xlabel('PE / (g dTh / Th0)')
    ax2.set_ylabel('rho0')
    ax2.set_title('cos(ir)')
    ax3.set_xlabel('PE / (g dTh / Th0)')
    ax3.set_title('cos(r)')
    fig.savefig('./rho_initial_envelope_zvar_rvar.png')
    plt.close(fig)

    k_max = kmax(i_arr,j_arr,ic,jc,1000,1000)
    print('shapes', k_max.shape, x_arr.shape)
    # plt.figure()
    # plt.plot(k_max)
    # plt.show()

    # def PE_density_const(r, z, rstar, zstar):
    def dPdV(PE, z):
        # zstar = 1000
        # aux = zstar**2
        aux = z
        return aux

    def test_ode(y,t):
        k = 0.3
        dydt = -k * y
        return dydt

    y0 = 5
    z = np.linspace(0,20)
    y = odeint(test_ode, y0, z)

    z0 = 0
    PE0 = 0.
    z_ = np.arange(0,5)
    print(z_.shape)
    V0 = odeint(dPdV, PE0, z_)


    # zstar_range = [500, 1000, 1500, 2000, 2500]
    # rstar_range = [250, 500, 1000, 1500, 2000, 2500]
    #
    # # varying height; constant radius
    # # (1) with constant rho
    # V_const = np.zeros(len(zstar_range))
    # # (2) with reference profile rho
    # V_rho = np.zeros(len(zstar_range))
    # k_max = np.zeros((len(zstar_range), nx_half, ny_half))
    # z_max = np.zeros((len(zstar_range), nx_half, ny_half))
    # rho0_mat = np.tile(rho0_stats, [nx_half, 1])
    #
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    # lvls = np.linspace(rho0_stats[0], rho0_stats[nz_half], 10)
    # cf = ax1.contourf(x_arr, z_arr, rho0_mat[:, :np.int(zstar / dz) + 10].T, alpha=0.7)  # , levels=lvls)
    # plt.colorbar(cf, ax=ax1)
    # for count,zstar_ in enumerate(zstar_range):
    #     for i in range(nx_half):
    #         for j in range(ny_half):
    #             # pass
    #             for k in range(nz_half):
    #                 k_max[count, i,j] = kmax(i, j, ic, jc, rstar, zstar_)
    #                 if k <= k_max[count,i,j]:
    #                     z_max[count, i,j] = zmax(i*dx, j*dy, xc, yc, rstar, zstar_)
    #                     V_const[count] += (k*dz)*dV
    #                     V_rho[count] += rho0_stats[k]*(k*dz)*dV
    #     min = rho0_stats[0]
    #     max = rho0_stats[np.amax(k_max)]
    #     rho_test = np.linspace(min, max, 10)
    #     cl = str(np.double(count)/len(zstar_range))
    #     ax1.plot(x_arr, z_max[count, :,jc], linewidth=3, label='zstar='+str(zstar_), color=cl)
    #     ax2.plot(rho_test*V_const[count], rho_test, '-x', color=cl, linewidth=2, label= 'V = rho0 * int(z dz)')
    #     ax2.plot([V_rho[count], V_rho[count]], [min, max], color=cl, linewidth=2, label='V = int( rho0(z)*z dz )')
    #     rho0_half = rho0_stats[np.int(np.round(np.amax(k_max[count, :, :]) / 2))]
    #     ax2.plot(rho0_half*V_const[count], rho0_half, 'or')
    # ax1.legend(loc='best', fontsize=8)
    # ax1.set_ylabel('z  [m]')
    # ax1.set_title('rstar = '+str(rstar)+'m')
    # ax2.legend(loc='best', fontsize=8)
    # ax2.set_xlabel('PE / (g dTh / Th0)')
    # ax2.set_ylabel('rho0')
    # fig.savefig('./rho_initial_envelope_zvar.png')
    # plt.close(fig)
    #
    #
    # # constant height; varying radius
    # zstar_ = 1000
    # V_const = np.zeros(len(rstar_range))
    # V_rho = np.zeros(len(rstar_range))
    # k_max = np.zeros((len(rstar_range), nx_half, ny_half))
    # z_max = np.zeros((len(rstar_range), nx_half, ny_half))
    #
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    # lvls = np.linspace(rho0_stats[0], rho0_stats[nz_half], 10)
    # cf = ax1.contourf(x_arr, z_arr, rho0_mat[:, :np.int(zstar / dz) + 10].T, alpha=0.7)  # , levels=lvls)
    # plt.colorbar(cf, ax=ax1)
    # for count,rstar_ in enumerate(rstar_range):
    #     for i in range(nx_half):
    #         for j in range(ny_half):
    #             # pass
    #             for k in range(nz_half):
    #                 k_max[count, i,j] = kmax(i, j, ic, jc, rstar_, zstar_)
    #                 if k <= k_max[count,i,j]:
    #                     z_max[count, i,j] = zmax(i*dx, j*dy, xc, yc, rstar_, zstar_)
    #                     V_const[count] += (k*dz)*dV
    #                     V_rho[count] += rho0_stats[k]*(k*dz)*dV
    #     min = rho0_stats[0]
    #     max = rho0_stats[np.amax(k_max)]
    #     rho_test = np.linspace(min, max, 10)
    #     cl = str(np.double(count)/len(rstar_range))
    #     ax1.plot(x_arr, z_max[count, :,jc], linewidth=3, label='rstar='+str(rstar_), color=cl)
    #     ax2.plot(rho_test*V_const[count], rho_test, '-x', color=cl, linewidth=2, label= 'V = rho0 * int(z dz)')
    #     ax2.plot([V_rho[count], V_rho[count]], [min, max], color=cl, linewidth=2, label='V = int( rho0(z)*z dz )')
    #     rho0_half = rho0_stats[np.int(np.round(np.amax(k_max[count, :, :]) / 2))]
    #     ax2.plot(rho0_half*V_const[count], rho0_half, 'or')
    # ax1.legend(loc='best', fontsize=8)
    # ax1.set_ylabel('z  [m]')
    # ax1.set_title('zstar = '+str(zstar_)+'m')
    # ax2.legend(loc='best', fontsize=8)
    # ax2.set_xlabel('PE / (g dTh / Th0)')
    # ax2.set_ylabel('rho0')
    # fig.savefig('./rho_initial_envelope_rvar.png')
    # plt.close(fig)


    return





if __name__ == '__main__':
    main()