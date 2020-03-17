import os
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import netCDF4 as nc
import json as simplejson

# import ../thermodynamic_profiles
# from thermodynamic_profiles import eos
from thermodynamic_functions import thetas_c
# from ../Initialization_PE.py import compute_PE, compute_envelope

execfile('settings.py')
label_size = 15
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 21
# plt.rcParams['xtick.direction']='out'
# plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['grid.linewidth'] = 20


def main():
    # # path_out_figs = '/nbi/ac/cond1/meyerbe/paper_CP_single'
    # path_out_figs = '/nbi/home/meyerbe/paper_CP'
    path_data_dx100m = '/nbi/ac/conv4/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/single_3D_noise/run5_PE_scaling_dx100m'
    path_data_ref_dx100m = '/nbi/ac/cond1/meyerbe/ColdPools/3D_sfc_fluxes_off/single_3D_noise/run2_dx100m/dTh3_z1000_r1000'
    path_data_dx50m = '/nbi/ac/conv4/rawdata/ColdPools_PyCLES/3D_sfc_fluxes_off/single_3D_noise/run6_PE_scaling_dx50m'


    'geometry'
    global nz, dz
    dz = 25
    nz = 2e3/dz
    # nz = 3

    ''' parameters '''
    global Rd, Rv, eps_v, eps_vi
    Rd = 287.1
    Rv = 461.5
    eps_v = 0.62210184182
    eps_vi = 1.60745384883
    global g, cpd, cpv, kappa, Lv
    g = 9.80665
    cpd = 1004.0
    cpv = 1859.0
    Lv = 2264.705e3       # sensible latent heat of water evaporation
    kappa = 0.285956175299
    global rho_w
    rho_w = 997 #kg/m3
    global T_tilde, p_tilde, sd_tilde
    T_tilde = 298.15
    p_tilde = 100000.0
    sd_tilde = 6864.8

    '''surface values'''
    global Pg, Tg, qtg
    Pg = 1.0e5
    Tg = 300.0
    qtg = 0.0
    # sg = entropy_dry(Pg, Tg, qtg, 0.0, 0.0)
    # sg = 6e3

    ''' reference / background profile '''
    # compute reference pressure and density profiles (assuming dry thermodynamics)
    p_ref, al_d = compute_pressure_profile_dry()
    rho_d = 1./al_d


    ''' evaporation parameters '''
    z0 = 800.       # height of evaporation
    theta0 = 298.   # temperature at level of evaporation
    q0 = 298.       # moisture at level of evaporatio
    evap = 0.1      # fraction of rain water that is evaporated
    ''' rain event '''
    z_BL = 1e3      # height of sub-cloud layer [m]
    A = 1e3**2      # area of precipitation cell
    intensity = 5   # [mm/h]
    tau = 1         # duration of rain event



    ''' compute PE from reference simulation (r*=z*=1km, dTh=5K) '''
    print('')
    print('------ compute PE for reference simulation -----')
    ''' (a) analytically '''
    ''' (b) numerical '''
    ''' (c) from output field '''
    rstar = 1000.
    zstar = 1000.
    # dx = 100
    nml = simplejson.loads(open(os.path.join(path_data_ref_dx100m, 'ColdPoolDry_single_3D.in')).read())
    nx = nml['grid']['nx']
    # ny = nml['grid']['ny']
    # nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dV = dx**3
    ic = nml['init']['ic']
    jc = nml['init']['jc']
    marg = nml['init']['marg']
    imin = ic - np.int((rstar + 1000)/dx)
    imax = nx - imin
    ni = imax - imin
    ni2 = ni**2
    kmax = np.int((zstar+1000.)/dx)
    root = nc.Dataset(os.path.join(path_data_ref_dx100m, 'fields', '0.nc'))
    temp_ = root.groups['fields']['temperature'][imin:imax,imin:imax,:kmax]
    s_ = root.groups['fields']['s'][imin:imax,imin:imax,:kmax]
    root.close()
    theta_ = thetas_c(s_, 0.0)
    root = nc.Dataset(os.path.join(path_data_ref_dx100m, 'Stats.ColdPoolDry_single_3D.nc'))
    rho0_stats = root.groups['reference'].variables['rho0'][:kmax]
    # alpha0_stats = root.groups['reference'].variables['alpha0'][:kmax]
    zhalf_stats = root.groups['reference'].variables['z'][:kmax]
    # rho_unit = root.groups['reference'].variables['rho0'].units
    root.close()
    Th_ref = Tg     # at surface, temperature equal to pot. temp. (Tg=Theta_g)
    print('imin, imax, kmax: ', imin, imax, kmax)
    print('shape: ', ni, 200, theta_.shape)
    print('reference temperature: ', Th_ref, Tg, theta_[10,10,10])
    print('reference density: ', np.amax(np.abs(rho_d[:kmax]-rho0_stats[:kmax])))
    theta = theta_.reshape(ni2,kmax)
    print('reshaping: ', theta.shape, theta_.shape, kmax)
    PE0 = 0.
    for i in range(ni2):
        for k in range(kmax):
            PE0 += (theta[i,k]-Th_ref)*zhalf_stats[k]*rho_d[k]*dV
    PE0 = g/Th_ref*PE0
    print('')

    print('POTENTIAL ENERGY')
    print('(c) from output field theta: ' + str(PE0))




    ''' (1) compute precipitation, given potential energy '''
    ''' dry '''
    PE_ref = -0.8e11
    # PE_ref = -0.5e11
    # evaporated water:     V * rho_w * evap
    # volume:               V = I * A * tau, A: area, I: intensity, tau: duration
    # released heat:        Q = Lv * V * rho_w * evap
    # temperature change:   Q = dT * cp <<>> dT = Q / cp
    # >> Lv * V * rho_w * evap = dT * cp >> V = dT * cp / (Lv * rho_w * evap)
    print('tau: ' + str(tau) + ' h, area: '+str(A) + 'm2')
    I_ref = compute_intensity_from_PE(PE_ref, rho_d, p_ref, tau, A, z0, z_BL, theta0)
    P = I_ref * A * tau
    height = P/A
    print('   Precip. intensity:  ' + str(np.round(I_ref*1e3,0)) + ' mm/h')
    print('   Precip. total:      ' + str(np.round(P,2)) + ' m^3')
    print('   height watercolumn: ' + str(np.round(1e3*height,0)) + ' mm')
    # This corresponds to the cooling generated by the fractional evaporation of rain from a cell of area 1km2 that
    # precipitates at moderate intensity of 5 mm/h for half an hour at an evaporation rate of 10%
    tau_ = height / (intensity*1e-3)
    intensity_ = height / tau * 1e3
    print('if intensity is '+str(intensity)+'mm/h, then duration is: tau='+str(np.round(tau_,2))+' h')
    print('if duration is '+str(tau)+'h, then intensity is: I='+str(np.round(intensity_,2))+' mm/h')

    ''' (2) plot histogram '''
    print('')
    print('PE range: ')
    PE_range = 2. ** np.arange(-1, 4)
    print(PE_range)
    print('tau: ' + str(tau) + ' h, area: ' + str(A) + 'm2')
    I = compute_intensity_from_PE(PE_range * PE_ref, rho_d, p_ref, tau, A, z0, z_BL, theta0)
    # P = -PE_ref*PE_range * theta0 / (g * z0 * rho_d[k0]) * cpd / (rho_w * evap) / exner_c(p[k0])
    P = I * A * tau
    height = P / A
    print('   Precip. intensity:  ' + str(np.round(I * 1e3, 0)) + ' mm/h')
    print('   Precip. total:      ' + str(np.round(P, 2)) + ' m^3')
    print('   height watercolumn: ' + str(np.round(1e3 * height, 0)) + ' mm')
    plot_histogram_PE_vs_Intensity(PE_range, P, PE_ref, I, I_ref, path_out_figs)


    ''' (3) plot PE vs. R (run5)'''
    dTh_ref = 3
    rstar_ref = 1000
    zstar_ref = 1000
    # run5
    dTh = 5
    r_params = [500, 1100, 1600, 2300]
    r_params_ = [500, 1000, 1100, 1600, 2300]
    z_params = [1000]
    PE_array = [0.5, 2, 4, 8]
    PE_array_log = 2. ** np.arange(-1, 4)
    # print('PE: ' + str(PE_array))
    n_params = len(r_params)
    plot_PE_vs_R(r_params, z_params, n_params, dTh,
                 rstar_ref, zstar_ref, dTh_ref,
                 PE_array, PE_array_log, I,
                 path_out_figs)


    # ''' moist '''
    #
    # T0 = theta
    # print('!!')
    # qt0 = 1e-3
    # qv0 = qt0
    # # evaporated rain
    # # dqv =
    #
    #
    # # rho = 1./alpha(T, p, qt, qv)
    # buoy = g*rho_d[k]*(alpha_tot(T0, p[k], qt0, qv0))

    return



# -----------------------------------------
def plot_histogram_PE_vs_Intensity(PE_range, P, PE_ref, I, I_ref, path_out_figs):

    fig, [ax0, ax1, ax2] = plt.subplots(1, 3, figsize=(20, 6))
    ax0.plot(PE_range, P, '-o')
    ax0.set_xlabel(r'PE / PE$_0$')
    ax1.plot(PE_ref * PE_range, P, '-o')
    ax1.set_xlabel('PE / PE_ref')
    ax2.bar(PE_range, I * 1e3, align='center', facecolor='lightblue')
    ax2.set_xlabel('PE / PE_ref')
    ax2.set_ylabel('Intensity [mm/h]')
    plt.subplots_adjust(bottom=0.2, right=.95, left=0.07, top=0.9, wspace=0.25)
    plt.savefig('./preciptation_run5.png')

    print(I * 1e3)
    aux = PE_range[0] * np.ones(np.int(I[0] * 1e3))
    aux = np.append(aux, np.ones(np.int(I_ref * 1e3)))
    aux = np.append(aux, PE_range[1] * np.ones(np.int(I[1] * 1e3)))
    aux = np.append(aux, PE_range[2] * np.ones(np.int(I[2] * 1e3)))
    aux = np.append(aux, PE_range[3] * np.ones(np.int(I[3] * 1e3)))
    aux = np.append(aux, PE_range[4] * np.ones(np.int(I[4] * 1e3)))
    print(I * 1e3)
    print(PE_range)
    print(aux)

    fig_name = 'precipitation_run5_hist.png'
    fig, [ax0, ax1, ax2] = plt.subplots(1, 3, figsize=(20, 6))
    ax0.plot(PE_range, np.round(I * 1e3, 0), '-o')
    ax0.set_xlabel(r'PE / PE$_0$')
    ax1.hist(aux)
    ax1.set_xlabel(r'PE / PE$_0$')
    ax1.set_ylabel('Intensity [mm/h]')
    # ax2.hist(np.round(I * 1e3, 0), '-o')
    ax2.bar(PE_range, I * 1e3, align='center', facecolor='lightblue')
    # plt.bar(range(len(D)), D.values(), align='center')
    ax2.set_xlabel(r'PE / PE$_0$')
    ax2.set_ylabel('Intensity [mm/h]')
    ax2.xaxis.get_ticklabels()[3].set_visible(False)
    ax2.xaxis.get_ticklabels()[5].set_visible(False)
    ax2.xaxis.get_ticklabels()[6].set_visible(False)
    ax2.xaxis.get_ticklabels()[7].set_visible(False)
    ax2.xaxis.get_ticklabels()[9].set_visible(False)
    for ax in [ax0, ax1, ax2]:
        ax.set_xticklabels([np.int(ti) for ti in ax.get_xticks()])
        ax.set_yticklabels([np.int(ti) for ti in ax.get_yticks()])
    plt.subplots_adjust(bottom=0.2, right=.95, left=0.07, top=0.9, wspace=0.25)
    #plt.savefig('./preciptation_run5_hist.png')
    plt.savefig(os.path.join(path_out_figs, fig_name))
    return
# -----------------------------------------
def plot_PE_vs_R(r_params, z_params, n_params, dTh, rstar_ref, zstar_ref, dTh_ref,
                 PE_array, PE_array_log, intensity,
                 path_out_figs):
    fig_name = 'PE_scaling_run5.png'

    r_params_ = list(r_params)
    i = 0
    while (rstar_ref>r_params[i]) and (i<len(r_params)):
        i+=1
    r_params_.insert(i, rstar_ref)
    print('TESTING: ', i, r_params[i])
    print(r_params)
    print(r_params_)

    ''' envelope '''
    Lx = 6e3
    dx_ = 10
    nx_ = Lx / dx_
    x_arr = np.arange(0, Lx, dx_)
    ic = np.int(nx_ / 2)
    xc = x_arr[ic]
    # print nx_, ic, xc
    #
    zmax = np.zeros((n_params + 1, nx_))
    for istar in range(n_params):
        zstar = z_params[0]
        rstar = r_params[istar]
        irstar = np.int(np.double(rstar) / dx_)
        x_arr_ = np.array(x_arr, copy=True)
        x_arr_[:ic - irstar] = rstar + xc
        x_arr_[ic + irstar:] = rstar + xc
        zmax[istar, :] = zstar * np.cos((x_arr_ - xc) / rstar * np.pi / 2) ** 2
    x_arr_ = np.array(x_arr, copy=True)
    irstar = np.int(np.double(rstar_ref) / dx_)
    x_arr_[:ic - irstar] = rstar_ref + xc
    x_arr_[ic + irstar:] = rstar_ref + xc
    zmax[-1, :] = zstar_ref * np.cos((x_arr_ - xc) / rstar_ref * np.pi / 2) ** 2

    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(20, 6))
    for istar in range(n_params):
        zstar = z_params[0]
        rstar = r_params[istar]
        ax0.plot(x_arr - xc, zmax[istar, :], label='dTh' + str(dTh) + ', z*' + str(zstar) + ', r*' + str(rstar))
    ax0.plot(x_arr - xc, zmax[-1], 'k',
             label='dTh' + str(dTh_ref) + ', z*' + str(zstar_ref) + ', r*' + str(rstar_ref))
    ax0.set_ylabel('Height z  [km]', fontsize=21)
    ax0.set_ylim(0, z_params[-1] + 100)
    ax0.legend(loc='upper right', bbox_to_anchor=(-0.2, 1.0),
               fancybox=True, ncol=1)

    ax1.plot(r_params_, PE_array_log, 'k', linewidth=0.5)
    ax1.plot(r_params, PE_array, '--k', linewidth=0.5)
    for istar in range(n_params):
        # #     axes[0].plot(np.log2(PE_array[istar]), r_params_[istar], 'o', markersize=10, markeredgecolor='w', )
        ax1.plot(r_params[istar], PE_array[istar], 'o', markersize=10, markeredgecolor='w', )
    ax1.plot(r_params_[1], PE_array_log[1], 'ko', markersize=10, markeredgecolor='w', )

    ax2.bar(PE_array_log, intensity * 1e3, align='center', facecolor='lightblue')
    ax2.set_xlabel(r'PE / PE$_0$')
    ax2.set_ylabel('Intensity [mm/h]')
    ax2.xaxis.get_ticklabels()[3].set_visible(False)
    ax2.xaxis.get_ticklabels()[5].set_visible(False)
    ax2.xaxis.get_ticklabels()[6].set_visible(False)
    ax2.xaxis.get_ticklabels()[7].set_visible(False)
    ax2.xaxis.get_ticklabels()[9].set_visible(False)
    # for ax in [ax0, ax1, ax2]:
    #     ax.set_xticklabels([np.int(ti) for ti in ax.get_xticks()])
    #     ax.set_yticklabels([np.int(ti) for ti in ax.get_yticks()])
    ax2.set_xticklabels([np.int(ti) for ti in ax2.get_xticks()])
    ax2.set_yticklabels([np.int(ti) for ti in ax2.get_yticks()])

    ax1.set_ylabel(r'PE / PE$_0$')
    ax1.set_xlim(400, 2500)
    ax1.set_ylim(0, 8.5)
    for ax in [ax0, ax1]:
        ax.set_xlabel('Radius r  [km]')
    ax0.set_xticklabels([np.int(ti*1e-3) for ti in ax0.get_xticks()])
    ax1.set_xticklabels([np.round(ti*1e-3,1) for ti in ax1.get_xticks()])
    ax0.set_yticklabels([np.round(ti*1e-3,0) for ti in ax0.get_yticks()])
    ax1.set_yticklabels([np.int(ti) for ti in ax.get_yticks()])
    for label in ax1.yaxis.get_ticklabels()[3::2]:
        label.set_visible(False)
    ax0.yaxis.get_ticklabels()[1].set_visible(False)
    ax0.yaxis.get_ticklabels()[2].set_visible(False)
    ax0.yaxis.get_ticklabels()[3].set_visible(False)
    ax0.yaxis.get_ticklabels()[4].set_visible(False)
    ax1.yaxis.get_ticklabels()[6].set_visible(False)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, right=.85, left=0.2, top=0.9, wspace=0.25)
    plt.savefig(os.path.join(path_out_figs, fig_name))
    plt.close(fig)
    return
# -----------------------------------------
def compute_intensity_from_PE(PE0, rho_d, p, tau, A, z0, z_BL, theta0):
    ''' evaporation parameters '''
    # tau         # duration of rain event
    # A           # area of precipitation cell
    # z0          # height of evaporation
    # z_BL        # height of sub-cloud layer [m]
    # theta0      # temperature at level of evaporation
    # p           # reference pressure
    # rho_d       # reference density profile
    k0 = np.int(z0 / dz)
    evap = 0.1  # fraction of rain water that is evaporated

    V = A*z_BL
    dTh = PE0 * theta0 / (g*z0*rho_d[k0]*V)
    dT = dTh / exner_c(p[k0])
    I = cpd * (-dT) * (rho_d[k0]*A*z_BL) / (evap*rho_w*tau*A*Lv)
    print('   PE ref:             ' + str(PE0))
    print('   dTh:                ' + str(np.round(dTh, 2)) + ' K')
    print('   dT:                 ' + str(np.round(dT, 2)) + ' K')
    return I
# -----------------------------------------
def compute_pressure_profile_dry():
    # compute entropy from pressure and temperature
    # sg = Thermodynamics.entropy(self.Pg, self.Tg, self.qtg, 0.0, 0.0)
    def entropy_dry(pd, T, qt, ql, qi):
        print 'entropy dry'
        sd_c = sd_tilde + cpd * np.log(T / T_tilde) - Rd * np.log(pd / p_tilde)
        return sd_c

    # compute Temperature from pressure and entropy (dry)
    def eos(pd, s, qt):
        ql = np.zeros(pd.shape)
        qi = np.zeros(pd.shape)
        eos_c = T_tilde * (np.exp((s - sd_tilde + Rd * np.log(pd / p_tilde)) / cpd))
        return eos_c, ql, qi

    def rhs(p, z):
        ql = 0.0
        qi = 0.0
        # # given sfc values for pressure, temperature and moisture
        # # >> compute sfc entropy (= constant value throught atmosphere for reference profile being defines as constant-entropy profile)
        # # compute temperature from pressure and entropy (at any given height)
        T, ql, qi = eos(np.exp(p), sg, qtg)
        rhs_ = -g / (Rd * T * (1.0 - qtg + eps_vi * (qtg - ql - qi)))
        return rhs_

    sg = entropy_dry(Pg, Tg, qtg, 0.0, 0.0)
    p0 = np.log(Pg)
    # Construct arrays for integration points
    z_ = np.array(np.arange(nz) * dz)
    # Perform the integration
    p = odeint(rhs, p0, z_, hmax=1.0)[:, 0]  # type: object
    # p_half = odeint(rhs, p0, z_half, hmax=1.0)[1:, 0]
    p = np.exp(p)
    # temperature[k], ql[k], qi[k] = Thermodynamics.eos(p_[k], self.sg, self.qtg)
    # print('...', sg, qtg)
    T, ql, qi = eos(p, sg, qtg)
    # qv[k] = self.qtg - (ql[k] + qi[k])
    qv = np.zeros(p.shape)
    # alpha[k] = Thermodynamics.alpha(p_[k], temperature[k], self.qtg, qv[k])
    al = alpha_tot(T, p, qtg, qv)
    # al = alpha_dry(T, p)

    fig, axes = plt.subplots(1, 3, figsize=(24, 10))
    ax = axes[0]
    ax.plot(p, z_, '-')
    ax.set_title('pressure')
    ax = axes[1]
    ax.plot(T, z_, '-')
    ax.set_title('temperature')
    ax = axes[2]
    ax.plot(1./al, z_, '-')
    ax.plot([np.amin(1./al), np.amax(1./al)], [800, 800], 'k')
    # # ax.plot([3.0, 3.0],[0, 5000], 'k')
    ax.plot([1./al[np.int(800. / dz)], 1./al[np.int(800. / dz)]], [0, 2e3], 'k')
    ax.set_title('rho')
    plt.savefig('./thermodynamic_profiles_pycles.png')
    plt.close()

    return p, al


# -----------------------------------------
def compute_envelope(dTh, rstar, zstar, marg, ic, jc, th_g, dx, nx, ny, kmax):
    z_max_arr = np.zeros((2, nx, ny), dtype=np.double)
    theta_z = th_g * np.ones(shape=(nx, ny, kmax))
    # theta_pert = np.random.random_sample(npg)

    x_half, y_half, z_half = compute_dimension_arrays(dx, dx, dx, nx, ny, kmax)
    xc = x_half[ic]  # center of cold-pool
    yc = y_half[jc]  # center of cold-pool


    for i in range(nx):
        # ishift = i * nlg[1] * nlg[2]
        for j in range(nx):
            # jshift = j * nlg[2]
            r = np.sqrt((x_half[i] - xc) ** 2 +
                        (y_half[j] - yc) ** 2)
            if r <= rstar:
                z_max = zstar * (np.cos(r / rstar * np.pi / 2) ** 2)
                z_max_arr[0, i, j] = z_max
            if r <= (rstar + marg):
                z_max = (zstar + marg) * (np.cos(r / (rstar + marg) * np.pi / 2) ** 2)
                z_max_arr[1, i, j] = z_max

            # kstar = np.int(np.round(zstar / dz))
            for k in range(kmax):
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


def compute_dimension_arrays(dx, dy, dz, nx, ny, nz):
    global x_half, y_half, z_half
    x_half = np.empty(nx, dtype=np.double, order='c')
    y_half = np.empty(ny, dtype=np.double, order='c')
    z_half = np.empty(nz, dtype=np.double, order='c')
    count = 0
    for i in range(nx):
        x_half[count] = (i + 0.5) * dx
        count += 1
    count = 0
    for j in range(ny):
        y_half[count] = (j + 0.5) * dy
        count += 1
    count = 0
    for i in range(nz):
        z_half[count] = (i + 0.5) * dz
        count += 1
    return x_half, y_half, z_half
# -----------------------------------------

def alpha_tot(T, p0, qt, qv):
    # qd*alpha_dry + qv*alpha_moist, with qd=1-qt
    return Rd*T/p0 * (1.-qt+qv*eps_vi)

def alpha_dry(T, p0):
    return Rd*T / p0

def alpha_moist(T, p0):
    return Rv*T/p0

def exner_c(p0):
    return pow((p0/p_tilde),kappa)

# -----------------------------------------

if __name__ == '__main__':
    main()
