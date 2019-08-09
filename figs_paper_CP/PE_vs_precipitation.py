import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


# import ../thermodynamic_profiles
# from thermodynamic_profiles import eos

def main():
    'geometry'
    global nz, dz
    dz = 25
    nz = 2e3/dz
    # nz = 3

    # define parameters
    global Rd, Rv, eps_v, eps_vi
    Rd = 287.1
    Rv = 461.5
    eps_v = 0.62210184182
    eps_vi = 1.60745384883
    global g, cpd, cpv, kappa
    g = 9.80665
    cpd = 1004.0
    cpv = 1859.0
    kappa = 0.285956175299

    global rho_w
    rho_w = 997 #kg/m3
    global T_tilde, p_tilde, sd_tilde
    T_tilde = 298.15
    p_tilde = 100000.0
    sd_tilde = 6864.8

    'surface values'
    global Pg, Tg, qtg
    Pg = 1.0e5
    Tg = 300.0
    qtg = 0.0
    # sg = entropy_dry(Pg, Tg, qtg, 0.0, 0.0)
    # sg = 6e3

    # compute reference pressure and density profiles (assuming dry thermodynamics)
    p, al_d = compute_pressure_profile_dry()
    rho_d = 1./al_d

    ''' evaporation parameters '''
    # height of evaporation
    z0 = 800.
    k = z0/dz
    # temperature and moisture at level of evaporation
    theta0 = 298.
    evap = 0.1      # fraction of rain water that is evaporated



    ''' dry '''
    PE0 = -0.8e11
    P = -PE0 * theta0 / (g*z0*rho_d[k]) * cpd / (rho_w * evap) / exner_c(p[k])
    print('PE ref: ' + str(PE0))
    print('Precipitation: ' + str(P) +'m^3')
    print rho_d[k]


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
    print('...', sg, qtg)
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