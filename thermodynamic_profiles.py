import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

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
    T_tilde = 298.15
    p_tilde = 100000.0
    pv_star_t = 611.7
    sd_tilde = 6864.8
    sv_tilde = 10513.6

    nz = 500
    dz = 10
    # nz = 10
    # dz = 500


    'surface values'
    global Pg, Tg, qtg, sg
    Pg = 1.0e5
    Tg = 300.0
    qtg = 0.0
    # sg = entropy_dry(Pg, Tg, qtg, 0.0, 0.0)
    sg = 6e3

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
    # def eos_c(pd, s):
    #     return T_tilde * (exp((s - sd_tilde + Rd * log(pd / p_tilde)) / cpd))
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
    z = np.array(np.arange(nz)*dz)
    # Perform the integration
    p = odeint(rhs, p0, z, hmax=1.0)[:, 0]
    # p_half = odeint(rhs, p0, z_half, hmax=1.0)[1:, 0]
    p = np.exp(p)
    # temperature[k], ql[k], qi[k] = Thermodynamics.eos(p_[k], self.sg, self.qtg)
    T, ql, qi = eos(p, sg, qtg)
    print('shapes', z.shape, T.shape)
    # qv[k] = self.qtg - (ql[k] + qi[k])
    qv = np.zeros(p.shape)
    # alpha[k] = Thermodynamics.alpha(p_[k], temperature[k], self.qtg, qv[k])
    alpha = alpha_c(p, T, qtg, qv)
    rho = 1./alpha



    fig, axes = plt.subplots(1,3)
    ax = axes[0]
    ax.plot(p, z, '-')
    ax.set_title('pressure')
    ax = axes[1]
    ax.plot(T, z, '-')
    ax.set_title('temperature')
    ax = axes[2]
    ax.plot(rho, z, '-')
    ax.set_title('rho')
    plt.savefig('./thermodynamic_profiles.png')
    plt.close()
    return




def eos_c(pd, s):
    T_tilde * (exp((s - sd_tilde + Rd * log(pd / p_tilde)) / cpd))
    return
if __name__ == '__main__':
    main()