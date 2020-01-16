import numpy as np

# ----------------------- PARAMETERS -----------------------------

#include "parameters.h"
g = 9.80665
Rd = 287.1
Rv = 461.5
eps_v = 0.62210184182
eps_vi = 1.60745384883
cpd = 1004.0
cpv = 1859.0
kappa = 0.285956175299

T_tilde = 298.15
p_tilde = 100000.0
sd_tilde = 6864.8
sv_tilde = 10513.6


# ----------------------- DIVERSE -----------------------------

def exner_c(p0):
    return np.power((p0/p_tilde), kappa)

def theta_c(p0, T):
    # Dry potential temperature
    return T / exner_c(p0)

# inline double thetali_c(const double p0, const double T, const double qt, const double ql, const double qi, const double L){
#     // Liquid ice potential temperature consistent with Tripoli and Cotton (1981)
#     return theta_c(p0, T) * exp(-L*(ql/(1.0 - qt) + qi/(1.0 - qt))/(T*cpd));
# }

def pd_c(p0, qt, qv):
    return p0*(1.0-qt)/(1.0 - qt + eps_vi * qv)

def pv_c(p0, qt, qv):
    return p0 * eps_vi * qv /(1.0 - qt + eps_vi * qv)

# inline double density_temperature_c(const double T, const double qt, const double qv){
#     return T * (1.0 - qt + eps_vi * qv);
# }

# inline double theta_rho_c(const double p0, const double T, const double qt, const double qv){
#     return density_temperature_c(T,qt,qv)/exner_c(p0);
# }

def cpm_c(qt):
    return (1.0-qt) * cpd + qt * cpv

def thetas_c(s, qt):
    return T_tilde*np.exp((s-(1.0-qt)*sd_tilde - qt*sv_tilde)/cpm_c(qt))

# inline double thetas_t_c(const double p0, const double T, const double qt, const double qv, const double qc, const double L){
#     const double qd = 1.0 - qt;
#     const double pd = pd_c(p0,qt,qt-qc);
#     const double pv = pv_c(p0,qt,qt-qc);
#     const double cpm = cpm_c(qt);
#     return T * pow(p_tilde/pd,qd * Rd/cpm)*pow(p_tilde/pv,qt*Rv/cpm)*exp(-L * qc/(cpm*T));
# }

def entropy_from_thetas_c(thetas, qt):
    return cpm_c(qt) * np.log(thetas/T_tilde) + (1.0 - qt)*sd_tilde + qt * sv_tilde

def buoyancy_c(alpha0, alpha):
    return g * (alpha - alpha0)/alpha0

def qv_star_c(p0, qt, pv):
    return eps_v * (1.0 - qt) * pv / (p0 - pv)

def alpha_c(p0, T, qt, qv):
    return (Rd * T)/p0 * (1.0 - qt + eps_vi * qv)



# ----------------------- ENTROPIES -----------------------------

def sd_c(pd, T):
    return sd_tilde + cpd*np.log(T/T_tilde) -Rd*np.log(pd/p_tilde)


def sv_c(pv, T):
    return sv_tilde + cpv*np.log(T/T_tilde) - Rv * np.log(pv/p_tilde)

def sc_c(L, T):
    return -L/T

def s_tendency_c(p0, qt, qv, T, qt_tendency, T_tendency):
    pv = pv_c(p0, qt, qv)
    pd = pd_c(p0, qt, qv)
    return cpm_c(qt) * T_tendency / T +  (sv_c(pv,T) - sd_c(pd,T)) * qt_tendency
