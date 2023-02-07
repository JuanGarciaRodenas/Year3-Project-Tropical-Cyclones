import numpy as np

P0 = 1000
R_air = 287
cp_air = 1005
cp_water = 4184

def saturation_mix_ratio(p, T):
    return 3.802/p * np.exp((17.67 * (T - 273.15)) /  (T - 29.65))

def mix_ratio(p, T, RH):
    return RH/100 * saturation_mix_ratio(p, T)

def potential_temp(p, T, q):
    # Potential temperature of moist air neglecting cp variation (Bolton 1980)
    return (T - 273.15) * (P0/p) ** ((R_air / cp_air) * (1 - 0.28) * q)

def equiv_potential_temp(p, T, q):
    # Empirical formula for equivalent potential temperature (Bolton 1980)
    e = (p * q)/(0.622 + q)
    T_L = 2840/(3.5 * np.log(T) - np.log(e) - 4.805)

    return potential_temp(p, T, q) * np.exp((1000 * q) * (1 + .81*q) * (3.376/T_L - 0.00254))


def eyewall_delta_T(theta_ES, p, T, T_env, q_star):
    return (theta_ES * (p/1000) ** (R_air/cp_air)) / np.exp((1000 * q_star) * (1 + .81*q_star) * (3.376/T - 0.00254)) - T_env

def eye_delta_T(theta_ES, p, T, T_env, q, q_star):
    return (theta_ES * (p/1000) ** (R_air/cp_air)) / np.exp((1000 * q_star) * (1 + .81*q) * (3.376/T - 0.00254)) - T_env


if __name__ == "__main__":
    P_env_willis = 1007
    P_env_barbados = 1010

    RH = 80
    SST = 27 + 273.15
    q = mix_ratio(P_env_willis, SST, RH)
    q_star = saturation_mix_ratio(P_env_willis, SST)
    theta_ES_start = equiv_potential_temp(P_env_willis, SST, q)

    deltaT = eyewall_delta_T(theta_ES_start, P_env_willis, SST, SST, q_star)a
    print(deltaT)
    