import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from enum import Enum

class constants(float, Enum):
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
    return T * (constants.P0/p) ** ((constants.R_air / constants.cp_air) * (1 - 0.28) * q)

def equiv_potential_temp(p, T, q):
    # Empirical formula for equivalent potential temperature (Bolton 1980)
    e = (p * q)/(0.622 + q)
    T_L = 2840/(3.5 * np.log(T) - np.log(e) - 4.805)

    return potential_temp(p, T, q) * np.exp((1000 * q) * (1 + .81*q) * (3.376/T_L - 0.00254))

def eyewall_delta_T(theta_ES, p, T, T_env, q_star):
    return (theta_ES * (p/1000) ** (constants.R_air/constants.cp_air)) / np.exp((1000 * q_star) * (1 + .81*q_star) * (3.376/T - 0.00254)) - T_env

def eye_delta_T(theta_ES, p, T, T_env, q, q_star):
    return (theta_ES * (p/1000) ** (constants.R_air/constants.cp_air)) / np.exp((1000 * q_star) * (1 + .81*q) * (3.376/T - 0.00254)) - T_env

def T_v(T, q):
    return T * (1 + 0.61*q)

def hydrostatic_pressure_change(P_s, Tv_s, p, delta_Tv):
    return P_s/Tv_s * np.trapz(delta_Tv, 1/p)

def eye_RH(P, P_s_E):
    RH = np.zeros(len(P))

    for i, p in enumerate(P):
        if p <= 200:
            RH[i] = 0
            continue
        
        C = -10
        X = C + ((P_s_E - 800) * (p - 200))/1000

        if p > 200 and p <= 450:
            RH[i] = max(X, 0)
            continue

        if p > 450 and p <= 700:
            RH[i] = min(X, 80)
            continue

        if p > 700:
            RH[i] = 0.075*p + 27.5

    return RH

def update_Ps(P_s, T_s, RH, p, T, T_env):
    q_s = mix_ratio(P_s, T_s, RH)
    q_star_s = saturation_mix_ratio(P_s, T_s)
    theta_ES = equiv_potential_temp(P_s, T_s, q_s)

    q = mix_ratio(p, T, RH)
    q_star = saturation_mix_ratio(p, T)

    delta_T = eyewall_delta_T(theta_ES, p, T, T_env, q_star)
    print(delta_T)
    delta_Tv = T_v(delta_T, q)
    delta_Ps = hydrostatic_pressure_change(P_s, T_v(T_s, q_s), p, delta_Tv)
    P_s += delta_Ps

    return P_s


if __name__ == "__main__":
    P_env = 1007
    SST = 27 + 273.15
    RH = 80
    df = pd.read_csv('data/holland_willis_island_january.csv', header=0)

    T_env = np.flip(df['Temperature'].to_numpy() + 273.15)
    p = np.flip(df['Pressure'].to_numpy())
    T = np.copy(T_env)
    P_s = P_env
    T_s = SST

    q_s = mix_ratio(P_env, T_s, RH)
    q_star_s = saturation_mix_ratio(P_env, T_s)
    theta_ES_start = equiv_potential_temp(P_env, T_s, q_s)
    print(theta_ES_start)

    # plt.ylim(1020, 0)
    # plt.scatter(df['Temperature'], df['Pressure'])    
    # plt.show()

    P_s_change = 1.1
    while abs(P_s_change) > 1:
        P_s_new = update_Ps(P_s, T_s, RH, p, T, T_env)
        P_s_change = P_s_new - P_s
        P_s = P_s_new
        print(f"Surface pressure: {P_s}, Pressure Change {P_s_change}")

    eye_RH = eye_RH(p, P_s)
    print(eye_RH)

