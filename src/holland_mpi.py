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

if __name__ == "__main__":
    P_env = 1007
    SST = 27 + 273.15
    RH = 80
    df = pd.read_csv('data/holland_willis_island_january.csv', header=0)

    T_env = np.flip(df['Temperature'].to_numpy() + 273.15)
    T = np.copy(T_env)
    p = np.flip(df['Pressure'].to_numpy())
    P_s = p[0]
    T_s = T[0]

    q = mix_ratio(P_env, T_s, RH)
    q_star = saturation_mix_ratio(P_env, T_s)
    theta_ES_start = equiv_potential_temp(P_env, T_s, q)
    print(theta_ES_start)

    delta_Tv = eyewall_delta_T(theta_ES_start, p, T, T_env, q_star)
    print(delta_Tv)

    delta_Ps = hydrostatic_pressure_change(P_s, T_v(T_s, q), p, delta_Tv)
    print(delta_Ps)

    P_s += delta_Ps
    while abs(delta_Ps) > 1:
        theta_ES = equiv_potential_temp(P_s, T_env, q)
        delta_Tv = eye_delta_T(theta_ES, p, T, T_env, q, q_star)
        delta_Ps = hydrostatic_pressure_change(P_s, T_v(T_s, q), p, delta_Tv)

        P_s += delta_Ps 

    print(P_s)

    # plt.ylim(1020, 0)
    # plt.scatter(df['Temperature'], df['Pressure'])    
    # plt.show()
