import numpy as np
import math
import globe
from scipy.optimize import fsolve
from equation_of_state import Equation_of_state, Integrate_DeltaV_BL, Integrate_DeltaV_BI
from occupancy import IntegrateIsotherm


def DissociationTemperature(P):
    ndata = len(P)
    T = np.zeros(ndata)
    T0 = 280 # K
    for i in range(0,ndata):
        p = P[i]
        t_l = fsolve(Get_DissT,T0,args=(p,0))
        t_i = fsolve(Get_DissT,T0,args=(p,1))
        t = min(t_i,t_l)
        T[i] = t
    
    return T

def DissociationPressure(T):
    ndata = len(T)
    P = np.zeros(ndata)
    P0 = 200e6 # Pa
    lnp0 = math.log(P0)
    for i in range(0,ndata):
        t = T[i]
        p_l = math.exp(fsolve(Get_DissP,lnp0,args=(t,0)))
        p_i = math.exp(fsolve(Get_DissP,lnp0,args=(t,1)))
        p = max(p_i,p_l)
        P[i] = p
    
    return P

def Get_DissT(T,*Params):
    P, aq_state = Params
    lnfwh_fwb = IntegrateIsotherm(T,P)
    if aq_state == 0:
        delv_bl = Integrate_DeltaV_BL(T,P)
        mu_bl = globe.mu_bi0+globe.delh_bl0*(1/T - 1/globe.t0)-(globe.cp0_bl-globe.b_bl*globe.t0)*math.log(T/globe.t0)-0.5*(2*globe.cp0_bl*globe.t0-globe.b_bl*(globe.t0**2))*(1/T-1/globe.t0)-globe.b_bl*0.5*(T-globe.t0)+delv_bl
        return mu_bl + lnfwh_fwb

    elif aq_state == 1:
        delv_bi = Integrate_DeltaV_BI(T,P)
        mu_bi = globe.mu_bi0+globe.delh_bi0*(1/T - 1/globe.t0)-(globe.cp0_bi-globe.b_bi*globe.t0)*math.log(T/globe.t0)-0.5*(2*globe.cp0_bi*globe.t0-globe.b_bi*(globe.t0**2))*(1/T-1/globe.t0)-globe.b_bi*0.5*(T-globe.t0)+delv_bi
        return mu_bi + lnfwh_fwb

def Get_DissP(lnP,*Params):
    T, aq_state = Params
    P = math.exp(lnP)
    lnfwh_fwb = IntegrateIsotherm(T,P)

    if aq_state == 0:
        delv_bl = Integrate_DeltaV_BL(T,P)
        mu_bl = globe.mu_bi0+globe.delh_bl0*(1/T - 1/globe.t0)-(globe.cp0_bl-globe.b_bl*globe.t0)*math.log(T/globe.t0)-0.5*(2*globe.cp0_bl*globe.t0-globe.b_bl*(globe.t0**2))*(1/T-1/globe.t0)-globe.b_bl*0.5*(T-globe.t0)+delv_bl
        return mu_bl + lnfwh_fwb

    elif aq_state == 1:
        delv_bi = Integrate_DeltaV_BI(T,P)
        mu_bi = globe.mu_bi0+globe.delh_bi0*(1/T - 1/globe.t0)-(globe.cp0_bi-globe.b_bi*globe.t0)*math.log(T/globe.t0)-0.5*(2*globe.cp0_bi*globe.t0-globe.b_bi*(globe.t0**2))*(1/T-1/globe.t0)-globe.b_bi*0.5*(T-globe.t0)+delv_bi
        return mu_bi + lnfwh_fwb