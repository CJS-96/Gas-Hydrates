import numpy as np
import math
import globe
from sympy import Symbol, solve
from utils import Simpsons_OneThird_Rule

def Equation_of_state(T,P):
    fugacity, Volume = SRK_eos(T,P)

    return fugacity, Volume

def SRK_eos(T,P):
    T = float(T)
    P = float(P)
    Tr = T/globe.Tc
    Pr = P/globe.Pc
    psi=0.42747
    omega=0.08664
    alpha = (1+(0.480+1.574*globe.w-0.176*(globe.w**2))*(1-Tr**0.5))**2
    A = psi*alpha*Pr/(Tr**2)
    B = omega*Pr/Tr
    q = psi*alpha/omega/Tr
    Z = Symbol('Z')
    Cubic_eos = Z**3-Z**2+Z*(A-B-B**2)-A*B
    Zs = solve(Cubic_eos)
    Z_com = [complex(z.evalf()) for z in Zs]
    Zs_real = [z.real for z in Z_com if abs(z.imag) < 1e-16]
    Z = max(Zs_real)                 # Max or Min will depend on whether we want a root for liquid or vapor phase.
    fug = P*math.exp(Z-1-math.log(Z-B)-(A/B)*math.log((Z+B)/Z))
    Vgas = Z*globe.R*T/P

    return fug, Vgas

def Integrate_DeltaV_BL(T,P):
    nsteps = 100
    npoints = nsteps+1
    Avag_Number = 6.022e23

    al1 = 10.9241; al2 = 2.5e-4; al3 = 3.532e-4; al4 = 1.559e-7
    ab1 = 17.13; ab2 = 2.429e-4; ab3 = 2.013e-6; ab4 = 1.009e-9; ab5 = 8.006e-9; ab6 = 5.448e-12

    Prange = np.linspace(globe.p0*1e-6,P*1e-6,npoints)
    step = (P-globe.p0)/nsteps
    dv = np.zeros(npoints)

    for i in range(0,npoints):
        v_l = math.exp(-al1+al2*(T-273.15)-al3*(Prange[i]-0.101325)+al4*((Prange[i]-0.101325)**2))
        v_b = ((ab1+ab2*T+ab3*(T**2)+ab4*(T**3))**3)*1e-30*Avag_Number/globe.Nwaters-ab5*Prange[i]+ab6*(Prange[i]**2)
        dv[i] = v_b-v_l

    delv_bl = Simpsons_OneThird_Rule(dv,step)/globe.R/T

    return delv_bl

def Integrate_DeltaV_BI(T,P):
    nsteps = 100
    npoints = nsteps+1
    Avag_Number = 6.022e23

    ai1 = 1.912e-5; ai2 = 8.387e-10; ai3 = 4.016e-12;
    ab1 = 17.13; ab2 = 2.429e-4; ab3 = 2.013e-6; ab4 = 1.009e-9; ab5 = 8.006e-9; ab6 = 5.448e-12

    Prange = np.linspace(globe.p0*1e-6,P*1e-6,npoints)
    step = (P-globe.p0)/nsteps
    dv = np.zeros(npoints)
    v_i = ai1 + ai2*T + ai3*(T**2);
    for i in range(0,npoints):
        v_b = ((ab1+ab2*T+ab3*(T**2)+ab4*(T**3))**3)*1e-30*Avag_Number/globe.Nwaters-ab5*Prange[i]+ab6*(Prange[i]**2)
        dv[i] = v_b-v_i

    delv_bi = Simpsons_OneThird_Rule(dv,step)/globe.R/T

    return delv_bi