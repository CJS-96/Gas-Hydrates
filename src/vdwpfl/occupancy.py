import numpy as np
import math
from . import globe
from scipy.optimize import fsolve
from .equation_of_state import Equation_of_state
from .utils import Simpsons_OneThird_Rule

def CalculateLangmuirConstants(T,P):
    T=float(T)
    P=float(P)
    LCform = np.array([1.0, math.log(P), 1/T, math.log(P)/T, P])
    ln_LCP = np.matmul(globe.LC_coeffs,LCform)
    LC = np.exp(ln_LCP)/P

    return LC

def CalculateTheta(theta,*args):
    fug, aij, LC = args
    return theta - np.array([(fug*LC[0]*math.exp(-(aij[0,0]*theta[0] + aij[0,1]*theta[1])))/(1 + fug*LC[0]*math.exp(-(aij[0,0]*theta[0] + aij[0,1]*theta[1]))), \
                   (fug*LC[1]*math.exp(-(aij[1,0]*theta[0] + aij[1,1]*theta[1])) + (fug**2)*LC[2]*math.exp(-2*(aij[1,0]*theta[0] + aij[1,1]*theta[1])))/  \
                   (1 + (math.exp(-(aij[1,0]*theta[0] + aij[1,1]*theta[1]))*fug*LC[1] + 0.5*(fug**2)*LC[2]*math.exp(-2*(aij[1,0]*theta[0] + aij[1,1]*theta[1]))))])

def CalculateIsotherm(fug_max,aij,LC):
    nsteps = 118
    npoints = nsteps+1
    ln_fug_max = np.log(fug_max)
    ln_fug = np.linspace(np.log(1),np.log(fug_max),npoints)
    fug = np.exp(ln_fug)
    frac_occ = np.zeros((npoints,2))
    theta0 = np.array([0,1])

    for i in range(0,npoints):
        frac_occ[i,:] = fsolve(CalculateTheta,theta0,args=(fug[i],aij,LC))

    return frac_occ, ln_fug

def Occupancies(T,P):
    beta = 1/T
    aij = beta*globe.aij
    fugacity, Volume = Equation_of_state(T,P)
    theta0 = np.array([0,1])

    LC = CalculateLangmuirConstants(T,P)

    if globe.Property == 0:
        frac_occ = fsolve(CalculateTheta,theta0,args=(fugacity,aij,LC))
        return frac_occ
    elif globe.Property == 1:
        frac_occ, ln_fug = CalculateIsotherm(fugacity,aij,LC)
        isotherm = np.array([frac_occ[:,0],frac_occ[:,1],ln_fug])
        return isotherm
    else:
        frac_occ, ln_fug = CalculateIsotherm(fugacity,aij,LC)
        return frac_occ, ln_fug

def IntegrateIsotherm(T,P):
    nsteps = 118
    frac_occ, ln_fug = Occupancies(T,P)
    step = ln_fug[-1]-ln_fug[-2]
    frac_occ_tot = frac_occ[:,0]*globe.R_small_cav+frac_occ[:,1]*globe.R_large_cav
    Integrands = frac_occ_tot*globe.R_cav_water
    lnfwh_fwb = -Simpsons_OneThird_Rule(Integrands,step)
    return lnfwh_fwb
