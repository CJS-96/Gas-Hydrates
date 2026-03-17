#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
from scipy.optimize import fsolve
import math
import pandas as pd


# In[ ]:


############### Virial EOS function for Nitrogen ####################
def N2gasprop(T,P):
    
    ## Coefficients of quadratic fit for the virial coefficients as a function of temperature
    IndB = np.array([-8.1451e-10, 5.6637e-07,-9.8201e-05])    
    IndC = np.array([-4.5881e-14,3.8305e-11,-8.3299e-09])
    IndD = np.array([1.3720e-18,-1.2292e-15,3.9277e-13])
    
    R = 8.314              ##(J/mol/K)
    T = np.atleast_1d(T)   ## (K)
    nT = np.size(T)
    P = np.atleast_1d(P)   ## (Pa)
    nP = np.size(P)
    dens = np.zeros((nT,nP))    ## density (mol/m^3)
    fug = np.zeros((nT,nP))     ## Fugacity (Pa)
    Vgas = np.zeros((nT,nP))    ## Volume (m^3/mol)
    HRgas = np.zeros((nT,nP))   ## Residual Enthalpy (J/mol)

    for i in range(0,nT):       ## Loop over all the input Temperatures
        B = IndB[0]*(T[i]**2) + IndB[1]*T[i] + IndB[2]
        C = IndC[0]*(T[i]**2) + IndC[1]*T[i] + IndC[2]
        D = IndD[0]*(T[i]**2) + IndD[1]*T[i] + IndD[2]
        Bdash = 2*IndB[0]*(T[i]**2) + IndB[1]*T[i]
        Cdash = 0.5*(2*IndC[0]*(T[i]**2) + IndC[1]*T[i])
        Ddash = (1/3)*(2*IndD[0]*(T[i]**2) + IndD[1]*T[i])
        for j in range(0,nP):                                ## Loop over all the input Pressures
            def virfun(rho):                                 ## Function to get the density at T and P
                return D*(rho**4) + C*(rho**3) + B*(rho**2) + rho - P[j]/(R*T[i])
            rho0 = 10000                                     ## Inital guess for the density
            dens[i,j] = fsolve(virfun,rho0)                  ## Using fsolve to solve for the density
            currden = dens[i,j]                              ## Storing the densities for all the (T,P) points
            Z = 1 + B*currden + C*(currden**2) + D*(currden**3)   ## Compressibility
            Vgas[i,j] = Z*R*T[i]/P[j]
            lnphi = 2*B*currden + 1.5*C*(currden**2) +(4/3)*D*(currden**3) - math.log(1+B*currden+C*(currden**2)+D*(currden**3))
            fug[i,j] = P[j]*math.exp(lnphi)
            HRgas[i,j] = R*T[i]*(-Bdash*currden - Cdash*(currden**2) -Ddash*(currden**3)    \
                              + B*currden + C*(currden**2) + D*(currden**3))
    
    ### Transposing the quantities to make the arrays compatible with Phase Equilibria codes.
    
    HRgas = np.transpose(HRgas)
    fug = np.transpose(fug)
    Vgas = np.transpose(Vgas)
    return HRgas,fug, Vgas

