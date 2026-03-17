#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from scipy.optimize import fsolve
import math
import pandas as pd


# In[3]:

## This code is written based on the EOS proposed by
## Kolafa and Nezbeda, Fluid Phase Equilib., 100 (1994), 1-34
######### EOS function for Argon ############################
def Argasprop(T,P):
    Temp=np.atleast_1d(T)
    eps = 119.8  ## (K)
    sigma = 3.405 ## (Angstorm)
    Kb=1.381e7 ## (JA^3/K)
    R=8.314  ## (J/mol/K)
    Pfac = sigma**3/(eps*Kb) ## (Pa^-1)
    Tr = Temp/eps                ## Reduced Temperature
    Pr = np.atleast_1d(P)*Pfac   ## Reduced Pressure
    nT = np.size(Tr)
    nP = np.size(Pr)
    
    fug = np.zeros((nT,nP))        ## Fugacity (Pa)
    vfugr = np.zeros((nT,nP))      ## reduced fugacity
    densities = np.zeros((nT,nP))  ## Density 
    Vgas = np.zeros((nT,nP))       ## Volume (m^3/mol)
    Hreal = np.zeros((nT,nP))      ## Residual Enthalpy (J/mol)
    for trloop in range(0,nT):           ## Loop over all input temperature values
        for prloop in range(0,nP):       ## Loop over all input pressure values 
            def vapor_compress(rho,*TP):    ## Function to get the density at T and P
                Tr, Pr = TP
                dhBH = 0.011117524*(Tr**(-2*0.5)) - 0.076383859*(Tr**(-1*0.5)) + 1.080142248*(Tr**(0*0.5)) + 0.000693129*(Tr**(1*0.5)) - 0.063920968*math.log(Tr)
                deltaB2 = -0.58544978*Tr**(-7*0.5) + 0.43102052*Tr**(-6*0.5) + 0.87361369*Tr**(-5*0.5) - 4.13749995*Tr**(-4*0.5) + 2.90616279*Tr**(-3*0.5)  \
                          -7.02181962*Tr**(-2*0.5) + 0.02459877*Tr**(0);

                term4 = 2*2.0154679*Tr**(-1)*rho**2 - 3*28.17881636*Tr**(-1)*rho**3 + 4*28.28313847*Tr**(-1)*rho**4 - 5*10.42402873*Tr**(-1)*rho**5  \
                        - 2*19.58371655*Tr**(-1.5)*rho**2 + 3*75.62340289*Tr**(-1.5)*rho**3 - 4*120.70586598*Tr**(-1.5)*rho**4 + 5*93.93740328*Tr**(-1.5)*rho**5 - 6*27.37737354*Tr**(-1.5)*rho**6  \
                        + 2*29.34470520*Tr**(-2)*rho**2 - 3*112.35356937*Tr**(-2)*rho**3 + 4*170.64908980*Tr**(-2)*rho**4 - 5*123.06669187*Tr**(-2)*rho**5 + 6*34.42288969*Tr**(-2)*rho**6         \
                        - 2*13.37031968*Tr**(-3)*rho**2 + 3*65.38059570*Tr**(-3)*rho**3 - 4*115.09233113*Tr**(-3)*rho**4 + 5*88.91973082*Tr**(-3)*rho**5 - 6*25.62099890*Tr**(-3)*rho**6
                y = 1.92907278
                eta = math.pi*rho*(dhBH**3)/6
                Zhs = (1+eta+eta**2-(2/3)*(eta**3)*(1+eta))/((1-eta)**3)
                return Zhs + rho*(1-2*y*(rho**2))*deltaB2 + term4 - Pr/(rho*Tr)

            rho0 = 0.005                                 ## initial guess for density
            TP = (Tr[trloop],Pr[prloop])                 ## T and P values to supply through fsolve
            den = fsolve(vapor_compress,rho0,args=TP)    ## Using fsolve to solve for density
            densities[trloop,prloop] = den               ## Storing the density at each (T,P) point
            
            ### Calculating the rest of the properties using the density calculated above
            
            y = 1.92907278
            Cdi = np.array([0.011117524,-0.076383859,1.080142248,0.00063920968])
            di = np.array([-2,-1,0,1])
            Cdln = -0.063920968
            poly1=0
            poly2=0
            for sloop in range (0,4):
                poly1 = poly1 + Cdi[sloop]*(Tr[trloop]**(di[sloop]*0.5))
                poly2 = poly2 + (di[sloop]/2)*Cdi[sloop]*(Tr[trloop]**(di[sloop]*0.5-1))

            dhBH = poly1 + Cdln*math.log(Tr[trloop])
            ddhBH = poly2 + Cdln/Tr[trloop]

            eta = (math.pi/6)*den*(dhBH**3)                               ## Hard-Sphere packing fraction
            Zhs = (1+eta+eta**2-(2/3)*(eta**3)*(1+eta))/((1-eta)**3)      ## Hard-Sphere compressibility
            I4 = Zhs - 1
            I1 = - (25/6 - (2*eta + 5*math.log(abs(eta-1)))/3 + 10/(3*(eta-1)) - 5/(6*((eta-1)**2)))

            CBi = np.array([-0.58544978,0.43102052,0.87361369,-4.13749995,2.90616279,-7.02181962,0.02459877])
            bi = np.array([-7,-6,-5,-4,-3,-2,0])
            deltaB2=0
            ddeltaB2=0
            for tloop in range(0,7):
                deltaB2 = deltaB2 + CBi[tloop]*Tr[trloop]**(bi[tloop]*0.5)
                ddeltaB2 = ddeltaB2 + (bi[tloop]/2)*CBi[tloop]*(Tr[trloop]**(bi[tloop]*0.5-1))

            I2 = 2*deltaB2*den*math.exp(-y*(den**2))*(1-y*(den**2))
            Z2 = deltaB2*den*math.exp(-y*den**2)*(1-2*y*den**2)

            Cij = np.array([[2.01546797,-28.17881636,28.28313847,-10.42402873, 0.0000000],
                    [-19.58371655,75.62340289,-120.70586598,93.92740328,-27.37737354],
                    [29.34470520,-112.35356937,170.64908980,-123.06669187,34.42288969],
                    [-13.37031968,65.38059570,-115.09233113,88.91973082,-25.62099890]])
            i = np.array([0,-1,-2,-4])
            j = np.array([2,3,4,5,6])
            I3 = 0
            Z3 = 0
            dI3 = 0
            for uloop in range(0,4):
                for vloop in range(0,5):
                    I3 = I3 + (j[vloop]+1)*Cij[uloop,vloop]*Tr[trloop]**(i[uloop]/2-1)*den**(j[vloop])
                    dI3 = dI3 + (i[uloop]*0.5 - 1)*Cij[uloop,vloop]*(Tr[trloop]**(i[uloop]/2-2))*den**(j[vloop])
                    Z3 = Z3 + (j[vloop])*Cij[uloop,vloop]*Tr[trloop]**(i[uloop]/2-1)*den**(j[vloop])

            Z = Z2 + Z3 + Zhs                       ## Compressibility
            Prcal = Z*den*Tr[trloop]                ## Pressure calculated back from Temperature and density
            I5 = -math.log(Z)
            lnphi = I1 + I2 + I3 + I4 + I5
            Press=Pr[prloop]/Pfac                                ## Pressure (Pa) from reduced pressure
            fug[trloop,prloop] = Press*math.exp(lnphi)           
            vfugr[trloop,prloop] = Pr[prloop]*math.exp(lnphi)
            Vgas[trloop,prloop] = Z*R*Temp[trloop]/Press
            IH1 = 6*eta - (11/(eta-1)) - (1/(eta-1)**2) - 5/(3*(eta-1)**3) + 16*math.log(abs(eta-1)) - 35/3
            H1 = 3*ddhBH*IH1/dhBH
            H2 = den*math.exp(-y*(den**2))*ddeltaB2
            H3 = dI3
            Hr = -Tr[trloop]*(H1+H2+H3)+Z-1
            Hreal[trloop,prloop] = R*Temp[trloop]*Hr
    
    ### Transposing the quantities to make the arrays compatible with Phase Equilibria codes.
    
    Hreal=np.transpose(Hreal)                           
    fug=np.transpose(fug)
    Vgas=np.transpose(Vgas)
    return Hreal,fug,Vgas

