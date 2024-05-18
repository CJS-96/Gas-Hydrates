#!/usr/bin/env python
# coding: utf-8

# In[12]:


import numpy as np
from scipy.optimize import fsolve
import math
import pandas as pd
from Nitrogen_EOS import N2gasprop


# In[20]:


######### Function to integrate the isotherm till vapor fugacity at current T,P ##########
def Th_Int_lnfg_N2(T,P,fug):
    #print('part3 ', T, P)
    nstep = 118
    lnfg = np.linspace(math.log(1),np.log(fug),nstep+1)
    lnfgm = np.log(fug)
    fg = np.exp(lnfg)
    fg = np.reshape(fg,-1)
    NCav = 24
    NWater = 136
    h = lnfgm/nstep
    npoints = nstep+1

    LC_Co = pd.read_csv("LCP1_N2_IW.dat",header=None);  ## Reading the fitted coefficients for LCs.
    LC_Co = LC_Co.to_numpy()

    ## LC_Co = load("LCP2_N2_IW.dat");

    beta = 1/T;
    aij = beta*np.array([[-105.864806350734,  -38.5156603319702],[-80.7268703391869, -25.629540534039]])

    LnCP = np.zeros((3))
    Lc = np.zeros((3))
    focci = np.zeros((npoints,2));
    Thetatot = np.zeros((npoints));
    Ints = np.zeros((npoints,1));
    for i in range(0,npoints):
        for j in range(0,3):
            a = LC_Co[j,0]
            b = LC_Co[j,1]*math.log(P)
            c = LC_Co[j,2]*(1/T)
            d = LC_Co[j,3]*(math.log(P)/T)
            LnCP[j] =  a + b + c + d 
            Lc[j]=math.exp(LnCP[j])/P;

        thetan = np.array([0,1])
        def Occu(theta):

            return theta - np.array([((math.exp(-(aij[0,0]*theta[0] + aij[0,1]*theta[1]))*fg[i]*Lc[0])/(1 + (math.exp(-(aij[0,0]*theta[0] + aij[0,1]*theta[1]))*fg[i]*Lc[0]))),   \
                (math.exp(-(aij[1,0]*theta[0] + aij[1,1]*theta[1]))*fg[i]*Lc[1] + (fg[i]**2)*Lc[2]*math.exp(-2*(aij[1,0]*theta[0] + aij[1,1]*theta[1])))/           \
                (1 + (math.exp(-(aij[1,0]*theta[0] + aij[1,1]*theta[1]))*fg[i]*Lc[1] + 0.5*(fg[i]**2)*Lc[2]*math.exp(-2*(aij[1,0]*theta[0] + aij[1,1]*theta[1]))))])
        
        focci[i,:] = fsolve(Occu,thetan)
        Thetatot[i] = (2/3)*(focci[i,0]) + (1/3)*(focci[i,1]);

        Ints[i] = Thetatot[i]*NCav/NWater;

    bunch1 = Ints[0]+Ints[npoints-1]
    bunch2 = 0
    bunch4 = 0
    for i in range(2,nstep-1,2):
        bunch2 = bunch2 + Ints[i]
    for i in range(1,nstep,2):
        bunch4 = bunch4 + Ints[i]

    lnfwh_fwb = -(h/3)*(bunch1+4*bunch4+2*bunch2)
    return lnfwh_fwb

######### Function for I-H-V Calculation ###############

def N2_IHV_TheoryFlFh(T):
    #print('par1 ', T, P)
    HRgas, fug, Vgas = N2gasprop(T,P)  ## Calculating vapor fugacity at T and P;

    lnfwh_fwb = Th_Int_lnfg_N2(T,P,fug) ## Calculating Thetatot and Integrating till vapor fugacity
    lnfwh_fwb = np.reshape(lnfwh_fwb,-1)

    R = 8.314
    t0 = 280
    p0 = 45.06e5
    ti0 = 260
    pi0 = 2.531e5

    fbw0 = 38.24   ## Pa
    fiw0 = 3.219   ## Pa

    Fbco = np.array([21.885e-06,38.276e-10,-93.671e-17,-54.191e-19,-68.826e3,23.389])*(1/R)
    Fico = np.array([1.9066e-05,3.03e-09,-5.6556e-16,-4.3158e-18,-69954,24.119])*(1/R)
    
    lnfwi = Fico[0]*(P/T - pi0/ti0) + Fico[1]*(P-pi0) + 0.5*Fico[2]*(P**2/T - pi0**2/ti0)   \
            + 0.5*Fico[3]*(P**2 - pi0**2) + Fico[4]*(1/T - 1/ti0) - Fico[5]*math.log(T/ti0) + math.log(fiw0)

    lnfwb = Fbco[0]*(P/T - p0/t0) + Fbco[1]*(P-p0) + 0.5*Fbco[2]*(P**2/T - p0**2/t0)    \
            + 0.5*Fbco[3]*(P**2 - p0**2) + Fbco[4]*(1/T - 1/t0) - Fbco[5]*math.log(T/t0) + math.log(fbw0)

    lnfwh = lnfwh_fwb + lnfwb

    return lnfwi - lnfwh


######### Function for L-H-V Calculation ###############

def N2_LHV_TheoryFlFh(T):
    #print('par1 ', T, P)
    HRgas, fug, Vgas = N2gasprop(T,P)  ## Calculating vapor fugacity at T and P;

    lnfwh_fwb = Th_Int_lnfg_N2(T,P,fug) ## Calculating Thetatot and Integrating till vapor fugacity
    lnfwh_fwb = np.reshape(lnfwh_fwb,-1)

    R = 8.314
    t0 = 280
    p0 = 45.06e5
    tw0 = 280
    pw0 = 20e5

    fbw0 = 38.24   ## Pa
    flw0 = 24.35   ## Pa

    Fbco = np.array([21.885e-06,38.276e-10,-93.671e-17,-54.191e-19,-68.826e3,23.389])*(1/R)
    Flco = np.array([19.074e-06,-30.365e-10,-23.990e-15,56.020e-18,-74.170e03,58.817])*(1/R)
    
    lnfwl = Flco[0]*(P/T - pw0/tw0) + Flco[1]*(P-pw0) + 0.5*Flco[2]*(P**2/T - pw0**2/tw0)   \
            + 0.5*Flco[3]*(P**2 - pw0**2) + Flco[4]*(1/T - 1/tw0) - Flco[5]*math.log(T/tw0) + math.log(flw0)

    lnfwb = Fbco[0]*(P/T - p0/t0) + Fbco[1]*(P-p0) + 0.5*Fbco[2]*(P**2/T - p0**2/t0)    \
            + 0.5*Fbco[3]*(P**2 - p0**2) + Fbco[4]*(1/T - 1/t0) - Fbco[5]*math.log(T/t0) + math.log(fbw0)

    lnfwh = lnfwh_fwb + lnfwb

    return lnfwl - lnfwh

######## Main program to calculate Phase Equilibria ##########
#Pressure = 1e5*np.array([51.2085596507515,60.535052084374,71.5601562678224,80.6444333101443,100,121.072352955239,143.122970888239]) ### N2_Ice Pressure points 
Pressure = 1e6*np.array([20,30,40,50,60,70,80]) ### N2_W Pressure points
#Pressure = np.exp(np.linspace(math.log(5.121e6),math.log(80e6),100)); ### For calculating Phase-Equilibria for the total Ice-Liquid range

#Tnot = np.linspace(200,270,100) ### Initial Temperature guesses
Tnot=280  ## Can give a single guess for multiple the eqm points
npE = np.size(Pressure)
TempTh = np.zeros((npE))

for i in range(0,npE):
    print(i)
    P = Pressure[i]
    #TempTh[i] = fsolve(N2_IHV_TheoryFlFh,Tnot)
    TempTh[i] = fsolve(N2_LHV_TheoryFlFh,Tnot)

PhasePts = pd.DataFrame({'Pressure (Pa)' : Pressure , 'Temperature (K)' : TempTh}, columns=['Pressure (Pa)', 'Temperature (K)'])
#PhasePts.to_csv('TempThI_N2.dat', index=False, sep='\t')
PhasePts.to_csv('TempThW_N2.dat', index=False, sep='\t')


# In[21]:


TempTh


# In[ ]:




