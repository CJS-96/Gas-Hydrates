#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
from scipy.optimize import fsolve
import math
import pandas as pd
from Nitrogen_EOS import N2gasprop


# In[6]:


######## Function to calculate I-H-V data ###############
def SimFlFh_Ice(T,*ZHPdata):
    Zterm, HconstP, P = ZHPdata
    
    R = 8.314
    t0 = 280
    
    ti0 = 260 
    pi0 = 2.531e5    
    Fico = np.array([1.9066e-05,3.03e-09,-5.6556e-16,-4.3158e-18,-69954,24.119])*(1/R)
    fiw0 = 3.219

    fhw0 = 33.5367                        ### Fugacity of water in a filled N2-hydrate
    Temperature = np.linspace(200,300,21)  ### Enthalpy temperature range for Nitrogen_Ice

    lnfhw0 = math.log(fhw0)
    OnebyT = 1/Temperature
    poly = np.polyfit(OnebyT,HconstP,3)    ### Fitting H integrand with 1/T as variable for N2 and Argon.
    poly = poly/R


    ###### Evaluating the integral of Hr term at instantaneous temperature T ######
    Hterm = poly[0]*(1/4)*((1/T)**4 - (1/t0)**4) + poly[1]*(1/3)*((1/T)**3 - (1/t0)**3)    \
            + poly[2]*(1/2)*((1/T)**2 - (1/t0)**2) + poly[3]*((1/T) - (1/t0))       ### N2 and Argon

    lnfwhS = Hterm + Zterm + lnfhw0  ## log(fug of water in hydrate phase at T and Pressure(a))

    ### EOS of fugacity of water in ice
    lnfwi = Fico[0]*(P/T - pi0/ti0) + Fico[1]*(P-pi0) + 0.5*Fico[2]*(P**2/T - pi0**2/ti0)   \
            + 0.5*Fico[3]*(P**2 - pi0**2) + Fico[4]*(1/T - 1/ti0) - Fico[5]*math.log(T/ti0) + math.log(fiw0)

    return lnfwhS - lnfwi   ## Function equating fug of water in ice and hyd phases, to be minimized.

######## Function to calculate L-H-V data ###############
def SimFlFh_Water(T,*ZHPdata):
    Zterm, HconstP, P = ZHPdata
    
    R = 8.314
    t0 = 280
    
    tw0 = 280;          
    pw0 = 20e5;
    Flco = np.array([19.074e-06,-30.365e-10,-23.990e-15,56.020e-18,-74.170e03,58.817])*(1/R);
    flw0 = 24.35;


    fhw0 = 33.5367          ### Fugacity of water in a real Ar-hydrate
    Temperature = np.linspace(270,300,7)  ### Enthalpy temperature range for Argon

    lnfhw0 = math.log(fhw0)
    OnebyT = 1/Temperature
    poly = np.polyfit(OnebyT,HconstP,3)    ### Fitting H integrand with 1/T as variable for N2 and Argon.
    poly = poly/R


    ###### Evaluating the integral of Hr term at instantaneous temperature T ######
    Hterm = poly[0]*(1/4)*((1/T)**4 - (1/t0)**4) + poly[1]*(1/3)*((1/T)**3 - (1/t0)**3)    \
            + poly[2]*(1/2)*((1/T)**2 - (1/t0)**2) + poly[3]*((1/T) - (1/t0))       ### N2 and Argon

    lnfwhS = Hterm + Zterm + lnfhw0  ## log(fug of water in hydrate phase at T and Pressure(a))

    ### EOS of fugacity of water in ice
    lnfwl = Flco[0]*(P/T - pw0/tw0) + Flco[1]*(P-pw0) + 0.5*Flco[2]*(P**2/T - pw0**2/tw0)   \
            + 0.5*Flco[3]*(P**2 - pw0**2) + Flco[4]*(1/T - 1/tw0) - Flco[5]*math.log(T/tw0) + math.log(flw0)

    return lnfwhS - lnfwl   ## Function equating fug of water in ice and hyd phases, to be minimized.


# In[9]:


#### Reading Data for I-H-V
def Read_IHV_data():

    ETPHgasdata = pd.read_csv('ETPH_gas_N2_Ice.dat',sep='\s+')       
    ETPHgas = ETPHgasdata.to_numpy()
    ETPHhyddata = pd.read_csv('ETPH_hyd_N2_Ice.dat',sep='\s+')
    ETPHhyd = ETPHhyddata.to_numpy()
    ETPVhyddata = pd.read_csv('EPZpoints_all_N2_Ice.dat', sep='\s+')
    ETPVhyd = ETPVhyddata.to_numpy()
    
    return ETPHgas, ETPHhyd, ETPVhyd

#### Reading Data for L-H-V
def Read_LHV_data():

    ETPHgasdata = pd.read_csv('ETPH_gas_N2_W.dat',sep='\s+')       
    ETPHgas = ETPHgasdata.to_numpy()
    ETPHhyddata = pd.read_csv('ETPH_hyd_N2_W.dat',sep='\s+')
    ETPHhyd = ETPHhyddata.to_numpy()
    ETPVhyddata = pd.read_csv('EPZpoints_all_N2_W.dat', sep='\s+')
    ETPVhyd = ETPVhyddata.to_numpy()

    return ETPHgas, ETPHhyd, ETPVhyd

##################################################################

#ETPHgas, ETPHhyd, ETPVhyd = Read_IHV_data()   ### For Ice
ETPHgas, ETPHhyd, ETPVhyd = Read_LHV_data()   ### For Liquid water    

Nwater=136
NA = 6.023E+023
R = 8.314

HRgas = ETPHgas[:,2]
Thyd = ETPHhyd[:,0]
Phyd = ETPHhyd[:,1]
Numhyd = ETPHhyd[:,2]
Ntothyd = Numhyd+Nwater   ### T and P and N details from hyd phase
IEhyd = ETPHhyd[:,3]
Volhyd = (ETPHhyd[:,4]*NA)/(Ntothyd)   ## Energies and Vol from hydrate simulations

xHgH = Numhyd/(Numhyd+Nwater)
xHwH = 1-xHgH                           ## hyd phase mol fractions
HRhyd = IEhyd + Phyd*Volhyd - R*Thyd   ## Hr of hyd phase = U + PV -RT
Hcomb = HRhyd/xHwH - HRgas*(xHgH/xHwH) ## H integrand term.
nd = np.size(Hcomb)
npH = 7                                  ## Number of Pressure points for each temperature in EPTempSim.
HcombMat = np.zeros((nd//npH,npH))

for a in range(0,(nd//npH)):
    A = Hcomb[npH*(a):npH*(a+1)]         ## Putting it into  T x P matrix.
    HcombMat[a,:] = np.transpose(A)

T0 = 220              ## Starting guess for finding T
t0 = 280              ## Reference temperature of Hydrate

#Pressure = 1e5*np.array([51.2085596507515,60.535052084374,71.5601562678224,80.6444333101443,100,121.072352955239,143.122970888239])   ##  Ice -- Pressure (N2) at which the eqm T is calculated
Pressure = 1e6*np.array([20,30,40,50,60,70,80]) ##  Water -- Pressure (N2) at which the eqm T is calculated
npE = np.size(Pressure)
TempSim = np.zeros((npE))
nz=11
ZcombMat=np.zeros((nz,npE))

for a in range(0,npE):
#### Integrating Z term ####
    xHgV = ETPVhyd[nz*a:nz*(a+1),2]      ## Taking the values needed to evaluate
    xHwV = 1 - xHgV                      ## Z integrand at cont T=t0, upto different
    P = ETPVhyd[nz*a:nz*(a+1),1]         ## P0 to Pressure(a).
    Vhyd = ETPVhyd[nz*a:nz*(a+1),4]
    HRgas, fug, Vgas = N2gasprop(t0,P)   ## For fugacity of Nitrogen
    Vgas = np.reshape(Vgas,-1)
    Vcomb = Vhyd/xHwV - Vgas*(xHgV/xHwV)
    Zcomb = Vcomb*P/(R*t0)               ## Evaluated Z int at Ps leading to Pressure(a)
    ZcombMat[:,a]=Zcomb
    bunch1 = Zcomb[0] + Zcomb[nz-1]      ## Integrating Z from P=P0 to P=Pressure(a) at T=t0
    bunch2 = 0
    bunch4 = 0
    for j in range(2,nz-2,2):            ## 2:2:nz-1
        bunch2 = bunch2 + Zcomb[j]
      
    for j in range(1,nz-1,2):            ## 3:2:nz-2
        bunch4 = bunch4 + Zcomb[j]

    Zh = (math.log(P[nz-1]) - math.log(P[0]))/(nz-1)
    Zterm = (Zh/3)*(bunch1 + 4*bunch4 + 2*bunch2)

    HconstP = HcombMat[:,a]

    ZHPdata = (Zterm,HconstP,Pressure[a])
    #TempSim[a] = fsolve(SimFlFh_Ice,T0,args=ZHPdata)
    TempSim[a] = fsolve(SimFlFh_Water,T0,args=ZHPdata)
PhasePts = pd.DataFrame({'Pressure (Pa)' : Pressure , 'Temperature (K)' : TempSim}, columns=['Pressure (Pa)', 'Temperature (K)'])
#PhasePts.to_csv('TempSimI_N2.dat', index=False, sep='\t')
PhasePts.to_csv('TempSimW_N2.dat', index=False, sep='\t')


# In[10]:


TempSim


# In[ ]:




