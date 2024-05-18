#!/usr/bin/env python
# coding: utf-8

import numpy as np
from scipy.optimize import fsolve
import math
import pandas as pd
from Argon_EOS import Argasprop


######## Function to calculate I-H-V data ###############
def SimFlFh_Ice(T,*ZHPdata):
    Zterm, HconstP, P = ZHPdata
    
    R = 8.314
    t0 = 280
    
    ti0 = 260                ## ref temperature for ice
    pi0 = 2.531e5            ## ref pressure for ice
    fiw0 = 3.219             ## fugacity of ice at ti0 and pi0
    Fico = np.array([1.8701e-05,4.4465e-09,-5.9784e-16,-4.28e-18,-68558,18.74])*(1/R)   ## Constants for Ice Equation state 

    fhw0 = 31.2272           ## Fugacity of water in a real Ar-hydrate at ref temperature and pressure
    Temperature = np.linspace(200,260,13)  ### Enthalpy temperature range
    lnfhw0 = math.log(fhw0)
    OnebyT = 1/Temperature
    poly = np.polyfit(OnebyT,HconstP,3)    ### Fitting H integrand with 1/T as variable
    poly = poly/R


    #### Obtaining fugacity of water in the hydrate phase at input pressure and current temperature ####
    Hterm = poly[0]*(1/4)*((1/T)**4 - (1/t0)**4) + poly[1]*(1/3)*((1/T)**3 - (1/t0)**3)    \
            + poly[2]*(1/2)*((1/T)**2 - (1/t0)**2) + poly[3]*((1/T) - (1/t0))

    lnfwhS = Hterm + Zterm + lnfhw0      ## log(fug of water in hydrate phase at T and Pressure(a))

    #### Using fitted EOS for fugacity of water in liquid water ####
    lnfwi = Fico[0]*(P/T - pi0/ti0) + Fico[1]*(P-pi0) + 0.5*Fico[2]*(P**2/T - pi0**2/ti0)   \
            + 0.5*Fico[3]*(P**2 - pi0**2) + Fico[4]*(1/T - 1/ti0) - Fico[5]*math.log(T/ti0) + math.log(fiw0)

    return lnfwhS - lnfwi    ## Function equating fug of water in ice and hyd phases, to be minimized.

######## Function to calculate L-H-V data ###############
def SimFlFh_Water(T,*ZHPdata):
    Zterm, HconstP, P = ZHPdata
    
    R = 8.314
    t0 = 280
    
    tw0 = 280;               ## ref temperature for liquid water
    pw0 = 20e5;              ## ref pressure for liquid water
    flw0 = 24.35;            ## fugacity of ice at tw0 and pw0
    Flco = np.array([19.074e-06,-30.365e-10,-23.990e-15,56.020e-18,-74.170e03,58.817])*(1/R);

    fhw0 = 31.2272           ## Fugacity of water in a real Ar-hydrate at ref temperature and pressure
    Temperature = np.linspace(260,310,51)  ### Enthalpy temperature range

    lnfhw0 = math.log(fhw0)
    OnebyT = 1/Temperature
    poly = np.polyfit(OnebyT,HconstP,2)    ### Fitting H integrand with 1/T as variable
    poly = poly/R


    #### Obtaining fugacity of water in the hydrate phase at input pressure and current temperature ####
    Hterm = poly[0]*(1/3)*((1/T)**3 - (1/t0)**3) + poly[1]*(1/2)*((1/T)**2 - (1/t0)**2) + poly[2]*((1/T) - (1/t0))

    lnfwhS = Hterm + Zterm + lnfhw0  ## log(fug of water in hydrate phase at T and Pressure(a))

    #### Using fitted EOS for fugacity of water in liquid water ####
    lnfwl = Flco[0]*(P/T - pw0/tw0) + Flco[1]*(P-pw0) + 0.5*Flco[2]*(P**2/T - pw0**2/tw0)   \
            + 0.5*Flco[3]*(P**2 - pw0**2) + Flco[4]*(1/T - 1/tw0) - Flco[5]*math.log(T/tw0) + math.log(flw0)

    return lnfwhS - lnfwl   ## Function equating fug of water in liquid water and hyd phases, to be minimized.


############# Main program to calculate I-H-V data ##################

#### Reading Data for I-H-V
def Read_IHV_data():
    
    ETPHgasdata = pd.read_csv('ETPH_gas_Ar_Ice.dat',sep='\s+')       
    ETPHgas = ETPHgasdata.to_numpy()
    ETPHhyddata = pd.read_csv('ETPH_hyd_Ar_Ice.dat',sep='\s+')
    ETPHhyd = ETPHhyddata.to_numpy()

    ETPVhyddata = pd.read_csv('EPZpoints_all_Ar_Ice.dat', sep='\s+')
    ETPVhyd = ETPVhyddata.to_numpy()
    
    return ETPHgas, ETPHhyd, ETPVhyd

#### Read Data for L-H-V
def Read_LHV_data():
    ETPHgasdata = pd.read_csv('ETPH_gas_Ar_W.dat',sep='\s+')       
    ETPHgas = ETPHgasdata.to_numpy()
    ETPHhyddata = pd.read_csv('ETPH_hyd_Ar_W.dat',sep='\s+')
    ETPHhyd = ETPHhyddata.to_numpy()

    ETPVhyddata = pd.read_csv('EPZpoints_all_Ar_W.dat', sep='\s+')
    ETPVhyd = ETPVhyddata.to_numpy()
    
    return ETPHgas, ETPHhyd, ETPVhyd

############### Main Program starts here ####################

ETPHgas, ETPHhyd, ETPVhyd = Read_IHV_data()   ### To collect data for Ice-Hydrate equilibria
#ETPHgas, ETPHhyd, ETPVhyd = Read_LHV_data()  ### To collect data for Liquid water-Hydrate equilibria

Nwater=136
NA = 6.023E+023
R = 8.314

HRgas = ETPHgas[:,2]
Thyd = ETPHhyd[:,0]
Phyd = ETPHhyd[:,1]
Numhyd = ETPHhyd[:,2]
Ntothyd = Numhyd+Nwater                ### T and P and N details from hyd phase
IEhyd = ETPHhyd[:,3]
Volhyd = (ETPHhyd[:,4]*NA)/(Ntothyd)   ## Energies and Vol from hydrate simulations

xHgH = Numhyd/(Numhyd+Nwater)
xHwH = 1-xHgH                          ## hyd phase mol fractions
HRhyd = IEhyd + Phyd*Volhyd - R*Thyd   ## Residual Enthalpy (Hr) of hyd phase = U + PV - RT
Hcomb = HRhyd/xHwH - HRgas*(xHgH/xHwH) 
nd = np.size(Hcomb)
npH = 6                                  ## Number of Pressure points for each temperature in EPTempSim.
HcombMat = np.zeros((nd//npH,npH))

for a in range(0,(nd//npH)):
    A = Hcomb[npH*(a):npH*(a+1)]         ## Putting it into  T x P matrix.
    HcombMat[a,:] = np.transpose(A)


T0 = 245              ## Starting guess for finding T. Give a suitable initial guess, or even a range
t0 = 280              ## Reference temperature of Hydrate

#### Specify Pressures at which dissociation temperature is to be calculated ####
Pressure = 1e6*np.array([0.8,1.6,2.4,3.2,4.0,5.0])               ## Pressure (Ar) at which the eqm T is calculated -- Ice-Hydrate
#Pressure = 1e6*np.array([7.01,7.18,7.37,7.79,8.03,8.21,8.49])   ## Pressure (Ar) at which the eqm T is calculated -- Liquid Water-Hydrate
npE = np.size(Pressure)           ## Number of Pressure points            
TempSim = np.zeros((npE))         ## Dissociation temperature array initialized

nz=11                             ## Number of points in the integration of occupancy isotherm at reference temperature.
ZcombMat=np.zeros((nz,npE))       ## The isotherm is integrated from reference pressure to the input pressure.

for a in range(0,npE):            ## Loop over the input pressure array
#### Integrating Z term ####
    xHgV = ETPVhyd[nz*a:nz*(a+1),2]      ## Taking the values needed to evaluate
    xHwV = 1.0 - xHgV                    ## Z integrand at cont T=t0, upto different
    P = ETPVhyd[nz*a:nz*(a+1),1]         ## P0 to Pressure(a).
    Vhyd = ETPVhyd[nz*a:nz*(a+1),4]
    HRgas, fug, Vgas = Argasprop(t0,P)   ## Calculating fugacity of Argon at t0 and P.
    Vgas = np.reshape(Vgas,-1)
    Vcomb = Vhyd/xHwV - Vgas*(xHgV/xHwV)
    Zcomb = Vcomb*P/(R*t0)               
    ZcombMat[:,a]=Zcomb
    
    #### Isotherm integration from P=P0 to P=Pressure(a) at T=t0 ####
    bunch1 = Zcomb[0] + Zcomb[nz-1]      
    bunch2 = 0
    bunch4 = 0
    for j in range(2,nz-2,2):            
        bunch2 = bunch2 + Zcomb[j]
      
    for j in range(1,nz-1,2):            
        bunch4 = bunch4 + Zcomb[j]

    Zh = (math.log(P[nz-1]) - math.log(P[0]))/(nz-1)
    Zterm = (Zh/3)*(bunch1 + 4*bunch4 + 2*bunch2)

    HconstP = HcombMat[:,a]

    ZHPdata = (Zterm,HconstP,Pressure[a])
    
#### Uncomment the respective function to get equilibria with either ice or liquid water ####
    TempSim[a] = fsolve(SimFlFh_Ice,T0,args=ZHPdata)
    #TempSim[a] = fsolve(SimFlFh_Water,T0,args=ZHPdata)

#### Printing temperature values to file ####
PhasePts = pd.DataFrame({'Pressure (Pa)' : Pressure , 'Temperature (K)' : TempSim}, columns=['Pressure (Pa)', 'Temperature (K)'])
PhasePts.to_csv('TempSimI_Ar.dat', index=False, sep='\t')
#PhasePts.to_csv('TempSimW_Ar.dat', index=False, sep='\t')





