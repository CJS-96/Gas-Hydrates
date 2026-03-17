#!/usr/bin/env python
# coding: utf-8

import numpy as np
from scipy.optimize import fsolve
import math
import pandas as pd

lnfug = np.linspace(math.log(10),math.log(2.5e9),100)       ## Range of fugacity to calculate occupancies.
fug = np.exp(lnfug)
len=np.size(fug)

T=271.15
beta = 1/T
focci = np.zeros((len,2))       ## individual occupancy array
focc = np.zeros((len))          ## Total occupancy array
fsmall1 = np.zeros((len))       ## Contribution of first occupancy of small cage
flarge1 = np.zeros((len))       ## Contribution of first occupancy of large cage
flarge2 = np.zeros((len))       ## Contribution of second occupancy of large cage

### Gas Specific Info
### Argon:
#Lc = np.array([5.8233249788031E-07,1.44062161093766E-06,5.80040718576832E-15])                        ## Langmuir Constants
#aij = beta*np.array([[-100.128361799375, -36.4780370459466],[-76.7174404618716,-24.3812364686072]])   ## Guest-Guest interactions

### Nitrogen
Lc = np.array([1.66741366798993E-07,1.18368970284654E-06,5.87532363872445E-16])                        ## Langmuir Constants
aij = beta*np.array([[-105.864806350734, -38.5156603319702],[-80.7268703391869, -25.629540534039]])    ## Guest-Guest interactions

for i in range(0,len):                          ## Loop over fugacity points
    thetan = np.array([0,1]);                   ## Initial guesses

    def Occu(theta):                            ## Function to calculate occupancy
        
        return theta - np.array([((math.exp(-(aij[0,0]*theta[0] + aij[0,1]*theta[1]))*fug[i]*Lc[0])/(1 + (math.exp(-(aij[0,0]*theta[0] + aij[0,1]*theta[1]))*fug[i]*Lc[0]))), \
                (math.exp(-(aij[1,0]*theta[0] + aij[1,1]*theta[1]))*fug[i]*Lc[1] + (fug[i]**2)*Lc[2]*math.exp(-2*(aij[1,0]*theta[0] + aij[1,1]*theta[1])))/           \
                (1 + (math.exp(-(aij[1,0]*theta[0] + aij[1,1]*theta[1]))*fug[i]*Lc[1] + 0.5*(fug[i]**2)*Lc[2]*math.exp(-2*(aij[1,0]*theta[0] + aij[1,1]*theta[1]))))])
    
    focci[i,:] = fsolve(Occu,thetan)
    focc[i] = (2.0/3.0)*focci[i,0] + (1.0/3.0)*focci[i,1]
    fsmall1[i] = (2/3)*((math.exp(-(aij[0,0]*focci[i,1] + aij[0,1]*focci[i,0]))*fug[i]*Lc[0])/(1 + (math.exp(-(aij[0,0]*focci[i,0] + aij[0,1]*focci[i,0]))*fug[i]*Lc[0])))
    flarge1[i] = (1/3)*(math.exp(-(aij[1,0]*focci[i,1] + aij[1,1]*focci[i,1]))*fug[i]*Lc[1])/ \
                (1 + (math.exp(-(aij[1,0]*focci[i,1] + aij[1,1]*focci[i,1]))*fug[i]*Lc[1] + 0.5*(fug[i]**2)*Lc[2]*math.exp(-2*(aij[1,0]*focci[i,1] + aij[1,1]*focci[i,1]))))
    flarge2[i] = (1/3)*((math.exp(-2*(aij[1,0]*focci[i,1] + aij[1,1]*focci[i,1]))*fug[i]**2)*Lc[2])/ \
                (1 + (math.exp(-(aij[1,0]*focci[i,1] + aij[1,1]*focci[i,1]))*fug[i]*Lc[1] + 0.5*(fug[i]**2)*Lc[2]*math.exp(-2*(aij[1,0]*focci[i,1] + aij[1,1]*focci[i,1]))));
    
Occupancy = pd.DataFrame({'Fugacity (Pa)' : fug,'Total Occupancy' : focc, 'Small Cage Occupancy' : fsmall1, '1$^st$ Large Cage Occupancy' : flarge1, '2$^nd$ Large Cage Occupancy' : flarge2}, 
                          columns=['Fugacity (Pa)','Total Occupancy','Small Cage Occupancy','1$^st$ Large Cage Occupancy','2$^nd$ Large Cage Occupancy'])

#### Gas Specific name of the file ####
#Occupancy.to_csv('Occupancy_Ar.dat', index=False, sep='\t')
Occupancy.to_csv('Occupancy_N2.dat', index=False, sep='\t')
