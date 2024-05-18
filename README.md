# Description for the codes to calculate gas hydrates properties

There are two directories in this repository: Simulation/ and vdWP-FL_Theory/

The codes provided in the vdWP-FL_Theory/ are useful to:
1. get the cage occupancies 
2. get the dissociation temperature of gas hydrates for a given Pressure.
using vdWP-Fl theory.

The codes provided in the Simulation/ are useful to:
1. get the dissociation temperature of gas hydrates for a given Pressure.
using simulation data.

## ########################################################
## Description about the files in vdWP-Fl_Theory/
## ########################################################

1. Occupancy_Calculations.py : Python script for obtaining cage occupancies for gas hydrates.
2. Ar_PhaseEqm_Theory.py, N2_PhaseEqm_Theory.py : These files are python scripts to calculate the dissociation temperatures for the input pressure values.
3. Argon_EOS.py, Nitrogen_EOS.py : These are python scripts to calculate argon gas and nitrogen gas properties from equations of state.
4. LCP1_Ar_IW.dat, LCP1_N2_IW.dat : These files contain the constants from Langmuir constant fit.

These individual files contain comments on the working of the codes which can be of further help.

## ########################################################
## Description about the files in Simulation/
## ########################################################
To use these codes to obtain dissociation temperature at any given pressure, we need to simulation data corresponding to that pressure as described in point 3 and 4 below.

1. Ar_PhaseEqm_Sim.py, N2_PhaseEqm_Sim.py : These files are python scripts to calculate dissociation temperatures for input pressure values but using the GCMC simulation data.
2. EPZpoints_all_* : These files contain simulation data - guest mole fraction, water mole fraction and volume of the system, at reference temperature, with pressures varying from reference pressure to the input pressure.
3. ETPH_gas_* : These files contain residual enthalpy data of gases calculated from equations of state.
4. Argon_EOS.py, Nitrogen_EOS.py : These are python scripts to calculate argon gas and nitrogen gas properties from equations of state.
