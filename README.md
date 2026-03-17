# Gas Hydrates — Thermodynamic Property Calculations

Welcome to this repository! This page contains Python codes for computing thermodynamic
properties of gas hydrates using the van der Waals-Platteeuw Flexible Lattice (vdWP-FL) theory.

The vdWP-FL theory accurately predicts phase equilibria as well as cage occupancy of gas
hydrates, including those of multiple occupancy hydrates. The current implementation supports
pure gases; inclusion of gas mixtures is planned for a future version.

---

## Repository Contents

### `vdwpfl` — Recommended Package

A clean, installable Python package for computing gas hydrate properties using the parameterized vdWP-FL theory.
This is the recommended way to use the code. It supports interactive and batch modes, auto-generated
output files, and works on any operating system.

**Supported gases:** Nitrogen, Methane, Ethane, CO2, Propane, Iso-butane

**Supported calculations:**
- Cage occupancy at a given temperature and pressure
- Occupancy isotherm
- Dissociation temperature as a function of pressure
- Dissociation pressure as a function of temperature

See [INSTALL.md](INSTALL.md) for installation and usage instructions.

## References:
- Shah C. J. and Punnathanam S. N., *J. Phys. Chem. C*, 2025, **129**, 16434–16444
- Veesam S. K. and Punnathanam S. N., *J. Phys. Chem. C*, 2019, **123**, 26406–26414

---

### `Shah_JPCC-2025_PaperCodes` — Original Paper Codes

Contains the original codes used for the calculations in:

> Shah C. J. and Punnathanam S. N., *J. Phys. Chem. C*, 2025, **129**, 16434–16444

There are two folders:
1. Simulation: Contains the code to calculate phase equilibria from simulation data.
2. vdWP-FL_Theory: Contains code to calculate the phase equilibria using the non-parameterized vdWP-FL theory.

These are made available for reproducibility of the results in fig. 1 and 2 of the paper.

---

## Contact

For questions or feedback, feel free to reach out:

- chaitanyasj@iisc.ac.in
- chaitanyashah19@gmail.com
