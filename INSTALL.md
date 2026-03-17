# vdWP-FL: van der Waals-Platteeuw Flexible Lattice Theory for Gas Hydrates

A Python package for computing thermodynamic properties of gas hydrates using the
van der Waals-Platteeuw Flexible Lattice (vdWP-FL) theory.

**Supported gases:** Nitrogen, Methane, Ethane, CO2, Propane, Iso-butane

**Supported calculations:**
- Cage occupancy at a given temperature and pressure
- Occupancy isotherm
- Dissociation temperature as a function of pressure
- Dissociation pressure as a function of temperature

---

## Installation

Requires Python 3.9 or higher.

```
pip install git+https://github.com/your-username/vdwpfl.git
```

---

## Usage

### Interactive mode

Simply run:

```
vdwpfl
```

The program will guide you through selecting the gas, calculation type, and input conditions using arrow-key menus.

### Batch mode (no prompts)

```
# Dissociation pressure for Methane over a temperature range
vdwpfl --gas Methane --property 3 --T-range 273,290,10

# Dissociation temperature for Nitrogen over a pressure range (log spacing)
vdwpfl --gas Nitrogen --property 2 --P-range 1e5,1e7,20

# Occupancy at a single point
vdwpfl --gas CO2 --property 0 --T 280 --P 5000000

# Occupancy with a custom output filename
vdwpfl --gas CO2 --property 0 --T 280 --P 5000000 --output my_results
```

### Command-line options

| Option | Description |
|---|---|
| `--gas GAS` | Gas species (Nitrogen, Methane, Ethane, CO2, Propane, Iso-butane) |
| `--property N` | 0=occupancy, 1=isotherm, 2=dissociation T, 3=dissociation P |
| `--T VALUE` | Temperature in K (for property 0 or 1) |
| `--P VALUE` | Pressure in Pa (for property 0 or 1) |
| `--T-values a,b,c` | Comma-separated temperatures in K (for property 3) |
| `--T-range s,e,n` | Temperature range: start, end, npoints — linear spacing (for property 3) |
| `--P-values a,b,c` | Comma-separated pressures in Pa (for property 2) |
| `--P-range s,e,n` | Pressure range: start, end, npoints — log spacing (for property 2) |
| `--output FILE` | Custom output filename (.dat added automatically if omitted) |

---

## Output

Results are written to space-aligned `.dat` files in the current directory.
The filename is auto-generated from the gas name, structure, property, and conditions, e.g.:

```
Methane_sI_dissociation_pressure_T273-290K.dat
Nitrogen_sII_occupancy_T280K_P5000000Pa.dat
```

Files can be loaded with:
```python
import numpy as np
data = np.loadtxt("Methane_sI_dissociation_pressure_T273-290K.dat", comments="#")
```

---

## Citation

If you use this package in your research, please cite:

- Shah C. J. and Punnathanam S. N., *J. Phys. Chem. C*, 2025, **129**, 16434–16444
- Veesam S. K. and Punnathanam S. N., *J. Phys. Chem. C*, 2019, **123**, 26406–26414
