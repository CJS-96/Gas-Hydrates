import numpy as np
import argparse
import globe
from utils import (
    select_gas, select_calculation, prompt_float, prompt_values,
    prompt_filename, prompt_continue,
    Get_Gas_name, Get_Structure_Info,
)
from read_files import Read_LC_Coeffs, Read_Ugg, Read_Gas_Props
from occupancy import Occupancies
from phase_equilibrium import DissociationTemperature, DissociationPressure
from output import (
    auto_filename,
    save_occupancy, save_isotherm,
    save_dissociation_temperature, save_dissociation_pressure,
)


# ── Banner messages ───────────────────────────────────────────────────────────

def Introduction():
    print("╔══════════════════════════════════════════════════════════╗")
    print("║       VdWP-FL: Gas Hydrate Thermodynamics                ║")
    print("╚══════════════════════════════════════════════════════════╝")
    print()


def ExitMessage():
    print()
    print("╔══════════════════════════════════════════════════════════╗")
    print("║  Thank you for using VdWP-FL!                            ║")
    print("║  If this work was useful, please cite:                   ║")
    print("║    Shah & Punnathanam, JPCC 2025, 129, 16434             ║")
    print("║    Veesam & Punnathanam, JPCC 2019, 123, 26406           ║")
    print("╚══════════════════════════════════════════════════════════╝")
    print()


# ── Gas setup ─────────────────────────────────────────────────────────────────

def _setup_gas(gas_index):
    """Load all parameters for the selected gas into the globe namespace."""
    globe.Gas = gas_index
    globe.name, globe.structure = Get_Gas_name(globe.Gas)
    struct_info = Get_Structure_Info(globe.structure)
    globe.LC_coeffs = Read_LC_Coeffs(globe.name)
    globe.aij = Read_Ugg(globe.name)
    GasProperties = Read_Gas_Props(globe.name)

    globe.Nwaters, Ncavities, Nsmall, Nlarge = struct_info
    globe.R_small_cav = Nsmall / Ncavities
    globe.R_large_cav = Nlarge / Ncavities
    globe.R_cav_water = Ncavities / globe.Nwaters

    globe.Mass, globe.w, globe.Tc, globe.Pc = GasProperties

    # Reference state parameters
    globe.R       = 8.314           # J/mol/K
    globe.t0      = 273.15          # K
    globe.p0      = 0               # Pa
    globe.cp0_bl  = -38.12 / globe.R    # unitless
    globe.cp0_bi  = 0.0
    globe.b_bl    = 0.141  / globe.R    # 1/K
    globe.b_bi    = 0.0
    globe.delh_il0 = 6009.5 / globe.R  # K

    if globe.structure == 'sI':
        globe.delh_bi0 = 1058 / globe.R
        globe.mu_bi0   = 834  / globe.R / globe.t0
    elif globe.structure == 'sII':
        globe.delh_bi0 = 1064 / globe.R
        globe.mu_bi0   = 709  / globe.R / globe.t0
    else:
        raise ValueError("Unknown structure. Cannot set reference parameters.")

    globe.delh_bl0 = globe.delh_bi0 - globe.delh_il0

    print(f"\n  Gas loaded: {globe.name}  |  Structure: {globe.structure}\n")


# ── Calculation runner ────────────────────────────────────────────────────────

def _run_property(prop, T=None, P=None, Pvalues=None, Tvalues=None, output_file=None):
    """Execute the chosen calculation, print results, and save to a CSV file."""
    globe.Property = prop

    if prop == 0:
        result = Occupancies(T, P)
        print(f"\nFractional Occupancies at T = {T} K, P = {P} Pa:")
        print(f"  Small cavity: {result[0]:.6f}")
        print(f"  Large cavity: {result[1]:.6f}")
        save_occupancy(T, P, result, globe.name, globe.structure, filename=output_file)

    elif prop == 1:
        result = Occupancies(T, P)
        print(f"\nIsotherm at T = {T} K, P = {P} Pa:")
        print(result)
        save_isotherm(T, P, result, globe.name, globe.structure, filename=output_file)

    elif prop == 2:
        Tresult = DissociationTemperature(Pvalues)
        print("\nDissociation Temperatures (K):")
        for p, t in zip(Pvalues, Tresult):
            print(f"  P = {p:.4e} Pa  →  T = {t:.4f} K")
        save_dissociation_temperature(Pvalues, Tresult, globe.name, globe.structure, filename=output_file)

    elif prop == 3:
        Presult = DissociationPressure(Tvalues)
        print("\nDissociation Pressures (Pa):")
        for t, p in zip(Tvalues, Presult):
            print(f"  T = {t:.4f} K  →  P = {p:.4e} Pa")
        save_dissociation_pressure(Tvalues, Presult, globe.name, globe.structure, filename=output_file)


# ── Interactive mode ──────────────────────────────────────────────────────────

def Theory():
    """Run the program interactively, prompting the user for all inputs."""
    while True:
        # --- Gas selection ---
        gas_index = select_gas()
        _setup_gas(gas_index)

        # --- Calculation type ---
        prop = select_calculation()
        globe.Property = prop

        # --- Inputs for the selected calculation ---
        if prop in (0, 1):
            T = prompt_float("  Enter temperature (K): ")
            P = prompt_float("  Enter pressure    (Pa): ")
            Pvalues = Tvalues = None

        elif prop == 2:
            result = prompt_values("pressure", "Pa")
            if isinstance(result, np.ndarray):
                Pvalues = result
            else:
                start, end, npoints = result
                Pvalues = np.exp(np.linspace(np.log(start), np.log(end), npoints))
            T = P = Tvalues = None

        elif prop == 3:
            result = prompt_values("temperature", "K")
            if isinstance(result, np.ndarray):
                Tvalues = result
            else:
                start, end, npoints = result
                Tvalues = np.linspace(start, end, npoints)
            T = P = Pvalues = None

        # --- Output filename ---
        fname = auto_filename(prop, globe.name, globe.structure,
                              T=T, P=P, Pvalues=Pvalues, Tvalues=Tvalues)
        output_file = prompt_filename(fname)

        _run_property(prop, T=T, P=P, Pvalues=Pvalues, Tvalues=Tvalues,
                      output_file=output_file)

        # --- Continue? ---
        if not prompt_continue():
            break


# ── Batch mode ────────────────────────────────────────────────────────────────

_GAS_NAMES = ["Nitrogen", "Methane", "Ethane", "CO2", "Propane", "Iso-butane"]


def run_batch(args):
    """Run a single calculation using command-line arguments (no prompts)."""
    gas_index = _GAS_NAMES.index(args.gas)
    _setup_gas(gas_index)
    prop = args.property
    globe.Property = prop

    if prop in (0, 1):
        if args.T is None or args.P is None:
            raise ValueError("--T and --P are required for --property 0 or 1.")
        _run_property(prop, T=args.T, P=args.P, output_file=args.output)

    elif prop == 2:
        if args.P_values:
            Pvalues = np.array([float(v) for v in args.P_values.split(",")])
        elif args.P_range:
            start, end, npoints = [float(v) for v in args.P_range.split(",")]
            Pvalues = np.exp(np.linspace(np.log(start), np.log(end), int(npoints)))
        else:
            raise ValueError("--P-values or --P-range is required for --property 2.")
        _run_property(prop, Pvalues=Pvalues, output_file=args.output)

    elif prop == 3:
        if args.T_values:
            Tvalues = np.array([float(v) for v in args.T_values.split(",")])
        elif args.T_range:
            start, end, npoints = [float(v) for v in args.T_range.split(",")]
            Tvalues = np.linspace(start, end, int(npoints))
        else:
            raise ValueError("--T-values or --T-range is required for --property 3.")
        _run_property(prop, Tvalues=Tvalues, output_file=args.output)


def parse_args():
    parser = argparse.ArgumentParser(
        description="VdWP-FL: Thermodynamic properties of gas hydrates.\n"
                    "Run without arguments for interactive mode.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode (default)
  python vdWP_FL.py

  # Batch: dissociation pressure for Methane over a temperature range
  python vdWP_FL.py --gas Methane --property 3 --T-range 273,290,10

  # Batch: occupancy at a single point, custom output filename
  python vdWP_FL.py --gas CO2 --property 0 --T 280 --P 5000000 --output my_results
""",
    )
    parser.add_argument(
        "--gas",
        choices=_GAS_NAMES,
        metavar="GAS",
        help=f"Gas species. Choices: {', '.join(_GAS_NAMES)}",
    )
    parser.add_argument(
        "--property", type=int, choices=[0, 1, 2, 3],
        help="Calculation type: 0=occupancy, 1=isotherm, 2=dissociation T, 3=dissociation P",
    )
    parser.add_argument("--T", type=float, help="Temperature in K (for --property 0 or 1)")
    parser.add_argument("--P", type=float, help="Pressure in Pa    (for --property 0 or 1)")
    parser.add_argument(
        "--T-values", dest="T_values",
        help="Comma-separated temperatures in K (for --property 3)",
    )
    parser.add_argument(
        "--T-range", dest="T_range",
        help="Temperature range as start,end,npoints in K — linear spacing (for --property 3)",
    )
    parser.add_argument(
        "--P-values", dest="P_values",
        help="Comma-separated pressures in Pa (for --property 2)",
    )
    parser.add_argument(
        "--P-range", dest="P_range",
        help="Pressure range as start,end,npoints in Pa — log spacing (for --property 2)",
    )
    parser.add_argument(
        "--output", metavar="FILE",
        help="Output filename (default: auto-generated from gas, property, and conditions). "
             "The .csv extension is added automatically if omitted.",
    )
    return parser.parse_args()


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    args = parse_args()
    Introduction()
    try:
        if args.gas is not None:
            run_batch(args)
        else:
            Theory()
    except KeyboardInterrupt:
        print("\n\n  Interrupted by user.")
    ExitMessage()
