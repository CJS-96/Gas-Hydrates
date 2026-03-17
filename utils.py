import numpy as np
import questionary

# ── Gas and calculation metadata ─────────────────────────────────────────────

_GAS_NAMES = ["Nitrogen", "Methane", "Ethane", "CO2", "Propane", "Iso-butane"]

_CALC_DESCRIPTIONS = [
    "Occupancy at a single (T, P) point",
    "Isotherm for a (T, P) point",
    "Dissociation temperatures for given pressures",
    "Dissociation pressures for given temperatures",
]

# ── Interactive input helpers ─────────────────────────────────────────────────

def select_gas():
    """Arrow-key menu for gas selection. Returns integer index (0–5)."""
    choices = [f"{i}:  {name}" for i, name in enumerate(_GAS_NAMES)]
    answer = questionary.select(
        "Select the hydrate former:",
        choices=choices,
    ).ask()
    if answer is None:
        raise KeyboardInterrupt
    return int(answer.split(":")[0])


def select_calculation():
    """Arrow-key menu for calculation type. Returns integer index (0–3)."""
    choices = [f"{i}:  {desc}" for i, desc in enumerate(_CALC_DESCRIPTIONS)]
    answer = questionary.select(
        "Select type of calculation:",
        choices=choices,
    ).ask()
    if answer is None:
        raise KeyboardInterrupt
    return int(answer.split(":")[0])


def prompt_float(message):
    """Prompt for a single float; re-prompts on invalid input."""
    while True:
        try:
            return float(input(message))
        except ValueError:
            print("  Invalid input — please enter a numeric value.")


def prompt_values(quantity, unit):
    """Prompt for an array or range of values; returns a numpy array.

    Parameters
    ----------
    quantity : str  e.g. "pressure" or "temperature"
    unit     : str  e.g. "Pa" or "K"

    Returns
    -------
    numpy.ndarray
    """
    choice = questionary.select(
        f"How would you like to enter {quantity} values?",
        choices=[
            "Array  —  enter comma-separated values",
            "Range  —  enter start, end, npoints",
        ],
    ).ask()
    if choice is None:
        raise KeyboardInterrupt

    if choice.startswith("Array"):
        while True:
            try:
                raw = input(f"  Enter {quantity} values ({unit}), separated by commas: ")
                values = np.array([float(v.strip()) for v in raw.split(",")])
                if len(values) == 0:
                    raise ValueError("No values entered.")
                return values
            except ValueError as e:
                print(f"  Invalid input: {e}  Please enter numbers separated by commas.")
    else:
        while True:
            try:
                raw = input(f"  Enter start, end, npoints for {quantity} ({unit}): ")
                parts = [v.strip() for v in raw.split(",")]
                if len(parts) != 3:
                    raise ValueError("Expected exactly 3 values.")
                start, end, npoints = float(parts[0]), float(parts[1]), int(float(parts[2]))
                if npoints < 2:
                    raise ValueError("npoints must be at least 2.")
                return start, end, npoints   # caller decides spacing (log vs linear)
            except ValueError as e:
                print(f"  Invalid input: {e}  Please enter: start, end, npoints.")


def prompt_filename(auto_name):
    """Show the auto-generated filename and let the user accept or override it.

    Returns the chosen filename (always ends with .csv).
    """
    choice = questionary.select(
        "Output filename:",
        choices=[
            f"Auto-generated: {auto_name}",
            "Enter a custom filename",
        ],
    ).ask()
    if choice is None:
        raise KeyboardInterrupt

    if choice.startswith("Auto"):
        return auto_name

    while True:
        name = input("  Enter filename (the .dat extension will be added if omitted): ").strip()
        if name:
            return name if "." in name.split("/")[-1] else name + ".dat"
        print("  Filename cannot be empty.")


def prompt_continue():
    """Ask whether to perform another calculation. Returns bool."""
    answer = questionary.confirm("Perform another calculation?", default=True).ask()
    if answer is None:
        return False
    return answer


# ── Legacy display functions (kept for reference) ────────────────────────────

def show_gases():
    print("\nSelect the hydrate former: \n")
    for i, name in enumerate(_GAS_NAMES):
        print(f"{i}: {name}")


def show_calculations():
    print("\n\nSelect type of Calculation: \n")
    for i, desc in enumerate(_CALC_DESCRIPTIONS):
        print(f"{i}: {desc}")
    print("Note: T - Temperature and P - Pressure")


# ── Data helpers ──────────────────────────────────────────────────────────────

def Get_Gas_name(Gas):
    available_gases = ("Nitrogen", "Methane", "Ethane", "CO2", "Propane", "Iso-butane")
    name = available_gases[Gas]

    if name in ("Nitrogen", "Propane", "Iso-butane"):
        structure = "sII"
    elif name in ("Methane", "Ethane", "CO2"):
        structure = "sI"
    else:
        raise ValueError("Unknown Gas")

    return name, structure


def Get_Structure_Info(structure):
    if structure == "sI":
        Nwaters, Ncavities, Nsmall, Nlarge = 46, 8, 6, 2
    elif structure == "sII":
        Nwaters, Ncavities, Nsmall, Nlarge = 136, 24, 16, 8
    else:
        raise ValueError("Unknown hydrate structure")

    return Nwaters, Ncavities, Nsmall, Nlarge


# ── Numerical utilities ───────────────────────────────────────────────────────

def Simpsons_OneThird_Rule(Integrands, h):
    nsteps = len(Integrands) - 1
    bunch1 = Integrands[0] + Integrands[-1]
    bunch2 = 0
    bunch4 = 0
    for i in range(2, nsteps - 1, 2):
        bunch2 += Integrands[i]
    for i in range(1, nsteps, 2):
        bunch4 += Integrands[i]
    return (h / 3) * (bunch1 + 4 * bunch4 + 2 * bunch2)
