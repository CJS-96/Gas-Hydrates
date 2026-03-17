"""
output.py — Write calculation results to space-padded .dat files.

Auto-generated filename format:
  {gas}_{structure}_{property}_{conditions}.dat

Examples:
  Methane_sI_dissociation_pressure_T273-290K.dat
  CO2_sI_occupancy_T280K_P5000000Pa.dat
  Nitrogen_sII_dissociation_temperature_P1e+05-1e+07Pa.dat

Files are plain-text, space-delimited, and can be loaded with:
  numpy.loadtxt(filename, comments="#")
  pandas.read_csv(filename, sep=r"\s+", comment="#")
"""

import datetime

_PROPERTY_LABELS = {
    0: "occupancy",
    1: "isotherm",
    2: "dissociation_temperature",
    3: "dissociation_pressure",
}

# Column gap (number of spaces between columns)
_GAP = "   "


# ── Number formatters ─────────────────────────────────────────────────────────

def _fmt_T(v):      return f"{v:.4f}"        # Temperature  e.g. 280.0000
def _fmt_P(v):      return f"{v:.6e}"        # Pressure     e.g. 5.000000e+06
def _fmt_theta(v):  return f"{v:.8f}"        # Occupancy    e.g. 0.48265048
def _fmt_lnfug(v):  return f"{v:.6f}"        # ln(fugacity) e.g. 12.345678


# ── Filename helpers ──────────────────────────────────────────────────────────

def _format_value(value):
    """Compact scalar string for use inside a filename."""
    if value == int(value):
        return str(int(value))
    return f"{value:.3g}"


def _ensure_dat(filename):
    """Append .dat extension if the filename has no extension."""
    if "." not in filename.split("/")[-1]:
        return filename + ".dat"
    return filename


def _make_filename(gas, structure, property_label, conditions):
    return f"{gas}_{structure}_{property_label}_{conditions}.dat"


def _conditions_point(T, P):
    return f"T{_format_value(T)}K_P{_format_value(P)}Pa"


def _conditions_array(values, unit, label):
    lo = f"{values[0]:.3g}"
    hi = f"{values[-1]:.3g}"
    return f"{label}{lo}-{hi}{unit}"


def auto_filename(prop, gas, structure, T=None, P=None, Pvalues=None, Tvalues=None):
    """Return the auto-generated filename for a calculation, without writing anything."""
    label = _PROPERTY_LABELS[prop]
    if prop in (0, 1):
        conditions = _conditions_point(T, P)
    elif prop == 2:
        conditions = _conditions_array(Pvalues, "Pa", "P")
    else:
        conditions = _conditions_array(Tvalues, "K", "T")
    return _make_filename(gas, structure, label, conditions)


# ── Table writer ──────────────────────────────────────────────────────────────

def _write_header(f, gas, structure, property_label, conditions):
    """Write metadata comment lines at the top of the file."""
    f.write("# VdWP-FL — Gas Hydrate Thermodynamics\n")
    f.write(f"# Generated:  {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"# Gas:        {gas}   Structure: {structure}\n")
    f.write(f"# Property:   {property_label}\n")
    f.write(f"# Conditions: {conditions}\n")
    f.write("#\n")


def _write_table(f, headers, rows):
    """Write an aligned, space-separated data table.

    Column widths are determined by the wider of the header label or the
    longest value in that column, so nothing ever overlaps.
    """
    # Compute column widths
    col_widths = [len(h) for h in headers]
    for row in rows:
        for i, val in enumerate(row):
            col_widths[i] = max(col_widths[i], len(val))

    # Header row
    f.write(_GAP.join(h.ljust(col_widths[i]) for i, h in enumerate(headers)) + "\n")

    # Separator line
    f.write(_GAP.join("-" * w for w in col_widths) + "\n")

    # Data rows
    for row in rows:
        f.write(_GAP.join(v.ljust(col_widths[i]) for i, v in enumerate(row)) + "\n")


# ── Public save functions ─────────────────────────────────────────────────────

def save_occupancy(T, P, frac_occ, gas, structure, filename=None):
    """Save single-point occupancy result.

    Parameters
    ----------
    T, P      : float   Temperature (K) and Pressure (Pa)
    frac_occ  : array   [theta_small, theta_large]
    filename  : str, optional   Custom filename; auto-generated if omitted.
    """
    conditions = _conditions_point(T, P)
    if filename is None:
        filename = _make_filename(gas, structure, "occupancy", conditions)
    else:
        filename = _ensure_dat(filename)

    headers = ["Temperature (K)", "Pressure (Pa)", "Theta_small (-)", "Theta_large (-)"]
    rows = [[_fmt_T(T), _fmt_P(P), _fmt_theta(frac_occ[0]), _fmt_theta(frac_occ[1])]]

    with open(filename, "w") as f:
        _write_header(f, gas, structure, "occupancy", conditions)
        _write_table(f, headers, rows)

    print(f"  Results saved to: {filename}")
    return filename


def save_isotherm(T, P, isotherm, gas, structure, filename=None):
    """Save occupancy isotherm (theta vs ln_fugacity).

    Parameters
    ----------
    isotherm : array  shape (3, N) — rows: theta_small, theta_large, ln_fug
    filename : str, optional   Custom filename; auto-generated if omitted.
    """
    conditions = _conditions_point(T, P)
    if filename is None:
        filename = _make_filename(gas, structure, "isotherm", conditions)
    else:
        filename = _ensure_dat(filename)

    theta_small = isotherm[0]
    theta_large = isotherm[1]
    ln_fug      = isotherm[2]

    headers = ["ln_fugacity (-)", "Theta_small (-)", "Theta_large (-)"]
    rows = [
        [_fmt_lnfug(lf), _fmt_theta(ts), _fmt_theta(tl)]
        for lf, ts, tl in zip(ln_fug, theta_small, theta_large)
    ]

    with open(filename, "w") as f:
        _write_header(f, gas, structure, "isotherm", conditions)
        _write_table(f, headers, rows)

    print(f"  Results saved to: {filename}")
    return filename


def save_dissociation_temperature(Pvalues, Tvalues, gas, structure, filename=None):
    """Save dissociation temperature vs pressure results.

    Parameters
    ----------
    filename : str, optional   Custom filename; auto-generated if omitted.
    """
    conditions = _conditions_array(Pvalues, "Pa", "P")
    if filename is None:
        filename = _make_filename(gas, structure, "dissociation_temperature", conditions)
    else:
        filename = _ensure_dat(filename)

    headers = ["Pressure (Pa)", "Dissociation Temperature (K)"]
    rows = [[_fmt_P(p), _fmt_T(t)] for p, t in zip(Pvalues, Tvalues)]

    with open(filename, "w") as f:
        _write_header(f, gas, structure, "dissociation_temperature", conditions)
        _write_table(f, headers, rows)

    print(f"  Results saved to: {filename}")
    return filename


def save_dissociation_pressure(Tvalues, Pvalues, gas, structure, filename=None):
    """Save dissociation pressure vs temperature results.

    Parameters
    ----------
    filename : str, optional   Custom filename; auto-generated if omitted.
    """
    conditions = _conditions_array(Tvalues, "K", "T")
    if filename is None:
        filename = _make_filename(gas, structure, "dissociation_pressure", conditions)
    else:
        filename = _ensure_dat(filename)

    headers = ["Temperature (K)", "Dissociation Pressure (Pa)"]
    rows = [[_fmt_T(t), _fmt_P(p)] for t, p in zip(Tvalues, Pvalues)]

    with open(filename, "w") as f:
        _write_header(f, gas, structure, "dissociation_pressure", conditions)
        _write_table(f, headers, rows)

    print(f"  Results saved to: {filename}")
    return filename
