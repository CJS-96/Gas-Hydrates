import numpy as np
from importlib.resources import files


def _open_data(filename):
    """Open a data file bundled inside the vdwpfl.data package."""
    return files("vdwpfl.data").joinpath(filename).open("r")


def Read_LC_Coeffs(name):
    with _open_data("Langmuir_Coefficients.dat") as f:
        lines = [line.rstrip('') for line in f]

    n = 0
    for line in lines:
        if line.strip() == name:
            start = n
            break
        n = n + 1

    coeffs_raw = lines[start+1]
    c11, c12, c13, c14, c15 = [float(c) for c in coeffs_raw.split(" ")]
    coeffs_raw = lines[start+2]
    c21, c22, c23, c24, c25 = [float(c) for c in coeffs_raw.split(" ")]
    coeffs_raw = lines[start+3]
    c31, c32, c33, c34, c35 = [float(c) for c in coeffs_raw.split(" ")]

    LC_coeffs = np.array([[c11, c12, c13, c14, c15],[c21, c22, c23, c24, c25],[c31, c32, c33, c34, c35]])

    return LC_coeffs


def Read_Ugg(name):
    with _open_data("Ugg.dat") as f:
        lines = f.readlines()

    n = 0
    for line in lines:
        if line.strip() == name:
            start = n
            break
        n = n + 1

    coeffs_raw = lines[start+1]
    aSS, aSL, aLS, aLL = [float(c) for c in coeffs_raw.split(" ")]

    aij = np.array([[aSS, aSL],[aLS, aLL]])

    return aij


def Read_Gas_Props(name):
    with _open_data("Gas_Properties.dat") as f:
        lines = f.readlines()

    n = 0
    for line in lines:
        if line.strip() == name:
            start = n
            break
        n = n + 1

    values = lines[start+1]
    Mass, w, Tc, Pc = [float(v) for v in values.split(" ")]

    return Mass, w, Tc, Pc
