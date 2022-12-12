import os as os
import numpy as np
import pandas as pd
from pysagas.wrappers import Cart3DWrapper
from pysagas.utilities import isentropic_dPdp


def run_main(data_path):
    """Analyse the Cart3D solution."""

    # Filepaths
    sensitivity_filepath = os.path.join(data_path, "combined.csv")
    celldata_filepath = os.path.join(data_path, "cells.csv")
    pointdata_filepath = os.path.join(data_path, "points.csv")

    wrapper = Cart3DWrapper(sensitivity_filepath, celldata_filepath, pointdata_filepath)
    F_sense = wrapper.calculate()
    # F_sense = wrapper.calculate(isentropic_dPdp, P_inf=-24)

    M_inf = 6
    A_ref = 1  # m2

    rho_lookup = {
        5: 0.0450814,
        6: 0.0308742,
        7: 0.0225510,
        7.5: 0.0195525,
        8: 0.0171295,
        9: 0.0134424,
        10: 0.0107105,
    }
    a_lookup = {
        5: 297.891,
        6: 299.499,
        7: 300.840,
        7.5: 301.575,
        8: 302.018,
        9: 303.061,
        10: 305.562,
    }
    a_inf = a_lookup[M_inf]
    rho_inf = rho_lookup[M_inf]
    V_inf = M_inf * a_inf

    # Non-dimensionalise
    coef_sens = F_sense / (0.5 * rho_inf * A_ref * V_inf**2)

    # Print results
    print("PySAGAS Result:")
    print("      dFx/dP           dFy/dP          dFz/dP")
    print(coef_sens)

    print("Cart 3D Finite Difference Result:")
    c3d_sens = np.array([[0.147517, 0.126153, 0]])
    print(c3d_sens)

    print("\nError (%):")
    errors = np.nan_to_num(100 * (coef_sens - c3d_sens) / c3d_sens, posinf=0, neginf=0)
    print(errors)

    return errors


def cart3d_fd(data_path):
    """Calculate differences from Cart3D result using finite
    differencing."""
    # Load data
    cols = {
        "thickness": "thickness",
        "C_D-entire": "CD",
        "C_L-entire": "CL",
        "C_m-entire": "Cm",
    }
    data = pd.read_csv(os.path.join(data_path, "sims.csv"), usecols=cols.keys())
    data.rename(columns=cols, inplace=True)

    # Transform into CX, CY and CZ
    aoa = np.deg2rad(3)
    data["CX"] = -data["CL"] * np.sin(aoa) + data["CD"] * np.cos(aoa)
    data["CY"] = data["CL"] * np.cos(aoa) + data["CD"] * np.sin(aoa)

    # Calculate sensitivities
    pc1 = (data.iloc[4] - data.iloc[3]) / (
        data.iloc[4]["thickness"] - data.iloc[3]["thickness"]
    )


def test_ramp():
    """Run the test."""
    file_dir = os.path.abspath(os.path.dirname(__file__))
    data_path = os.path.join(file_dir, "data")
    errors = run_main(data_path)

    assert np.max(np.abs(errors)) < 40, "Adjoints inaccurately calculated"


if __name__ == "__main__":
    np.seterr(all="ignore")
    data_path = "tests/data"
    run_main(data_path)
