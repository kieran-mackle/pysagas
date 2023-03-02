import os as os
import numpy as np
import pandas as pd
from pysagas.geometry import Vector
from pysagas.wrappers import Cart3DWrapper


def run_main(data_path):
    """Analyse the Cart3D solution."""
    # Define flow conditions
    M_inf = 6
    L_ref = 1  # m
    A_ref = 1  # m2
    a_inf = 299.499  # m/s
    rho_inf = 0.0308742  # kg/m3
    V_inf = M_inf * a_inf  # m/s

    # Filepaths
    sensitivity_filepath = os.path.join(data_path, "combined.csv")

    # Data
    pointdata = pd.read_csv(os.path.join(data_path, "points.csv"))
    celldata = pd.read_csv(os.path.join(data_path, "cells.csv"))

    # Create PySAGAS wrapper
    wrapper = Cart3DWrapper(
        a_inf=a_inf,
        rho_inf=rho_inf,
        sensitivity_filepath=sensitivity_filepath,
        pointdata=pointdata,
        celldata=celldata,
    )

    # Calculate sensitivities
    cog = Vector(2.25, 0, 0)
    F_sense, M_sense = wrapper.calculate(cog=cog)

    # Non-dimensionalise
    coef_sens = F_sense / (0.5 * rho_inf * A_ref * V_inf**2)
    M_coef_sens = M_sense / (0.5 * rho_inf * A_ref * L_ref * V_inf**2)

    print("\nCart 3D Finite Difference Result:")
    c3d_sens = np.array([[0.147517, 0.126153, 0]])
    print(c3d_sens)

    # Print results
    print("\nPySAGAS Result:")
    print(coef_sens)
    print("")
    print(M_coef_sens)

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

    return pc1


def test_cart3d_wedge():
    """Run the test."""
    file_dir = os.path.abspath(os.path.dirname(__file__))
    data_path = os.path.join(file_dir, "data")
    errors = run_main(data_path)

    assert np.max(np.abs(errors)) < 5.6, "Adjoints inaccurately calculated"


if __name__ == "__main__":
    # cart3d_fd("tests/data")
    np.seterr(all="ignore")
    data_path = "tests/data"
    run_main(data_path)
