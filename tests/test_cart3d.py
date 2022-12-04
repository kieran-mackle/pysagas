import os as os
import numpy as np
import pandas as pd
from adjoint.flow import FlowState
from adjoint.utilities import all_dfdp
from adjoint.geometry import Vector, Cell


def run_main(data_path):
    # Configuration
    parameters = ["P"]  # TODO - automate extraction of parameters

    # Filepaths
    sensitivity_filepath = os.path.join(data_path, "combined.csv")
    celldata_filepath = os.path.join(data_path, "cells.csv")
    pointdata_filepath = os.path.join(data_path, "points.csv")

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

    # Load data
    sensdata = pd.read_csv(sensitivity_filepath)
    celldata = pd.read_csv(celldata_filepath)
    pointdata = pd.read_csv(pointdata_filepath)

    coordinates = pointdata[["Points_0", "Points_1", "Points_2"]]
    cell_vertex_ids = celldata[["Point Index 0", "Point Index 1", "Point Index 2"]]

    params_sens_cols = []
    for p in parameters:
        for d in ["x", "y", "z"]:
            params_sens_cols.append(f"d{d}d_{p}")

    # Construct cells
    cells = []
    for cell in celldata.index:
        vertex_ids = cell_vertex_ids.loc[cell].values
        cell_v_coords = [coordinates.loc[v_id].values for v_id in vertex_ids]
        vertices = [Vector.from_coordinates(c) for c in cell_v_coords]

        # Extract sensitivity information for the cell
        vertex_sensitivity_info = [sensdata.loc[v_id] for v_id in vertex_ids]

        dvdp = np.empty((9, len(parameters)))
        for i, p in enumerate(parameters):
            r = 0
            for v in vertex_sensitivity_info:
                for c in ["x", "y", "z"]:
                    dvdp[r, i] = v[f"d{c}d{p}"]
                    r += 1

        # Create Cell
        newcell = Cell.from_points(vertices)

        # Add geometry sensitivity information
        newcell._add_sensitivities(np.array(dvdp))

        # Add flow state at cell
        cell_velocity = Vector.from_coordinates(
            celldata.loc[cell][["U", "V", "W"]].values
        )
        cell_flowstate = FlowState(
            mach=celldata.loc[cell]["M"],
            pressure=celldata.loc[cell]["p"],
            temperature=celldata.loc[cell]["T"],
            direction=cell_velocity,
        )
        newcell.flowstate = cell_flowstate

        # Append new Cell
        cells.append(newcell)

    # Calculate force sensitivity
    F_sense = all_dfdp(cells=cells)

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

    assert np.max(np.abs(errors)) < 40, "Adjoints inaccurately calculated"


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
    file_dir = os.path.abspath(os.path.dirname(__file__))
    data_path = os.path.join(file_dir, "data")
    run_main(data_path)


if __name__ == "__main__":
    np.seterr(all="ignore")
    path = "tests/data"
    run_main(path)
