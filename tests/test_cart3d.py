import os as os
import numpy as np
import pandas as pd
import gdtk.ideal_gas_flow as igf
import xml.etree.ElementTree as ET

from adjoint.flow import FlowState
from adjoint.utilities import (
    calculate_pressures,
    calculate_force_vector,
    all_dfdp,
)
from adjoint.geometry import (
    Vector,
    Cell,
)


# ------------------------------------------------------
# Function definition
# ------------------------------------------------------

# create matrix linking vertex positons to parameter
def create_positon_sense_matrix(A, B, C, N_parameters):
    # [[dA.x/dp0, dA.x/dp1 ....],
    #  [dA.y/dp0, .....
    #  ...
    # [dC.z/dp0, .... dC.z/dpN]]
    M_dp = np.zeros((9, N_parameters))

    for i in range(N_parameters):
        M_dp[0:3, i] = [A[12 + i], A[13 + i], A[14 + i]]
        M_dp[3:6, i] = [B[12 + i], B[13 + i], B[14 + i]]
        M_dp[6:9, i] = [C[12 + i], C[13 + i], C[14 + i]]

    return M_dp


def parse_tri_file(filepath: str) -> tuple[list, list]:
    """Parses a .tri file."""

    def parse_rows(data: str, dtype: type):
        rows = [r.split() for r in data.strip().splitlines()]
        return [[dtype(c) for c in r] for r in rows]

    tree = ET.parse(filepath)
    root = tree.getroot()
    points = root[0][0][0]
    cells = root[0][0][1]

    vertices = parse_rows(data=points[0].text, dtype=float)
    triangles = parse_rows(data=cells[0].text, dtype=int)

    return vertices, triangles


# ------------------------------------------------------
# Main script
# ------------------------------------------------------

# The objective is to get the cells making up the surface
# of the geometry from a Cart3D simulation, including the
# flow properties on those cells.

# Hypervehicle allows us to obtain the vertex-parameter
# sensitivity information (dvdp). This is, how each vertex
# changes position with a change in the parameter.

# The cells are defined from vertices. So one cell needs to
# know which vertices make it up.

# The sensitivity files outputted by hypervehicle have the
# point coordinates, but the ordering in the file may not
# necessarily match the .tri file indexing.
# In the future, the .tri index could be added to this file.


# Configuration
verbosity = 2
directory = "/home/kieran/Documents/PoorMansAdjoint/coarse/testpv"
components_filepath = os.path.join(directory, "Components.i.tri")

# Filepaths
# sim_result_filepath = os.path.join(directory, "simdata.csv")
sensitivity_filepath = os.path.join(directory, "combined.csv")
celldata_filepath = os.path.join(directory, "cells.csv")
pointdata_filepath = os.path.join(directory, "points.csv")

# Define freestream flow proerties
dynamic_pressure_inf = 50e3  # kPa
rho_inf = 0.0308742  # kg/m3
a_inf = 299.499  # m/s
L_ref = 1  # m
Area_ref = 1  # m2


# Load data
# simdata = pd.read_csv(sim_result_filepath)
sensdata = pd.read_csv(sensitivity_filepath)
celldata = pd.read_csv(celldata_filepath)
pointdata = pd.read_csv(pointdata_filepath)

coordinates = pointdata[["Points_0", "Points_1", "Points_2"]]
cell_vertex_ids = celldata[["Point Index 0", "Point Index 1", "Point Index 2"]]

parameters = ["chord", "thickness", "wingspan"]

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

    dvdp = []
    for v in vertex_sensitivity_info:
        for p in parameters:
            dvdp.append([v[f"dxd_{p}"], v[f"dyd_{p}"], v[f"dzd_{p}"]])

    # Create Cell
    newcell = Cell.from_points(vertices)

    # Add geometry sensitivity information
    newcell._add_sensitivities(np.array(dvdp))

    # Add flow state at cell
    cell_velocity = Vector.from_coordinates(celldata.loc[cell][["U", "V", "W"]].values)
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
