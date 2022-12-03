import os as os
import numpy as np
import pandas as pd
from adjoint.flow import FlowState
from adjoint.utilities import all_dfdp
from adjoint.geometry import Vector, Cell


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
directory = "/mnt/c/Users/kemac/Documents/pysagas"
parameters = ["P"]  # TODO - automate extraction of parameters

# Filepaths
sensitivity_filepath = os.path.join(directory, "combined.csv")
celldata_filepath = os.path.join(directory, "cells.csv")
pointdata_filepath = os.path.join(directory, "points.csv")

L_ref = 1  # m
Area_ref = 1  # m2

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
