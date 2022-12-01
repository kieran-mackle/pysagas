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

    # Append new Cell
    cells.append(newcell)


# # Calculate force sensitivity
# F_sense = all_dfdp(
#     cells=cells, P=1000, rho=rho_inf, a=a_inf, vel_vector=vel_ramp_vector
# )


# # Load and process Components.i.tri geometry file
# vertices, triangles = parse_tri_file(components_filepath)
# vertices_df = pd.DataFrame(vertices, columns=["x", "y", "z"])
# triangles_df = pd.DataFrame(triangles, columns=["p0", "p1", "p2"])

# # Load simdata.csv and combined.csv
# sim_result_df = pd.read_csv(sim_result_filepath)
# sensitivity_df = pd.read_csv(sensitivity_filepath)

# # Construct vertices and cells
# vector_vertices = [Vector.from_coordinates(c) for c in vertices]
# cells = [Cell.from_points([vector_vertices[pi] for pi in t]) for t in triangles]


# # Add vertex-parameter sensitivities


# # print(sim_result_df)  # [x,y,z, U, V, W, p, rho, T, a, V_mag, Mach]
# # [index, x, y, z, dx/dp0, dy/dp0, dz/dp0, dx/dp1, dy/dp1, dz/dp1, ...]

# # convert to numpy array so that we have matching indexing.
# points = points_df.to_numpy()
# connect = connectivity_df.to_numpy()
# sim_result = sim_result_df.to_numpy()
# sensitivity = sensitivity_df.to_numpy()

# N_points = np.shape(points)[0]
# print(f"Number of Points in points list: {N_points}")

# # check that points matches combined_sensitivity.csv and simdata.csv
# # extra points in combined_sensitivities will be dropped.
# data = np.zeros((N_points, np.shape(sim_result)[1] + np.shape(sensitivity)[1] - 4))
# # [x,y,z, U, V, W, p, rho, T, a, V_mag, Mach, dx/dp0, dy/dp0, dz/dp0, dx/dp1, dy/dp1, dz/dp1, ...]

# i_points = 0
# i_sense = 0
# tolerance = 1e-5
# while i_points < N_points:
#     # print(i_sense, i_points)
#     xp, yp, zp = points[i_points, :]
#     xr, yr, zr = sim_result[i_points, 0:3]
#     xs, ys, zs = sensitivity[i_sense, 1:4]
#     # Note there is an off-set in index.

#     # print(xp, yp, zp)
#     flag = 0
#     if abs(xp - xr) > tolerance or (xp - xs) > tolerance:
#         if verbosity > 1:
#             print(f"Sense_entry {i_sense}: xp={xp}, xr={xr}, xs={xs}")
#         flag = flag + 1
#         # raise Exception(f"X-coordinates for point {i} don't match")
#     if abs(yp - yr) > tolerance or (yp - ys) > tolerance:
#         if verbosity > 1:
#             print(f"Sense_entry {i_sense}: yp={yp}, yr={yr}, ys={ys}")
#         flag = flag + 1
#         # raise Exception(f"Y-coordinates for point {i} don't match")
#     if abs(zp - zr) > tolerance or (zp - zs) > tolerance:
#         if verbosity > 1:
#             print(f"Sense_entry {i_sense}: zp={zp}, zr={zr}, zs={zs}")
#         flag = flag + 1
#         # raise Exception(f"Z-coordinates for point {i} don't match")

#     if flag == 0:
#         # points match.
#         data[i_points, 0 : np.shape(sim_result)[1]] = sim_result[i_points, :]
#         data[i_points, np.shape(sim_result)[1] :] = sensitivity[i_sense, 4:]
#         i_points = i_points + 1
#         i_sense = i_sense + 1
#     else:
#         # no mathcin point
#         i_sense = i_sense + 1
# print(
#     f"Number of lines in data = {i_points}, from possible {i_sense} available in sensitivity.csv"
# )


# print("\nstarting to evaluate mean properties for each triangle")
# N_triangles = np.shape(connect)[0]
# N_parameters = int((np.shape(sensitivity)[1] - 4) / 3)
# # N_parameters = 2
# print(f"Number of triangles = {N_triangles}")
# print(f"Number of design parameter = {N_parameters}")
# flow = np.zeros((N_triangles, 9 + 1))
# # [t_index, U, V, W, pressure, density, T, speed of sound, V_mag, Mach]
# n_sense = np.zeros((N_triangles, N_parameters * 3 + 1))
# # [t_index, n.x/dp0, n.y/dp0, n.z/dp0, ...]
# fx_list = np.zeros((N_triangles, N_parameters + 1))  # [t_index, f_x/dp, f_x/dp1, ...]
# fy_list = np.zeros((N_triangles, N_parameters + 1))  # [t_index, f_x/dp, f_x/dp1, ...]
# fz_list = np.zeros((N_triangles, N_parameters + 1))  # [t_index, f_x/dp, f_x/dp1, ...]


# for i in range(N_triangles):
#     a, b, c = connect[i, :]
#     A = data[int(a), :]
#     B = data[int(b), :]
#     C = data[int(c), :]

#     # Calculate mean properties for each triangle
#     u_mean = 1 / 3 * (A[3] + B[3] + C[3])
#     v_mean = 1 / 3 * (A[4] + B[4] + C[4])
#     w_mean = 1 / 3 * (A[5] + B[5] + C[5])
#     p_mean = 1 / 3 * (A[6] + B[6] + C[6])
#     rho_mean = 1 / 3 * (A[7] + B[7] + C[7])
#     T_mean = 1 / 3 * (A[8] + B[8] + C[8])
#     a_mean = 1 / 3 * (A[9] + B[9] + C[9])
#     V_mag_mean = 1 / 3 * (A[10] + B[10] + C[10])
#     M_mean = 1 / 3 * (A[11] + B[11] + C[11])

#     flow[i, :] = [
#         i,
#         u_mean,
#         v_mean,
#         w_mean,
#         p_mean,
#         rho_mean,
#         T_mean,
#         a_mean,
#         V_mag_mean,
#         M_mean,
#     ]

#     # evaluate sensitivities
#     # create matrix linling n_x, n_y, n_z to vertex positons
#     a_vector = Vector(A[0], A[1], A[2])
#     b_vector = Vector(B[0], B[1], B[2])
#     c_vector = Vector(C[0], C[1], C[2])

#     n0, M_sense = create_sensitivity_matrix(a_vector, b_vector, c_vector)

#     # create matrix linking vtx positons to parameters
#     M_dp = create_positon_sense_matrix(A, B, C, N_parameters)

#     # create matrix linking n_x, n_y, n_z to parameters
#     N_sense = np.dot(M_sense, M_dp)
#     if verbosity > 1:
#         print(f"M_sense={M_sense}")
#         print(f"M_dp={M_dp}")
#         print(f"N_sense = {N_sense}")

#     # evaluate dP/dn
#     vel_vector = flow[i, 1:4]
#     pressure_sense = flow[i, 5] * flow[i, 7] * np.dot(vel_vector, -N_sense)
#     if verbosity > 1:
#         print(f"A: x={A[0]}, y={A[1]}, z={A[2]}")
#         print(f"B: x={B[0]}, y={B[1]}, z={B[2]}")
#         print(f"C: x={C[0]}, y={C[1]}, z={C[2]}")
#         print(f"vel_vector = {vel_vector}")
#         print(f"pressure_sense = {pressure_sense}")

#     # evaluate force contribution

#     Area, Area_sense = create_area_sensitivity_matrix(a_vector, b_vector, c_vector)
#     A_sense = np.dot(Area_sense, M_dp)
#     if verbosity > 1:
#         print(f"Area = {Area}")
#         print(f"Area_sense = {Area_sense}")
#         print(f"A_sense = {A_sense}")  ## NEED TO CHECK THIS SEEMS WRONG

#     pressure = flow[i, 4]

#     def calculate_force_sensitivity(n_direction=[1.0, 0.0, 0.0], verbosity=0):
#         sense = (
#             pressure_sense * Area * np.dot(n0, n_direction)
#             + pressure * A_sense * np.dot(n0, n_direction)
#             + pressure * Area * np.dot(-N_sense.T, n_direction)
#         )
#         # sense = (pressure_sense*Area*np.dot(-n0,n_direction)
#         #        + pressure*A_sense*np.dot(-n0,n_direction)
#         #        + pressure*Area*np.dot(-pressure_sense.T,n_direction))
#         return sense

#     Fx_sense = calculate_force_sensitivity(n_direction=[1.0, 0.0, 0.0])
#     Fy_sense = calculate_force_sensitivity(n_direction=[0.0, 1.0, 0.0])
#     Fz_sense = calculate_force_sensitivity(n_direction=[0.0, 0.0, 1.0])

#     if verbosity > 1:
#         print("Force Sensitivities from Local Piston Theory")
#         print(f"Fx_sense={Fx_sense}")
#         print(f"Fy_sense={Fy_sense}")
#         print(f"Fz_sense={Fz_sense}")

#     fx_list[i, 0] = i
#     fx_list[i, 1:] = Fx_sense
#     fy_list[i, 0] = i
#     fy_list[i, 1:] = Fy_sense
#     fz_list[i, 0] = i
#     fz_list[i, 1:] = Fz_sense

# print("")
# print("RESULTS")
# print(f"fx_list={fx_list}")
# print(f"fy_list={fy_list}")
# print(f"fz_list={fz_list}")
# print("")
# Fx = np.sum(fx_list[:, 1:], axis=0) / (dynamic_pressure_inf * Area_ref)
# Fy = np.sum(fy_list[:, 1:], axis=0) / (dynamic_pressure_inf * Area_ref)
# Fz = np.sum(fz_list[:, 1:], axis=0) / (dynamic_pressure_inf * Area_ref)

# for i, (x, y, z) in enumerate(zip(Fx, Fy, Fz)):
#     print(f"dFx/dp{i} = {x}; dFy/dp{i} = {y}; dFz/dp{i} = {z}")
