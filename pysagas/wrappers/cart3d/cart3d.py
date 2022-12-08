import os as os
import numpy as np
import pandas as pd
from pysagas.flow import FlowState
from pysagas.utilities import all_dfdp
from pysagas.geometry import Vector, Cell
from pysagas.wrappers.wrapper import AbstractWrapper
from pysagas.wrappers.cart3d.utilities import process_components_file


class Cart3DWrapper(AbstractWrapper):
    """PySAGAS Cart3D wrapper."""

    def __init__(
        self,
        sensitivity_filepath: str,
        celldata_filepath: str = None,
        pointdata_filepath: str = None,
        components_filepath: str = None,
        **kwargs,
    ) -> None:
        self.sensitivity_filepath = sensitivity_filepath

        if celldata_filepath is None and pointdata_filepath is None:
            if components_filepath is not None:
                # Generate cell and point data from Components.i.plt
                process_components_file()
            else:
                raise Exception(
                    "Please provide the filepaths to the "
                    + "cell data file and points data file. Alternatively, "
                    + "provide the filepath to the Components.i.plt file."
                )

        parameters = ["P"]  # TODO - automate extraction of parameters

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
