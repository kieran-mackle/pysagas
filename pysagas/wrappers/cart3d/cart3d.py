import sys
import numpy as np
import pandas as pd
from typing import Union, List
from pysagas.flow import FlowState
from pysagas.geometry import Vector, Cell
from pysagas.wrappers.wrapper import Wrapper
from pysagas.wrappers.cart3d.utilities import process_components_file


class Cart3DWrapper(Wrapper):
    """PySAGAS Cart3D wrapper."""

    solver = "Cart3D"

    def __init__(
        self,
        sensitivity_filepath: str,
        celldata_filepath: str = None,
        pointdata_filepath: str = None,
        components_filepath: str = None,
        **kwargs,
    ) -> None:

        if celldata_filepath is None and pointdata_filepath is None:
            # Generate cell and point data from Components.i.plt
            try:
                pointdata_filepath, celldata_filepath = process_components_file(
                    components_filepath
                )
            except FileNotFoundError as e:
                if components_filepath is not None:
                    # A filepath was provided, but not valid
                    print(
                        f"Cannot find Components.i.plt file at {components_filepath}."
                    )
                else:
                    # No filepath was provided, and could not find the default
                    print(
                        "Please provide the filepaths to the "
                        + "cell data file and points data file. Alternatively, "
                        + "provide the filepath to the Components.i.plt file."
                    )
                sys.exit()

        # Load data
        self.sensdata = pd.read_csv(sensitivity_filepath)
        self.celldata = pd.read_csv(celldata_filepath)
        self.pointdata = pd.read_csv(pointdata_filepath)

    def _create_cells(self, parameters: list):
        coordinates = self.pointdata[["Points_0", "Points_1", "Points_2"]]
        cell_vertex_ids = self.celldata[
            ["Point Index 0", "Point Index 1", "Point Index 2"]
        ]

        # Construct cells
        cells = []
        for cell in self.celldata.index:
            vertex_ids = cell_vertex_ids.loc[cell].values
            cell_v_coords = [coordinates.loc[v_id].values for v_id in vertex_ids]
            vertices = [Vector.from_coordinates(c) for c in cell_v_coords]

            # Extract sensitivity information for the cell
            vertex_sensitivity_info = [self.sensdata.loc[v_id] for v_id in vertex_ids]

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
                self.celldata.loc[cell][["U", "V", "W"]].values
            )
            cell_flowstate = FlowState(
                mach=self.celldata.loc[cell]["M"],
                pressure=self.celldata.loc[cell]["p"],
                temperature=self.celldata.loc[cell]["T"],
                direction=cell_velocity,
            )
            newcell.flowstate = cell_flowstate

            # Append new Cell
            cells.append(newcell)

        return cells

    def _extract_parameters(self):
        parameters = set()
        for e in self.sensdata.columns:
            e: str
            if e.startswith("dxd") or e.startswith("dyd") or e.startswith("dzd"):
                # Sensitivity coluns
                parameters.add("".join(e.split("d")[2:]))
        return list(parameters)
