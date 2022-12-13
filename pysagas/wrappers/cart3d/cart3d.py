import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
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
        a_inf: float,
        rho_inf: float,
        sensitivity_filepath: str,
        components_filepath: str = None,
        write_data: bool = False,
        **kwargs,
    ) -> None:
        # Load data
        self.sensdata = pd.read_csv(sensitivity_filepath)
        self.pointdata, self.celldata = process_components_file(
            a_inf, rho_inf, components_filepath, write_data=write_data
        )

    def _create_cells(self, parameters: list):
        coordinates = self.pointdata[["Points_0", "Points_1", "Points_2"]]
        cell_vertex_ids = self.celldata[
            ["Point Index 0", "Point Index 1", "Point Index 2"]
        ]

        # Construct cells
        print("Transcribing cells from Components.i.plt...")
        pbar = tqdm(
            total=len(self.celldata.index),
            position=0,
            leave=True,
        )
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

            # Update progress bar
            pbar.update(1)
        pbar.close()
        print("Done.")

        return cells

    def _extract_parameters(self):
        parameters = set()
        for e in self.sensdata.columns:
            e: str
            if e.startswith("dxd") or e.startswith("dyd") or e.startswith("dzd"):
                # Sensitivity coluns
                parameters.add("".join(e.split("d")[2:]))
        return list(parameters)
