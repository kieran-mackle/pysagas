import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import Union, List
from pysagas.flow import FlowState
from pysagas.wrappers.wrapper import Wrapper
from pysagas.geometry import Vector, Cell, DegenerateCell
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
        pointdata: pd.DataFrame = None,
        celldata: pd.DataFrame = None,
        write_data: bool = False,
        verbosity: int = 1,
        **kwargs,
    ) -> None:
        """A PySAGAS wrapper for Cart3D.

        Parameters
        ----------
        a_inf : float
            The freestream speed of sound (m/s).
        rho_inf : float
            The freestream density (kg/m^3).
        sensitivity_filepath : str
            The filepath to the geometry sensitivities.
        components_filepath : str, optional
            The filepath to the Components.i.plt file to be processed.
            The default is None.
        pointdata : pd.DataFrame, optional
            The point data. Must be supplied with celldata.
        celldata : pd.DataFrame, optional
            The cell data. Must be supplied with pointdata.
        write_data : bool, optional
            Write the flow data to CSV files. The default is True.
        verbosity : int, optional
            The verbosity of the code. The defualt is 1.
        """
        # Load data
        self.sensdata = pd.read_csv(sensitivity_filepath)
        if pointdata is not None and celldata is not None:
            self.pointdata = pointdata
            self.celldata = celldata
        else:
            self.pointdata, self.celldata = process_components_file(
                a_inf,
                rho_inf,
                components_filepath,
                write_data=write_data,
                verbosity=verbosity,
            )
        self.verbosity = verbosity

        # Check dimensionality
        if len(self.sensdata) != len(self.pointdata):
            raise ValueError("The sensitivity data does not match the point data.")

        super().__init__(**kwargs)

    def _transcribe_cells(self, parameters: List[str]) -> List[Cell]:
        """Transcribes the cells from Components.i.plt files into
        PySAGAS Cell objects.

        Parameters
        -----------
        parameters : List[str]
            A list of the geometric design parameters.

        Returns
        --------
        cells : List[Cell]
            A list of all transcribed cells.

        See Also
        --------
        pysagas.geometry.Cell
        """
        coordinates = self.pointdata[["Points_0", "Points_1", "Points_2"]]
        cell_vertex_ids = self.celldata[
            ["Point Index 0", "Point Index 1", "Point Index 2"]
        ]

        # Construct cells
        if self.verbosity > 0:
            print("\nTranscribing cells from Components.i.plt...")
            pbar = tqdm(
                total=len(self.celldata.index),
                position=0,
                leave=True,
            )
        cells = []
        degen_cells = 0
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
            try:
                newcell = Cell.from_points(vertices)
            except DegenerateCell:
                # Bad cell, skip it
                degen_cells += 1
                if self.verbosity > 0:
                    pbar.update(1)

                    if self.verbosity > 2:
                        print("\033[1mWarning\033[0m: Degenerate cell.")
                continue

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
            if self.verbosity > 0:
                pbar.update(1)

        if self.verbosity > 0:
            pbar.close()
            print("Done.")

            if self.verbosity > 1:
                if degen_cells > 0:
                    print(f"{100*degen_cells/len(self.celldata):.2f}% degenerate cells")

        return cells

    def _extract_parameters(self):
        parameters = set()
        for e in self.sensdata.columns:
            e: str
            if e.startswith("dxd") or e.startswith("dyd") or e.startswith("dzd"):
                # Sensitivity coluns
                parameters.add(e[3:])
        return list(parameters)
