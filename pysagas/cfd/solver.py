from __future__ import annotations
import os
import copy
import meshio
import numpy as np
import pandas as pd
from abc import ABC, abstractmethod
from typing import List, Optional, Dict
from pysagas import Cell, FlowState, Vector


class AbstractFlowSolver(ABC):
    """Abstract interface for flow solvers."""

    method = None

    @abstractmethod
    def __init__(self, **kwargs) -> None:
        pass

    def __repr__(self) -> str:
        return f"PySAGAS {self.method} flow solver"

    def __str__(self) -> str:
        return f"PySAGAS {self.method} flow solver"

    @property
    @abstractmethod
    def method(self):
        # This is a placeholder for a class variable defining the solver method
        pass

    @abstractmethod
    def solve(self, *args, **kwargs):
        pass


class FlowSolver(AbstractFlowSolver):
    """Base class for CFD flow solver."""

    def __init__(
        self,
        cells: List[Cell],
        freestream: Optional[FlowState] = None,
        verbosity: Optional[int] = 1,
    ) -> None:
        """Instantiates the flow solver.

        Parameters
        ----------
        cells : list[Cell]
            A list of the cells describing the geometry.

        freestream : Flowstate, optional
            The free-stream flow state. The default is the freestream provided
            upon instantiation of the solver.

        verbosity : int, optional
            The verbosity of the solver. The default is 1.
        """

        # TODO - allow providing a geometry parser, eg. tri file parser for Cart3D,
        # or a (yet to be implemented) STL parser.

        self.cells = cells
        self.freestream = freestream
        self.verbosity = verbosity
        self._last_solve_freestream: FlowState = None
        self._last_sens_freestream: FlowState = None

        # Results
        self.flow_result: FlowResults = None
        self.flow_sensitivity: SensitivityResults = None

    def solve(
        self,
        freestream: Optional[FlowState] = None,
        Mach: Optional[float] = None,
        aoa: Optional[float] = None,
    ) -> bool:
        """Run the flow solver.

        Parameters
        ----------
        freestream : Flowstate, optional
            The free-stream flow state. The default is the freestream provided
            upon instantiation of the solver.

        Mach : float, optional
            The free-stream Mach number. The default is that specified in
            the freestream flow state.

        aoa : float, optional
            The free-stream angle of attack. The default is that specified in
            the freestream flow state.

        Returns
        -------
        result : FlowResults
            The flow results.

        Raises
        ------
        Exception : when no freestream can be found.
        """
        if not freestream:
            # No freestream conditions specified
            if not self.freestream:
                # No conditions provided on instantiation either
                raise Exception("Please specify freestream conditions.")
            else:
                # Use nominal freestream as base
                freestream = copy.copy(self.freestream)
                if Mach:
                    # Update mach number
                    freestream._M = Mach

                if aoa:
                    # Update flow direction
                    flow_direction = Vector(1, 1 * np.tan(np.deg2rad(aoa)), 0).unit
                    freestream.direction = flow_direction

        else:
            # Solve-specific freestream conditions provided
            if Mach or aoa:
                print("Using freestream provided; ignoring Mach/aoa provided.")

        # Check if already solved
        if self._last_solve_freestream and freestream == self._last_solve_freestream:
            return True
        else:
            # Update last solve freestream and continue
            self._last_solve_freestream = freestream
            return False

    def solve_sens(
        self,
        freestream: Optional[FlowState] = None,
        Mach: Optional[float] = None,
        aoa: Optional[float] = None,
    ) -> SensitivityResults:
        """Run the flow solver to obtain sensitivity information.

        Parameters
        ----------
        freestream : Flowstate, optional
            The free-stream flow state. The default is the freestream provided
            upon instantiation of the solver.

        Mach : float, optional
            The free-stream Mach number. The default is that specified in
            the freestream flow state.

        aoa : float, optional
            The free-stream angle of attack. The default is that specified in
            the freestream flow state.

        Returns
        -------
        SensitivityResults

        Raises
        ------
        Exception : when no freestream can be found.
        """
        if not freestream:
            # No freestream conditions specified
            if not self.freestream:
                # No conditions provided on instantiation either
                raise Exception("Please specify freestream conditions.")
            else:
                # Use nominal freestream as base
                freestream = copy.copy(self.freestream)
                if Mach:
                    # Update mach number
                    freestream._M = Mach

                if aoa:
                    # Update flow direction
                    flow_direction = Vector(1, 1 * np.tan(np.deg2rad(aoa)), 0).unit
                    freestream.direction = flow_direction

        else:
            # Solve-specific freestream conditions provided
            if Mach or aoa:
                print("Using freestream provided; ignoring Mach/aoa provided.")

        # Check if already solved
        if self._last_sens_freestream and freestream == self._last_sens_freestream:
            return True
        else:
            # Update last solve freestream and continue
            self._last_sens_freestream = freestream
            return False

    def save(self, name: str, attributes: Dict[str, list]):
        """Save the solution to VTK file format. Note that the
        geometry mesh must have been loaded from a parser which
        includes connectivity information (eg. PyMesh or TRI).

        Parameters
        ----------
        name : str
            The filename prefix of the desired output file.

        attributes : Dict[str, list]
            The attributes dictionary used to initialise the file writer.
            Note that this is not required when calling `solve` from a
            child FlowSolver.
        """
        # Import PyMesh
        try:
            import pymesh
        except ModuleNotFoundError:
            raise Exception(
                "Could not find pymesh. Please follow the "
                + "installation instructions at "
                + "https://pymesh.readthedocs.io/en/latest/installation.html"
            )

        # Construct vertices and faces
        faces = []
        vertex_faces = {}
        for cell in self.cells:
            # Get faces
            faces.append(cell._face_ids)

            # Get attributes
            for a in attributes:
                attributes[a].append(cell.attributes.get(a))

            # Get vertices
            for i, vid in enumerate(cell._face_ids):
                vertex_faces[vid] = cell.vertices[i]

        # Sort vertices
        sorted_vertices = dict(sorted(vertex_faces.items()))

        # Note: vertex attributes can be treated the same way as
        # the vertices - they must be sorted.

        # Finally
        vertices = np.array([v for v in sorted_vertices.values()])
        faces = np.array(faces)

        # Convert back to mesh object
        mesh_obj = pymesh.form_mesh(vertices, faces)

        # Try add attributes to the cells
        for a, v in attributes.items():
            mesh_obj.add_attribute(a)
            mesh_obj.set_attribute(a, np.array(v))

        # Save mesh
        pymesh.save_mesh(
            f"{name}.ply", mesh_obj, *mesh_obj.get_attribute_names(), ascii=True
        )

        # Convert to VTK format
        m = meshio.read(f"{name}.ply")
        m.write(f"{name}.vtk")

        # Remove .ply file
        os.remove(f"{name}.ply")

    @staticmethod
    def body_to_wind(v: Vector, aoa: float):
        """Converts a vector from body axes to wind axes.

        Parameters
        ----------
        v : Vector
            The result to transform.

        aoa : float
            The angle of attack, specified in degrees.

        Returns
        --------
        r : Vector
            The transformed result.
        """
        # TODO - use euler matrices
        A = v.x  # axial
        N = v.y  # normal

        L = N * np.cos(np.deg2rad(aoa)) - A * np.sin(np.deg2rad(aoa))
        D = N * np.sin(np.deg2rad(aoa)) + A * np.cos(np.deg2rad(aoa))

        # Construct new result
        r = Vector(x=D, y=L, z=v.z)

        return r


class FlowResults:
    def __init__(
        self, freestream: FlowState, net_force: Vector, net_moment: Vector
    ) -> None:
        self.freestream = freestream
        self.net_force = net_force
        self.net_moment = net_moment

        # Calculate angle of attack
        self.aoa = np.rad2deg(
            np.arctan(freestream.direction.y / freestream.direction.x)
        )

    def __str__(self) -> str:
        return f"Net force = {self.net_force} N"

    def __repr__(self) -> str:
        return self.__str__()

    def coefficients(self, A_ref: Optional[float] = 1.0, c_ref: Optional[float] = 1.0):
        """Calculate the aerodynamic coefficients CL, CD and Cm.

        Parameters
        -----------
        A_ref : float, optional
            The reference area (m^2). The default is 1 m^2.

        c_ref : float, optional
            The reference length (m). The default is 1 m.
        """
        w = FlowSolver.body_to_wind(v=self.net_force, aoa=self.aoa)
        C_L = w.y / (self.freestream.q * A_ref)
        C_D = w.x / (self.freestream.q * A_ref)

        mw = FlowSolver.body_to_wind(v=self.net_moment, aoa=self.aoa)
        C_m = mw.z / (self.freestream.q * A_ref * c_ref)
        return C_L, C_D, C_m


class SensitivityResults:
    def __init__(self, f_sens: pd.DataFrame, m_sens: pd.DataFrame) -> None:
        self.f_sens = f_sens
        self.m_sens = m_sens

    def __str__(self) -> str:
        s1 = f"Force sensitivties:\n{self.f_sens}"
        s2 = f"Moment sensitivties:\n{self.m_sens}"

        return f"{s1}\n\n{s2}"

    def __repr__(self) -> str:
        return self.__str__()
