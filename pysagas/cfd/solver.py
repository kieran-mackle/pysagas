import os
import meshio
import numpy as np
from pysagas import Cell, FlowState
from abc import ABC, abstractmethod
from typing import List, Optional, Dict


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
        self, cells: List[Cell], freestream: Optional[FlowState] = None
    ) -> None:
        """Instantiates the flow solver."""

        # TODO - allow providing a geometry parser, eg. tri file parser for Cart3D,
        # or a (yet to be implemented) STL parser.

        self.cells = cells
        self.freestream = freestream
        self._last_solve_freestream: FlowState = None

    def solve(self, freestream: Optional[FlowState] = None):
        """Run the flow solver."""
        # TODO - what is the return of this method?

        if not freestream:
            # No freestream conditions specified
            if not self.freestream:
                # No conditions provided on instantiation either
                raise Exception("Please specify freestream conditions.")
            else:
                # Assign last solve freestream for solver to use
                self._last_solve_freestream = self.freestream
        else:
            # Solve-specific freestream conditions provided
            self._last_solve_freestream = freestream

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
