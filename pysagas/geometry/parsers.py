from stl import mesh
from tqdm import tqdm
from typing import List
import xml.etree.ElementTree as ET
from abc import ABC, abstractmethod
from pysagas.geometry import Cell, Vector


class AbstractParser(ABC):
    """Interface for a geometry parser."""

    filetype = None

    @abstractmethod
    def __init__(self, **kwargs) -> None:
        pass

    def __repr__(self) -> str:
        return f"PySAGAS {self.filetype} parser"

    def __str__(self) -> str:
        return f"PySAGAS {self.filetype} parser"

    @property
    @abstractmethod
    def filetype(self):
        # This is a placeholder for a class variable defining the parser file type
        pass

    @abstractmethod
    def load(self) -> List[Cell]:
        """Parses the file."""


class Parser(AbstractParser):
    def __init__(self, filepath: str, verbosity: int = 1) -> None:
        self.filepath = filepath
        self.verbosity = verbosity


class STL(Parser):
    filetype = "STL"

    def load(self) -> List[Cell]:
        # Load the STL
        mesh_obj = mesh.Mesh.from_file(self.filepath)

        cells = []
        if self.verbosity > 0:
            print("Transcribing cells:")
            pbar = tqdm(
                total=len(mesh_obj.vectors),
                position=0,
                leave=True,
                desc="  Cell transcription progress",
            )

        for vector_triple in mesh_obj.vectors:
            vertices = [Vector.from_coordinates(v) for v in vector_triple]
            try:
                cell = Cell.from_points(vertices)
                cells.append(cell)
            except:
                pass

            # Update progress bar
            if self.verbosity > 0:
                pbar.update(1)

        if self.verbosity > 0:
            pbar.close()
            print("Done.")

        return cells


class TRI(Parser):
    filetype = ".tri"

    def load(self) -> List[Cell]:
        # Parse .tri file
        tree = ET.parse(self.filepath)
        root = tree.getroot()
        grid = root[0]
        piece = grid[0]
        points = piece[0]
        cells = piece[1]

        points_data = points[0].text
        cells_data = cells[0].text

        points_data_list = [el.split() for el in points_data.splitlines()[1:]]
        points_data_list = [[float(j) for j in i] for i in points_data_list]

        cells_data_list = [el.split() for el in cells_data.splitlines()[1:]]
        cells_data_list = [[int(j) for j in i] for i in cells_data_list]

        cells = []
        if self.verbosity > 0:
            print("Transcribing cells:")
            pbar = tqdm(
                total=len(cells_data_list),
                position=0,
                leave=True,
                desc="  Cell transcription progress",
            )
        for vertex_idxs in cells_data_list:
            vertices = [
                Vector.from_coordinates(points_data_list[i]) for i in vertex_idxs
            ]
            cell = Cell.from_points(vertices)
            cells.append(cell)

            # Update progress bar
            if self.verbosity > 0:
                pbar.update(1)

        if self.verbosity > 0:
            pbar.close()
            print("Done.")

        return cells
