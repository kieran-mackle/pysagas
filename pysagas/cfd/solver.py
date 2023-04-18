from abc import ABC, abstractmethod
from typing import Tuple, List, Optional
from pysagas import Cell, Vector, FlowState


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
