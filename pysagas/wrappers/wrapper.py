import numpy as np
import pandas as pd
from abc import ABC, abstractmethod
from typing import List, Callable, Tuple
from pysagas.utilities import all_dfdp, panel_dPdp


class AbstractWrapper(ABC):
    """Abstract Wrapper class defining the wrapper interface."""

    solver = None

    @abstractmethod
    def __init__(self, **kwargs) -> None:
        pass

    def __repr__(self) -> str:
        return f"PySAGAS Solver Wrapper for {self.solver}"

    def __str__(self) -> str:
        return f"PySAGAS Solver Wrapper for {self.solver}"

    @property
    @abstractmethod
    def solver(self):
        # This is a placeholder for a class variable defining the wrapper's solver
        pass

    @abstractmethod
    def _transcribe_cells(self, parameters: List):
        """Transcribes the unified Cells from the solver data.

        Returns
        --------
        cells : List[Cell]
            A list of the Cells from the flow solution.
        """
        # This method should be overloaded with the solver-specific method
        pass

    @abstractmethod
    def _extract_parameters(self):
        """Extract the geometry parameters.

        Returns
        -------
        parameters : List[str]
            A list of the parameter names.
        """
        # This method should be overloaded with the solver-specific method
        pass

    @abstractmethod
    def to_csv(self):
        """Dumps the sensitivity data to CSV file."""
        pass


class Wrapper(AbstractWrapper):
    """Wrapper base class."""

    def __init__(self, **kwargs) -> None:
        # Cell data
        self.cells = None

    def calculate(
        self, dPdp_method: Callable = panel_dPdp, **kwargs
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Calculate the force sensitivities of the surface to the
        parameters.
        """
        parameters = self._extract_parameters()
        if self.cells is None:
            self.cells = self._transcribe_cells(parameters)

        params_sens_cols = []
        for p in parameters:
            for d in ["x", "y", "z"]:
                params_sens_cols.append(f"d{d}d_{p}")

        # Calculate force sensitivity
        F_sense, M_sense = all_dfdp(cells=self.cells, dPdp_method=dPdp_method, **kwargs)

        # Construct dataframes to return
        df_f = pd.DataFrame(
            F_sense, columns=["dFx/dp", "dFy/dp", "dFz/dp"], index=parameters
        )
        df_m = pd.DataFrame(
            M_sense, columns=["dMx/dp", "dMy/dp", "dMz/dp"], index=parameters
        )

        return df_f, df_m

    def to_csv(self):
        """Dumps the sensitivity data to CSV file."""
        if self.cells is not None:
            pass

        else:
            raise Exception("No cells have been transcribed.")
