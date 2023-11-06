import numpy as np
import pandas as pd
from pysagas.flow import FlowState
from abc import ABC, abstractmethod
from pysagas.geometry import Vector, Cell
from pysagas.cfd.solver import SensitivityResults
from typing import List, Callable, Tuple, Optional
from pysagas.utilities import all_dfdp, piston_dPdp, add_sens_data


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
    def _extract_parameters(self) -> List[str]:
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
        self.cells: list[Cell] = None
        self.verbosity = 1

    def calculate(
        self,
        dPdp_method: Optional[Callable] = piston_dPdp,
        cog: Optional[Vector] = Vector(0, 0, 0),
        flowstate: Optional[FlowState] = None,
        **kwargs,
    ) -> SensitivityResults:
        """Calculate the force sensitivities of the surface to the
        parameters.

        Parameters
        ----------
        dPdp_method : callable, optional
            The sensitivity method to use - either piston_dPdp or van_dyke_dPdp. The
            default is piston_dPdp.

        cog : Vector, optional
            The centre of gravity. The default is Vector(0, 0, 0).

        flowstate : FlowState, optional
            The flowstate associated with this result. The default is None.

        Returns
        -------
        SensitivityResults : the sensitivity results object.
        """
        if self.verbosity > 0:
            print("\nCalculating aerodynamic sensitivities.")

        parameters = self._extract_parameters()

        # if parameters != ["theta_2", "theta_1"]:
        #     print("swapped.")

        if self.cells is None:
            self.cells = self._transcribe_cells(parameters=parameters)

        params_sens_cols = []
        for p in parameters:
            for d in ["x", "y", "z"]:
                params_sens_cols.append(f"d{d}d_{p}")

        # Calculate force sensitivity
        F_sense, M_sense = all_dfdp(
            cells=self.cells, dPdp_method=dPdp_method, cog=cog, **kwargs
        )

        # hm = np.array([[ 3.59853397e+04,  1.85845803e+04,  2.15575270e+03],
        #     [ 5.79130106e+02,  1.81481290e+03, -1.05987335e-06]])

        # if not np.all(np.isclose(hm, F_sense)):
        #     print("\n\n\n\nbad")

        # if parameters[0] != "theta_2":
        #     print("\n\n\n\nwtf")

        # Construct dataframes to return
        df_f = pd.DataFrame(
            F_sense, columns=["dFx/dp", "dFy/dp", "dFz/dp"], index=parameters
        )
        df_m = pd.DataFrame(
            M_sense, columns=["dMx/dp", "dMy/dp", "dMz/dp"], index=parameters
        )

        if self.verbosity > 0:
            print("Done.")

        # Construct results
        result = SensitivityResults(
            f_sens=df_f,
            m_sens=df_m,
            freestream=flowstate,
        )

        return result

    def to_csv(self):
        """Dumps the sensitivity data to CSV file."""
        if self.cells is not None:
            pass

        else:
            raise Exception("No cells have been transcribed.")

    def _extract_parameters(self):
        parameters = []
        for e in self.sensdata.columns:
            e: str
            if e.startswith("dxd") or e.startswith("dyd") or e.startswith("dzd"):
                # Sensitivity coluns
                if e[3:] not in parameters:
                    parameters.append(e[3:])
        return parameters


class GenericWrapper(Wrapper):
    solver = "Generic Flow Solver"

    def __init__(
        self,
        cells: List[Cell],
        sensitivity_filepath: Optional[str] = None,
        cells_have_sens_data: Optional[bool] = False,
        verbosity: Optional[int] = 1,
        **kwargs,
    ) -> None:
        """Instantiate a generic PySAGAS sensitivity wrapper.

        Parameters
        ----------
        cells : list[Cell]
            A list of transcribed Cell objects, containing the nominal flow solution.

        sensitivity_filepath : str
            The filepath to the geometry sensitivities. This must be provided, if the
            cells do not have geometric sensitivity information attached. This must
            also be provided anyway, to determine the geometric design parameters.

        cells_have_sens_data : bool, optional
            The cells already have geometric sensitivity data matched to them. When this
            is True, the sensitivity_filepath argument (if provided) is ignored. The
            default is False.

        verbosity : int, optional
            The verbosity of the code. The defualt is 1.
        """
        super().__init__(**kwargs)

        if cells_have_sens_data:
            self.cells = cells
        else:
            self._pre_transcribed_cells = cells

        self.sensdata = pd.read_csv(sensitivity_filepath)
        self.verbosity = verbosity

    def _transcribe_cells(self, **kwargs) -> List[Cell]:
        """This is a dummy method to satisfy the abstract base class. Transcribed cells
        are provided upon instantiation of the wrapper.
        """
        # Add sensitivity data to _pre_transcribed_cells
        add_sens_data(
            cells=self._pre_transcribed_cells,
            data=self.sensdata,
            verbosity=self.verbosity,
        )

        return self._pre_transcribed_cells
