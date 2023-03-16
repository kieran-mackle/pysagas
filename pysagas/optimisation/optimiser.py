import os
from pysagas import banner
from numpy.typing import ArrayLike
from abc import abstractmethod, ABC
from hypervehicle.generator import Generator
from pyoptsparse import Optimizer, Optimization


class ShapeOpt(ABC):
    """PySAGAS Shape Optimisation via a gradient descent algorithm."""

    def __init__(
        self,
        optimiser: Optimizer,
        generator: Generator,
        working_dir: str,
    ) -> None:
        """Initialise PySAGAS Shape Optimiser.

        Parameters
        ----------
        optimiser : Optimizer
            The pyoptsparse Optimizer object of choice.
        """
        # TODO - allow optimiser options
        self.optimiser = optimiser()
        self.generator = generator

        # Prepare working directory
        if not os.path.exists(working_dir):
            os.mkdir(working_dir)

    def add_variables(
        self,
        parameters_dict: dict = None,
        name: str = None,
        n: int = None,
        vartype: str = "c",
        initial: ArrayLike = None,
        lower_bounds: ArrayLike = None,
        upper_bounds: ArrayLike = None,
        **kwargs,
    ):
        """Adds variables to the optimisation problem."""
        if parameters_dict:
            # Unpack parameters dict
            for param, value in parameters_dict.items():
                # TODO - handle bounds
                self.opt_problem.addVarGroup(
                    name=param,
                    nVars=1,
                    varType=vartype,
                    value=value,
                    lower=lower_bounds,
                    upper=upper_bounds,
                    **kwargs,
                )

        else:
            self.opt_problem.addVarGroup(
                name=name,
                nVars=n,
                varType=vartype,
                value=initial,
                lower=lower_bounds,
                upper=upper_bounds,
                **kwargs,
            )

    def add_constriants(self):
        """Adds constraints to the optimisation problem."""
        pass

    def add_objective(self):
        """Adds an objective to the optimisation problem."""
        self.opt_problem.addObj("objective")

    def run(self):
        # Print banner
        banner()
        print("\033[4mPySAGAS Shape Optimisation\033[0m".center(50, " "))

        # Run optimiser
        self.sol = self.optimiser(
            self.opt_problem,
            storeHistory="history.hst",
        )

    # @abstractmethod
    # def evaluate_objective(self, x: dict) -> dict:
    #     """Evaluates the objective function to be minimised at x."""
    #     pass

    # @abstractmethod
    # def evaluate_gradient(self, x: dict, objective: dict) -> dict:
    #     """Evaluates the Jacobian (objective gradient) at x."""
    #     pass


def _unwrap_x(x: dict) -> dict:
    """Unwraps an ordered dictionary."""
    unwrapped = {}
    for key, val in x.items():
        if len(val) == 1:
            unwrapped[key] = val[0]
        else:
            unwrapped[key] = val
    return unwrapped
