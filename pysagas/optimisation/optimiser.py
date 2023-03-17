import os
from pysagas import banner
from numpy.typing import ArrayLike
from typing import Union, Optional, Dict
from pyoptsparse import Optimizer, Optimization


class ShapeOpt:
    """PySAGAS Shape Optimisation via a gradient descent algorithm."""

    def __init__(
        self,
        optimiser: Optimizer,
        working_dir: str,
        optimiser_options: dict = None,
    ) -> None:
        """Initialise PySAGAS Shape Optimiser.

        Parameters
        ----------
        optimiser : Optimizer
            The pyoptsparse Optimizer object of choice.
        """
        self.opt_problem: Optimization
        self.optimiser = optimiser(optimiser_options)

        # Prepare working directory
        if not os.path.exists(working_dir):
            os.mkdir(working_dir)

    def add_variables(
        self,
        parameters_dict: Optional[dict] = None,
        name: Optional[str] = None,
        n: Optional[int] = None,
        vartype: Optional[str] = "c",
        initial: Optional[ArrayLike] = None,
        lower_bounds: Optional[Union[ArrayLike, Dict]] = None,
        upper_bounds: Optional[Union[ArrayLike, Dict]] = None,
        **kwargs,
    ):
        """Adds variables to the optimisation problem."""
        if parameters_dict:
            # Check bounds
            if lower_bounds is None:
                lower_bounds = {}
            if upper_bounds is None:
                upper_bounds = {}

            # Unpack parameters dict
            for param, value in parameters_dict.items():
                self.opt_problem.addVarGroup(
                    name=param,
                    nVars=1,
                    varType=vartype,
                    value=value,
                    lower=lower_bounds.get(param),
                    upper=upper_bounds.get(param),
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

    def add_constriants(
        self,
        name: str,
        n: Optional[int] = 1,
        lower: Optional[Union[ArrayLike, Dict]] = None,
        upper: Optional[Union[ArrayLike, Dict]] = None,
        scale: Optional[float] = 1,
    ):
        """Adds constraints to the optimisation problem."""
        self.opt_problem.addConGroup(
            name=name,
            nCon=n,
            lower=lower,
            upper=upper,
            scale=scale,
            linear=None,
            wrt=None,
            jac=None,
        )

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
