from pysagas import banner
from pyoptsparse import Optimizer
from abc import abstractmethod, ABC


class ShapeOpt(ABC):
    """PySAGAS Shape Optimisation via a gradient descent algorithm."""

    def __init__(self, optimiser: Optimizer) -> None:
        """Initialise PySAGAS Shape Optimiser.

        Parameters
        ----------
        optimiser : Optimizer
            The pyoptsparse Optimizer object of choice.
        """
        self.optimiser = optimiser

    def run(self):
        # Print banner
        banner()
        print("\033[4mCart3D Shape Optimisation\033[0m".center(50, " "))

        # Construct optimisation problem
        opt_problem = None

        # Run optimiser
        sol = self.optimiser(opt_problem, sens="FD")

    @abstractmethod
    def evaluate_objective(self, x: dict) -> dict:
        """Evaluates the objective function to be minimised at x."""
        pass

    @abstractmethod
    def evaluate_gradient(self, x: dict) -> dict:
        """Evaluates the Jacobian (objective gradient) at x."""
        pass
