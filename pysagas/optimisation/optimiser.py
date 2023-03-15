from pysagas import banner
from abc import abstractmethod, ABC
from hypervehicle.generator import Generator
from pyoptsparse import Optimizer, Optimization


class ShapeOpt(ABC):
    """PySAGAS Shape Optimisation via a gradient descent algorithm."""

    def __init__(self, optimiser: Optimizer, generator: Generator) -> None:
        """Initialise PySAGAS Shape Optimiser.

        Parameters
        ----------
        optimiser : Optimizer
            The pyoptsparse Optimizer object of choice.
        """
        self.optimiser = optimiser
        self.generator = generator

        # Construct optimisation problem
        self.opt_problem = Optimization(
            name="PySAGAS Shape Optimisation",
            objFun=self.evaluate_objective,
            sens=self.evaluate_gradient,
        )

    def add_variables(self):
        """Adds variables to the optimisation problem."""
        pass

    def add_constriants(self):
        """Adds constraints to the optimisation problem."""
        pass

    def add_objective(self):
        """Adds an objective to the optimisation problem."""
        pass

    def run(self):
        # Print banner
        banner()
        print("\033[4mCart3D Shape Optimisation\033[0m".center(50, " "))

        # Run optimiser
        self.sol = self.optimiser(self.opt_problem, sens="FD")

    @abstractmethod
    def evaluate_objective(self, x: dict) -> dict:
        """Evaluates the objective function to be minimised at x."""
        pass

    @abstractmethod
    def evaluate_gradient(self, x: dict) -> dict:
        """Evaluates the Jacobian (objective gradient) at x."""
        pass
