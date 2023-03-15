import os
import numpy as np
import pandas as pd
from pysagas import banner
from abc import abstractmethod
import matplotlib.pyplot as plt


class ShapeOpt:
    """PySAGAS Shape Optimisation via a gradient descent algorithm."""

    def __init__(self) -> None:
        pass

    def run(self, convergence_tolerance: float = 1e-6, max_iterations: int = 10):
        # Print banner
        banner()
        print("\033[4mCart3D Shape Optimisation\033[0m".center(50, " "))

        # Iterate
        iteration = 0
        while True:
            # Evaluate gradient
            gradient = 0

            iteration += 1

            if iteration >= max_iterations:
                break

            if np.linalg.norm(gradient) <= convergence_tolerance:
                break

    @abstractmethod
    def objective_function(self, x):
        """Evaluates the objective function to be minimised at x."""

    @abstractmethod
    def jacobian(self, x):
        """Evaluates the Jacobian at x."""

    def post_process(
        self, plot_convergence: bool = True, theoretical_convergence: float = None
    ) -> pd.DataFrame:
        """Crawls through iteration directories to compile results and
        optionally plot convergence of the objective function."""
        iteration_dirs = [
            i
            for i in os.listdir(self.working_dir)
            if os.path.isdir(os.path.join(self.working_dir, i))
        ]

        results = []
        for directory in iteration_dirs:
            iteration = int(directory)

            completion_file = os.path.join(
                self.working_dir, directory, self.completion_filename
            )
            if os.path.exists(completion_file):
                # This iteration completed, load the results
                iter_result = pd.read_csv(completion_file)

                # Add iteration number
                iter_result.loc[len(iter_result)] = ["iteration", iteration]

                # Append
                results.append(iter_result.set_index("Unnamed: 0").to_dict()["0"])

        df = pd.DataFrame(results).set_index("iteration").sort_index()

        if plot_convergence:
            plt.plot(df.index, df["objective"], label="PySAGAS ShapeOpt convergence")
            plt.title("Convergence of Objective Function")
            plt.xlabel("Iteration")
            plt.ylabel("Objective Function")

            if theoretical_convergence:
                plt.axhline(
                    theoretical_convergence,
                    c="k",
                    ls="--",
                    label="Theoretical convergence",
                )

            plt.legend()
            plt.grid()
            plt.show()

        return df
