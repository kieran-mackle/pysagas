import numpy as np
from scipy.optimize import bisect
from pysagas.flow import FlowState
from typing import List, Tuple, Optional
from pysagas.geometry import Vector, Cell
from pysagas.cfd.solver import FlowSolver


class OPM(FlowSolver):
    """Oblique shock thoery and Prandtl-Meyer expansion theory flow solver.

    This implementation uses oblique shock theory for flow-facing cell
    elements, and Prandtl-Meyer expansion theory for rearward-facing elements.
    """

    method = "Oblique/Prandtl-Meyer combination"

    def solve(self, freestream: Optional[FlowState] = None):
        super().solve(freestream)

        raise NotImplementedError("Coming soon!")

    @staticmethod
    def pm(M: float, gamma: float = 1.4):
        """Solves the Prandtl-Meyer function and returns the Prandtl-Meyer angle
        in degrees."""
        v = ((gamma + 1) / (gamma - 1)) ** 0.5 * np.arctan(
            ((M**2 - 1) * (gamma - 1) / (gamma + 1)) ** 0.5
        ) - np.arctan((M**2 - 1) ** 0.5)
        return np.rad2deg(v)

    @staticmethod
    def inv_pm(angle: float, gamma: float = 1.4):
        """Solves the inverse Prandtl-Meyer function using a bisection algorithm."""
        func = lambda M: pm(M, gamma) - angle
        return bisect(func, 1.0, 42.0)

    @staticmethod
    def solve_pm(
        theta: float, M1: float, p1: float = 1.0, T1: float = 1.0, gamma: float = 1.4
    ):
        """Solves for the Mach number, pressure and temperature after a Prandtl-Meyer
        expansion fan.

        Parameters
        ----------
        theta : float
            The deflection angle, specified in degrees.

        p1 : float
            The pre-expansion pressure.

        M1 : float
            The pre-expansion Mach number.

        gamma : float, optional
            The ratio of specific heats. The default is 1.4.

        Returns
        --------
        p2 : float
            The post-expansion pressure.
        """
        # Solve for M2
        v_M1 = pm(M=M1, gamma=gamma)
        v_M2 = theta + v_M1
        M2 = inv_pm(v_M2)

        a = (gamma - 1) / 2
        n = 1 + a * M1**2
        d = 1 + a * M2**2

        p2 = p1 * (n / d) ** (gamma / (gamma - 1))

        T2 = T1 * (n / d)

        return M2, p2, T2


if __name__ == "__main__":
    # Example 9.9, page 654
    theta = 1
    M1 = 1.5
    p1 = 1
    T1 = 288

    M2, p2, T2 = solve_pm(theta, M1, p1, T1)

    print(f"M2 = {M2}")
    print(f"p2 = {p2}")
    print(f"T2 = {T2}")
