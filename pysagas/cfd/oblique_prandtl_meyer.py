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

        # Get flow
        flow = self._last_solve_freestream

        for cell in self.cells:
            # Calculate orientation of cell to flow
            theta = np.arccos(
                np.dot(cell.n.vec, flow.direction.vec)
                / np.cross(cell.n.vec, flow.direction.vec)
            )

            # Solve flow for this cell
            if theta < 0:
                # Rearward facing cell, use Prandtl-Meyer expansion theory
                M2, p2, T2 = self.solve_pm(theta, flow.M, flow.P, flow.T, flow.gamma)
            elif theta > 0:
                # Use oblique shock theory
                M2, p2, T2 = self.solve_oblique(
                    theta, flow.M, flow.P, flow.T, flow.gamma
                )
            else:
                # Cell is parallel to flow, assume no change
                M2, p2, T2 = (flow.M, flow.P, flow.T)

            # Save results for this cell
            # not sure what this should look like yet

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
        func = lambda M: OPM.pm(M, gamma) - angle
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

        M1 : float
            The pre-expansion Mach number.

        p1 : float, optional
            The pre-expansion pressure (Pa). The default is 1.0.

        T1 : float, optional
            The pre-expansion temperature (K). The default is 1.0.

        gamma : float, optional
            The ratio of specific heats. The default is 1.4.

        Returns
        --------
        M2 : float
            The post-expansion Mach number.

        p2 : float
            The post-expansion pressure (Pa).

        T2 : float
            The post-expansion temperature (K).
        """
        # Solve for M2
        v_M1 = OPM.pm(M=M1, gamma=gamma)
        v_M2 = theta + v_M1
        M2 = OPM.inv_pm(v_M2)

        a = (gamma - 1) / 2
        n = 1 + a * M1**2
        d = 1 + a * M2**2

        p2 = p1 * (n / d) ** (gamma / (gamma - 1))

        T2 = T1 * (n / d)

        return M2, p2, T2

    @staticmethod
    def solve_oblique(
        theta: float, M1: float, p1: float = 1.0, T1: float = 1.0, gamma: float = 1.4
    ):
        """Solves the flow using oblique shock theory.

        Parameters
        ----------
        theta : float
            The flow deflection angle, specified in radians.

        M1 : float
            The pre-shock Mach number.

        p1 : float, optional
            The pre-shock pressure (Pa). The default is 1.0.

        T1 : float, optional
            The pre-expansion temperature (K). The default is 1.0.

        gamma : float, optional
            The ratio of specific heats. The default is 1.4.

        Returns
        --------
        M2 : float
            The post-shock Mach number.

        p2 : float
            The post-shock pressure (Pa).

        T2 : float
            The post-shock temperature (K).
        """
        # Calculate angle and ratios
        beta = OPM.oblique_beta(M1, theta, gamma)
        p2_p1 = OPM.oblique_p2_p1(M1, beta, gamma)
        rho2_rho1 = OPM.oblique_rho2_rho1(M1, beta, gamma)
        T2_T1 = p2_p1 / rho2_rho1

        # Calculate properties
        M2 = OPM.oblique_M2(M1, beta, theta, gamma)
        p2 = p2_p1 * p1
        T2 = T1 * T2_T1

        return M2, p2, T2

    @staticmethod
    def oblique_p2_p1(M1: float, beta: float, gamma: float = 1.4):
        """Returns the pressure ratio p2/p1 across an oblique shock.

        Parameters
        -----------
        M1 : float
            The pre-shock Mach number.

        beta : float
            The shock angle specified in radians.

        Returns
        --------
        p2/p1 : float
            The pressure ratio across the shock.

        References
        ----------
        Peter Jacobs
        """
        M1n = M1 * abs(np.sin(beta))
        p2p1 = 1.0 + 2.0 * gamma / (gamma + 1.0) * (M1n**2 - 1.0)
        return p2p1

    @staticmethod
    def oblique_beta(
        M1: float, theta: float, gamma: float = 1.4, tolerance: float = 1.0e-6
    ):
        """Calculates the oblique shock angle using a recursive algorithm.

        Parameters
        ----------
        M1 : float
            The pre-shock Mach number.

        theta : float
            The flow deflection angle, specified in radians.

        gamma : float, optional
            The ratio of specific heats. The default is 1.4.

        References
        ----------
        Peter Jacobs
        """

        def func(beta):
            """Recursive function to solve."""
            return OPM.oblique_beta(M1, beta, gamma) - theta

        # Initialise
        sign_beta = np.sign(theta)
        theta = abs(theta)
        b1 = np.arcsin(1.0 / M1)
        b2 = b1 * 1.05

        # Check f1
        f1 = func(b1)
        if abs(f1) < tolerance:
            return sign_beta * b1

        # Check f2
        f2 = func(b2)
        if abs(f2) < tolerance:
            return sign_beta * b2

        # Finally, use bisection
        beta = sign_beta * bisect(func, b1, b2)
        return beta

    @staticmethod
    def oblique_M2(M1: float, beta: float, theta: float, gamma: float = 1.4):
        """Calculates the Mach number following an oblique shock.

        Parameters
        ----------
        M1 : float
            The pre-expansion Mach number.

        beta : float
            The shock angle specified in radians.

        theta : float
            The deflection angle, specified in degrees.

        gamma : float, optional
            The ratio of specific heats. The default is 1.4.

        Returns
        --------
        M2 : float
            The post-shock Mach number.

        References
        ----------
        Peter Jacobs
        """
        M1n = M1 * abs(np.sin(beta))
        a = 1 + (gamma - 1) * 0.5 * M1n**2
        b = gamma * M1n**2 - (gamma - 1) * 0.5
        M2 = (a / b / (np.sin(beta - theta)) ** 2) ** 0.5
        return M2

    def oblique_T2_T1(M1: float, beta: float, gamma: float = 1.4):
        """Returns the temperature ratio T2/T1 across an oblique shock.

        Parameters
        -----------
        M1 : float
            The pre-shock Mach number.

        beta : float
            The shock angle specified in radians.

        gamma : float, optional
            The ratio of specific heats. The default is 1.4.

        Returns
        --------
        T2/T1 : float
            The temperature ratio across the shock.

        References
        ----------
        Peter Jacobs
        """
        T2_T1 = OPM.oblique_p2_p1(M1, beta, gamma) / OPM.oblique_rho2_rho1(
            M1, beta, gamma
        )
        return T2_T1

    @staticmethod
    def oblique_rho2_rho1(M1: float, beta: float, gamma: float = 1.4):
        """Returns the density ratio rho2/rho1 across an oblique shock.

        Parameters
        -----------
        M1 : float
            The pre-shock Mach number.

        beta : float
            The shock angle specified in radians.

        gamma : float, optional
            The ratio of specific heats. The default is 1.4.

        Returns
        --------
        T2/T1 : float
            The temperature ratio across the shock.

        References
        ----------
        Peter Jacobs
        """
        M1n = M1 * abs(np.sin(beta))
        rho2_rho1 = (gamma + 1) * M1n**2 / 2 + (gamma - 1) * M1n**2
        return rho2_rho1


if __name__ == "__main__":
    # Example 9.9, page 654
    theta = 1
    M1 = 1.5
    p1 = 1
    T1 = 288

    M2, p2, T2 = OPM.solve_pm(theta, M1, p1, T1)

    print(f"M2 = {M2}")
    print(f"p2 = {p2}")
    print(f"T2 = {T2}")
