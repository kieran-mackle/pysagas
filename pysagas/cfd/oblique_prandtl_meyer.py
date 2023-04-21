import numpy as np
from typing import Optional
from scipy.optimize import bisect
from pysagas.flow import FlowState
from pysagas.geometry import Vector
from scipy.optimize import root_scalar
from pysagas.cfd.solver import FlowSolver


class OPM(FlowSolver):
    """Oblique shock thoery and Prandtl-Meyer expansion theory flow solver.

    This implementation uses oblique shock theory for flow-facing cell
    elements, and Prandtl-Meyer expansion theory for rearward-facing elements.

    Extended Summary
    ----------------
    Data attribute 'method' refers to which method was used for a particular
    cell, according to:
        -1 : invalid / skipped (eg. 90 degree face)
        0 : parallel face, do nothing
        1 : Prandlt-Meyer
        2 : normal shock
        3 : oblique shock
    """

    method = "Oblique/Prandtl-Meyer combination"

    def solve(
        self,
        freestream: Optional[FlowState] = None,
        Mach: Optional[float] = None,
        aoa: Optional[float] = None,
    ):
        super().solve(freestream=freestream, Mach=Mach, aoa=aoa)

        # Get flow
        flow = self._last_solve_freestream

        # Iterate over all cells
        # TODO - add progress bar
        net_force = Vector(0, 0, 0)
        for cell in self.cells:
            # Calculate orientation of cell to flow
            theta = np.pi / 2 - np.arccos(
                np.dot(flow.direction.vec, cell.n.vec)
                / (cell.n.norm * flow.direction.norm)
            )

            # Solve flow for this cell
            if theta < 0:
                # Rearward facing cell
                if theta == np.deg2rad(-90):
                    # TODO: what do do here?
                    M2, p2, T2 = (flow.M, 0.0, flow.T)
                    method = -1

                else:
                    # Use Prandtl-Meyer expansion theory
                    M2, p2, T2 = self.solve_pm(
                        abs(theta), flow.M, flow.P, flow.T, flow.gamma
                    )
                    method = 1

            elif theta > 0:
                # Use shock theory
                beta_max = OPM.beta_max(M=flow.M, gamma=flow.gamma)
                theta_max = OPM.theta_from_beta(
                    M1=flow.M, beta=beta_max, gamma=flow.gamma
                )
                if theta > theta_max:
                    # Detached shock
                    method = 2
                    M2, p2, T2 = self.solve_normal(flow.M, flow.P, flow.T, flow.gamma)

                else:
                    # Use oblique shock theory
                    method = 3
                    M2, p2, T2 = self.solve_oblique(
                        theta, flow.M, flow.P, flow.T, flow.gamma
                    )

            else:
                # Cell is parallel to flow, assume no change
                method = 0
                M2, p2, T2 = (flow.M, flow.P, flow.T)

            # Save results for this cell
            cell.attributes["pressure"] = p2
            cell.attributes["Mach"] = M2
            cell.attributes["temperature"] = T2
            cell.attributes["method"] = method

            # Calculate force vector
            net_force += cell.n * p2 * cell.A

        # Construct results
        # TODO - improve return results, perhaps create FlowResults class

        return net_force

    @staticmethod
    def pm(M: float, gamma: float = 1.4):
        """Solves the Prandtl-Meyer function and returns the Prandtl-Meyer angle
        in radians."""
        v = ((gamma + 1) / (gamma - 1)) ** 0.5 * np.arctan(
            ((M**2 - 1) * (gamma - 1) / (gamma + 1)) ** 0.5
        ) - np.arctan((M**2 - 1) ** 0.5)
        return v

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
            The deflection angle, specified in radians.

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
    def solve_normal(M1: float, p1: float = 1.0, T1: float = 1.0, gamma: float = 1.4):
        """Solves the flow using normal shock theory.

        Parameters
        ----------
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
        rho2_rho1 = (gamma + 1) * M1**2 / (2 + (gamma - 1) * M1**2)
        p2_p1 = 1 + 2 * gamma * (M1**2 - 1) / (gamma + 1)

        M2 = ((1 + M1 * (gamma - 1) / 2) / (gamma * M1**2 - (gamma - 1) / 2)) ** 0.5
        p2 = p1 * p2_p1
        T2 = T1 * p2_p1 / rho2_rho1

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
        func = lambda beta: OPM.theta_from_beta(M1, beta, gamma) - theta

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

        # Finally, solve with secant method
        # beta = sign_beta * secant(func, b1, b2, tol=tolerance)
        root_result = root_scalar(f=func, x0=b1, x1=b2, xtol=tolerance, method="secant")
        if root_result.converged:
            beta = sign_beta * root_result.root
        else:
            raise Exception(
                f"Cannot converge beta (theta={np.rad2deg(theta):.2f} degrees)."
            )

        return beta

    @staticmethod
    def theta_from_beta(M1: float, beta: float, gamma: float = 1.4):
        """Calculates the flow deflection angle from the oblique shock angle.

        Parameters
        ----------
        M1 : float
            The pre-expansion Mach number.

        beta : float
            The shock angle specified in radians.

        gamma : float, optional
            The ratio of specific heats. The default is 1.4.

        Returns
        --------
        theta : float
            The deflection angle, specified in radians.

        References
        ----------
        Peter Jacobs
        """
        M1n = M1 * abs(np.sin(beta))
        t1 = 2.0 / np.tan(beta) * (M1n**2 - 1)
        t2 = M1**2 * (gamma + np.cos(2 * beta)) + 2
        theta = np.arctan(t1 / t2)
        return theta

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
            The deflection angle, specified in radians.

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

    @staticmethod
    def beta_max(M: float, gamma: float = 1.4):
        """Returns the maximum shock angle for a given
        Mach number.
        """
        beta_max = np.arcsin(
            np.sqrt(
                (1 / (gamma * M**2))
                * (
                    (gamma + 1) * M**2 / 4
                    - 1
                    + np.sqrt(
                        (gamma + 1)
                        * ((gamma + 1) * M**4 / 16 + (gamma - 1) * M**2 / 2 + 1)
                    )
                )
            )
        )
        return beta_max

    @staticmethod
    def dp_dtheta_obl(
        theta: float, M1: float, p1: float = 1.0, T1: float = 1.0, gamma: float = 1.4
    ):
        """Solves for the sensitivities at the given point using oblique shock
        theory.

        Parameters
        ----------
        theta : float
            The deflection angle, specified in radians.

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
        dM2 : float
            The sensitivity of the post-expansion Mach number.

        dp2 : float
            The sensitivity of the post-expansion pressure (Pa).

        dT2 : float
            The sensitivity of the post-expansion temperature (K).
        """
        # Create function handle
        func = lambda theta: OPM.solve_oblique(theta, M1, p1, T1, gamma)

        # Calculate derivitive by finite differencing
        g = OPM._findiff(func, theta)

        return g

    @staticmethod
    def dp_dtheta_pm(
        theta: float, M1: float, p1: float = 1.0, T1: float = 1.0, gamma: float = 1.4
    ):
        """Solves for the sensitivities at the given point using Prandtl-Meyer
        theory.

        Parameters
        ----------
        theta : float
            The deflection angle, specified in radians.

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
        dM2 : float
            The sensitivity of the post-expansion Mach number.

        dp2 : float
            The sensitivity of the post-expansion pressure (Pa).

        dT2 : float
            The sensitivity of the post-expansion temperature (K).
        """
        # Create function handle
        func = lambda theta: OPM.solve_pm(theta, M1, p1, T1, gamma)

        # Calculate derivitive by finite differencing
        g = OPM._findiff(func, theta)

        return g

    @staticmethod
    def _findiff(func: callable, point: any, perturbation: float = 1e-3):
        # vals = [func(p)[1] for p in [point*(1-perturbation), point, point*(1+perturbation)]]
        vals = [
            func(p)
            for p in [point * (1 - perturbation), point, point * (1 + perturbation)]
        ]
        Mvals = [v[0] for v in vals]
        pvals = [v[1] for v in vals]
        Tvals = [v[2] for v in vals]

        dM = np.gradient(Mvals)
        dp = np.gradient(pvals)
        dT = np.gradient(Tvals)

        return dM[1], dp[1], dT[1]

    def save(self, name: str):
        # Initialise attributes dictionary
        attributes = {
            "pressure": [],
            "Mach": [],
            "temperature": [],
            "method": [],
        }

        # Call super method
        super().save(name, attributes)
