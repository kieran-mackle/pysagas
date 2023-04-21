import numpy as np
from pysagas.flow import FlowState
from typing import List, Tuple, Optional
from pysagas.geometry import Vector, Cell
from pysagas.cfd.solver import FlowSolver


class NewtonianImpact(FlowSolver):
    """Newtonian impact flow solver."""

    method = "Newtonian Impact"

    def solve(
        self,
        freestream: Optional[FlowState] = None,
        Mach: Optional[float] = None,
        aoa: Optional[float] = None,
    ):
        # Initial checks
        super().solve(freestream=freestream, Mach=Mach, aoa=aoa)

        raise NotImplementedError("Coming soon!")

        # Get flow
        flow = self._last_solve_freestream

        # TODO - review this code; what about ref area and length?

        # Prepare analysis
        a_inf = 299.499  # m/s
        V_inf = flow.M * a_inf
        p_inf = flow.rho * a_inf**2 / flow.gamma
        q_inf = 0.5 * flow.rho * V_inf**2

        # Run analysis
        F, M = self.newtonian_impact_solver(
            cells=self.cells, v=freestream.direction, p_inf=p_inf, q_inf=q_inf
        )

        C_f = F / q_inf

    @staticmethod
    def newtonian_cp(cell: Cell, v: Vector):
        """Calculates all direction force sensitivities.

        Parameters
        ----------
        cell : Cell
            The cell.
        flow : FlowState
            The freestream flow condition.

        Returns
        --------
        Cp : float
            The non-dimensional pressure coefficient for the cell.
        """
        return 2 * np.dot(cell.n.vec, v.unit.vec) ** 2

    @staticmethod
    def newtonian_impact_solver(
        cells: List[Cell], v: Vector, p_inf: float, q_inf: float
    ) -> Tuple[Vector, np.array]:
        c_ps = [NewtonianImpact.newtonian_cp(c, v) for c in cells]
        ps = [c_p * q_inf + p_inf for c_p in c_ps]

        # Calculate forces
        fs = [c.n * ps[i] * c.A for i, c in enumerate(cells)]
        net_force = Vector(0, 0, 0)
        for f in fs:
            net_force += f

        return net_force, None
