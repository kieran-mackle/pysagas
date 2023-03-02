import numpy as np
import gdtk.ideal_gas_flow as igf
from pysagas.flow import FlowState
from typing import List, Callable, Tuple
from pysagas.geometry import Vector, Cell


def calculate_pressures(flow: FlowState, theta: float) -> float:
    """Calculates the pressure from a flow state and deflecion angle
    using ideal gas oblique shock theory.

    Parameters
    ----------
    flow : FlowState
        The flow state.
    theta : float
        The deflection angle (radians).

    Returns
    --------
    P2 : float
        The pressure behind the oblique shock.
    """
    beta = igf.beta_obl(M1=flow.M, theta=abs(theta), g=flow.gamma, tol=1.0e-6)
    P2_P1 = igf.p2_p1_obl(flow.M, beta, g=flow.gamma)
    P2 = P2_P1 * flow.P
    return P2


def calculate_force_vector(P: float, n: np.array, A: float) -> np.array:
    """Calculates the force vector components, acting on a
    surface defined by its area and normal vector.

    Parameters
    ----------
    P : float
        The pressure.
    n : np.array
        The normal vector.
    A : float
        The reference area (m^2).

    Returns
    --------
    forces : np.array
        The force components.
    """
    F_x = A * P * np.dot(n, np.array([-1, 0, 0]))
    F_y = A * P * np.dot(n, np.array([0, -1, 0]))
    F_z = A * P * np.dot(n, np.array([0, 0, -1]))

    return [F_x, F_y, F_z]


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


def newtonian_impact_solver(
    cells: List[Cell], v: Vector, p_inf: float, q_inf: float
) -> (Vector, np.array):
    c_ps = [newtonian_cp(c, v) for c in cells]
    ps = [c_p * q_inf + p_inf for c_p in c_ps]

    # Calculate forces
    fs = [c.n * ps[i] * c.A for i, c in enumerate(cells)]
    net_force = Vector(0, 0, 0)
    for f in fs:
        net_force += f

    return net_force, None


def cell_dfdp(
    cell: Cell, dPdp_method: Callable, cog: Vector = Vector(0, 0, 0), **kwargs
) -> Tuple[np.array, np.array]:
    """Calculates all direction force sensitivities.

    Parameters
    ----------
    cell : Cell
        The cell.
    dPdp_method : Callable
        The method to use when calculating the pressure sensitivities.
    cog : Vector, optiona
        The reference centre of gravity, used in calculating the moment
        sensitivities. The defualt is Vector(0,0,0).

    Returns
    --------
    sensitivities : np.array
        An array of shape n x 3, for a 3-dimensional cell with
        n parameters.

    See Also
    --------
    all_dfdp : a wrapper to calculate force sensitivities for many cells
    """
    # Initialisation
    all_directions = [Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)]
    sensitivities = np.empty(shape=(cell.dndp.shape[1], 3))
    moment_sensitivities = np.empty(shape=(cell.dndp.shape[1], 3))

    # Calculate moment dependencies
    r = cell.c - cog
    F = cell.flowstate.P * cell.A * cell.n.vec

    # For each parameter
    for p_i in range(cell.dndp.shape[1]):
        # Calculate pressure sensitivity
        dPdp = dPdp_method(cell=cell, p_i=p_i, **kwargs)

        # Evaluate for sensitivities for each direction
        for i, direction in enumerate(all_directions):
            dF = (
                dPdp * cell.A * np.dot(cell.n.vec, direction.vec)
                + cell.flowstate.P * cell.dAdp[p_i] * np.dot(cell.n.vec, direction.vec)
                + cell.flowstate.P * cell.A * np.dot(-cell.dndp[:, p_i], direction.vec)
            )
            sensitivities[p_i, i] = dF

        # Now evaluate moment sensitivities
        moment_sensitivities[p_i, :] = np.cross(
            r.vec, sensitivities[p_i, :]
        ) + np.cross(cell.dcdp[:, p_i], F)

    # Append to cell
    cell.sensitivities = sensitivities
    cell.moment_sensitivities = moment_sensitivities

    return sensitivities, moment_sensitivities


def panel_dPdp(cell: Cell, p_i, **kwargs):
    """Calculates the pressure-parameter sensitivity using
    the Panel method approximation."""
    dPdp = (
        cell.flowstate.rho
        * cell.flowstate.a
        * np.dot(cell.flowstate.vec, -cell.dndp[:, p_i])
    )
    return dPdp


def isentropic_dPdp(cell: Cell, p_i: int, **kwargs):
    """Calculates the pressure-parameter sensitivity using
    the isentropic flow relation directly."""
    gamma = cell.flowstate.gamma
    power = (gamma + 1) / (gamma - 1)
    dPdW = (cell.flowstate.P * gamma / cell.flowstate.a) * (
        1 + cell.flowstate.v * (gamma - 1) / (2 * cell.flowstate.a)
    ) ** power
    dWdn = -cell.flowstate.vec
    dndp = cell.dndp[:, p_i]
    dPdp = dPdW * np.dot(dWdn, dndp)
    return dPdp


def all_dfdp(
    cells: List[Cell],
    dPdp_method: Callable = panel_dPdp,
    cog: Vector = Vector(0, 0, 0),
    **kwargs
) -> Tuple[np.array, np.array]:
    """Calcualtes the force sensitivities for a list of Cells.

    Parameters
    ----------
    cells : list[Cell]
        The cells to be analysed.
    dPdp_method : Callable, optional
        The method used to calculate the pressure/parameter sensitivities.
        The default is the Panel method approximation panel_dPdp (see below).
    cog : Vector, optiona
        The reference centre of gravity, used in calculating the moment
        sensitivities. The defualt is Vector(0,0,0).

    Returns
    --------
    dFdp : np.array
        The force sensitivity matrix with respect to the parameters.
    dMdp : np.array
        The moment sensitivity matrix with respect to the parameters.

    See Also
    --------
    cell_dfdp : the force sensitivity per cell
    panel_dPdp : pressure sensitivities calculated using Panel method
        approximations
    """
    dFdp = 0
    dMdp = 0
    for cell in cells:
        # Calculate force sensitivity
        dFdp_c, dMdp_c = cell_dfdp(
            cell=cell, dPdp_method=dPdp_method, cog=cog, **kwargs
        )

        dFdp += dFdp_c
        dMdp += dMdp_c

    return dFdp, dMdp
