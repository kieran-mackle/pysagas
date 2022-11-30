import numpy as np
from typing import List
import gdtk.ideal_gas_flow as igf
from adjoint.flow import FlowState
from adjoint.geometry import Vector, Cell


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


def cell_dfdp(
    P: float,
    pressure_sense: np.array,
    cell: Cell,
) -> np.array:
    """Calculates all direction force sensitivities.

    Parameters
    ----------
    P : float
        The flow pressure.
    pressure_sense : np.array
        The pressure sensitivity array.
    cell : Cell
        The cell.

    Returns
    --------
    sensitivities : np.array
        An array of shape n x m, for a 3-dimensional cell with
        n parameters. The value of m is derived from the length
        of area_dp.

    See Also
    --------
    all_dfdp : a wrapper to calculate force sensitivities for many cells
    """
    sensitivities = np.zeros(shape=(len(cell.dAdp), 3))
    all_directions = [Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)]
    for i, direction in enumerate(all_directions):
        dir_sens = (
            pressure_sense * cell.A * np.dot(cell.n.vec, direction.vec)
            + P * cell.dAdp * np.dot(cell.n.vec, direction.vec)
            + P * cell.A * np.dot(-cell.dndp, direction.vec)
        )
        sensitivities[:, i] = dir_sens

    return sensitivities


def all_dfdp(
    cells: List[Cell], P: float, rho: float, a: float, vel_vector: np.array
) -> np.array:
    """Calcualtes the force sensitivities for a list of Cells.

    Parameters
    ----------
    cells : list[Cell]
        The cells to be analysed.
    P : float
        The flow pressure.
    rho : float
        The flow density.
    a : float
        The flow speed of sound.
    vel_vector : np.array
        The flow velocity vector.

    Returns
    --------
    dFdp : np.array
        The force sensitivity matrix with respect to the parameters.

    See Also
    --------
    cell_dfdp : the force sensitivity per cell
    """
    dFdp = 0
    for cell in cells:
        # Calculate pressure sensitivity
        dPdp = rho * a * np.dot(vel_vector, -cell.dndp)

        # Calculate force sensitivity
        cell_dFdp = cell_dfdp(
            P=P,
            pressure_sense=dPdp,
            cell=cell,
        )

        dFdp += cell_dFdp

    return dFdp
