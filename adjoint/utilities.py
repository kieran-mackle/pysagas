import numpy as np
import gdtk.ideal_gas_flow as igf
from adjoint.flow import FlowState
from adjoint.geometry import Vector, Cell


def calculate_pressures(flow: FlowState, theta: float) -> float:
    beta = igf.beta_obl(M1=flow.M, theta=abs(theta), g=flow.gamma, tol=1.0e-6)
    P2_P1 = igf.p2_p1_obl(flow.M, beta, g=flow.gamma)
    P2 = P2_P1 * flow.P
    return P2


def calculate_force_vector(P: float, n: np.array, A: float):
    """Calculates the force vector components, acting on a
    surface defined by its area and normal vector.
    """

    F_x = A * P * np.dot(n, np.array([-1, 0, 0]))
    F_y = A * P * np.dot(n, np.array([0, -1, 0]))
    F_z = A * P * np.dot(n, np.array([0, 0, -1]))

    return [F_x, F_y, F_z]


def force_sensitivity(
    P: float,
    area: float,
    area_dp,
    pressure_sense,
    combined_sense,
    n: np.array,
    direction: np.array,
) -> np.array:
    """Calculates the force sensitivity."""
    sensitivity = (
        pressure_sense * area * np.dot(n, direction)
        + P * area_dp * np.dot(n, direction)
        + P * area * np.dot(-combined_sense, direction)
    )
    return sensitivity


def all_force_sensitivities(
    P: float,
    area: float,
    area_dp,
    pressure_sense,
    combined_sense,
    n: np.array,
) -> np.array:
    """Calculates all direction force sensitivities.

    Parameters
    ----------
    To be documented.

    Returns
    --------
    sensitivities : np.array
        An array of shape n x m, for a 3-dimensional cell with
        n parameters. The value of m is derived from the length
        of area_dp.
    """
    sensitivities = np.zeros(shape=(len(area_dp), 3))
    all_directions = [Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)]
    for i, direction in enumerate(all_directions):
        dir_sens = force_sensitivity(
            P_ramp=P,
            area=area,
            area_dp=area_dp,
            pressure_sense=pressure_sense,
            combined_sense=combined_sense,
            n=n,
            direction=direction.vec,
        )
        sensitivities[:, i] = dir_sens

    return sensitivities


def dfdp(
    cells: list[Cell], P: float, rho: float, a: float, vel_vector: np.array
) -> np.array:
    """Calcualtes the force sensitivities for a list of Cells."""
    dFdp = 0
    for cell in cells:
        # Calculate pressure sensitivity
        dPdp = rho * a * np.dot(vel_vector, -cell.dndp)

        # Calculate force sensitivity
        cell_dFdp = all_force_sensitivities(
            P_ramp=P,
            area=cell.A,
            area_dp=cell.dAdp,
            pressure_sense=dPdp,
            combined_sense=cell.dndp,
            n=cell.n.vec,
        )

        dFdp += cell_dFdp

    return dFdp
