import numpy as np
from pysagas import FlowState
from pysagas.geometry import Cell


def piston_sensitivity(cell: Cell, p_i: int, **kwargs):
    """Calculates the pressure-parameter sensitivity using
    local piston theory.

    Parameters
    ----------
    cell : Cell
        The cell object.

    p_i : int
        The index of the parameter to find the sensitivity for. This is used to
        index cell.dndp.
    """
    M_l = cell.flowstate.M
    if M_l < 1.0:
        # Subsonic cell, skip
        return 0

    dPdp = (
        cell.flowstate.rho
        * cell.flowstate.a
        * np.dot(cell.flowstate.vec, -cell.dndp[:, p_i])
    )
    return dPdp


def van_dyke_sensitivity(
    cell: Cell,
    p_i,
    **kwargs,
):
    """
    Calculates the pressure-parameter sensitivity using
    Van Dyke second-order theory.

     Parameters
    ----------
    cell : Cell
        The cell object.

    p_i : int
        The index of the parameter to find the sensitivity for. This is used to
        index cell.dndp.
    """
    M_l = cell.flowstate.M
    if M_l < 1.0:
        # Subsonic cell, skip
        return 0

    piston = piston_sensitivity(cell=cell, p_i=p_i)
    dPdp = piston * M_l / (M_l**2 - 1) ** 0.5
    return dPdp


def isentropic_sensitivity(cell: Cell, p_i: int, **kwargs):
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


def freestream_isentropic_sensitivity(cell: Cell, p_i: int, eng_outflow: FlowState, eng_sens, **kwargs):
    """Calculates the pressure-parameter sensitivity, including
    the sensitivity to the incoming flow state (for use on nozzle cells
    where the engine outflow changes due to parameter change"""

    gamma1 = eng_outflow.gamma
    gamma2 = cell.flowstate.gamma
    beta1 = np.sqrt(eng_outflow.M ** 2 - 1)
    beta2 = np.sqrt(cell.flowstate.M ** 2 - 1)
    fun1 = 1 + (gamma1 - 1) / 2 * eng_outflow.M ** 2
    fun2 = 1 + (gamma2 - 1) / 2 * cell.flowstate.M ** 2

    # Calculate sens of P to grid changes using the isentropic method
    dP2_dgeom = isentropic_sensitivity(cell=cell, p_i=p_i, **kwargs)

    # Calculate sens to inflow Mach number
    dP2_dM1 = (-gamma2 * cell.flowstate.P * cell.flowstate.M ** 2 / eng_outflow.M
              * beta1 / (beta2 * fun1))

    # Calculate sens to inflow pressure
    dP2_dP1 = (fun1 / fun2) ** (gamma2 / (gamma2 - 1))

    # Calculate sens to inflow aoa
    dP2_daoa = -cell.flowstate.M ** 2 / beta2 * gamma2 * cell.flowstate.P

    # sum contributions
    dPdp = (dP2_dgeom
            + dP2_dM1 * eng_sens.flow_sens.loc['dMout'][p_i]
            + dP2_dP1 * eng_sens.flow_sens.loc['dPout'][p_i]
            + dP2_daoa * eng_sens.flow_sens.loc['daoa'][p_i])

    return dPdp











