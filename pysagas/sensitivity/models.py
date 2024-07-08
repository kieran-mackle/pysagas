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
    piston = piston_sensitivity(cell=cell, p_i=p_i)
    M_l = cell.flowstate.M
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


def freestream_isentropic_sensitivity(cell: Cell, p_i: int, inflow: FlowState, inflow_sens, **kwargs):
    """Calculates the pressure-parameter sensitivity, including
    the sensitivity to the incoming flow state (for use on nozzle cells
    where the engine outflow changes due to parameter change)

    Parameters
        ----------
        cell : Cell
            The cell to be analysed.

        p_i : int
            Index of the design variable (in inflow_sens) to differeciate with respect to

        inflow : FlowState
            "Free stream" of the infoming flow (sould be the engine combustor outflow for nozzle claculation)

        inflow_sens :
            HyperPro engine outflow sensitivities

        Returns
        --------
        dPdp : np.array
            The pressure sensitivity matrix with respect to the parameter.
    """

    gamma1 = inflow.gamma
    gamma2 = cell.flowstate.gamma
    M1 = inflow.M
    M2 = cell.flowstate.M
    P1 = inflow.P
    P2 = cell.flowstate.P

    beta1 = np.sqrt(M1 ** 2 - 1)
    beta2 = np.sqrt(M2 ** 2 - 1)
    fun1 = 1 + (gamma1 - 1) / 2 * M1 ** 2
    fun2 = 1 + (gamma2 - 1) / 2 * M2 ** 2

    # Calculate sens to inflow Mach number
    dM2_dM1 = (M2 / M1) * (beta1 / beta2) * (fun2 / fun1)
    t1 = M1 * fun1 ** (1/(gamma1-1)) * fun2 ** ((-gamma1)/(gamma1-1))
    t2 = M2 * fun1 ** (gamma1 / (gamma1 - 1)) * fun2 ** ((1 - 2 * gamma1) / (gamma1 - 1)) * dM2_dM1
    dP2_dM1 = P1 * gamma1 * (t1 - t2)


    # Calculate sens to inflow pressure
    dP2_dP1 = (fun1 / fun2) ** (gamma2 / (gamma2 - 1))

    # Calculate sens to inflow aoa
    dP2_daoa = M2 ** 2 / beta2 * gamma2 * P2

    # Calculate sens to inflow gamma
    gp1 = gamma1 + 1
    gm1 = gamma1 - 1
    fg = gp1 / gm1
    fM = beta1**2 / fg
    num1 = np.sqrt(fg) *(beta1 / gp1) ** 2
    den1 = np.sqrt(fM) * (fM + 1)
    num2 = 1 / gm1 * np.arctan(np.sqrt(fM))
    den2 = np.sqrt(gm1 * gp1)
    dnu_dg1 = num1 / den1 - num2 / den2

    q = (fun1/fun2)**(gamma1/gm1)
    r1 = gamma1 * (M1 ** 2 - (M2 ** 2 * fun1) / fun2)
    r2 = 2*gm1*fun1
    s1 = 1/gm1 - gamma1/(gm1**2)
    s2 = np.log(fun1/fun2)
    df_dg = q*(r1/r2 + s1*s2)

    dP2_dg1 = P1*df_dg - dP2_daoa*dnu_dg1

    # sum contributions
    dPdp = (  dP2_dM1 * inflow_sens.loc['M'][p_i]
            + dP2_dP1 * inflow_sens.loc['P'][p_i]
            + dP2_daoa * inflow_sens.loc['flow_angle'][p_i]
            + dP2_dg1 * inflow_sens.loc['gamma'][p_i])

    return dPdp











