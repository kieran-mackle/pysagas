import numpy as np
import pandas as pd
from tqdm import tqdm
import gdtk.ideal_gas_flow as igf
from numpy.typing import ArrayLike
from pysagas.flow import FlowState
from pysagas.geometry import Vector, Cell
from typing import List, Callable, Tuple, Union, Optional


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
        The pressure (Pa).

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
    cell: Cell, dPdp_method: Callable, cog: Vector = Vector(0, 0, 0), **kwargs
) -> Tuple[np.array, np.array]:
    """Calculates all direction force sensitivities.

    Parameters
    ----------
    cell : Cell
        The cell.

    dPdp_method : Callable
        The method to use when calculating the pressure sensitivities.

    cog : Vector, optional
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


def piston_dPdp(cell: Cell, p_i: int, **kwargs):
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
    dPdp = (
        cell.flowstate.rho
        * cell.flowstate.a
        * np.dot(cell.flowstate.vec, -cell.dndp[:, p_i])
    )
    return dPdp


def van_dyke_dPdp(
    cell: Cell,
    p_i,
    freestream: FlowState,
    dp: Union[List, ArrayLike],
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
    mach_inf = freestream.M
    a_inf = freestream.a
    gamma = freestream.gamma

    mach = cell.flowstate.M

    if mach < 1.0:
        # Subsonic cell, skip
        return 0

    beta = np.sqrt(mach**2 - 1)

    # Calculate normal velocity
    v_n = np.dot(cell.flowstate.vec, -cell.dndp[:, p_i] * dp[p_i])

    a = (
        -cell.flowstate.vec
        * (2 / (mach_inf**2))
        * (
            mach_inf / (a_inf * beta)
            + (v_n / a_inf)
            * ((gamma + 1) * mach_inf**4 - 4 * (mach**2 - 1))
            / (2 * a_inf * (mach**2 - 1))
        )
    )

    dCPdp = np.dot(a, cell.dndp[:, p_i])

    # Normalise to pressure sensitivity
    dPdp = dCPdp * freestream.q

    return dPdp


def van_dyke_dPdp_ingo(
    cell: Cell,
    p_i,
    **kwargs,
):
    """
    Calculates the pressure-parameter sensitivity using
    Van Dyke second-order theory.
    """
    M_l = cell.flowstate.M
    if M_l < 1.0:
        # Subsonic cell, skip
        return 0

    piston = piston_dPdp(cell=cell, p_i=p_i)
    dPdp = piston * M_l / (M_l**2 - 1) ** 0.5

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
    dPdp_method: Callable = piston_dPdp,
    cog: Vector = Vector(0, 0, 0),
    **kwargs,
) -> Tuple[np.array, np.array]:
    """Calcualtes the force sensitivities for a list of Cells.

    Parameters
    ----------
    cells : list[Cell]
        The cells to be analysed.

    dPdp_method : Callable, optional
        The method used to calculate the pressure/parameter sensitivities.
        The default is the Panel method approximation panel_dPdp (see below).

    cog : Vector, optional
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


def add_sens_data(
    cells: List[Cell],
    data: pd.DataFrame,
    verbosity: Optional[int] = 1,
    match_tolerance: Optional[float] = 1e-5,
    rounding_tolerance: Optional[float] = 1e-8,
    force: Optional[bool] = False,
) -> float:
    """Appends shape sensitivity data to a list of Cell objects.

    Parameters
    ----------
    cells : List[Cell]
        A list of cell objects.

    data : pd.DataFrame
        The sensitivity data to be transcribed onto the cells.

    verbosity : int, optional
        The verbosity. The default is 1.

    match_tolerance : float, optional
        The precision tolerance for matching point coordinates. The
        default is 1e-5.

    rounding_tolerance : float, optional
        The tolerance to round data off to. The default is 1e-8.

    force : bool, optional
        Force the sensitivity data to be added, even if a cell
        already has sensitivity data. This can be used if new
        data is being used, or if the append is being repeated with
        a different matching tolerance. The default is False.

    Returns
    --------
    match_fraction : float
        The fraction of points matched.
    """
    # Extract parameters
    parameters = []
    param_cols = data.columns[3:]
    for i in range(int(len(param_cols) / 3)):
        parameters.append(param_cols[int(i * 3)].split("dxd")[-1])

    # Construct progress bar
    if verbosity > 0:
        print()
        desc = "Adding geometry-parameter sensitivity data to cells"
        pbar = tqdm(
            total=len(cells),
            desc=desc,
            position=0,
            leave=True,
        )

    matched_points = 0
    total_points = 0
    skipped_cells = 0
    total_cells = 0
    for cell in cells:
        # Check if cell already has sensitivity data
        total_cells += 1
        if cell.dvdp is None or force:
            # Initialise sensitivity
            dvdp = np.zeros((9, len(parameters)))
            for i, point in enumerate([cell.p0, cell.p1, cell.p2]):
                match_x = (point.x - data["x"]).abs() < match_tolerance
                match_y = (point.y - data["y"]).abs() < match_tolerance
                match_z = (point.z - data["z"]).abs() < match_tolerance

                match = match_x & match_y & match_z
                try:
                    # Get match
                    matched_data = data[match].iloc[0][param_cols]

                    # Round off infinitesimally small values
                    matched_data[abs(matched_data) < rounding_tolerance] = 0

                    for j, p in enumerate(parameters):
                        # For each parameter (column)
                        for k, c in enumerate(["x", "y", "z"]):
                            dvdp[3 * i + k, j] = matched_data[f"d{c}d{p}"]

                    # Update match count
                    matched_points += 1

                except IndexError:
                    # No match found, leave as zero sensitivity
                    pass

                # Update total_points
                total_points += 1

            cell._add_sensitivities(np.array(dvdp))

        else:
            # Skip this cell
            skipped_cells += 1

        # Update progress bar
        if verbosity > 0:
            pbar.update(1)

    try:
        match_fraction = matched_points / total_points
    except:
        match_fraction = 1

    if verbosity > 0:
        pbar.close()
        print(
            f"Done - matched {100*match_fraction:.2f}% of points "
            + f"(skipped {100*skipped_cells/total_cells:.2f}% of cells)."
        )

    # TODO - allow dumping data to file

    return match_fraction
