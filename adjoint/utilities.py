import sys
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


def cell_dfdp(cell: Cell) -> np.array:
    """Calculates all direction force sensitivities.

    Parameters
    ----------
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
    # Calculate pressure sensitivity
    dPdp = (
        cell.flowstate.rho * cell.flowstate.a * np.dot(cell.flowstate.vec, -cell.dndp)
    )

    sensitivities = np.zeros(shape=(len(cell.dAdp), 3))
    all_directions = [Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)]
    for i, direction in enumerate(all_directions):
        dir_sens = (
            dPdp * cell.A * np.dot(cell.n.vec, direction.vec)
            + cell.flowstate.P * cell.dAdp * np.dot(cell.n.vec, direction.vec)
            + cell.flowstate.P * cell.A * np.dot(-cell.dndp, direction.vec)
        )
        sensitivities[:, i] = dir_sens

    return sensitivities


def all_dfdp(cells: List[Cell]) -> np.array:
    """Calcualtes the force sensitivities for a list of Cells.

    Parameters
    ----------
    cells : list[Cell]
        The cells to be analysed.

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
        # Calculate force sensitivity
        cell_dFdp = cell_dfdp(
            cell=cell,
        )

        dFdp += cell_dFdp

    return dFdp


def process_components_file():
    """A ParaView script to process Components.i.plt."""

    try:
        # import the simple module from the paraview
        from paraview.simple import (
            VisItTecplotBinaryReader,
            _DisableFirstRenderCameraReset,
            GetActiveViewOrCreate,
            Show,
            ProgrammableFilter,
            Hide,
            PointDatatoCellData,
            ExportView,
        )
    except ModuleNotFoundError:
        # Cannot find paraview python package
        print(
            "Cannot find ParaView Python package. If ParaView is already "
            + "installed, please append the bin/ directory to the Python path. "
            + "If it is not installed, please do so. If you are using "
            + "an Anaconda environment, you can install using "
            + "'conda install -c conda-forge paraview'."
        )
        sys.exit()

    # disable automatic camera reset on 'Show'
    _DisableFirstRenderCameraReset()

    # create a new 'VisItTecplotBinaryReader'
    print("Loading Components.i.plt.")
    componentsiplt = VisItTecplotBinaryReader(FileName=["Components.i.plt"])
    componentsiplt.MeshStatus = ["Surface"]
    componentsiplt.PointArrayStatus = []

    # Properties modified on componentsiplt
    componentsiplt.PointArrayStatus = [
        "Cp",
        "Pressure",
        "Rho",
        "U",
        "V",
        "W",
        "x",
        "y",
        "z",
    ]

    # get active view
    spreadSheetView1 = GetActiveViewOrCreate("SpreadSheetView")
    # uncomment following to set a specific view size
    # spreadSheetView1.ViewSize = [400, 400]

    # show data in view
    componentsipltDisplay = Show(componentsiplt, spreadSheetView1)

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    # create a new 'Programmable Filter'
    programmableFilter1 = ProgrammableFilter(Input=componentsiplt)
    programmableFilter1.Script = ""
    programmableFilter1.RequestInformationScript = ""
    programmableFilter1.RequestUpdateExtentScript = ""
    programmableFilter1.PythonPath = ""

    # Properties modified on programmableFilter1
    print("  Dimensionalising attributes.")
    programmableFilter1.Script = """# Paraview script to normalise Cart3D variables.

    ######## Define freestream conditions here ###########
    M_inf = 6
    ######################################################

    R_gas = 287.058
    gamma = 1.4

    # Look up the freestream and sound speed for the given Mach number and q = 50 kPa
    rho_lookup = {5:0.0450814, 6:0.0308742, 7:0.0225510, 7.5:0.0195525, 8:0.0171295, 9:0.0134424, 10:0.0107105}
    a_lookup = {5:297.891, 6:299.499, 7:300.840, 7.5:301.575, 8:302.018, 9:303.061, 10:305.562}
    a_inf = a_lookup[M_inf]
    rho_inf = rho_lookup[M_inf]

    # Redefine the Cart3D variables - normalised values are *_tilde
    output.PointData.append(inputs[0].PointData["U"], "U_tilde")
    output.PointData.append(inputs[0].PointData["V"], "V_tilde")
    output.PointData.append(inputs[0].PointData["W"], "W_tilde")
    output.PointData.append(inputs[0].PointData["Pressure"], "p_tilde")
    output.PointData.append(inputs[0].PointData["Rho"], "rho_tilde")

    # Define the dimensional flow properties
    output.PointData.append(a_inf * inputs[0].PointData["U"], "U")
    output.PointData.append(a_inf * inputs[0].PointData["V"], "V")
    output.PointData.append(a_inf * inputs[0].PointData["W"], "W")

    output.PointData.append(a_inf * a_inf * rho_inf * output.PointData["p_tilde"], "p")
    output.PointData.append(rho_inf * output.PointData["rho_tilde"], "rho")
    output.PointData.append(a_inf*a_inf*rho_inf*output.PointData["p_tilde"] / (rho_inf*output.PointData["rho_tilde"]*R_gas), "T")
    output.PointData.append((abs(gamma * R_gas * output.PointData["T"]))**0.5, "a")
    output.PointData.append((output.PointData["U"]**2 + output.PointData["V"]**2 + output.PointData["W"]**2)**0.5, "V_mag")
    output.PointData.append(output.PointData["V_mag"]/output.PointData["a"], "M")"""
    programmableFilter1.RequestInformationScript = ""
    programmableFilter1.RequestUpdateExtentScript = ""
    programmableFilter1.PythonPath = ""

    # show data in view
    programmableFilter1Display = Show(programmableFilter1, spreadSheetView1)

    # hide data in view
    Hide(componentsiplt, spreadSheetView1)

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    # create a new 'Point Data to Cell Data'
    pointDatatoCellData1 = PointDatatoCellData(Input=programmableFilter1)

    # show data in view
    pointDatatoCellData1Display = Show(pointDatatoCellData1, spreadSheetView1)

    # hide data in view
    Hide(programmableFilter1, spreadSheetView1)

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    # Properties modified on pointDatatoCellData1
    pointDatatoCellData1.PassPointData = 1

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    # Properties modified on spreadSheetView1
    spreadSheetView1.GenerateCellConnectivity = 1

    # export view
    print("  Saving point data.")
    ExportView("points.csv", view=spreadSheetView1)

    # Properties modified on spreadSheetView1
    spreadSheetView1.FieldAssociation = "Cell Data"

    # export view
    print("  Saving cell data.")
    ExportView("cells.csv", view=spreadSheetView1)

    print("Complete.")
