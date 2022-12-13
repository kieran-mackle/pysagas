import os
import sys
import pandas as pd
from typing import Tuple


def process_components_file(
    a_inf: float,
    rho_inf: float,
    filepath: str = "Components.i.plt",
    write_data: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A ParaView script to process Components.i.plt.

    Parameters
    ----------
    a_inf : float
        The freestream speed of sound (m/s).
    rho_inf : float
        The freestream density (kg/m^3).
    filepath : str, optional
        The filepath to the Components.i.plt file to be processed.
        The default is Components.i.plt.
    write_data : bool, optional
        Write the flow data to CSV files. The default is True.

    Returns
    -------
    points : pd.DataFrame
        A DataFrame of the point data.
    cells : pd.DataFrame
        A DataFrame of the cell data.
    """
    # Check file exists
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File '{filepath}' does not exist.")

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
    print(f"Loading {filepath}")
    componentsiplt = VisItTecplotBinaryReader(FileName=[filepath])
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
    programmableFilter1.Script = f"""
    # Paraview script to normalise Cart3D variables.
    R_gas = 287.058
    gamma = 1.4

    # Redefine the Cart3D variables - normalised values are *_tilde
    output.PointData.append(inputs[0].PointData["U"], "U_tilde")
    output.PointData.append(inputs[0].PointData["V"], "V_tilde")
    output.PointData.append(inputs[0].PointData["W"], "W_tilde")
    output.PointData.append(inputs[0].PointData["Pressure"], "p_tilde")
    output.PointData.append(inputs[0].PointData["Rho"], "rho_tilde")

    # Define the dimensional flow properties
    output.PointData.append({a_inf} * inputs[0].PointData["U"], "U")
    output.PointData.append({a_inf} * inputs[0].PointData["V"], "V")
    output.PointData.append({a_inf} * inputs[0].PointData["W"], "W")

    output.PointData.append({a_inf} * {a_inf} * {rho_inf} * output.PointData["p_tilde"], "p")
    output.PointData.append({rho_inf} * output.PointData["rho_tilde"], "rho")
    output.PointData.append({a_inf}*{a_inf}*{rho_inf}*output.PointData["p_tilde"] / ({rho_inf}*output.PointData["rho_tilde"]*R_gas), "T")
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
    # TODO - get parent dir from components file
    points_filename = "points.csv"
    ExportView(points_filename, view=spreadSheetView1)

    # Properties modified on spreadSheetView1
    spreadSheetView1.FieldAssociation = "Cell Data"

    # export view
    print("  Saving cell data.")
    # TODO - get parent dir from components file
    cells_filename = "cells.csv"
    ExportView(cells_filename, view=spreadSheetView1)

    print("Complete.")

    # Read the files to Pandas DataFrames
    points = pd.read_csv(points_filename)
    cells = pd.read_csv(cells_filename)

    if not write_data:
        # Delete the files
        os.remove(points_filename)
        os.remove(cells_filename)

    return points, cells
