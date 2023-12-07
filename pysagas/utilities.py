import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import List, Optional
from pysagas.geometry import Cell


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
