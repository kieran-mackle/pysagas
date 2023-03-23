import os
import time
import glob
import pickle
import shutil
import subprocess
import numpy as np
import pandas as pd
from pysagas.optimisation import ShapeOpt
from pysagas.wrappers import Cart3DWrapper
from hypervehicle.generator import Generator
from pyoptsparse import Optimizer, Optimization
from typing import List, Dict, Optional, Optional
from pysagas.optimisation.cart3d.utilities import C3DPrep
from hypervehicle.utilities import SensitivityStudy, append_sensitivities_to_tri


np.seterr(all="ignore")


class Cart3DShapeOpt(ShapeOpt):
    """A wrapper to perform shape optimisation with Cart3D."""

    C3D_errors = [
        "==> ADAPT failed",
        "Check cart3d.out in AD_A_J for more clues",
        "==> adjointErrorEst_quad failed again, status = 1",
        "ERROR: CUBES failed",
        "ERROR: ADAPT failed with status = 1",
        "ERROR",
    ]

    def __init__(
        self,
        a_inf: float,
        rho_inf: float,
        V_inf: float,
        A_ref: float,
        optimiser: Optimizer,
        vehicle_generator: Generator,
        objective_callback: callable,
        jacobian_callback: callable,
        sensitivity_filename: str = "all_components_sensitivity.csv",
        working_dir_name: str = "working_dir",
        sim_directory_name: str = "simulation",
        basefiles_dir_name: str = "basefiles",
        c3d_log_name: str = "C3D_log",
        c3d_info_file: str = None,
        matching_tolerance: float = 1e-5,
        optimiser_options: dict = None,
        save_evolution: bool = True,
    ) -> None:
        """Initialise Cart3D Shape Optimiser.

        Parameters
        ----------
        optimiser : Optimizer
            The pyoptsparse Optimizer object of choice.

        vehicle_generator : Generator
            The vehicle generator object.

        objective_callback : callable
            The callback function to compute and return the objective function
            value and constraint violation values.

        jacobian_callback : callable
            The callback function to compute and return the Jacobian of the
            objective function value and constraint violation values.

        working_dir_name : str
            The name of the working directory.

        optimiser_options : dict, optional
            The options to pass to the optimiser.
        """
        # Define global variable so that functions can access them
        global c3d_logname, _matching_tolerance, _max_matching_tol, _matching_target
        global sens_filename, basefiles_dir, _c3dprepper, sim_dir_name
        global _rho_inf, _V_inf, _a_inf, _A_ref
        global working_dir, f_sense_filename
        global obj_cb, jac_cb
        global generator
        global save_comptri, evo_dir

        save_comptri = save_evolution

        generator = vehicle_generator

        # Construct paths
        home_dir = os.getcwd()
        basefiles_dir = os.path.join(home_dir, basefiles_dir_name)
        working_dir = os.path.join(home_dir, working_dir_name)
        sim_dir_name = sim_directory_name
        sens_filename = sensitivity_filename
        f_sense_filename = "F_sensitivities.csv"
        c3d_logname = c3d_log_name

        if save_comptri:
            # Create evolution history directory
            evo_dir = os.path.join(home_dir, "evolution")
            if not os.path.exists(evo_dir):
                os.mkdir(evo_dir)

        # Save callback functions
        obj_cb = objective_callback
        jac_cb = jacobian_callback

        # TODO - pass as flowstate
        _rho_inf = rho_inf
        _V_inf = V_inf
        _a_inf = a_inf
        _A_ref = A_ref

        # Create instance of Cart3D prepper
        _c3dprepper = C3DPrep(logfile=c3d_logname, info_file=c3d_info_file)

        # Other settings
        _matching_tolerance = matching_tolerance
        _max_matching_tol = 0.1
        _matching_target = 0.9

        # Construct optimisation problem
        self.opt_problem = Optimization(
            name="Cart3D-PySAGAS Shape Optimisation",
            objFun=evaluate_objective,
            sens=evaluate_gradient,
        )

        # Complete super initialisation
        super().__init__(
            optimiser=optimiser,
            working_dir=working_dir,
            optimiser_options=optimiser_options,
        )


def evaluate_objective(x: dict) -> dict:
    """Evaluates the objective function at the parameter set `x`."""
    # Pre-process parameters
    _process_parameters(x)

    # Run Cart3D simulation
    sim_success, loads_dict, _ = _run_simulation()

    if sim_success:
        # Load properties
        properties_dir = glob.glob("*_properties")
        if properties_dir:
            volmass = pd.read_csv(
                glob.glob(os.path.join(properties_dir[0], "*volmass.csv"))[0],
                index_col=0,
            )["0"]

            # Fetch user-defined properties
            properties_file = glob.glob(
                os.path.join(properties_dir[0], "*properties.csv")
            )
            if len(properties_file) > 0:
                # File exists, load it
                properties = pd.read_csv(
                    properties_file[0],
                    index_col=0,
                )["0"]

        else:
            # No properties data found
            volmass = None
            properties = None

        # Evaluate objective function
        funcs = obj_cb(loads_dict=loads_dict, volmass=volmass, properties=properties)
        failed = False

    else:
        # Simulation failed
        funcs = {}
        failed = True

    return funcs, failed


def evaluate_gradient(x: dict, objective: dict) -> dict:
    """Evaluates the gradient function at the parameter set `x`."""
    # Pre-process parameters
    _process_parameters(x)

    # Run Cart3D simulation
    _, loads_dict, components_plt_filepath = _run_simulation()

    # Initialise filepaths
    components_filepath = components_plt_filepath
    sensitivity_filepath = os.path.join(working_dir, sens_filename)

    # Create PySAGAS wrapper and run
    try:
        wrapper = Cart3DWrapper(
            a_inf=_a_inf,
            rho_inf=_rho_inf,
            sensitivity_filepath=sensitivity_filepath,
            components_filepath=components_filepath,
            verbosity=0,
        )

    except ValueError:
        # The sensitivity data does not match the point data, regenerate it
        tri_components_filepath = os.path.join(
            working_dir, sim_dir_name, "Components.i.tri"
        )
        sensitivity_files = glob.glob(os.path.join("*sensitivity*"))
        _combine_sense_data(
            tri_components_filepath,
            sensitivity_files=sensitivity_files,
            match_target=_matching_target,
            tol_0=_matching_tolerance,
            max_tol=_max_matching_tol,
            outdir=working_dir,
        )

        # Re-instantiate the wrapper
        wrapper = Cart3DWrapper(
            a_inf=_a_inf,
            rho_inf=_rho_inf,
            sensitivity_filepath=sensitivity_filepath,
            components_filepath=components_filepath,
            verbosity=0,
        )

    F_sense, _ = wrapper.calculate()

    # Non-dimensionalise
    coef_sens = F_sense / (0.5 * _rho_inf * _A_ref * _V_inf**2)

    # Calculate Jacobian
    properties_dir = glob.glob("*_properties")
    if properties_dir:
        scalar_sens_dir = os.path.join("scalar_sensitivities")
        vm = pd.read_csv(
            glob.glob(os.path.join(properties_dir[0], "*volmass.csv"))[0],
            index_col=0,
        )["0"]
        vm_sens = pd.read_csv(
            os.path.join(scalar_sens_dir, "volmass_sensitivity.csv"),
            index_col=0,
        )[x.keys()]

        # Fetch user-defined properties and sensitivities
        properties_file = glob.glob(os.path.join(properties_dir[0], "*properties.csv"))
        if len(properties_file) > 0:
            # File exists, load it
            properties = pd.read_csv(
                properties_file[0],
                index_col=0,
            )["0"]

            # Also load sensitivity file
            property_sens = pd.read_csv(
                os.path.join(scalar_sens_dir, "property_sensitivity.csv"),
                index_col=0,
            )[x.keys()]

    else:
        # No properties data found
        vm = None
        vm_sens = None
        properties = None
        property_sens = None

    # Call function
    jac = jac_cb(
        parameters=x,
        coef_sens=coef_sens,
        loads_dict=loads_dict,
        volmass=vm,
        volmass_sens=vm_sens,
        properties=properties,
        property_sens=property_sens,
    )

    return jac


def _process_parameters(x):
    # Move into working directory
    os.chdir(working_dir)

    # Load existing parameters to compare
    already_started = _compare_parameters(x)

    if not already_started:
        # These parameters haven't run yet, prepare the working directory
        if save_comptri:
            # Move components file into evolution directory
            comp_filepath = os.path.join(working_dir, sim_dir_name, "Components.i.tri")

            # Check if file exists
            if os.path.exists(comp_filepath):
                # Determine new name
                new_filename = f"{len(os.listdir(evo_dir)):04d}.tri"

                # Copy file over
                shutil.copyfile(
                    comp_filepath,
                    os.path.join(evo_dir, new_filename),
                )

        # Delete any files from previous run
        _clean_dir(working_dir)

    # Dump design parameters to file
    with open("parameters.pkl", "wb") as f:
        pickle.dump(x, f)

    # Generate vehicle and geometry sensitivities
    if len(glob.glob("*sensitivity*")) == 0 or not already_started:
        # No sensitivity files generated yet, or this is new geometry
        parameters = _unwrap_x(x)
        ss = SensitivityStudy(vehicle_constructor=generator)
        ss.dvdp(parameter_dict=parameters, perturbation=2, write_nominal_stl=True)
        ss.to_csv()


def _run_simulation(no_attempts: int = 3):
    """Prepare and run the CFD simulation with Cart3D. The simulation will be
    run in the 'simulation' subdirectory of the iteration directory.
    """
    # Make simulation directory
    sim_dir = os.path.join(working_dir, sim_dir_name)
    run_intersect = False
    components_filepath = os.path.join(sim_dir, "Components.i.tri")
    if not os.path.exists(sim_dir):
        os.mkdir(sim_dir)
        run_intersect = True
    else:
        # Check for intersected file
        run_intersect = not os.path.exists(components_filepath)
        intersected = True

    # Attempt component intersection
    sim_success = False
    for attempt in range(no_attempts):
        # Run intersect
        if run_intersect:
            _c3dprepper._log(f"SHAPEOPT INTERSECT ATTEMPT {attempt+1}")
            intersected = _c3dprepper.intersect_stls()

        # Check for intersection
        if intersected:
            # Prepare rest of simulation directory
            if not os.path.exists(os.path.join(sim_dir, "input.cntl")):
                # Move files to simulation directory (including Components.i.tri)
                _c3dprepper.run_autoinputs()
                os.system(
                    f"mv *.tri Config.xml input.c3d preSpec.c3d.cntl {sim_dir} >> {c3d_logname} 2>&1"
                )

                # Copy sim files and permissions
                for filename in ["input.cntl", "aero.csh"]:
                    shutil.copyfile(
                        os.path.join(basefiles_dir, filename),
                        os.path.join(sim_dir, filename),
                    )
                    shutil.copymode(
                        os.path.join(basefiles_dir, filename),
                        os.path.join(sim_dir, filename),
                    )

            # Create all_components_sensitivity.csv
            if not os.path.exists(sens_filename):
                _combine_sense_data(
                    components_filepath,
                    sensitivity_files=glob.glob("*sensitivity*"),
                    match_target=_matching_target,
                    tol_0=_matching_tolerance,
                    max_tol=_max_matching_tol,
                )

            # Run Cart3D and await result
            os.chdir(sim_dir)
            target_adapt = _infer_adapt(sim_dir)
            c3d_donefile = os.path.join(sim_dir, target_adapt, "FLOW", "DONE")
            run_cmd = "./aero.csh restart"
            _restarts = 0
            if not os.path.exists(c3d_donefile):
                with open(c3d_logname, "a") as f:
                    subprocess.run(
                        run_cmd, shell=True, stdout=f, stderr=subprocess.STDOUT
                    )
                while not os.path.exists(c3d_donefile):
                    # Wait...
                    time.sleep(5)

                    # Check for C3D failure
                    running, e = _c3d_running()

                    if not running:
                        # C3D failed, try restart it
                        if _restarts > 3:
                            return False

                        f = open(c3d_logname, "a")
                        subprocess.run(
                            run_cmd, shell=True, stdout=f, stderr=subprocess.STDOUT
                        )
                        f.close()
                        _restarts += 1

            sim_success = True
            break

    if sim_success:
        # Sim finished successfully, read loads file
        loads_dict = _read_c3d_loads(
            os.path.join(working_dir, sim_dir_name, "BEST/FLOW/loadsCC.dat")
        )
        components_plt_filepath = os.path.join(
            working_dir, sim_dir_name, "BEST/FLOW/Components.i.plt"
        )

    else:
        loads_dict = None
        components_plt_filepath = None

    # Change back to working directory
    os.chdir(working_dir)

    return sim_success, loads_dict, components_plt_filepath


def _read_c3d_loads(
    loadsCC_filepath: str,
    b_frame: bool = True,
    v_frame: bool = True,
    moments: bool = True,
) -> dict:
    load_dict = {}
    with open(loadsCC_filepath, "r") as file:
        for line in file:
            if line[0] == "#":
                # Commented line, skip
                continue

            # Remove linebreaks and multiple spaces
            line = " ".join(line.split())
            words = line.split(":")

            if len(words) == 0 or len(words) == 1:  # skip if empty line
                continue

            text = words[0]
            number = float(words[1])
            word = text.split(" ")
            tag = word[0]
            coeff = word[-1]
            coeff = coeff[1:4]
            if b_frame is True:
                if coeff in ["C_A", "C_Y", "C_N"]:  # get force in b_frame
                    load_dict["{0}-{1}".format(coeff, tag)] = number
            if v_frame is True:
                if coeff in ["C_D", "C_S", "C_L"]:  # get force in v_frame
                    load_dict["{0}-{1}".format(coeff, tag)] = number
            if moments is True:
                if coeff in ["C_l", "C_m", "C_n", "C_M"]:  # get moment coeff
                    load_dict["{0}-{1}".format(coeff, tag)] = number

    return load_dict


def _combine_sense_data(
    components_filepath: str,
    sensitivity_files: List[str],
    match_target: float = 0.9,
    tol_0: float = 1e-5,
    max_tol: float = 1e-1,
    outdir: Optional[str] = None,
):
    """Combine the component sensitivity data for intersected geometry."""
    match_frac = 0
    tol = tol_0
    while match_frac < match_target:
        # Check tolerance
        if tol > max_tol:
            raise Exception(
                "Cannot combine sensitivity data (match fraction: "
                + f"{match_frac}, tolerance: {tol}, max tolerance: {max_tol})."
            )

        # Run matching algorithm
        match_frac = append_sensitivities_to_tri(
            dp_filenames=sensitivity_files,
            components_filepath=components_filepath,
            match_tolerance=tol,
            verbosity=0,
            outdir=outdir,
        )

        if match_frac < match_target:
            print(
                "Failed to combine sensitivity data "
                f"({100*match_frac:.02f}% match rate)."
            )
            print("  Increasing matching tolerance and trying again.")

        # Increase matching tolerance
        tol *= 10


def _c3d_running() -> bool:
    """Watches the Cart3D log file to check for errors and return False
    if Cart3D has stopped running."""
    with open(c3d_logname) as f:
        # Get last line in log file
        for line in f:
            pass

        # Check if if it is in the known errors
        for e in Cart3DShapeOpt.C3D_errors:
            if e in line:
                return False, e

    # No errors
    return True, None


def _infer_adapt(sim_dir) -> str:
    with open(f"{sim_dir}/aero.csh", "r") as f:
        lines = f.readlines()

        for line in lines:
            if line.find("set n_adapt_cycles") != -1:
                return f"adapt{int(line.split('=')[-1]):02d}"


def _compare_parameters(x):
    """Compares the current parameters x to the last run simulation
    parameters."""
    try:
        with open("parameters.pkl", "rb") as f:
            xp = pickle.load(f)

        # Compare to current parameters
        already_run = x == xp

    except FileNotFoundError:
        # Simulation not run yet
        already_run = False

    return already_run


def _clean_dir(directory: str, keep: list = None):
    """Deletes everything in a directory except for what is specified
    in keep."""
    all_files = os.listdir(directory)
    if keep is None:
        # Convert to empty list
        keep = []
    rm_files = set(all_files) - set(keep)
    for f in rm_files:
        if os.path.isdir(f):
            shutil.rmtree(f)
        else:
            os.remove(f)


def _unwrap_x(x: dict) -> dict:
    """Unwraps an ordered dictionary."""
    unwrapped = {}
    for key, val in x.items():
        if len(val) == 1:
            unwrapped[key] = val[0]
        else:
            unwrapped[key] = val
    return unwrapped
