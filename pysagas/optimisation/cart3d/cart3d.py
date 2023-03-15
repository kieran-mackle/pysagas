import os
import time
import glob
import shutil
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from pysagas import banner
from pysagas.optimisation import ShapeOpt
from pysagas.wrappers import Cart3DWrapper
from hypervehicle.generator import Generator
from pysagas.optimisation.cart3d.utilities import C3DPrep
from typing import List, Dict, Optional, Optional, Callable
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
        home_dir: str,
        a_inf: float,
        rho_inf: float,
        V_inf: float,
        A_ref: float,
        generator: Generator,
        sensitivity_filename: str = "all_components_sensitivity.csv",
        working_dir_name: str = "working_dir",
        sim_dir_name: str = "simulation",
        basefiles_dir_name: str = "basefiles",
        c3d_logname: str = "C3D_log",
        c3d_info_file: str = None,
        matching_tolerance: float = 1e-5,
    ) -> None:
        # Construct paths
        self.root_dir = home_dir
        self.basefiles_dir = os.path.join(home_dir, basefiles_dir_name)
        self.working_dir = os.path.join(home_dir, working_dir_name)
        self.sim_dir_name = sim_dir_name
        self.sensitivity_filename = sensitivity_filename
        self.completion_filename = "ITERATION_COMPLETE"
        self.parameters_filename = "parameters.csv"
        self.f_sense_filename = "F_sensitivities.csv"
        self.jacobian_filename = "jacobian.csv"
        self.objective_filename = "objective.txt"
        self._c3d_checkpoint_rename = "ref_checkpoint"

        self.generator = generator

        self.c3d_logname = c3d_logname

        # Save callback function
        # TODO - implement standard functions (eg. min drag, L/D, etc.)
        self._obj_jac_cb = None

        # TODO - pass as flowstate
        self.rho_inf = rho_inf
        self.V_inf = V_inf
        self.a_inf = a_inf
        self.A_ref = A_ref

        # Create instance of Cart3D prepper
        self._c3dprepper = C3DPrep(logfile=c3d_logname, info_file=c3d_info_file)

        # Other settings
        self._matching_tolerance = matching_tolerance
        self._max_matching_tol = 0.1
        self._matching_target = 0.9

    def evaluate_objective(self, x: dict) -> dict:
        pass

    def evaluate_gradient(self, x: dict) -> dict:
        pass

    def _run_sensitivity_study(
        self, iter_dir: str, param_names: List[str], x: List[float]
    ):
        """Runs the geometric sensitivity study for the specified iteration
        directory iter_dir. This method will produce the base geometry STL
        files, and the sensitivity files.
        """
        # Change into iteration directory
        os.chdir(iter_dir)

        # Check if sensitivity study has been run
        sens_files = glob.glob(f"{iter_dir}{os.sep}*sensitivity*")
        if len(sens_files) == 0:
            print("Running sensitivity study.")

            # Run the sensitivity study and generate the nominal geometry
            parameters = dict(zip(param_names, x))

            # Save these parameters for future reference
            pd.Series(parameters).to_csv(
                os.path.join(iter_dir, self.parameters_filename)
            )

            # Run sensitivity study
            ss = SensitivityStudy(vehicle_constructor=self.generator)
            ss.dvdp(parameter_dict=parameters, perturbation=2, write_nominal_stl=True)
            ss.to_csv()

        else:
            print("Sensitivity study already run.")

    def _run_simulation(
        self,
        basefiles_dir: str,
        iter_dir: str,
        max_adapt: int = None,
        warmstart: bool = False,
    ):
        """Prepare and run the CFD simulation with Cart3D. The simulation will be
        run in the 'simulation' subdirectory of the iteration directory. If warmstart
        is True, the previous iteration will be used to warm-start the solution.
        """
        # Make simulation directory
        sim_dir = os.path.join(iter_dir, self.sim_dir_name)
        run_intersect = False
        components_filepath = os.path.join(sim_dir, "Components.i.tri")
        if not os.path.exists(sim_dir):
            os.mkdir(sim_dir)
            run_intersect = True
        else:
            # Check for intersected file
            run_intersect = not os.path.exists(components_filepath)
            intersected = True
            print("Intersected components located.")

        if warmstart:
            # Create warm_iter_dir
            warm_iter_dir = (
                Path(iter_dir)
                .parent.joinpath(f"{int(Path(iter_dir).name)-1:04d}")
                .as_posix()
            )

        # Attempt component intersection
        N = 3
        for attempt in range(N):
            # Run intersect
            if run_intersect:
                self._c3dprepper._log(f"SHAPEOPT INTERSECT ATTEMPT {attempt+1}")
                intersected = self._c3dprepper.intersect_stls()

            # Check for intersection
            if intersected:
                # Prepare rest of simulation directory
                if not os.path.exists(os.path.join(sim_dir, "input.cntl")):
                    # Prepare remaining C3D files
                    if warmstart:
                        # Copy necessary files from warm-start directory
                        self._copy_warmstart_files(warm_iter_dir, iter_dir)

                    else:
                        # Move files to simulation directory (including Components.i.tri)
                        self._c3dprepper.run_autoinputs()
                        os.system(
                            f"mv *.tri Config.xml input.c3d preSpec.c3d.cntl {sim_dir} >> {self.c3d_logname} 2>&1"
                        )

                        # Copy sim files
                        os.system(
                            f"cp {basefiles_dir}/input.cntl {basefiles_dir}/aero.csh {sim_dir} >> {self.c3d_logname} 2>&1"
                        )

                        # Modify iteration aero.csh using max_adapt
                        if max_adapt:
                            self._overwrite_adapt(sim_dir, max_adapt)

                # Create all_components_sensitivity.csv
                if not os.path.exists(self.sensitivity_filename):
                    self._combine_sense_data(
                        components_filepath,
                        sensitivity_files=glob.glob("*sensitivity*"),
                        match_target=self._matching_target,
                        tol_0=self._matching_tolerance,
                        max_tol=self._max_matching_tol,
                    )

                # Override warmstart flag
                if os.path.exists(os.path.join(sim_dir, "aero.csh")) and warmstart:
                    # Warmstart set to True, but aero.csh file present
                    warmstart = False

                # Run Cart3D and await result
                os.chdir(sim_dir)
                if warmstart:
                    # Get checkpoint filename and Cart3D commands
                    warm_sim_dir = os.path.join(warm_iter_dir, self.sim_dir_name)
                    commands = self._get_cart_commands(warm_sim_dir)

                    # Define wait files and run command
                    c3d_donefile = os.path.join(sim_dir, "loadsCC.dat")
                    run_cmd = commands["flowCart"]

                else:
                    target_adapt = self._infer_adapt(sim_dir)
                    c3d_donefile = os.path.join(sim_dir, target_adapt, "FLOW", "DONE")
                    run_cmd = "./aero.csh restart"

                _restarts = 0
                if not os.path.exists(c3d_donefile):
                    # Cart3D has not started / didn't finish
                    print(
                        "\nStarting Cart3D, awaiting",
                        os.sep.join(c3d_donefile.split(os.sep)[-6:]),
                    )

                    if warmstart:
                        # Prepare for warm-start
                        with open(self.c3d_logname, "a") as f:
                            f.write(f"\nRunning cubes: {commands['cubes']}")
                            subprocess.run(
                                f"{commands['cubes']} -remesh",
                                shell=True,
                                stdout=f,
                                stderr=subprocess.STDOUT,
                            )
                            f.write(f"\n{commands['mgPrep']}")
                            subprocess.run(
                                f"{commands['mgPrep']}",
                                shell=True,
                                stdout=f,
                                stderr=subprocess.STDOUT,
                            )
                            f.write(f"\nRunning mesh2mesh...")
                            subprocess.run(
                                f"mesh2mesh -v -m1 refMesh.mg.c3d -m2 Mesh.mg.c3d -q1 {self._c3d_checkpoint_rename} -q2 Restart.file",
                                shell=True,
                                stdout=f,
                                stderr=subprocess.STDOUT,
                            )

                    _start = time.time()
                    with open(self.c3d_logname, "a") as f:
                        subprocess.run(
                            run_cmd, shell=True, stdout=f, stderr=subprocess.STDOUT
                        )
                    while not os.path.exists(c3d_donefile):
                        # Wait...
                        time.sleep(5)

                        # Check for C3D failure
                        running, e = self._c3d_running()

                        if not running:
                            # C3D failed, try restart it
                            if _restarts > 3:
                                print("Too many Cart3D failures... Something is wrong.")
                                return False

                            print(f"\033[1mERROR\033[0m: Cart3D failed with error {e}")
                            print("  Restarting Cart3D.")
                            f = open(self.c3d_logname, "a")
                            subprocess.run(
                                run_cmd, shell=True, stdout=f, stderr=subprocess.STDOUT
                            )
                            f.close()
                            _restarts += 1

                    _end = time.time()
                    print(f"Cart3D simulations complete in {(_end-_start):.2f} s.")

                else:
                    # Cart3D already finished for this iteration
                    print("Cart3D DONE file located.")

                return True

            else:
                if attempt < N - 1:
                    print("Could not intersect components. Trying again.")
                else:
                    print("Could not intersect components. Exiting.")
                    return False

    def _process_results(self, parameters: Dict[str, float], iter_dir: str):
        # Construct simulation directory
        sim_dir = os.path.join(iter_dir, self.sim_dir_name)
        target_adapt = self._infer_adapt(sim_dir)

        # Extract drag coefficient of geometry
        if target_adapt is not None:
            # Adaptation simulation
            loads_filepath = os.path.join(sim_dir, target_adapt, "FLOW", "loadsCC.dat")
            components_filepath = os.path.join(
                sim_dir, target_adapt, "FLOW/Components.i.plt"
            )
        else:
            # Possibly warm-started simulation
            loads_filepath = os.path.join(sim_dir, "loadsCC.dat")
            components_filepath = os.path.join(sim_dir, "Components.i.plt")
        loads_dict = self._read_c3d_loads(loads_filepath)

        # Approximate flow sensitivities
        f_sense_filename = os.path.join(iter_dir, self.f_sense_filename)
        if not os.path.exists(f_sense_filename):
            # Filepaths
            sensitivity_filepath = os.path.join(iter_dir, self.sensitivity_filename)

            # Create PySAGAS wrapper and run
            print("\nEvaluating sensitivities.")
            try:
                wrapper = Cart3DWrapper(
                    a_inf=self.a_inf,
                    rho_inf=self.rho_inf,
                    sensitivity_filepath=sensitivity_filepath,
                    components_filepath=components_filepath,
                    verbosity=0,
                )

            except ValueError:
                # The sensitivity data does not match the point data, regenerate it
                tri_components_filepath = os.path.join(sim_dir, "Components.i.tri")
                sensitivity_files = glob.glob(os.path.join(iter_dir, "*sensitivity*"))
                self._combine_sense_data(
                    tri_components_filepath,
                    sensitivity_files=sensitivity_files,
                    match_target=self._matching_target,
                    tol_0=self._matching_tolerance,
                    max_tol=self._max_matching_tol,
                    outdir=iter_dir,
                )

                # Re-instantiate the wrapper
                wrapper = Cart3DWrapper(
                    a_inf=self.a_inf,
                    rho_inf=self.rho_inf,
                    sensitivity_filepath=sensitivity_filepath,
                    components_filepath=components_filepath,
                    verbosity=0,
                )

            F_sense, _ = wrapper.calculate()
            print("  Done.")

            # Save F_sense
            F_sense.to_csv(f_sense_filename)

        else:
            # Load F_sense
            F_sense = pd.read_csv(f_sense_filename, index_col=0)
            print("Force sensitivities loaded from file.")

        # Non-dimensionalise
        coef_sens = F_sense / (0.5 * self.rho_inf * self.A_ref * self.V_inf**2)

        # Load parameter values
        x_df = pd.read_csv(
            os.path.join(iter_dir, self.parameters_filename), index_col=0
        )
        x = x_df.loc[parameters.keys()]["0"].values

        # Get objective function and Jacobian
        jacobian_filepath = os.path.join(iter_dir, self.jacobian_filename)
        if not os.path.exists(jacobian_filepath):
            # Load data
            properties_dir = glob.glob(f"{iter_dir}{os.sep}*_properties")
            if properties_dir:
                scalar_sens_dir = os.path.join(iter_dir, "scalar_sensitivities")
                vm = pd.read_csv(
                    glob.glob(os.path.join(properties_dir[0], "*volmass.csv"))[0],
                    index_col=0,
                )["0"]
                vm_sens = pd.read_csv(
                    os.path.join(scalar_sens_dir, "volmass_sensitivity.csv"),
                    index_col=0,
                )[parameters.keys()]
            else:
                # No properties data found
                vm = None
                vm_sens = None

            # Call function
            obj, jac_df = self._obj_jac_cb(
                parameters=parameters,
                coef_sens=coef_sens,
                loads_dict=loads_dict,
                volmass=vm,
                volmass_sens=vm_sens,
            )

            # TODO - dimensionality check?

            # Save
            with open(os.path.join(iter_dir, self.objective_filename), "w") as f:
                f.write(f"objective: {obj}\n")
            jac_df.to_csv(jacobian_filepath)

        else:
            # Load from file
            with open(os.path.join(iter_dir, self.objective_filename), "r") as f:
                lines = f.readlines()
                obj = float(lines[0].strip().split(":")[-1])
            jac_df = pd.read_csv(jacobian_filepath, index_col=0)["0"]

        # Extract ordered jacobian values
        jac = jac_df.loc[parameters.keys()].values

        # TODO - also return step size?

        return obj, jac, x

    def _read_c3d_loads(
        self,
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

    def _infer_adapt(self, sim_dir) -> str:
        try:
            with open(f"{sim_dir}/aero.csh", "r") as f:
                lines = f.readlines()

                for line in lines:
                    if line.find("set n_adapt_cycles") != -1:
                        return f"adapt{int(line.split('=')[-1]):02d}"
        except FileNotFoundError:
            return None

    def _overwrite_adapt(self, sim_dir: str, max_adapt: int) -> str:
        """Overwrites the adapt cycle in the iteration directory."""
        original = os.path.join(sim_dir, "aero.csh")
        new = os.path.join(sim_dir, "temp_aero.csh")
        with open(new, "w+") as new_file:
            with open(original, "r") as original_file:
                for line in original_file:
                    if line.startswith("set n_adapt_cycles"):
                        # This is the line setting the adapt number
                        line = f"{line.split('=')[0]}= {max_adapt}"

                    # Write line to new file
                    new_file.write(line)

        # Replace original aero.csh file with updated file
        shutil.copymode(original, new)
        os.remove(original)
        os.rename(new, original)

    @staticmethod
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

        print("Component sensitivity data combined successfully.")

    def _c3d_running(self) -> bool:
        with open(self.c3d_logname) as f:
            # Get last line in log file
            for line in f:
                pass

            # Check if if it is in the known errors
            for e in self.C3D_errors:
                if e in line:
                    return False, e

        # No errors
        return True, None
