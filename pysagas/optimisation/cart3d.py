import os
import time
import glob
import numpy as np
import pandas as pd
from random import random
from pysagas import banner
from typing import List, Dict
import matplotlib.pyplot as plt
from pysagas.wrappers import Cart3DWrapper
from hypervehicle.generator import Generator
from hypervehicle.utilities import SensitivityStudy, append_sensitivities_to_tri


np.seterr(all="ignore")


class ShapeOpt:
    """A wrapper to perform shape optimisation with Cart3D."""

    C3D_errors = ["==> ADAPT failed"]

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
    ) -> None:

        # Construct paths
        self.root_dir = home_dir
        self.basefiles_dir = os.path.join(home_dir, basefiles_dir_name)
        self.working_dir = os.path.join(home_dir, working_dir_name)
        self.sim_dir_name = sim_dir_name
        self.sensitivity_filename = sensitivity_filename
        self.completion_filename = "ITERATION_COMPLETE"
        self.parameters_filename = "parameters.csv"
        self.jacobian_filename = "jacobian.csv"

        self.generator = generator

        self.c3d_logname = c3d_logname
        self.loads_key = None

        # TODO - pass as flowstate
        self.rho_inf = rho_inf
        self.V_inf = V_inf
        self.a_inf = a_inf
        self.A_ref = A_ref

        # Create instance of Cart3D prepper
        self._c3dprepper = _C3DPrep(logfile=c3d_logname)

    def _prepare(self, warmstart: bool, param_names: List[str]):
        """Prepares the working directory for the optimisation
        problem.
        """
        # Initialise 'older' results
        x_older = None
        jac_older = None

        # Check base files directory exist
        if not os.path.exists(self.basefiles_dir):
            raise Exception("Cart3D base file directory does not exist!")

        # Create working directory
        if not os.path.exists(self.working_dir):
            os.mkdir(self.working_dir)

        # Create heirarchy for this run
        iteration_dirs = [
            int(i)
            for i in os.listdir(self.working_dir)
            if os.path.isdir(os.path.join(self.working_dir, i))
        ]
        if len(iteration_dirs) > 0:
            # A previous iteration has run, check for completion file
            if os.path.exists(
                os.path.join(
                    self.working_dir,
                    f"{int(max(iteration_dirs)):04d}",
                    self.completion_filename,
                )
            ):
                # This iteration ran to completion
                if warmstart:
                    # Warmstarting, load results from this iteration
                    current_iter = max(iteration_dirs)
                    print(f"\n\x1B[3mWarmstarting from iteration {current_iter}\x1B[0m")

                    # Look for x_older
                    prev_iter_dir = os.path.join(
                        self.working_dir, f"{int(current_iter-1):04d}"
                    )
                    x_older_path = os.path.join(prev_iter_dir, self.parameters_filename)
                    jac_older_path = os.path.join(prev_iter_dir, self.jacobian_filename)
                    if os.path.exists(x_older_path):
                        # Older iteration directory exists, pick up x_older
                        x_older = pd.read_csv(x_older_path, index_col=0).values

                        # Also load previous jacobian
                        F_sense = pd.read_csv(jac_older_path, index_col=0)
                        coef_sens = F_sense / (
                            0.5 * self.rho_inf * self.A_ref * self.V_inf**2
                        )
                        jac_older = coef_sens.loc[param_names]["dFx/dP"].values

                else:
                    # Not warmstarting, and this iteration completed: move to next
                    # But if we are here, we should be warmstarting...
                    # Not unless we've just picked up the previous...
                    # raise Exception("Old iterations detected. Either delete them, "+\
                    #     "or resume by setting warmstart to True.")
                    current_iter = max(iteration_dirs) + 1
                    print(f"\n\x1B[3mMoving onto iteration {current_iter}\x1B[0m")

            else:
                # This iteration did not complete, try resume it
                current_iter = max(iteration_dirs)
                print(f"\n\x1B[3mResuming from iteration {current_iter}\x1B[0m")

        else:
            # First iteration
            current_iter = 0

        # Print iteration header
        its = f"Iteration {int(current_iter)}".center(43, " ")
        print(f"{'':=>43}\n{its}\n{'':=>43}")

        # Define the current iteration directory
        iter_dir = os.path.join(self.working_dir, f"{int(current_iter):04d}")
        if not os.path.exists(iter_dir):
            # Create the directory
            os.mkdir(iter_dir)

        return iter_dir, x_older, jac_older

    def _run_sensitivity_study(
        self, iter_dir: str, param_names: List[str], x: List[float]
    ):
        # Change into iteration directory
        os.chdir(iter_dir)

        # Check if sensitivity study has been run
        sens_files = glob.glob(f"{iter_dir}/*sensitivity*")
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

    def _run_simulation(self, basefiles_dir: str, iter_dir: str):

        target_adapt = self._infer_adapt()

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

        # Run intersect
        if run_intersect:
            intersected = self._c3dprepper.intersect_stls()

        # Check for intersection
        if intersected:
            # Prepare rest of sim
            if not os.path.exists(os.path.join(sim_dir, "aero.csh")):
                # Prepare remaining C3D files
                os.system(f"autoInputs -r 2 >> {self.c3d_logname} 2>&1")

                # Move files to simulation directory
                os.system(
                    f"mv *.tri Config.xml input.c3d preSpec.c3d.cntl {sim_dir} >> {self.c3d_logname} 2>&1"
                )

                # Copy sim files
                os.system(
                    f"cp {basefiles_dir}/input.cntl {basefiles_dir}/aero.csh {sim_dir} >> {self.c3d_logname} 2>&1"
                )

            # Create all_components_sensitivity.csv
            if not os.path.exists(self.sensitivity_filename):
                self._combine_sense_data(components_filepath)

            # Run Cart3D and await result
            os.chdir(sim_dir)
            c3d_donefile = os.path.join(sim_dir, target_adapt, "FLOW", "DONE")
            if not os.path.exists(c3d_donefile):
                print(
                    "\nStarting Cart3D, awaiting",
                    os.sep.join(c3d_donefile.split(os.sep)[-6:]),
                )

                os.system(f"./aero.csh restart >> {self.c3d_logname} 2>&1")
                while not os.path.exists(c3d_donefile):
                    # Wait...
                    time.sleep(5)

                    # Check for C3D failure
                    running, e = self._c3d_running()

                    if not running:
                        # C3D failed, try restart it
                        print(f"\033[1mERROR\033[0m: Cart3D failed with error {e}")
                        os.system(f"./aero.csh restart >> {self.c3d_logname} 2>&1")

            print("Cart3D simulations complete.")

            complete = True

        else:
            print("Could not intersect components.")
            complete = False

        return complete

    def _process_results(self, param_names: List[str], iter_dir: str):

        target_adapt = self._infer_adapt()

        # Construct simulation directory
        sim_dir = os.path.join(iter_dir, self.sim_dir_name)

        # Extract drag coefficient of geometry
        loads_filepath = os.path.join(sim_dir, target_adapt, "FLOW", "loadsCC.dat")
        loads_dict = self._read_c3d_loads(loads_filepath)

        # Approximate flow sensitivities
        jacobian_filepath = os.path.join(iter_dir, self.jacobian_filename)
        if not os.path.exists(jacobian_filepath):

            # Filepaths
            sensitivity_filepath = os.path.join(iter_dir, self.sensitivity_filename)
            components_filepath = os.path.join(
                sim_dir, target_adapt, "FLOW/Components.i.plt"
            )

            # Create PySAGAS wrapper and run
            wrapper = Cart3DWrapper(
                a_inf=self.a_inf,
                rho_inf=self.rho_inf,
                sensitivity_filepath=sensitivity_filepath,
                components_filepath=components_filepath,
                verbosity=0,
            )
            F_sense = wrapper.calculate()

            # Save Jacobian
            F_sense.to_csv(jacobian_filepath)

        else:
            # Load Jacobian
            F_sense = pd.read_csv(jacobian_filepath, index_col=0)
            print("Jacobian loaded from file.")

        # Non-dimensionalise
        coef_sens = F_sense / (0.5 * self.rho_inf * self.A_ref * self.V_inf**2)

        # Load parameter values
        x_df = pd.read_csv(
            os.path.join(iter_dir, self.parameters_filename), index_col=0
        )

        # Construct output (note sorting of coef_sens!)
        # TODO - search for "dFx/dP" (jac_older) - this needs to be dynamic
        obj = loads_dict[self.loads_key]
        jac = coef_sens.loc[param_names]["dFx/dP"].values
        x = x_df.loc[param_names]["0"].values

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

    def _iterate(
        self, x: List[float], param_names: List[str], warmstart: bool, gamma: float
    ):
        """Wrapper function to perform an iteration of the
        shape optimisation problem.

        Steps involved:
            Generating geometry and geometric sensitivities
            Running Cart3D
            Approximating the gradient information
            Return the objective function and gradient.
        """
        # Prepare this iteration
        iter_dir, x_older, jac_older = self._prepare(warmstart, param_names)

        # Run the sensitivity study
        self._run_sensitivity_study(iter_dir, param_names, x)

        # Run simulation
        success = self._run_simulation(self.basefiles_dir, iter_dir)

        if success:
            # Simulation completed successfully
            obj, jac, x = self._process_results(param_names, iter_dir)

            # Create completion file
            pd.Series(
                {"objective": obj, "gamma": gamma, **dict(zip(param_names, x))}
            ).to_csv(os.path.join(iter_dir, self.completion_filename))

        else:
            raise Exception("Simulation failed.")

        # Change back to root dir
        os.chdir(self.root_dir)

        return obj, jac, x, x_older, jac_older

    def _gradient_search(self, parameters: Dict[str, float], warmstart: bool = False):
        """Performs a steepest descent search.

        Parameters
        ----------
        parameters: Dict[str, float]
            A dictionary of geometric parameters to pass to the vehicle generator.
        warmstart : bool, optional
            If you are resuming a previous run, set to True. This will accelerate
            convergence by improving the step size. The default is Fale.
        """

        param_names = list(parameters.keys())
        x0 = list(parameters.values())

        # Define initial step size
        gamma = 0.05

        # Constrain iterations
        max_iterations = 10

        # Iteration parameters
        i = 0
        tolerance = 1e-3
        change = 2 * tolerance
        obj_prev = 10 * tolerance
        bailout = False

        while change > tolerance:
            if i + 1 > max_iterations:
                # Exit now
                bailout = True
                break

            # Start timer
            _start = time.time()

            # Get objective and jacobian
            obj, jac, x_old, x_older, jac_older = self._iterate(
                x=x0,
                param_names=param_names,
                warmstart=warmstart,
                gamma=gamma,
            )

            # Check for non-zero Jacobian
            if np.linalg.norm(jac) == 0:
                print("\033[1mERROR\033[0m: Exiting due to zero-norm Jacobian.")
                bailout = True
                break

            # Calculate step size
            if x_older is not None:
                _gamma = (
                    np.linalg.norm((x_old - x_older) * (jac - jac_older))
                    / np.linalg.norm(jac - jac_older) ** 2
                )

                # Correct for nan
                if not np.isnan(_gamma):
                    # Update gamma
                    gamma = _gamma

            # Update x0
            x0 = x_old - gamma * jac

            # Calculate change
            change = abs((obj - obj_prev) / obj_prev)
            # TODO - detect divergence?

            # Break out of warmstart routine
            warmstart = False

            # Update iteration
            _end = time.time()
            i += 1
            obj_prev = obj
            jac_older = jac
            x_older = x_old

            # Print Information
            print("\nIteration complete:")
            print(f"Time to complete: {(_end-_start):.2f} s")
            print("Objective function:", obj)
            print("Step size:", gamma)
            print("New guess for next iteration:")
            print(pd.Series(x0, index=param_names, name="Parameters").to_string())
            print("")

        # Finished
        if not bailout:
            print(f"\nExited with change = {change}")

    def optimise(
        self,
        parameters: Dict[str, float],
        loads_key: str = "C_D-entire",
        warmstart: bool = True,
    ):
        """Performs a steepest descent search.

        Parameters
        ----------
        parameters: Dict[str, float]
            A dictionary of geometric parameters to pass to the vehicle generator.
        loads_key : str, optional
            The key to use for extracting the objective from the Cart3D loads file.
            The default is 'C_D-entire'.
        warmstart : bool, optional
            If you are resuming a previous run, set to True. This will accelerate
            convergence by improving the step size. The default is True.
        """
        # Print banner
        banner()
        print("\033[4mCart3D Shape Optimisation\033[0m".center(50, " "))

        # Save loadsCC key for objective function
        self.loads_key = loads_key

        # Run
        _opt_start = time.time()
        try:
            self._gradient_search(parameters=parameters, warmstart=warmstart)
        except KeyboardInterrupt:
            # Change back to root dir and exit
            os.chdir(self.root_dir)
        except Exception as e:
            os.chdir(self.root_dir)
            raise Exception(e)

        _opt_end = time.time()
        print(f"\nTotal run time: {(_opt_end-_opt_start):.2f}")

    def post_process(self, plot_convergence: bool = True) -> pd.DataFrame:
        """Crawls through iteration directories to compile results."""
        iteration_dirs = [
            i
            for i in os.listdir(self.working_dir)
            if os.path.isdir(os.path.join(self.working_dir, i))
        ]

        results = []
        for directory in iteration_dirs:
            iteration = int(directory)

            completion_file = os.path.join(
                self.working_dir, directory, self.completion_filename
            )
            if os.path.exists(completion_file):
                # This iteration completed, load the results
                iter_result = pd.read_csv(completion_file)

                # Add iteration number
                iter_result.loc[len(iter_result)] = ["iteration", iteration]

                # Append
                results.append(iter_result.set_index("Unnamed: 0").to_dict()["0"])

        df = pd.DataFrame(results).set_index("iteration").sort_index()

        if plot_convergence:
            plt.plot(df.index, df["objective"])
            plt.title("Convergence of Objective Function")
            plt.xlabel("Iteration")
            plt.ylabel("Objective Function")

            plt.grid()
            plt.show()

        return df

    def _infer_adapt(self) -> str:
        with open(f"{self.basefiles_dir}/aero.csh", "r") as f:
            lines = f.readlines()

            for line in lines:
                if line.find("set n_adapt_cycles") != -1:
                    return f"adapt{int(line.split('=')[-1]):02d}"

    def _combine_sense_data(
        self, components_filepath: str, match_target: float = 0.9, max_tol: float = 1e-1
    ):
        """Combine the component sensitivity data for intersected geometry."""
        match_frac = 0
        tol = 1e-5
        while match_frac < match_target:
            # Run matching algorithm
            match_frac = append_sensitivities_to_tri(
                dp_filenames=glob.glob("*sensitivity*"),
                components_filepath=components_filepath,
                match_tolerance=tol,
                verbosity=0,
            )

            # Reduce matching tolerance
            tol *= 10

            # Check new tolerance
            if tol > max_tol:
                raise Exception("Cannot combine sensitivity data.")

            if match_frac < match_target:
                print("Failed to combine sensitivity data.")
                print("  Reducing matching tolerance and trying again.")

    def _c3d_running(self) -> bool:
        with open(self.c3d_logname) as f:
            # Get last line in log file
            for line in f:
                pass

            # Check if if it is in the known errors
            for e in ShapeOpt.C3D_errors:
                if e in line:
                    return False, e

        # No errors
        return True, None


class _C3DPrep:
    def __init__(self, logfile) -> None:
        self.logfile = logfile

    def _run_stl2tri(self, stl_files: list):
        tri_files = []
        for file in stl_files:
            prefix = file.split(".")[0]
            tri_file = prefix + ".tri"
            os.system(f"stl2tri.pl {file} {tri_file} >> {self.logfile} 2>&1")
            tri_files.append(tri_file)

        os.system(f"rm *.comp.tri *.off >> {self.logfile} 2>&1")
        return tri_files

    def _jitter_tri_files(self, tri_files):
        for file in tri_files:
            prefix = file.split(".")[0]
            denom = 1000  # Max of 0.0001, min of 0
            x_pert = random() / denom
            y_pert = random() / denom
            z_pert = random() / denom
            os.system(
                f"trix -x {x_pert} -y {y_pert} -z {z_pert} -o {prefix} {file} >> {self.logfile} 2>&1"
            )

    def _shift_all(
        self,
        tri_files,
        x_shift,
        y_shift,
        z_shift,
        component: str = None,
        reverse: bool = False,
    ):
        if component is None:
            transform_files = tri_files
        else:
            transform_files = [component]

        for file in transform_files:
            prefix = file.split(".")[0]
            if reverse:
                os.system(
                    f"trix -x {-x_shift} -y {-y_shift} -z {-z_shift} -o {prefix} {file} >> {self.logfile} 2>&1"
                )
            else:
                os.system(
                    f"trix -x {x_shift} -y {y_shift} -z {z_shift} -o {prefix} {file} >> {self.logfile} 2>&1"
                )

    def _rotate_all(
        self,
        tri_files,
        x_rot,
        y_rot,
        z_rot,
        component: str = None,
        reverse: bool = False,
    ):
        if component is None:
            transform_files = tri_files
        else:
            transform_files = [component]

        for file in transform_files:
            prefix = file.split(".")[0]
            if reverse:
                order = ["z", "y", "x"]
            else:
                order = ["x", "y", "z"]

            for axis in order:
                rotation = vars()[f"{axis}_rot"]
                os.system(f"trix -r{axis} {rotation} -o {prefix} {file}")

    def _run_comp2tri(self, tri_files):
        tri_files_str = " ".join(tri_files)
        os.system(
            f"comp2tri -makeGMPtags {tri_files_str} -config >> {self.logfile} 2>&1"
        )

    def _run_intersect(self):
        os.system(f"intersect >> {self.logfile} 2>&1")

    @staticmethod
    def _get_stl_files():
        path = os.getcwd()
        all_files = [
            f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))
        ]
        stl_files = []
        for file in all_files:
            if file.split(".")[-1] == "stl":
                stl_files.append(file)

        # Sort files
        stl_files.sort()

        return stl_files

    @staticmethod
    def _check_for_success():
        success = False
        if os.path.isfile("Components.i.tri"):
            success = True
        return success

    def intersect_stls(self) -> bool:
        # Check for existing intersected file
        if self._check_for_success():
            return True

        # Continue
        stl_files = self._get_stl_files()
        tri_files = self._run_stl2tri(stl_files)

        # First try intersect original files
        self._run_comp2tri(tri_files)
        self._run_intersect()
        successful = self._check_for_success()
        if successful:
            return True

        # First attempt failed, try jittering components
        self._jitter_tri_files(tri_files)
        self._run_comp2tri(tri_files)
        self._run_intersect()
        successful = self._check_for_success()
        if successful:
            return True

        # That also failed, try arbitrary shift away
        for attempt in range(3):
            # Define shifts
            x_shift = random() * 10  # Max of 10, min of 0
            y_shift = random() * 10  # Max of 10, min of 0
            z_shift = random() * 10  # Max of 10, min of 0

            # Define rotations
            x_rot = random() * 10  # Max of 10, min of 0
            y_rot = random() * 10  # Max of 10, min of 0
            z_rot = random() * 10  # Max of 10, min of 0

            # Apply transformations
            self._shift_all(tri_files, x_shift, y_shift, z_shift)
            self._rotate_all(tri_files, x_rot, y_rot, z_rot)

            # Make attempt
            self._run_comp2tri(tri_files)
            self._run_intersect()
            successful = self._check_for_success()

            if successful:
                # Move configuration back to original location
                self._shift_all(
                    tri_files, x_shift, y_shift, z_shift, "Components.i.tri", True
                )
                self._rotate_all(
                    tri_files, x_rot, y_rot, z_rot, "Components.i.tri", True
                )
                if successful:
                    return True
            else:
                # Need to reset tri files
                tri_files = self._run_stl2tri(stl_files)

        return False
