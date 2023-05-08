import os
import time
import shutil
import subprocess
import numpy as np
from typing import Dict, List, Optional
from pysagas import Cell, FlowState, Vector
from pysagas.optimisation.cart3d.utilities import C3DPrep
from pysagas.cfd.solver import FlowSolver, FlowResults, SensitivityResults


class Cart3D(FlowSolver):
    method = "Cart3D"

    _C3D_errors = [
        "==> ADAPT failed",
        "Check cart3d.out in AD_A_J for more clues",
        "==> adjointErrorEst_quad failed again, status = 1",
        "ERROR: CUBES failed",
        "ERROR: ADAPT failed with status = 1",
        "ERROR",
    ]

    def __init__(
        self,
        stl_files: List[str],
        aero_csh: str = "aero.csh",
        input_cntl: str = "input.cntl",
        freestream: Optional[FlowState] = None,
        verbosity: Optional[int] = 1,
        c3d_log_name: str = "C3D_log",
        c3d_info_file: str = None,
        write_config_xml: bool = True,
        ref_area: Optional[float] = 1.0,
        ref_length: Optional[float] = 1.0,
    ) -> None:
        """Instantiate the Cart3D solver wrapper.

        Parameters
        -----------
        stl_files : list[str]
            A list of filepaths to the STL files to intersect and
            simulate.
        aero_csh : str
            The filepath to the reference aero.csh file. The default
            option will look in the current working directory.

        c3d_log_name : str, optional
            The name to use for the Cart3D logfile. The default is C3D_log.

        write_config_xml : bool, optional
            A boolean flag to write the Cart3D Config.xml file when running
            comp2tri. If the geometry may be perturbed to a state where one
            component is no longer part of the wetted surface, writing the
            comp2tri file can cause set-up issues, and so it can be beneficial
            to turn it off. The default is True.
        """
        # Save STL files for processing
        self.stl_files = stl_files
        self._root_dir = os.getcwd()
        self._aerocsh_path = os.path.join(self._root_dir, aero_csh)
        self._inputcntl_path = os.path.join(self._root_dir, input_cntl)

        # Instantiate helper
        self.c3d_log_name = c3d_log_name
        self._prepper = C3DPrep(
            logfile=c3d_log_name,
            info_file=c3d_info_file,
            write_config=write_config_xml,
        )

        # Dimensional properties
        self.ref_area = ref_area
        self.ref_length = ref_length

        # Complete instantiation
        super().__init__(None, freestream, verbosity)

    def solve(
        self,
        freestream: Optional[FlowState] = None,
        mach: Optional[float] = None,
        aoa: Optional[float] = None,
    ) -> FlowResults:
        already_run = super().solve(freestream=freestream, mach=mach, aoa=aoa)
        if already_run:
            # Already have a result
            if self.verbosity > 1:
                print("Attention: this flowstate has already been solved.")
            result = self.flow_result

        else:
            # Get flow
            flow = self._last_solve_freestream
            loads = self._run_sim(mach=flow.M, aoa=flow.aoa)

            # TODO - need to be careful and aware of the reference coordinates here
            C_F = np.array(
                [loads["C_A-entire"], loads["C_N-entire"], loads["C_Y-entire"]]
            )
            C_M = np.array(
                [loads["C_M_x-entire"], loads["C_M_z-entire"], -loads["C_M_y-entire"]]
            )

            # Dimensionalise results
            net_force = C_F * flow.q * self.ref_area
            net_moment = C_M * flow.q * self.ref_area * self.ref_length

            # Convert to vectors
            net_force = Vector.from_coordinates(net_force)
            net_moment = Vector.from_coordinates(net_moment)

            # Construct results
            result = FlowResults(
                freestream=flow, net_force=net_force, net_moment=net_moment
            )

            # Save
            self.flow_result = result

        return result

    def solve_sens(
        self,
        freestream: Optional[FlowState] = None,
        Mach: Optional[float] = None,
        aoa: Optional[float] = None,
    ) -> SensitivityResults:
        raise NotImplementedError("Coming soon.")

    def save(self, name: str, attributes: Dict[str, list]):
        raise NotImplementedError("Coming soon.")

    def _prepare_sim_dir(self, sim_dir_name: dir):
        """Prepares a new directory to run the simulation in."""
        # Make simulation directory
        sim_dir = os.path.join(self._root_dir, sim_dir_name)
        run_intersect = False
        components_filepath = os.path.join(sim_dir, "Components.i.tri")
        if not os.path.exists(sim_dir):
            os.mkdir(sim_dir)
            run_intersect = True
            intersected = False
        else:
            # Check for intersected file
            run_intersect = not os.path.exists(components_filepath)
            intersected = True

        if self.verbosity > 0:
            print(f"Running simulation in {sim_dir}.")

        return run_intersect, intersected

    def _run_sim(self, mach: float, aoa: float):
        """Run the simulation."""
        # Prepare sim directory
        sim_dir = os.path.join(self._root_dir, f"M{mach}A{aoa}")
        run_intersect, intersected = self._prepare_sim_dir(sim_dir)

        # Attempt simulation
        sim_success = False
        no_attempts = 3
        for attempt in range(no_attempts):
            # Run intersect
            if run_intersect:
                self._prepper._log(f"SHAPEOPT INTERSECT ATTEMPT {attempt+1}")
                intersected = self._prepper.intersect_stls(stl_files=self.stl_files)

            # Check for intersection
            if intersected:
                # Prepare rest of simulation directory
                if not os.path.exists(os.path.join(sim_dir, "input.cntl")):
                    # Move files to simulation directory (including Components.i.tri)
                    self._prepper.run_autoinputs()
                    os.system(
                        f"mv *.tri Config.xml input.c3d preSpec.c3d.cntl {sim_dir} >> {self.c3d_log_name} 2>&1"
                    )

                    # Copy sim files and permissions
                    for filename, fp in {
                        "input.cntl": self._inputcntl_path,
                        "aero.csh": self._aerocsh_path,
                    }.items():
                        shutil.copyfile(
                            fp,
                            os.path.join(sim_dir, filename),
                        )
                        shutil.copymode(
                            fp,
                            os.path.join(sim_dir, filename),
                        )

                    # Update flow conditions
                    new_input_cntl = os.path.join(sim_dir, "input.cntl")
                    self._edit_input_cntl(
                        mach=mach,
                        aoa=aoa,
                        original=self._inputcntl_path,
                        new=new_input_cntl,
                    )

                # Run Cart3D and await result
                os.chdir(sim_dir)
                target_adapt = self._infer_adapt(self._aerocsh_path)
                c3d_donefile = os.path.join(sim_dir, target_adapt, "FLOW", "DONE")
                run_cmd = "./aero.csh restart"
                _restarts = 0
                if not os.path.exists(c3d_donefile):
                    with open(self.c3d_log_name, "a") as f:
                        subprocess.run(
                            run_cmd, shell=True, stdout=f, stderr=subprocess.STDOUT
                        )
                    while not os.path.exists(c3d_donefile):
                        # Wait...
                        time.sleep(5)

                        # Check for C3D failure
                        running, e = self._c3d_running(c3d_logname=self.c3d_log_name)

                        if not running:
                            # C3D failed, try restart it
                            if _restarts > 3:
                                return False

                            f = open(self.c3d_log_name, "a")
                            subprocess.run(
                                run_cmd, shell=True, stdout=f, stderr=subprocess.STDOUT
                            )
                            f.close()
                            _restarts += 1

                sim_success = True
                break

        if sim_success:
            # Sim finished successfully, read loads file
            loads_dict = self._read_c3d_loads(
                os.path.join(sim_dir, "BEST/FLOW/loadsCC.dat"),
                tag="entire",
            )

        else:
            loads_dict = None

        # Change back to working directory
        os.chdir(self._root_dir)

        return loads_dict

    @staticmethod
    def _c3d_running(c3d_log_name) -> bool:
        """Watches the Cart3D log file to check for errors and return False
        if Cart3D has stopped running.

        Returns
        -------
        running : bool
            Whether Cart3D is still running.

        error : str
            The error message, if an error has occured.
        """
        with open(c3d_log_name) as f:
            # Get last line in log file
            for line in f:
                pass

            # Check if if it is in the known errors
            for e in Cart3D._C3D_errors:
                if e in line:
                    return False, e

        # No errors
        return True, None

    @staticmethod
    def _infer_adapt(aero_csh_fp: str) -> str:
        """Returns the target adapt number by reading the aero.csh file."""
        with open(aero_csh_fp, "r") as f:
            lines = f.readlines()

            for line in lines:
                if line.find("set n_adapt_cycles") != -1:
                    return f"adapt{int(line.split('=')[-1]):02d}"

    @staticmethod
    def _edit_input_cntl(mach: float, aoa: float, original: str, new: str):
        """Edits the Mach number and angle of attack in the input.cntl file."""
        with open(new, "w+") as new_file:
            with open(original, "r") as old_file:
                for line in old_file:
                    if line.startswith("Mach"):
                        line = f"Mach  {mach}\n"
                    elif line.startswith("alpha"):
                        line = f"alpha  {aoa}\n"

                    # Write line
                    new_file.write(line)

    @staticmethod
    def _read_c3d_loads(
        loadsCC_filepath: str,
        b_frame: bool = True,
        v_frame: bool = True,
        moments: bool = True,
        tag: str = None,
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
                _tag = word[0]

                if (not tag) or (tag is not None and tag == _tag):
                    coeff = word[-1]
                    short_coeff = coeff[1:4]
                    long_coeff = coeff[1:6]
                    if b_frame is True:
                        if short_coeff in ["C_A", "C_Y", "C_N"]:  # get force in b_frame
                            load_dict["{0}-{1}".format(short_coeff, tag)] = number
                    if v_frame is True:
                        if short_coeff in ["C_D", "C_S", "C_L"]:  # get force in v_frame
                            load_dict["{0}-{1}".format(short_coeff, tag)] = number
                    if moments is True:
                        if short_coeff in [
                            "C_l",
                            "C_m",
                            "C_n",
                            "C_M",
                        ]:  # get moment coeff
                            load_dict["{0}-{1}".format(short_coeff, tag)] = number
                        if b_frame is True:
                            if long_coeff in [
                                "C_M_x",
                                "C_M_y",
                                "C_M_z",
                            ]:  # get force in b_frame
                                load_dict["{0}-{1}".format(long_coeff, tag)] = number

        return load_dict
