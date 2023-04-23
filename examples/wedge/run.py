import numpy as np
from pysagas.cfd import OPM
from pysagas.flow import FlowState
from pysagas.geometry import Vector
from pysagas.geometry.parsers import PyMesh, TRI


# Load STL file
cells = PyMesh.load_from_file("wedge.stl")
A_ref = 1

# Instantiate solver
freestream = FlowState(mach=8, pressure=101e3, temperature=288)
solver = OPM(cells, freestream)

# Run solver
result = solver.solve()

# Can also evaluate at other flow conditions, with:
# result = solver.solve(aoa=5, Mach=6) # Modified AoA/Mach
# result = solver.solve(freestream=new_freestream) # Updated freestream

# Save results to VTK file for visualisation of solution
solver.save("wedge")
