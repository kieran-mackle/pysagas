<a name="readme-top"></a>

# PySAGAS

<!-- start intro -->
**PySAGAS** is a Python package for the generation of **S**ensitivity **A**pproximations
for **G**eometric and **A**erodynamic **Surface** properties.
It provides a computationally-efficient method for generating 
surface sensitivity approximations from existing flow solutions,
to use in aerodynamic shape optimisation studies.
<!-- end intro -->

<a><img src="https://github.com/0x6080604052/analytics/actions/workflows/tests.yml/badge.svg" alt="Test Status"></a>



<!-- TABLE OF CONTENTS -->
<details>
  <summary><h2>Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
  </ol>
</details>



## Getting Started
<!-- start getting started -->

### Prerequisites

#### GDTK
*PySAGAS* depends on the [GDTK](https://github.com/gdtk-uq/gdtk) python 
package for ideal gas dynamics. Note that a full install is 
not required. Instead, do a 
[sparse checkout](https://stackoverflow.com/questions/600079/how-do-i-clone-a-subdirectory-only-of-a-git-repository)
of the relevant files, using the commands below.

```
mkdir gdtk
cd gdtk/
git init
git remote add -f origin https://github.com/gdtk-uq/gdtk.git
git config core.sparseCheckout true
echo "src/lib/" >> .git/info/sparse-checkout
git pull origin master
cd src/lib
python3 -m pip install .
cd ../../../
```

#### ParaView
If you plan on using *PySAGAS* with Cart3D solutions, you should
also install ParaView, or at least the ParaView Python bindings.
If you are using an Anaconda environment, you can install the 
ParaView Python packages via
[Conda Forge](https://anaconda.org/conda-forge/paraview) using 
the command below:

```
conda install -c conda-forge paraview
```

If you already have ParaView installed, you can append the path
to the binaries to the Python path using the snippet below.

```python
import sys
sys.path.insert(0, "/opt/ParaView-5.6.2-MPI-Linux-64bit/bin")

# Now the import should work
import paraview
```

For more information on ParaView's Python packages, see the 
[ParaView Wiki](https://www.paraview.org/Wiki/PvPython_and_PvBatch).


### Installation
After installing the dependencies above, clone this repo to your 
machine.

```
git clone https://github.com/kieran-mackle/pysagas
```

Next, use pip to install the `pysagas` package from the repo you 
just cloned.

```
python3 -m pip install pysagas
```

<!-- end getting started -->

<p align="right">[<a href="#readme-top">back to top</a>]</p>




## Usage

<!-- start usage -->

*PySAGAS* uses low-order methods to approximate sensitivities
on the surface of aerodynamic geometries. The user must provide
a nominal condition of the flow properties on the surface, along
with the sensitivity of the geometric vertices to the design 
parameters. With these things, a *wrapper* can be conveniently 
used.

### Cart3D Wrapper
To generate surface sensitivities from a Cart3D solution, 
you are required to provide 2 things to *PySAGAS*:
1. The `Components.i.plt` file, containing the properties
on the geometry surface.
2. A parameter sensitivities file, for example generated
by [hypervehicle](https://github.com/kieran-mackle/hypervehicle).

With these two things, simply use the Cart3D wrapper, as 
shown in the snippet below. Note you must specify the 
freestream speed of sound and density, along with the path
to both the sensitivity file, and the components file.

```python
from pysagas.wrappers import Cart3DWrapper

wrapper = Cart3DWrapper(
    a_inf=a_inf,
    rho_inf=rho_inf,
    sensitivity_filepath=sensitivity_filepath,
    components_filepath=components_filepath,
)
F_sense = wrapper.calculate()
```


<!-- end usage -->

<p align="right">[<a href="#readme-top">back to top</a>]</p>




## Roadmap

PySAGAS is being developed along the following roadmap.

* [x] Modularisation of models
* [x] API for convenient, generalised usage
* [x] Integration with [hypervehicle](https://github.com/kieran-mackle/hypervehicle) for shape optimisation
* [x] Moment sensitivities
* [x] Control over surface tags / faces being analysed
* [x] Visualisation of sensitivities on mesh (cell-wise)
* [ ] Implementation of higher-fidelity correction models
* [ ] CLI
* [ ] Testing with flow solutions from [Eilmer](https://github.com/gdtk-uq/gdtk)


<p align="right">[<a href="#readme-top">back to top</a>]</p>




## Contributing 

<!-- start contribution guidelines -->

To contribute to `pysagas`, please read the instructions below,
and stick to the styling of the code.

### Setting up for Development

1. Create a new Python virtual environment to isolate the package. You 
can do so using [`venv`](https://docs.python.org/3/library/venv.html) or
[anaconda](https://www.anaconda.com/).

2. Install the code in editable mode using the command below (run from
inside the `pysagas` root directory). Also install all dependencies 
using the `[all]` command, which includes the developer dependencies.

```
pip install -e .[all]
```

3. Install the [pre-commit](https://pre-commit.com/) hooks.

```
pre-commit install
```

4. Start developing! After following the steps above, you are ready
to start developing the code. Make sure to follow the guidelines 
below.


### Developing PySAGAS

- Before making any changes, create a new branch to develop on using 
`git checkout -b new-branch-name`.

- Run [black](https://black.readthedocs.io/en/stable/index.html) on any
code you modify. This formats it according to 
[PEP8](https://peps.python.org/pep-0008/) standards.

- Document as you go: use 
[numpy style](https://numpydoc.readthedocs.io/en/latest/format.html) 
docstrings, and add to the docs where relevant.

- Write unit tests for the code you add, and include them in `tests/`. 
This project uses [pytest](https://docs.pytest.org/en/7.2.x/).

- Commit code regularly to avoid large commits with many changes. 

- Write meaningful commit messages, following the 
[Conventional Commits standard](https://www.conventionalcommits.org/en/v1.0.0/).
The python package [commitizen](https://commitizen-tools.github.io/commitizen/)
is a great tool to help with this, and is already configured for this
repo. Simply stage changed code, then use the `cz c` command to make a 
commit.

- Open a [Pull Request](https://github.com/kieran-mackle/pysagas/pulls) 
when your code is complete and ready to be merged.


### Building the Docs
To build the documentation, run the commands below. 

```
cd docs/
make html
xdg-open build/html/index.html
```

If you are actively developing the docs, consider using
[sphinx-autobuild](https://pypi.org/project/sphinx-autobuild/).
This will continuosly update the docs for you to see any changes
live, rather than re-building repeatadly. From the `docs/` 
directory, run the following command:

```
sphinx-autobuild source/ build/html --open-browser
```

<!-- end contribution guidelines -->

<p align="right">[<a href="#readme-top">back to top</a>]</p>



## License
To be confirmed.


<p align="right">[<a href="#readme-top">back to top</a>]</p>
