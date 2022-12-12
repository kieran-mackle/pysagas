<a name="readme-top"></a>

# PySAGAS

<!-- start intro -->
**PySAGAS** is a Python package for the generation of **S**ensitivity **A**pproximations
for **G**eometric and **A**erodynamic **Surface** properties.
It provides a computationally-efficient method for generating 
surface sensitivity approximations from existing flow solutions,
to use in geometric parameter optimisation studies.
<!-- end intro -->

<a><img src="https://github.com/0x6080604052/analytics/actions/workflows/tests.yml/badge.svg" alt="Test Status"></a>



<!-- TABLE OF CONTENTS -->
<details>
  <summary><h2>Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
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
This code depends on the [Eilmer](https://github.com/gdtk-uq/gdtk) python 
package. Note that a full Eilmer install is not required. Instead, do a 
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

### Installation
After installing the dependencies above, clone this repo to your machine.

```
git clone https://github.com/kieran-mackle/pysagas
```

Next, use pip to install the `pysagas` from repo you just cloned.

```
python3 -m pip install pysagas
```

<!-- end getting started -->

<p align="right">[<a href="#readme-top">back to top</a>]</p>




## Usage
Detailed usage instructions coming soon.


<p align="right">[<a href="#readme-top">back to top</a>]</p>




## Roadmap

PySAGAS is being developed along the following roadmap.

* [ ] Modularisation of models
* [ ] API for convenient, generalised usage
* [ ] Implementation of higher-fidelity correction models
* [ ] Control over surface tags / faces being analysed
* [ ] CLI
* [ ] Testing with flow solutions from [Eilmer](https://github.com/gdtk-uq/gdtk)
* [ ] Integration with [hypervehicle](https://github.com/kieran-mackle/hypervehicle)


<p align="right">[<a href="#readme-top">back to top</a>]</p>




## Contributing 

<!-- start contribution guidelines -->

To contribute to `pysagas`, please read the instructions below,
and stick to the styling of the code.

### Setting up for Development

1. Create a new Python virtual environment to isolate the package. You 
can do so using [`venv`](https://docs.python.org/3/library/venv.html) or
[anaconda](https://www.anaconda.com/).

2. Install the code in editable mode using the command below. Also install
all dependencies using the `[all]` command, which includes developer 
dependencies.

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


### Building the docs
To build the documentation, run the commands below. 

```
cd docs/
make html
xdg-open build/html/index.html
```

If you are actively developing the docs, consider using
[sphinx-autobuild](https://pypi.org/project/sphinx-autobuild/).
This will continuosly update the docs for you to see any changes
live, rather than re-building repeatadly. 

```
sphinx-autobuild source/ build/ --open-browser
```

<!-- end contribution guidelines -->

<p align="right">[<a href="#readme-top">back to top</a>]</p>



## License
To be confirmed.


<p align="right">[<a href="#readme-top">back to top</a>]</p>
