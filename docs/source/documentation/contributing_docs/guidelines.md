# Contribution Guidelines
To contribute to `pysagas`, please read the instructions below,
and stick to the styling of the code.

## Setting up for Development

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


## Developing PySAGAS

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


## Building the docs
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

