<a name="readme-top"></a>

# Poor Mans Adjoint

A cheap alternative to adjoint flow solvers.


<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
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
Clone this repo to your machine.

```
git clone https://github.com/kieran-mackle/py-adjoint
```

Use pip to install the repo you just cloned.

```
python3 -m pip install py-adjoint
```

<p align="right">[<a href="#readme-top">back to top</a>]</p>



## Usage
Coming soon.


## Roadmap
To be determined.


## Contributing 

More details coming soon.

1. Create a new virtual environment.

2. Install the code in editable mode.

```
pip install -e .[all]
```

3. Install the [pre-commit](https://pre-commit.com/) hooks, and run against 
all files.

```
pre-commit install
pre-commit run --all-files
```

### Building the docs

```
cd docs/
make html
xo build/html/index.html
```

If you are actively developing the docs, consider using
[sphinx-autobuild](https://pypi.org/project/sphinx-autobuild/).

```
sphinx-autobuild source/ build/ --open-browser
```

<p align="right">[<a href="#readme-top">back to top</a>]</p>



## License
To be confirmed.

