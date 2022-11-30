# Getting Started

## Prerequisites
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

## Installation
After installing the dependencies above, clone this repo to your machine.

```
git clone https://github.com/kieran-mackle/py-adjoint
```

Next, use pip to install the `py-adjoint` from repo you just cloned.

```
python3 -m pip install py-adjoint
```


