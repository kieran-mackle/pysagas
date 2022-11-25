import codecs
import os.path
import setuptools


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


# Define extra dependencies
dev_dep = ["pytest >= 7.1.1", "black >= 22.10.0", "commitizen >= 2.35.0", "pre-commit"]

all_dep = dev_dep

# Load README
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="adjoint",
    version=get_version("adjoint/__init__.py"),
    author="Kieran Mackle",
    author_email="kemackle98@gmail.com",
    description="The poor man's adjoint solver.",
    long_description=long_description,
    project_urls={
        "Bug Tracker": "https://github.com/kieran-mackle/py-adjoint/issues",
        "Source Code": "https://github.com/kieran-mackle/py-adjoint",
    },
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    package_dir={"": "."},
    packages=setuptools.find_packages(where="."),
    python_requires=">=3.7",
    install_requires=[],
    extras_require={
        "dev": dev_dep,
        "all": all_dep,
    },
    setup_requires=[
        "setuptools_git",
        "setuptools_scm",
    ],
)
