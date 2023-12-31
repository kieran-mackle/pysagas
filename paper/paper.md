---
title: 'PySAGAS: Aerodynamic Sensitivities for Parametric Geometries'
tags:
  - Python
  - engineering
  - cfd
  - hypersonics
  - optimisation
authors:
  - name: Kieran Mackle
    orcid: 0000-0001-8245-250X
    affiliation: 1
  - name: Ingo Jahn
    orcid: 0000-0001-6962-8187
    affiliation: 2
affiliations:
 - name: The University of Queensland, Australia
   index: 1
 - name: The University of Southern Queensland, Australia
   index: 2
date: 01 January 2024
bibliography: paper.bib
---


# Summary

PySAGAS is a **Py**thon package for **S**ensitivity **A**pproximations of **G**eometric and **A**erodynamic **S**urface properties. 
It provides a computationally-efficient method for approximating the sensitivities of aerodynamic forces and moments with respect to design parameters, from nominal surface properties provided by a single CFD flow solution.
The sensitivities outputted by PySAGAS can be used for sensitivity studies, or for aerodynamic shape optimisation.
PySAGAS also provides a conventient API to multiple senstivity model implementations, as well as optimisation wrappers and geometry utilities.

The Pyhon package has been developed according to best practices, and extensive documentation is provided.

# Statement of need

Gradient-based optimization techniques possess many desireable characteristics for the application of aerodynamic shape optimisation.
Specifically, these techniques are known for being efficient and scalable with the number of design parameters.
As the name implies, they require information about the gradient of the objective function with respect to the design parameters.
For analytical or computationally-inexpensive problems, this is a relatively straightforward hurdle to overcome.
Analytical expressions can easily be differentiated, and gradient information can be extracted from inexpensive problems by use of numerical techniques, such as finite-differencing.
However, as the system being optimised becomes more computationally demanding, either by virtue of complexity or dimensionality, these approaches of calculating gradient information are either impossible or intractable.

PySAGAS implements a novel solution to this problem in the context of hypersonic vehicle design.
A detailed description of the methodology employed can be found in [@MackleShapeOpt].
This package has played a key role in research conducted for the efficient optimisation of hypersonic vehicle configurations [@MackleCoDesign].


# Illustrative Example: Hypersonic Waverider Optimisation

![Shape optimised hypersonic waverider using *PySAGAS*.](vd-optimised-outlined.png){ width=70% }

+-------------------+-------------------------------+
| Sensitivity model | Time to converge [hh:mm:ss]   |
|                   |                               |
+:=================:+:=============================:+
| Van Dyke's theory |           4:02:55             |
+-------------------+-------------------------------+
| Piston theory     |           6:15:18             |
+-------------------+-------------------------------+
| Adjoint CFD       |           6:40:00             |
+-------------------+-------------------------------+
| Finite difference |       2 days, 8:40:00         |
+===================+===============================+


# Installation and Documentation

PySAGAS is extensively documented in-code through the use of docstrings.
Comprehensive documentation is also provided online at Read the Docs.


# References
