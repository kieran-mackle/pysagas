# Inclined Ramp Example

This page covers a simple test case of an inclined ramp.


## Problem Definition

The parameters of this problem are 
$\mathbf{p} = [\theta, \, L, \, W]$, shown
diagramatically in the figure below. The geometry is
built from two [cells](cell-definition), sharing a common edge.


![Inclined Ramp Example Study](../_static/ramp.png)

### Free Stream Conditions
The freestream conditions are defined below.

$$\gamma = 1.4$$

$$ Mach_{freestream} = 6$$

$$ P_{freestream} = 700 Pa$$

$$ T_{freestream} = 70 K$$


## Analytical Solution
This problem can be solved analytically using simple
isentropic flow relations. This solution will later
serve as validation for the results obtained using
*PySAGAS*.




## PySAGAS Solution
Given the surface properties on the ramp calculated using
the analytical solution, the parameter sensitivities can
be approximated using *PySAGAS*.






