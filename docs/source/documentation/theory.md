(theory-docs)=

# Approximating Sensitivities Using Low-Order Models
The core idea underpinning *PySAGAS* is to use low-order, 
differentiable fluid models to evaluate the sensitivities 
of surface pressures to geometric parameters.


## High-Level Overview of Approach

```{tip}
Refer to the [nomenclature](nomenclature) for definitions of 
symbols used here.
```

(vertex-param-sens)=
### Vertex-Parameter Sensitivities
To make the link between [cell](cell-definition) data 
and geometric parameters, 
*PySAGAS* requires the user to provide the 
[vertex-parameter sensitivities](geom-param-sens), 
$d\mathbf{v}/d\mathbf{p}$. In simple cases,
such as the [inclined ramp example](example-inclined-ramp), 
$d\mathbf{v}/d\mathbf{p}$ can be found manually. However,
for more [complex geometries](example-nose-opt), 
this is not possible. For this
reason, a parametric geometry generation tool offering 
geometric differentiation is recommended. One such tool
is [*hypervehicle*](https://github.com/kieran-mackle/hypervehicle/).


(nominal-surface-properties)=
### Nominal Surface Properties
*PySAGAS* also requires the nominal surface properties of the
geometry being considered. Namely, the local Mach number, 
pressure and temperature. This information can come from
any source, making *PySAGAS* CFD-solver agnostic; whether it
be [Oblique Shock Theory](oblique-shock-relations), or 
[Cart3D](nose-cart3d-modelling).


### Pressure Sensitivities
Given the nominal surface properties, *PySAGAS* can calculate the
sensitivity of local pressure to changes in the normal vector of a 
[cell](cell-definition). Once these sensitivities are known, they 
can be linked to the geometric design parameters via application 
of the chain rule.

Let's start with *piston theory*. Surface pressure is given by 
the following expression, relative to the freestream properties,
denoted by the $\infty$ subscript.

$$ P = P_\infty \, \left( 1 + \frac{\gamma-1}{2} 
\frac{W}{a_{\infty}} \right) 
^{\frac{2 \, \gamma}{\gamma-1}}
$$

In this expression, $W$ is the downwash speed, $a_{\infty}$ is the 
freestream speed of sound, and $\gamma$ is the ratio of specific 
heats. By following the approach presented in Zhan (2009)[^1], keeping
only first-order terms and decomposing the downwash velocity, we 
obtain "local piston theory":

$$ P = P_l + \rho_l \, a_l \, W$$

$$ W = \mathbf{V}_l \cdot \delta \mathbf{n} + \mathbf{V}_b \cdot \mathbf{n} $$

$$ \delta \mathbf{n} = \mathbf{n}_0 - \mathbf{n} $$



[^1]: [Zhan, W.-W., Ye, Z.-Y., and Zhang, C.-A., “Supersonic Flutter 
Analysis Based on a Local Piston Theory,” AIAA Journal, Vol. 47, 
No. 10, 2009.](https://arc.aiaa.org/doi/10.2514/1.37750)


Here, the subscript $l$ refers to the local surface 
conditions - i.e. the 
[nominal surface properties](nominal-surface-properties)
as described above.
In this model, we are only interested in the first term, 
$\mathbf{V}_l \cdot \delta \mathbf{n}$, as this relates 
to geometric deformation. The second term, 
$\mathbf{V}_b \cdot \mathbf{n}$, relates to vibration, 
and can therefore be discarded. We therefore reach an
expression for the local pressure, relative to a cell's
normal vector.

$$ P = P_l + \rho_l \, a_l \, \mathbf{V}_l \cdot \left( \mathbf{n}_0 - \mathbf{n} \right)$$

This expression can be differentiated with respect to 
the geometric parameters $\mathbf{p}$, yielding the 
sensitivity of pressure to parameters $\mathbf{p}$:

$$ 
\frac{d \, P}{d \, \mathbf{p}} = 
\rho_l \, a_l \, \mathbf{V}_l \cdot 
\left(- \frac{d \, \mathbf{n}}{d \, \mathbf{p}} \right)
$$

Of course, this expression requires $d\mathbf{n}/d\mathbf{p}$. 
Conveniently, this can be evaluated using the 
[vertex parameter sensitivities](vertex-param-sens) and the 
[normal vector vertex sensitivities](normal-area-vertex-sens)
as follows:

$$ 
\frac{d \, \mathbf{n}}{d \, \mathbf{p}} = \frac{d \, \mathbf{n}}{d \, 
\mathbf{v}}  \frac{d \, \mathbf{v}}{d \, \mathbf{p}} 
$$



### Aerodynamic Sensitivities
With the sensitivity of pressure to parameters, $dP/d\mathbf{p}$,
*PySAGAS* is able to calculate the sensitivity of the aerodynamic
surface properties to the geometric design parameters
(i.e. $\frac{d \, C_X}{d \, \mathbf{p}}$). First, note that the 
force acting on a given cell is the product of the pressure acting 
on that cell, and the cell's area. The total force acting on the 
geometry is therefore the sum of cell forces:

<!-- the sum notation here can be improved for clarity -->

$$ F_x = \sum_{i=0, \, N} P_i(p) \, A_i(p) \, \mathbf{n}(p) \cdot 
[1, \, 0, \, 0]^T $$


This expression can be differentiated with respect to the geometric
parameters to yield the following.

$$ \frac{d}{dp}(F_x) = $$

$$
\sum_{i=0, \, N} \left[ \frac{d \, P_i}{dp} \, A_i(p) \, 
\mathbf{n}(p) \cdot [1, \, 0, \, 0]^T 
+ P_i(p) \, \frac{d \, A_i}{d \, p} \, \mathbf{n}(p) \cdot [1, \, 0, \, 0]^T
+ P_i(p) \, A_i(p) \, \frac{d \, \mathbf{n}}{dp} \cdot [1, \, 0, \, 0]^T
\right] 
$$


These sensitivities can also be converted into non-dimensional 
force coefficients by dividing through by 
$\frac{1}{2} \, \rho_\infty \, V_\infty^2 \, A_{ref} $.



## Avenues for Improvement
The critical component in this approach is that pressure variations have 
been linked to changes in the local surface inclintion, that is we have 
an approximate solution for $\frac{d \, P}{d \, \mathbf{n}}$. 

In the current formulation we use 1st order Piston Theory to find this 
relationship. Increasing the accuracy of this calculation will deliver 
more accuracte sensitivities. 

Alternate methods are:
- Second Order Piston Theory
- Shock Expansion Theory
- Van Dyke Second Order Theory

The merit (increased computational cost vs accuracy) of these should be 
explored and compared to actual adjoint calculations. 


Extra references to review:

- [Lighthill, M. J., “Oscillating Airfoils at High Mach Numbers,” 
Journal of the Aeronautical Sciences, Vol. 20, No. 6, June 1953, 
pp. 402–406](https://arc.aiaa.org/doi/abs/10.2514/8.2657?journalCode=jans)

- [Ashley, H., and Zartarian, G., “Piston Theory: A New Aerodynamic 
Tool for the Aeroelastician,” Journal of the Aeronautical Sciences, 
Vol. 23, No. 12, Dec. 1956, 
pp. 1109–1118.](https://arc.aiaa.org/doi/abs/10.2514/8.3740)

- [Van Dyke, M., “A Study of Second-Order Supersonic Flow Theory,” 
NACA TR 1081, 1951.](https://thesis.library.caltech.edu/10587/1/van-dyke-milton-1949-thesis.pdf)

- [Zartarian, G., Hsu, P. T., and Ashley, H., “Dynamic Airloads 
and Aeroelastic Problems at Entry Mach Numbers,” Journal of the 
Aeronautical Sciences, Vol. 28, No. 3, 
March 1961, pp. 209–222.](https://arc.aiaa.org/doi/10.2514/8.8927)


### More complex correction equations

One caveat of usng higher order correction equations is that these may 
be better/worse suited for specific Mach number ranges. That is different 
approaches may be required for different parts of the flight envelope. 

$$ q^2 = (U+u)^2 + v^2 + w^2 $$

For isentropic flow the pressure coefficient is given by:

$$C_p = \frac{P - P_\infty}{\frac{1}{2} \, \rho_\infty \, U^2} = 
\frac{2}{\gamma \, M^2} \, \left[ \left( 1 + \frac{\gamma-1}{2} \, M^2 \, \left( 1 - 
\frac{q^2}{U^2} \right) \right)^{\frac{\gamma}{\gamma-1}} - 1 \right] $$

Expanding the equation yields:

$$
C_p = -2 \, \frac{u}{U} - \frac{v^2 + w^2}{U^2} + \beta^2 \, \frac{u^2}{U^2} + 
M^2 \, \frac{u}{U} \, \frac{v^2+w^2}{U^2} + \frac{M^2}{4} \, 
\left( \frac{v^2+w^2}{U^2} \right)^2 
\\
+\left[ \frac{u^3}{U^3}, \, \frac{u^2}{U^2} \, 
\frac{v^2+w^2}{U^2}, \, \frac{u}{U} \, \left( \frac{v^2+w^2}{U^2} \right)^2, \, 
\left(\frac{v^2 +w^2}{U^2} \right)^3 \right] 
$$

$$
\beta = \sqrt{M^2-1}
$$

All the explicitly listed terms will contribute to second order solutions. 


### Directly using the Pressure Equation

First, take the derivative of $P$ with respect to $W$.

$$ 
\frac{dP}{dW} = 
\frac{P\gamma}{a} \cdot
\left(
    1 + W \frac{\gamma - 1}{2a}
\right) ^
\frac{\gamma + 1}{\gamma - 1}
$$


$$ W = \mathbf{V}_l \cdot (\mathbf{n}_0 - \mathbf{n}) + \mathbf{V}_b \cdot \mathbf{n}$$

$$ \frac{dW}{d\mathbf{n}} = -\mathbf{V}_l + \mathbf{V}_b $$


Next, calculate $\frac{dP}{d\mathbf{p}}$:


$$ \frac{dP}{d\mathbf{p}} = \frac{dP}{dW} \cdot \frac{dW}{d\mathbf{n}} \cdot \frac{d\mathbf{n}}{d\mathbf{p}} $$

$$ = \frac{P_{\infty} \gamma \left(2a+W\left(\gamma-1\right)\right)^
{\frac{\gamma+1}{\gamma-1}}}
{2^{\frac{\gamma+1}{\gamma-1}}a^{\frac{2\gamma}{\gamma-1}}}
\cdot
(\mathbf{V}_l + \mathbf{V}_b) 
\cdot 
\frac{d\mathbf{n}}{d\mathbf{p}}
$$


### Conditional Models
Different models for forward facing [cells](cell-definition) 
versus backwards facing [cells](cell-definition).

