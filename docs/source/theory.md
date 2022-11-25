# Poor Man's Adjoint Solver
The proposed idea is to use local piston theory to evaluate the senitivity between surface pressure, and design and flow parameters. 

From Piston Theory:

$$ P = P_\infty \, \left( 1 + \frac{\gamma-1}{2} \frac{W}{a} \right) ^{\frac{2 \, \gamma}{\gamma-1}}$$

This can be converted to local pressure variations as per Zhan (2009).

$$ P = P_l + \rho_l \, a_l \, W$$

$$ W = \mathbf{V}_l \cdot \delta \mathbf{n} + \mathbf{V}_b \cdot \mathbf{n}$$

$$ \delta \mathbf{n} = \mathbf{n}_0 - \mathbf{n} $$

where subscript $l$ refers to the local conditions at the vehicle surface 
(e.g. results from a CFD solver) and $W$ is the local downwash velocity.

In their formulation, they consider deformation (or change in geometry) away 
from the shape when the local conditions were evaluated using the 
$\mathbf{V}_l \cdot \delta \mathbf{n}$ and the contribution due to surface 
normal velocity $\mathbf{V}_b \cdot \mathbf{n}$.

We are only interested in the first term, as this relates to vehicle deformation. 

Thus:

$$ P = P_l + \rho_l \, a_l \, \mathbf{V}_l \cdot \left( \mathbf{n}_0 - \mathbf{n} \right)$$


Ref: 
Zhan, W.-W., Ye, Z.-Y., and Zhang, C.-A., “Supersonic Flutter Analysis Based on a Local Piston Theory,” AIAA Journal,
Vol. 47, No. 10, 2009.



## Approach
To get the geometry sensitivities (e.g. $\frac{d \, C_X}{d \, p}$) we will apply the following process:

1. Find effect of design parameter, $p$ on surface normal. 

$$ \rightarrow \frac{d \, \mathbf{n}}{d \, p} $$

2. Using local piston theory get sensitivty between $P$ and $\mathbf{n}$.

$$ \rightarrow \frac{d \, P}{d \, \mathbf{n}} $$

3. Get total force 

$$ F_x = \sum_{i=0, \, N} P_i(p) \, A_i(p) \, \mathbf{n}(p) \cdot [1, \, 0, \, 0]^T $$

$$ \frac{d}{dp}(F_x) = \sum_{i=0, \, N} \left[ \frac{d \, P_i}{dp} \, A_i(p) \, \mathbf{n}(p) \cdot [1, \, 0, \, 0]^T 
+ P_i(p) \, \frac{d \, A_i}{d \, p} \, \mathbf{n}(p) \cdot [1, \, 0, \, 0]^T
+ P_i(p) \, A_i(p) \, \frac{d \, \mathbf{n}}{dp} \cdot [1, \, 0, \, 0]^T
\right] $$


4. Convert to coefficient. 
Just divide by $\frac{1}{2} \, \rho_\infty \, V_\infty^2 \, A_{ref} $ 


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
- Lighthill, M. J., “Oscillating Airfoils at High Mach Numbers,” Journal of the Aeronautical Sciences, Vol. 20, No. 6, June 1953, pp. 402–406

- Ashley, H., and Zartarian, G., “Piston Theory: A New Aerodynamic Tool for the Aeroelastician,” Journal of the Aeronautical Sciences, Vol. 23, No. 12, Dec. 1956, pp. 1109–1118.

- Van Dyke, M., “A Study of Second-Order Supersonic Flow Theory,” NACA TR 1081, 1951.

- Zartarian, G., Hsu, P. T., and Ashley, H., “Dynamic Airloads and Aeroelastic Problems at Entry Mach Numbers,” Journal of the Aeronautical Sciences, Vol. 28, No. 3, March 1961, pp. 209–222.


## More complex correction equations

One caveat of usng higher order correction equations is that these may be better/worse suited for specific Mach number ranges. That is different approaches may be required for different parts of the flight envelope. 

$$ q^2 = (U+u)^2 + v^2 + w^2 $$

For isentropic flow the pressure coefficient is given by:

$$C_p = \frac{P - P_\infty}{\frac{1}{2} \, \rho_\infty \, U^2} = \frac{2}{\gamma \, M^2} \, \left[ \left( 1 + \frac{\gamma-1}{2} \, M^2 \, \left( 1 - \frac{q^2}{U^2} \right) \right)^{\frac{\gamma}{\gamma-1}} - 1 \right] $$

Expanding the equation yields:

$$C_p = -2 \, \frac{u}{U} - \frac{v^2 + w^2}{U^2} + \beta^2 \, \frac{u^2}{U^2} + M^2 \, \frac{u}{U} \, \frac{v^2+w^2}{U^2} + \frac{M^2}{4} \, \left( \frac{v^2+w^2}{U^2} \right)^2 + 0 \, \left[ \frac{u^3}{U^3}, \, \frac{u^2}{U^2} \, \frac{v^2+w^2}{U^2}, \, \frac{u}{U} \, \left( \frac{v^2+w^2}{U^2} \right)^2, \, \left(\frac{v^2 +w^2}{U^2} \right)^3 \right] $$

$$\beta = \sqrt{M^2-1}$$

All the explicitly listed terms will contribute to second order solutions. 



