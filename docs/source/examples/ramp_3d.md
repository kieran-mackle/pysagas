# Poor Man's Adjoint Solver - 3D

Same as Poor Man's Adjoint, but the ramp has been split into two triangles. 


## Free Stream Conditions - used by all examples

Use ideal gas air and a flow condition that approximates TUSQ
$$\gamma = 1.4$$
$$ Mach_{freestream} = 6$$
$$ P_{freestream} = 700 Pa$$
$$ T_{freestream} = 70 K$$


## Reference Case Study - Ramp (shock expansion)

Simplest case ramp at angle, $\theta$, with length, $L$ and width, $W$.
Parameters $\mathbf{p} = [\theta, \, L, \, W]$

Calculate nominal pressure and forces and normal and then also sensitivites. 


## Calculate sensitivities using local piston theory

### Step 1: Work out $\frac{d \, \mathbf{n}}{d \, \mathbf{p}}$

This is done in a two stage process
$$ \frac{d \, \mathbf{n}}{d \, \mathbf{p}} = \frac{d \, \mathbf{n}}{d \, \mathbf{pos}} \cdot \frac{d \, \mathbf{pos}}{d \, \mathbf{p}}$$ 



## Step 2: Do pressure Caculation and differentiate $\frac{d \, P}{d \, \mathbf{p}}$

This is a straightforward calculation
$$ P = P_l + \rho_l \, a_l \, \mathbf{V}_l \cdot \left( \mathbf{n}_0 - \mathbf{n} \right)$$

Hence
$$ \frac{d \, P}{d \, \mathbf{p}} = \rho_l \, a_l \, \mathbf{V}_l \cdot  \left(- \frac{d \, \mathbf{n}}{d \, \mathbf{p}} \right)$$



## Step 3: Integrate to get force contribution (actually summation) $\frac{d \, F_i}{d \, \mathbf{p}}$

$$ F_x = P(p) \, A(p) \, \mathbf{n}(p) \cdot [1, \, 0, \, 0]^T $$

$$ \frac{d}{d \, \mathbf{p}}(F_x) = \frac{d \, P}{d \, \mathbf{p}} \, A(p) \, \mathbf{n}(p) \cdot [1, \, 0, \, 0]^T 
                      + P(p) \, \frac{d \, A}{d \, \mathbf{p}} \, \mathbf{n}(p) \cdot [1, \, 0, \, 0]^T
                      + P(p) \, A(p) \, \frac{d \, \mathbf{n}}{d \, \mathbf{p}} \cdot [1, \, 0, \, 0]^T $$


