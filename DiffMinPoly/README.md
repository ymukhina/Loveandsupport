# DiffMinPoly.jl


## About

`DiffMinPoly.jl` is a Julia package for the elimination problem for polynomial dynamical systems.

## How to install

The package can be installed from this repository by

```julia
using Pkg
Pkg.add(url="https://github.com/ymukhina/Loveandsupport.git", subdir="DiffMinPoly", rev="y-input")
```

## How to use

The package can be loaded by `using DiffMinPoly`.

For the ODE system 
``` math 
\begin{cases} 
    x_1' = a_1 x_2,\\
    x_2' = a_2 x_1.
    \end{cases} 
```
and 

``` math 
    y = x_1 + x_2
```

to perform the elimination for the chosen function of the coordinates we use the function `eliminate`. 
For instance:

```julia
using DiffMinPoly
using StructuralIdentifiability


ode = @ODEmodel(
                  x1'(t) = a1 * x2(t),
                  x2'(t) = a2 * x1(t),
                  y(t) = x1(t) + x2(t)
              )

eliminate(ode)
```
will return

```
a1*a2*y(t) - y(t)^(2)
```



## Contacts

Maintained by Yulia Mukhina (yulia.mukhina@lix.polytechnique.fr) and Gleb Pogudin (gleb.pogudin@polytechnique.edu).

## References

Based on
 * [Projecting dynamical systems via a support bound](https://arxiv.org/abs/2501.13680), preprint, 2025
 * *Support bound for differential elimination in polynomial dynamical systems*, preprint, 2025

