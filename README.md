# DiffMinPoly.jl

## About

`DiffMinPoly.jl` is a Julia package for the elimination problem for polynomial dynamical systems.

## How to install

The package can be installed from this repository by

```julia
using Pkg
Pkg.add(url = "https://github.com/ymukhina/Loveandsupport.git", subdir="DiffMinPoly")
```

## How to use

The package can be loaded by `using DiffMinPoly`.

For the ODE system 
``` math 
\begin{cases} 
    x_1' = x_2^2,\\
    x_2' = x_1.
    \end{cases} 
```
to perform the elimination for variable $x_1$ we use the function `eliminate`. 
For instance:

```julia
using StructuralIdentifiability


ode = @ODEmodel(
    x1'(t) = x2(t)^2,
    x2'(t) = x1(t),
    y(t) = x1(t)
)

eliminate(ode, x1)
```
will return

```
x1(t)^2*x1(t)^(1) - 1//4*x1(t)^(2)^2
```



## Contacts

Maintained by Yulia Mukhina (yulia.mukhina@lix.polytechnique.fr) and Gleb Pogudin (gleb.pogudin@polytechnique.edu).