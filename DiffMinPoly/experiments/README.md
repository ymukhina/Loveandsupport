## Experiments

For all examples of Section 4 of the paper we provide corresponding files, which can be ran to reproduce our timings above,
after installing our package. 


To produce the tables for Section 4.1 of the paper aiming at experimental exploration of the accuracy of the produced bound
use
```julia
include("DiffMinPoly/experiments/bound_accuracy.jl")
```


To produce the table for Section 4.2 of the paper aiming at experimental exploration of the potential support-based improvements for the bound use
```julia
include("DiffMinPoly/experiments/optimal_support_based_bound.jl")
```


To produce the table for Section 4.3 of the paper aiming at experimental exploration of an alternative approach via tropical implicitization use
```julia
include("DiffMinPoly/experiments/tropical_comparison.jl")
```
