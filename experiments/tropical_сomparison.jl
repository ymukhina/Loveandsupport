# This script produces the table for Section 4.2 of the paper aiming
# at experimental exploration of an alternative approach via tropical implicitization

include("../solver_love_and_support.jl")
include("../test_utils.jl")

cases = [
    rand_ode([2, 1]),
    rand_ode([2, 1, 1]),
    rand_ode([3, 1, 1])
]
     
push!(
    cases,
    @ODEmodel(
        x1'(t) = x2(t)^2 + x1(t) * x2(t) + x1(t)^2 + 1,
        x2'(t) = x2(t),
        y(t) = x1(t)
    )
)            
                

for c in cases
    ode = c
    x = first(sort(ode.x_vars, rev = true))
    minpoly_ord = minpoly_order(ode, x)                
                    
    supp = f_min_support(ode, x, minpoly_ord)
                
    newton_polytope = my_convex_hull(supp) 
    np_size = length(lattice_points(newton_polytope))              
    np_equations = facets(newton_polytope)
    @info "Newton polytope for the system $c is"
    @info "$np_equations, size $np_size"
end          
