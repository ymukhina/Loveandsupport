# This script produces the table for Section 4.2 of the paper aiming
# at experimental exploration of the accuracy of the produced bound

include("../solver_love_and_support.jl")
include("../test_utils.jl")

cases = [
       rand_ode([2,1])
       rand_ode([2,1,1])         
    ]
     
                
    push!(
    cases,
    @ODEmodel(
        x1'(t) = x2(t)^2 + x1(t)*x2(t) + x1(t)^2 + 1,
        x2'(t) = x2(t),
        y(t) = x1(t)
    )
)            
                

     for c in cases
        ode = c
        x = first(sort(ode.x_vars, rev = true))
        jac = minpoly_order(ode, x)                
                        
        supp = f_min_support(ode, x, jac)
                    
                    
        #h = eliminate(ode, x, 0.99)
        #m = collect(exponent_vectors(h)) 
        
        l = my_convex_hull(supp) 
        s = length(lattice_points(l))              
        F = facets(l)          
        @info "Newton polytope for the system $c is"
        @info "$F, size $s"            
     end          
