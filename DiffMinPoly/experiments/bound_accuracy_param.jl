# This script produces the tables for Section 4 of the paper aiming
# at experimental exploration of the accuracy of the produced bound

include("../src/solver_love_and_support.jl")
include("../src/utils.jl")

using Oscar
using Nemo
using StructuralIdentifiability
using IterTools
import StructuralIdentifiability: _reduce_mod_p, reduce_ode_mod_p, power_series_solution, ps_diff, var_to_str, switch_ring 
using Random

using Polyhedra

const Ptype = QQMPolyRingElem

function my_convex_hull(points)
    pts = permutedims(hcat([Oscar.homogenize(v, 1) for v in points]...))
    ptype = Oscar._scalar_type_to_polymake(QQFieldElem)
    lin = zero_matrix(QQ, 0, size(pts, 2))
    return Oscar.Polyhedron{QQFieldElem}(Polymake.polytope.Polytope{ptype}(; VERTICES = pts, LINEALITY_SPACE=lin))
end


NUM_RUNS = 5

#non-parametric case
DEGREES = [
    [1,1,1],
    [2,2,1],
    [2,2,2],
    [2,2,3],
    [2,2,4],
    [2,2,5],
    [3,3,1],
    [3,3,2],
    [3,3,3],
    [3,3,4],
    [3,3,5],
    [4,4,1],
    [4,4,2],
    [1,1,1,1],
    [1,1,1,2],
    [1,1,1,3],
    [2,2,2,1],
    [3,3,3,1],
] 


DEGREES_PARAM = [
    # case d_{\mu} = 0,  D_{\mu} = 1
    [(1, 1), (1, 1), (1, 0)],
    [(1, 1), (1, 1), (2, 0)],
    [(1, 1), (1, 1), (3, 0)],
    [(2, 1), (2, 1), (1, 0)],
    [(2, 1), (2, 1), (2, 0)],
    [(3, 1), (3, 1), (1, 0)],
    [(3, 1), (3, 1), (2, 0)],
    [(4, 1), (4, 1), (1, 0)],
    # case d_{\mu} = D_{\mu} = 1
    [(1, 1), (1, 1), (1, 1)],
    [(1, 1), (1, 1), (2, 1)],
    [(1, 1), (1, 1), (3, 1)], 
    [(2, 1), (2, 1), (1, 1)],
    [(2, 1), (2, 1), (2, 1)],  
    [(3, 1), (3, 1), (1, 1)],
     [(3, 1), (3, 1), (2, 1)],
    [(4, 1), (4, 1), (1, 1)]
] 

println("Degrees| # terms in the bound | # terms in NP of f_min | # terms in f_min | %")


for ds in DEGREES
    for i in 1:NUM_RUNS
            ode = rand_ode(ds)        
            bound_size = size(f_min_support(ode, minpoly_order(ode)))[1]

            minpoly_exponents = collect(exponent_vectors(
                eliminate_with_love_and_support_modp(ode, 2^31 - 1)))
            minpoly_points_inside = length(lattice_points(my_convex_hull(minpoly_exponents)))
                
            accuracy = length(minpoly_exponents) * 100 / bound_size     

            println("$ds | $bound_size | $(size(minpoly_exponents)[1])| $minpoly_points_inside | $accuracy ") 
    end
end  
        

println("Degrees | # num of parameters | # terms in the bound | # terms in NP of f_min | # terms in f_min | %")


for ds in DEGREES_PARAM
    for i in 1:NUM_RUNS
        for k in 1:2
            ode = rand_ode(ds, num_params=k)        
            bound_size = size(f_min_support(ode, minpoly_order(ode)))[1]

            minpoly_exponents = collect(exponent_vectors(
                eliminate_with_love_and_support_modp(ode, 2^31 - 1)))
            minpoly_points_inside = length(lattice_points(my_convex_hull(minpoly_exponents)))
                
            accuracy = length(minpoly_exponents) * 100 / bound_size     

            println("$ds | $k | $bound_size | $(size(minpoly_exponents)[1])| $minpoly_points_inside | $accuracy ") 
        end    
    end
end         





