# This script produces the table for Section 4.1 of the paper aiming
# at experimental exploration of the accuracy of the produced bound

using DiffMinPoly
using Oscar
using Nemo
using StructuralIdentifiability
using IterTools
import StructuralIdentifiability: _reduce_mod_p, reduce_ode_mod_p, power_series_solution, ps_diff, var_to_str, switch_ring 
using Random



const Ptype = QQMPolyRingElem

function my_convex_hull(points)
    pts = permutedims(hcat([Oscar.homogenize(v, 1) for v in points]...))
    ptype = Oscar._scalar_type_to_polymake(QQFieldElem)
    lin = zero_matrix(QQ, 0, size(pts, 2))
    return Oscar.Polyhedron{QQFieldElem}(Polymake.polytope.Polytope{ptype}(; VERTICES = pts, LINEALITY_SPACE=lin))
end

NUM_RUNS = 1
DEGREES = [
    [2, 1, 1],
    [2, 2, 2],
    [2, 3, 3],
    [2, 4, 4],
    [2, 5, 5], 
    [3, 1, 1],         
    [3, 2, 2],
    [3, 3, 3],
    [1, 2, 2, 2],
    [2, 1, 1, 1],
] 
 
println("Degrees | # terms in the bound | # terms in NP of f_min | # terms in f_min | %")
       
for ds in DEGREES
    for i in 1:NUM_RUNS
        ode = rand_ode(ds) 
        x = first(sort(ode.x_vars, rev = true))            
        bound_size = size(DiffMinPoly.f_min_support(ode, x, DiffMinPoly.minpoly_order(ode, x)))[1]

        minpoly_exponents = collect(exponent_vectors(
            eliminate_with_love_and_support_modp(ode, x, Int(rand_bits_prime(ZZ, 32)), info = false)
        ))
        minpoly_points_inside = length(lattice_points(my_convex_hull(minpoly_exponents)))
                
        accuracy = length(minpoly_exponents) * 100 / bound_size     

        println("$ds | $bound_size | $(length(minpoly_exponents)) | $minpoly_points_inside | $accuracy ") 
      
    end
end         

