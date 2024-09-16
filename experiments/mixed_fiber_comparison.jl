# This script produces the table for Section 4.3 of the paper aiming
# at experimental exploration of an alternative approach via mixed fiber polytope

include("../solver_love_and_support.jl")
include("../test_utils.jl")

using Nemo
using Groebner
using Oscar
using Polyhedra
import GLPK


    DEGREES = [
        [2, 1, 1],
        [3, 1, 1],
    ] 

    R, (x3, x2, x1, x11, x111, x1111) = polynomial_ring(GF(65521), ["x3","x2", "x1", "x11", "x111", "x1111"])


for d in DEGREES

    g1 = rand_poly(d[1], [x1, x2, x3])
    g2 = rand_poly(d[2], [x1, x2, x3])
    g3 = rand_poly(d[3], [x1, x2, x3])

    #computing the derivatives
    polys = [x11 - g1, x111 - (derivative(g1, x1) * x11 + derivative(g1, x2) * g2 + derivative(g1, x3) * g3)]
    push!(polys, x1111 - (derivative(polys[2], x1) * x11 + derivative(polys[2], x2) * g2 + derivative(polys[2], x3) * g3 + derivative(polys[2], x11) * x111))

    #computing the generic system for polys
    polys_new = [sum([rand(1:1000)*m for m in Oscar.monomials(p)]) for p in polys]    


    #computing mixed fiber polytope
    S = parent(first(polys_new))
    gb = groebner_basis_f4(Oscar.ideal(S, polys_new), eliminate = 2, complete_reduction=true, info_level=2)
    supp = collect(exponent_vectors(gb[1]))

    newton_polytope = my_convex_hull(supp)
    np_size = length(lattice_points(newton_polytope))      



    ode = rand_ode(d) 
    x = first(sort(ode.x_vars, rev = true))            
    bound_size = size(f_min_support(ode, x, minpoly_order(ode, x)))[1]

    println("$d | $bound_size | $np_size ") 

end