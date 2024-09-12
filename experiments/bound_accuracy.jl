# This script produces the table for Section 4.1 of the paper aiming
# at experimental exploration of the accuracy of the produced bound

include("../solver_love_and_support.jl")
include("../test_utils.jl")

NUM_RUNS = 5
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
        bound_size = size(f_min_support(ode, x, minpoly_order(ode, x)))[1]

        minpoly_exponents = collect(exponent_vectors(
            eliminate_with_love_and_support_modp(ode, x, Int(rand_bits_prime(ZZ, 32)), info = false)
        ))
        minpoly_points_inside = length(lattice_points(my_convex_hull(minpoly_exponents)))
                
        accuracy = length(minpoly_exponents) * 100 / bound_size     

        println("$ds | $bound_size | $(length(minpoly_exponents)) | $minpoly_points_inside | $accuracy ") 
    end
end         

