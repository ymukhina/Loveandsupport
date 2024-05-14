using Logging
using Oscar
using StructuralIdentifiability
import StructuralIdentifiability: reduce_ode_mod_p, power_series_solution, ps_diff 

function f_min_support(ode::ODE)
    x1, x2 = ode.x_vars
    g1 = ode.x_equations[x1]
    g2 = ode.x_equations[x2]
    d1 = total_degree(g1)
    d2 = total_degree(g2)
    if d1 <= d2
        A = [1 d1 (d1 + d2 - 1); -1 0 0; 0 -1 0; 0 0 -1]
        b = [d1*(d1 + d2 - 1), 0, 0, 0]
    else
        A = [1 d1 (2*d1 - 1); 1 d2 (d1 + d2 - 1); -1 0 0; 0 -1 0; 0 0 -1]
        b = [d1*(2*d1 - 1), d1*(d1 + d2 - 1), 0, 0, 0]
    end
    return lattice_points(polyhedron(A, b))
end
    
#Randomize general polynomial
function generic_poly(R, d)
    (t1, t2) = gens(R)
f = zero(R)
    for i = 0:d
        for j = 0:(d-i)
               
           f += rand(1:10) * t1^i * t2^j 
       
                end
    end
return f 
end



function solve_with_love_and_support(ode::ODE, p::Int)
    @assert is_probable_prime(p) "This is not a prime number, Yulia!"

    ode_mod_p = reduce_ode_mod_p(ode, p)
    x1, x2 = ode_mod_p.x_vars
    y = first(ode_mod_p.y_vars)
    F = Nemo.Native.GF(p)

    # compute Newton polytope of f_min
    possible_supp = f_min_support(ode)
    @info "The size of the estimates support is $(length(possible_supp))"
    nterms = length(possible_supp) + 2

    # random initial conditions
    ic = Dict(x1 => rand(1:p-1), x2 => rand(1:p-1))
    # no parameters, no inputs
    par = empty(ic)
    inp = empty(Dict(x1 => [1]))

    ps_sol = power_series_solution(ode_mod_p, par, ic, inp, nterms)
    pss = [ps_sol[y]]
    for i in 1:2
        push!(pss, ps_diff(pss[end]))
    end
    prods = [prod(pss .^ exp) for exp in possible_supp]

    ls = matrix([coeff(pr, j) for j in 0:(nterms - 3), pr in prods])
    ker = kernel(ls, side = :right)
    @info "The dimension of the solution space $(size(ker)[2])"
    sol = ker[:, 1]

    R, _ = polynomial_ring(F, ["x1", "x1'", "x1''"])

    mons = [prod(gens(R) .^ exp) for exp in possible_supp]

    return sum([s * m for (s, m) in zip(sol, mons)])
end

