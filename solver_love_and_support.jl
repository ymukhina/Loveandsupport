using Oscar
using Nemo
using StructuralIdentifiability
using IterTools
import StructuralIdentifiability: _reduce_mod_p, reduce_ode_mod_p, power_series_solution, ps_diff 

const Ptype = QQMPolyRingElem

# -------- estimate support for f_min based on Theorems 1 and 2 -------- #

function f_min_support_bivariate(ode::ODE)
    @assert length(ode.x_vars) == 2 "system has more than two variables."
    x1, x2 = sort(ode.x_vars, rev = true)
    g1 = ode.x_equations[x1]
    g2 = ode.x_equations[x2]
    d1 = total_degree(g1)
    d2 = total_degree(g2)
    if d1 <= d2
        A = [1 d1 (d1 + d2 - 1); -1 0 0; 0 -1 0; 0 0 -1]
        b = [d1 * (d1 + d2 - 1), 0, 0, 0]
    else
        A = [1 d1 (2 * d1 - 1); 1 d2 (d1 + d2 - 1); -1 0 0; 0 -1 0; 0 0 -1]
        b = [d1*(2*d1 - 1), d1*(d1 + d2 - 1), 0, 0, 0]
    end
    return collect(lattice_points(polyhedron(A, b)))
end

function f_min_support(ode::ODE; info = false)
    # I guess by sorting you try to ensure that the first x is the output, i am not sure
    # this would always work. You could just choose x such that the only y_equation is equal to x
    x = sort(ode.x_vars, rev = true)
    info && @info "Inferred order $x"
    n = length(x)
    gs = [ode.x_equations[xi] for xi in x]
    d1 = total_degree(gs[1])
    @assert d1 > 0 "d1 = 0"
    D = maximum(total_degree, gs[2:end])
    info && @info "We have d1 = $d1 and D = $D"

    if d1 <= D
        ineq_lhs = reshape([1, [d1 + (k - 1) * (D - 1) for k in 1:n]...], 1, n + 1)
        ineq_rhs = [prod([d1 + (k - 1) * (D - 1) for k in 1:n])]
        A = vcat(matrix(QQ, ineq_lhs), -identity_matrix(QQ, n + 1))
        b = vcat(ineq_rhs, zeros(QQ, n + 1))
    else
        ineq_lhs1 = [k <= l ? k * (D - 1) + 1 : 0 for l in 0:(n - 1), k in 0:n]
        ineq_lhs2 = zeros(Int, n, n + 1)
        for l in 0:(n - 1)
            for i in 1:(n - l)
                ineq_lhs2[l + 1, i + l + 1] = i * (d1 - 1) + l * (D - 1) + 1
            end
        end
        ineq_rhs = [prod([d1 + (k - 1) * (D - 1) for k in 1:l]) * prod([i * (d1 - 1) + l * (D - 1) + 1 for i in 1:(n - l)])
                    for l in 0:(n - 1)]
        A = vcat(matrix(QQ, ineq_lhs1 + ineq_lhs2), -identity_matrix(QQ, n + 1))
        b = vcat(ineq_rhs, zeros(QQ, n + 1))
    end
    return collect(lattice_points(polyhedron(A, b)))
end

# ---------------------------------------------------------------------- #


# -------- Compute f_min using an ansatz equation -------- #

function solve_with_love_and_support(ode::ODE, p::Int; info = true)
    @assert is_probable_prime(p) "This is not a prime number, Yulia!"

    ode_mod_p = reduce_ode_mod_p(ode, p)
    x = sort(ode_mod_p.x_vars, rev = true)
    n = length(x)
    y = first(ode_mod_p.y_vars)
    F = Nemo.Native.GF(p)

    # compute Newton polytope of f_min
    possible_supp = f_min_support(ode)
                l = length(possible_supp)
    info && @info "The size of the estimates support is $(length(possible_supp))"
    nterms = length(possible_supp) + n

    # random initial conditions
    ic = Dict([x[i] => rand(1:p - 1) for i in 1:n]...)
    # no parameters, no inputs
    par = empty(ic)
    inp = empty(Dict(x[1] => [1]))

    ps_soltime = @elapsed ps_sol = power_series_solution(ode_mod_p, par, ic, inp, nterms)
    info && @info "Power series solution computed in $ps_soltime"

    start_system_time = time()
    pss = [ps_sol[y]]
    for i in 1:n
        push!(pss, ps_diff(pss[end]))
    end
        
    sort!(possible_supp, by = sum)
    prods = eval_at_support(possible_supp, pss)
                    
    ls = matrix([coeff(pr, j) for j in 0:(nterms - n - 1), pr in prods])
    info && @info "System created in $(time() - start_system_time)"

    info && @info "linear system dims $(size(ls))"
    
    system_soltime = @elapsed ker = kernel(ls, side=:right)
    info && @info "Linear system solved in $system_soltime"
            dim = size(kernel(ls, side=:right))[2]
    info && @info "The dimension of the solution space is $(dim)"

    start_constructing_time = time()

    R, _ = polynomial_ring(F, ["x1", ["x1^($i)" for i in 1:n]...])
            
    mons = [prod([gens(R)[k]^exp[k] for k in 1:ngens(R)]) for exp in possible_supp]    
    
    g = gcd([sum([s * m for (s, m) in zip(ker[:, i], mons)]) for i in 1:dim])
       
    info && @info "The resulting polynomial computes in $(time() - start_constructing_time)"

    return g
end
        

# -------------------------------------------------------- #

# -------- Faster evaluation of many monomials in power series -------- #

#friendly helper
function one_thing(e, vals, cacher)
    if !haskey(cacher, e)
        for i in 1:length(e)
            if e[i] != 0 #[1,?]
                    e_love = copy(e)
                    e_love[i] -= 1
                    res = vals[i] * one_thing(e_love, vals, cacher)
                    cacher[e] = res
            end
         end  
    end
    return cacher[e]
end    

function eval_at_support(supp, vals)
    s = sort(supp, by = sum) 
    cacher = Dict(s[1] => one(parent(first(vals))))
        for e in s[2:end] #index
            one_thing(e, vals, cacher)                 
        end
    return return [cacher[s] for s in supp]
end                    
                    

# ---------------------------------- #
