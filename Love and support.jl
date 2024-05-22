using Oscar
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
    return lattice_points(polyhedron(A, b))
end

function f_min_support(ode::ODE)
    x = sort(ode.x_vars, rev = true)
    n = length(x)
    gs = [ode.x_equations[xi] for xi in x]
    d1 = total_degree(gs[1])
    @assert d1 > 0 "d1 = 0"
    D = maximum(total_degree, gs[2:end])

    if d1 <= D
        ineq_lhs = reshape([1, [d1 + (k-1)*(D-1) for k in 1:n]...], 1, n+1)
        ineq_rhs = [prod([d1 + (k-1)*(D-1) for k in 1:n])]
        A = vcat(matrix(QQ, ineq_lhs), -identity_matrix(QQ, n+1))
        b = vcat(ineq_rhs, zeros(QQ, n+1))
    else
        ineq_lhs1 = [k <= l ? k*(D-1) + 1 : 0 for l in 0:(n-1), k in 0:n]
        ineq_lhs2 = zeros(Int, n, n+1)
        for l in 0:(n-1)
            for i in 1:(n-l)
                ineq_lhs2[l+1, i+l+1] = i*(d1-1) + l*(D-1) + 1
            end
        end
        ineq_rhs = [prod([d1 + (k-1)*(D-1) for k in 1:l])*prod([i*(d1-1) + l*(D-1) + 1 for i in 1:(n-l)])
                    for l in 0:(n-1)]
        A = vcat(matrix(QQ, ineq_lhs1 + ineq_lhs2), -identity_matrix(QQ, n+1))
        b = vcat(ineq_rhs, zeros(QQ, n+1))
    end
    return lattice_points(polyhedron(A, b))
end

# ---------------------------------------------------------------------- #


# -------- Compute f_min using an ansatz equation -------- #

function solve_with_love_and_support(ode::ODE, p::Int)
    @assert is_probable_prime(p) "This is not a prime number, Yulia!"

    ode_mod_p = reduce_ode_mod_p(ode, p)
    x = sort(ode_mod_p.x_vars, rev = true)
    n = length(x)
    y = first(ode_mod_p.y_vars)
    F = Nemo.Native.GF(p)

    # compute Newton polytope of f_min
    possible_supp = f_min_support(ode)
    @info "The size of the estimates support is $(length(possible_supp))"
    nterms = length(possible_supp) + n

    # random initial conditions
    ic = Dict([x[i] => rand(1:p-1) for i in 1:n]...)
    # no parameters, no inputs
    par = empty(ic)
    inp = empty(Dict(x[1] => [1]))

    ps_sol = power_series_solution(ode_mod_p, par, ic, inp, nterms)
    pss = [ps_sol[y]]
    for i in 1:n
        push!(pss, ps_diff(pss[end]))
    end
    prods = [prod([pss[k]^exp[k] for k in 1:length(pss)]) for exp in possible_supp]

    ls = matrix([coeff(pr, j) for j in 0:(nterms - n - 1), pr in prods])
    @info "linear system dims $(size(ls))"
    #ker = kernel(ls, side = :right)
                ker =    kernel(ls, side=:right)[2]
    @info "The dimension of the solution space is $(size(ker))"
    sol = ker[:, 1]

    R, _ = polynomial_ring(F, ["x1", ["x1^($i)" for i in 1:n]...])

    mons = [prod([gens(R)[k]^exp[k] for k in 1:ngens(R)]) for exp in possible_supp]

    return sum([s * m for (s, m) in zip(sol, mons)])
end

# -------------------------------------------------------- #


# -------- Test function -------- #

# test ansatz against IO enemy equation 
function check_ansatz(ode::ODE, p::Int)

    @info "Solving with love and support!"
    tim = @elapsed io_tocheck = solve_with_love_and_support(ode, p)
    @info "time: $(tim) seconds"

    @info "Solving without love and support :("
    tim = @elapsed io_correct = first(values(find_ioequations(ode)))
    io_correct = _reduce_mod_p(io_correct, p)
    @info "time: $(tim) seconds"

    R = parent(io_tocheck)
    S = parent(io_correct)

    # the variables of io_correct
    n = ngens(R)
    ys = gens(S)[2*(n-1)+2:3*(n-1)+2]
    @info "IO variables $(ys)"

    phi = hom(R, S, ys)

    quot, rem = divrem(phi(io_tocheck), io_correct)
    return iszero(rem) && iszero(total_degree(quot))
end

# ------------------------------- #


# -------- Helper functions -------- #

#Randomize general polynomial
function rand_poly(deg, vars)
    result = 0
    degs = [collect(0:deg) for v in vars]

    for m in IterTools.product(degs...)
        if sum(m) <= deg
            monom = rand(1:5)
            for i in 1:length(vars)
                monom *= vars[i]^m[i]
            end
            result += rand(1:1000)*monom
        end
    end

    return result
end

function rand_ode(degs::Vector{Int})
    n = length(degs)
    R, vars = polynomial_ring(QQ, vcat(["x$i(t)" for i in 1:n], ["y(t)"]))
    return StructuralIdentifiability.ODE{Ptype}(
        Dict(vars[i] => rand_poly(degs[i], vars[1:n]) for i in 1:n),
        Dict(vars[end] => vars[1]),
        Ptype[]
    )
end

# ---------------------------------- #





