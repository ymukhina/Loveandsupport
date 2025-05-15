using Oscar
using Nemo
using StructuralIdentifiability
using IterTools
import StructuralIdentifiability: reduce_ode_mod_p, var_to_str, switch_ring
using Random

const Ptype = QQMPolyRingElem

# -------- exported functions -------- #

"""
    eliminate(ode, prob = 0.99)

Computes the minimal polynomial for the only output variable of a polynomial ODE model `ode` (without inputs and parameters)
using evaluation-interpolation approach. The result is guaranteed to be correct with probability at least `prob`.
If `prob` is set to 1, the result is guaranteed to be correct.
"""
function eliminate(ode::ODE, prob = 0.99)
                                                                           
    @assert length(ode.y_vars) == 1
    minimal_poly, starting_prime = eliminate_with_love_and_support(ode, rand(2^25:2^32 - 1))

    # add a comment
                                                                                                                
    check = min_pol -> begin
        if isone(prob)
            is_zero_mod_ode(min_pol, ode)
        else
            is_zero_mod_ode_prob(min_pol, ode, prob)
        end
    end
    
    while check(minimal_poly) == false
        starting_prime = Hecke.next_prime(starting_prime)
        minimal_poly, starting_prime = eliminate_with_love_and_support(ode, starting_prime)
    end

    return minimal_poly                                                                                    
end 


"""
    eliminate_with_love_and_support_modp(ode, p, ord, possible_supp)

Computes a polynomial of the order `ord` over a finite field F_p for the only output of a polynomial ODE model `ode` (without inputs and parameters) with support `possible_supp`.

"""
                                                                                                                    
function eliminate_with_love_and_support_modp(
    ode::ODE,
    p::Int, 
    ord::Int = minpoly_order(ode),
    possible_supp::Vector{PointVector{ZZRingElem}} = f_min_support(ode, ord);
    info = true
)
                                                                    
    @assert is_probable_prime(p) "This is not a prime number, Yulia!"

    ode_mod_p = reduce_ode_mod_p(ode, p)
    n = length(ode_mod_p.x_vars)
    F = Nemo.Native.GF(p)     

    l = length(possible_supp)
    info && @info "The size of the estimates support is $(length(possible_supp))"

    tim2 = @elapsed ls = build_matrix_multipoint(F, ode_mod_p, ord, possible_supp, info = info)
                                                                    
    info && @info "eval method $(tim2)"

    info && @info "linear system dims $(size(ls))"
    
    system_soltime = @elapsed ker = kernel(ls, side=:right)
    info && @info "Linear system solved in $system_soltime"
    dim = size(ker)[2]
    info && @info "The dimension of the solution space is $(dim)"

    start_constructing_time = time()

    y_var = only(ode.y_vars)
    # TODO: this lines is duplicated - should be removed somewhere
    R, _ = polynomial_ring(
        F,
        vcat(
            [var_to_str(p) for p in ode.parameters],
            [var_to_str(y_var)],
            [var_to_str(y_var) * "^($i)" for i in 1:ord],
        )
    )
            
    mons = [prod([gens(R)[k]^exp[k] for k in 1:ngens(R)]) for exp in possible_supp]
    
    g = gcd([sum([s * m for (s, m) in zip(ker[:, i], mons)]) for i in 1:dim])
       
    info && @info "The resulting polynomial computes in $(time() - start_constructing_time)"

    return g * (1 // Oscar.leading_coefficient(g))
end


"""
    eliminate_with_love_and_support(ode, p)

Computes the minimal polynomial for the only output of a polynomial ODE model `ode` (without inputs and parameters)
using evaluation-interpolation approach over a finite field F_p.

"""

function eliminate_with_love_and_support(ode::ODE, starting_prime::Int)
    minpoly_ord = minpoly_order(ode) 
    possible_supp = f_min_support(ode, minpoly_ord)
    l_supp = length(possible_supp)
    y_var = only(ode.y_vars)
    R, _ = polynomial_ring(
        QQ,
        vcat(
            [var_to_str(p) for p in ode.parameters],
            [var_to_str(y_var)],
            [var_to_str(y_var) * "^($i)" for i in 1:minpoly_ord],
        )
    )
 
    prod_of_done_primes = one(ZZ)
    prim_cnt = 0
 
    sol_vector = zeros(QQ, l_supp)
    crts = zeros(ZZ, l_supp)
    found_cand = falses(l_supp)
    is_stable = falses(l_supp)
    
    is_first_prime = true
    while !all(is_stable)
        @label nxt_prm 
        

        starting_prime = Hecke.next_prime(starting_prime)
        p = ZZ(starting_prime)                                        
        prim_cnt += 1
        
        @info "Chose $prim_cnt th prime $p, $(length(findall(is_stable))) stable coefficients"
        sol_mod_p = eliminate_with_love_and_support_modp(ode, Int(p), minpoly_ord, possible_supp)  

        if is_first_prime
            filter!(exp -> !iszero(coeff(sol_mod_p, Vector{Int}(exp))), possible_supp)
            add_unit!(possible_supp)
            l_supp = length(possible_supp)
            resize!(sol_vector, l_supp)
            resize!(crts, l_supp)
            resize!(found_cand, l_supp)
            resize!(is_stable, l_supp)
            @info "updated support, new size is $(length(possible_supp))"

            is_first_prime = false
            starting_prime = Hecke.next_prime(rand(2^32:2^62))
        end
    
        sol_vector_mod_p = [coeff(sol_mod_p, Vector{Int}(exp)) for exp in possible_supp]
        for (i, a) in enumerate(sol_vector_mod_p)
            is_stable[i] && continue

            if found_cand[i]
                if Oscar.divides(denominator(sol_vector[i]), p)[1]
                    @info "bad prime, restarting"
                    @goto nxt_prm
                end
                sol_i_mod_p = qq_to_mod(sol_vector[i], p)
                if sol_i_mod_p == sol_vector_mod_p[i]
                    is_stable[i] = true
                    continue
                end
            end

            crts[i] = crt(Oscar.lift(ZZ, a), p, ZZ(crts[i]), prod_of_done_primes)
           
            succ, r, s = rational_reconstruction(
                crts[i],
                ZZ(p * prod_of_done_primes)
            )
            if succ
                sol_vector[i] = r//s
                found_cand[i] = true
            end
        end

        prod_of_done_primes *= p
    end 
    
    mons = [prod([gens(R)[k]^exp[k] for k in 1:ngens(R)]) for exp in possible_supp]
    g = sum([s * m for (s, m) in zip(sol_vector, mons)])
    return g, starting_prime
end 

# -------- estimate support for f_min  -------- #

function f_min_support(ode::ODE, jacobian_rank::Int; info = true)
    n = jacobian_rank
    m = length(ode.parameters)
    
    y = first(values(ode.y_equations))

    if (y in ode.x_vars) && m == 0
        @info "The output is a single variable and there are no parameters, using the refined bound"
        # Bound using Theorem 1 from https://arxiv.org/abs/2501.13680
        d = total_degree(ode.x_equations[y])
        @assert d > 0 "d = 0"
        D = maximum(total_degree, [eq for (v, eq) in ode.x_equations if v != y])
        D = max(D, 0)
        info && @info "We have d = $d and D = $D"
        if d <= D
            ineq_lhs = reshape([1, [d + (k - 1) * (D - 1) for k in 1:n]...], 1, n + 1)
            ineq_rhs = [prod([d + (k - 1) * (D - 1) for k in 1:n])]
            A_final = vcat(matrix(QQ, ineq_lhs), -identity_matrix(QQ, n + 1))
            b_final = vcat(ineq_rhs, zeros(QQ, n + 1))
        else
            ineq_lhs1 = [k <= l ? k * (D - 1) + 1 : 0 for l in 0:(n - 1), k in 0:n]
            ineq_lhs2 = zeros(Int, n, n + 1)
            for l in 0:(n - 1)
                for i in 1:(n - l)
                    ineq_lhs2[l + 1, i + l + 1] = i * (d - 1) + l * (D - 1) + 1
                end
            end
            ineq_rhs = Vector{Int}(undef, n)
            for l in 0:(n-1)
                fac1 = prod(Vector{Int}([d + (k - 1) * (D - 1) for k in 1:l]))
                fac2 = prod(Vector{Int}([i * (d - 1) + l * (D - 1) + 1 for i in 1:(n - l)]))
                ineq_rhs[l+1] = fac1 * fac2
            end
            A_final = vcat(matrix(QQ, ineq_lhs1 + ineq_lhs2), -identity_matrix(QQ, n + 1))
            b_final = vcat(ineq_rhs, zeros(QQ, n + 1))
        end
    else
        # bound from the new paper (todo: precise reference)
        @info "The output is not a single variable or there are parameters, using the general bound"
        dx, dp = subtotal_degree(y, ode.x_vars), subtotal_degree(y, ode.parameters)
        Dx = maximum(subtotal_degree(ode.x_equations[x], ode.x_vars) for x in ode.x_vars)
        Dp = maximum(subtotal_degree(ode.x_equations[x], ode.parameters) for x in ode.x_vars)
        Dx = max(Dx, 0)
        Dp = max(Dp, 0)
    
        info && @info "We have d = ($dx, $dp) and D = ($Dx, $Dp)"
    
        d = dx + dp
        D = Dx + Dp

        # parameters first
        ineq_lhs_bezout = reshape(vcat([1 for _ in ode.parameters], [d + (k - 1) * (D - 1) for k in 1:(n + 1)]), 1, n + m + 1)
        ineq_rhs_bezout = [prod([d + (k - 1) * (D - 1) for k in 1:(n + 1)])]

        ineq_lhs_param = reshape(vcat( [1 for _ in ode.parameters], [dp + (k - 1) * Dp  for k in 1:(n + 1)]), 1, n + m + 1)
        ineq_rhs_param = [sum((max(1,dp) + (i - 1) * Dp) * prod(dx + (j - 1) * (Dx - 1) for j in 1:(n + 1) if j != i) for i in 1:( n + 1))]
    
        ineq_lhs = reshape(vcat([0 for _ in ode.parameters], [dx + (k - 1) * (Dx - 1) for k in 1:(n + 1)]), 1, n + m + 1)
        ineq_rhs = [prod([dx + (k - 1) * (Dx - 1) for k in 1:(n + 1)])]
    
    
        A_bezout = vcat(matrix(QQ, ineq_lhs_bezout), -identity_matrix(QQ, n + m + 1))
        b_bezout = vcat(ineq_rhs_bezout, zeros(QQ, n + m + 1))

        A_param = vcat(matrix(QQ, ineq_lhs_param), -identity_matrix(QQ, n + m + 1))
        b_param = vcat(ineq_rhs_param, zeros(QQ, n + m + 1))

        A = vcat(matrix(QQ, ineq_lhs), -identity_matrix(QQ, n + m + 1))
        b = vcat(ineq_rhs, zeros(QQ, n + m + 1))

        A_final = vcat(A_bezout, A_param, A)
        b_final = vcat(b_bezout, b_param, b)

    end

    return sort_gleb!(collect(lattice_points(Oscar.polyhedron(A_final, b_final))))
end

# -------- Functions for test of correctness -------- #
                    
function is_zero_mod_ode(pol, ode::ODE)
    start_time = time()                    
    n = length(ode.x_vars)
     
    dervs = lie_derivatives(first(values(ode.y_equations)), ode, n)
 
    res = pol(dervs...)
                      
    @info "Checked membership deterministaically in $(time() - start_time) seconds"
    return iszero(res)
end
                        
function is_zero_mod_ode_prob(pol, ode::ODE, prob = 0.99) 
    start_time = time()
    n = length(ode.x_vars)
    m = length(ode.parameters)
    y = first(values(ode.y_equations))
    ord = minpoly_order(ode)
                                                                             
    lie_derivs = lie_derivatives(y, ode, ord) 
    D = vcat([1 for _ in ode.parameters], [total_degree(d) for d in lie_derivs])
    deg_bnd = findmax([sum(m .* D) for m in Oscar.exponents(pol)])[1]
    
    N = Int(1 + ceil(deg_bnd / (1 - prob)))
                                
    vec = [rand(1:N) for _ in 1:(n + m + 1)]
     
    evals = vcat([vec[findfirst(isequal(p), gens(parent(ode)))] for p in ode.parameters], [deriv(vec...) for deriv in lie_derivs])

    res = pol(evals...)
    @info "Checked membership probabilistically in $(time() - start_time) seconds"     
    return iszero(res)
end                            

# -------- Function for matrix construction for Ansatz -------- #

# This function assumes that support contains the unit vectors and is sorted by `sort_gleb!`
function build_matrix_multipoint(F, ode, minpoly_ord, support; info = true)
    y = first(values(ode.y_equations))
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1:(minpoly_ord + m + 1) ]                                           
    n = length(ode.x_vars)
    m = length(ode.parameters)
    dervs = lie_derivatives(y, ode, minpoly_ord)

    support = [Vector{Int64}(p) for p in support]
    
    lsup = length(support)                                                    
    S = matrix_space(F, lsup, lsup)
    M = zero(S)
    supp_to_index = Dict(s => i for (i, s) in enumerate(support))

    # filling the columns corresponding to the derivatives
    for i in 1:lsup
        M[i, 1] = 1
        vec = [rand(F) for _ in 1:(n + m + 1)] 
        evals = [derv(vec...) for derv in dervs]
               
        for j in 1:(minpoly_ord + m + 1)
            supp = var_to_sup(j)
            ind = supp_to_index[supp]
            if j > m
                M[i, ind] = evals[j - m]
            else
                M[i, ind] = vec[findfirst(isequal(ode.parameters[j]), gens(parent(ode)))]
            end
        end
    end

    # filling the rest of the columns
    for i in (minpoly_ord + m + 3):lsup
        supp = support[i]
        supp_divisor = copy(supp)
        nonzero_ind = findfirst(e -> e > 0, supp_divisor)
        supp_divisor[nonzero_ind] -= 1                                                 
        multiplier = zeros(Int, minpoly_ord + m + 1)
        multiplier[nonzero_ind] += 1
        while !haskey(supp_to_index, supp_divisor)
            nonzero_ind = findfirst(e -> e > 0, supp_divisor)
            supp_divisor[nonzero_ind] -= 1
            multiplier[nonzero_ind] += 1
        end                                                    
        
        supp_div_ind = supp_to_index[supp_divisor]
        mult_ind = get(supp_to_index, multiplier, -1)
        for j in 1:lsup
            if mult_ind == -1
                multiplier_eval = prod(M[j, 2:(minpoly_ord + m + 2)] .^ multiplier)
            else
                multiplier_eval = M[j, mult_ind]
            end
            M[j, i] = M[j, supp_div_ind] * multiplier_eval           
        end
    end
    return M
end

# -------- Auxiliary Functions -------- #

# return the order of the largest upper left non vanishing minor of the jacobian matrix                          
function minpoly_order(ode)
    n = length(ode.x_vars)
    y = first(values(ode.y_equations))
                    
    dervs = lie_derivatives(y, ode, n - 1)
    J = jacobian_matrix(dervs)[1:n, :]
    res = rank(J)                                        
 
    return res                                                                                              
end                                            

function qq_to_mod(a::QQFieldElem, p)
    return numerator(a) * invmod(denominator(a), ZZ(p))
end

function add_unit!(supp)
    l_supp = length(supp)
    dim = length(first(supp))
    for j in 1:(dim + 1)         
        unit = [i == j ? one(ZZ) : zero(ZZ) for i in 1:dim] 
        !(unit in supp) && push!(supp, point_vector(ZZ, unit))  
    end                                                                                                                           
    l_supp < length(supp) && sort_gleb!(supp)
    return supp 
end                                                                                                                        

function lie_derivative(poly, ode)
    result = zero(poly)
    for v in ode.x_vars
        result += derivative(poly, v) * ode.x_equations[v]
    end
    return result
end
      
function lie_derivatives(poly, ode, ord)
    result = [poly]
    for i in 1:ord
        push!(result, lie_derivative(last(result), ode))
    end
    return result
end
            
function sort_gleb!(exp_vectors::Vector{PointVector{ZZRingElem}})
    sort!(exp_vectors, by = s -> [sum(s), s[end:-1:1]...])
end

# ————————————————— #
