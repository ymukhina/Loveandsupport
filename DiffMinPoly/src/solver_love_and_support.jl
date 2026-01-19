using Oscar
using Nemo
using StructuralIdentifiability
using IterTools
import StructuralIdentifiability: reduce_ode_mod_p, var_to_str, switch_ring
using Random

const Ptype = QQMPolyRingElem

# -------- exported functions -------- #

"""
    eliminate(ode, x, prob = 0.99)

Computes the minimal polynomial for the `x` variable of a polynomial ODE model `ode` (without inputs and parameters)
using evaluation-interpolation approach. The result is guaranteed to be correct with probability at least `prob`.
If `prob` is set to 1, the result is guaranteed to be correct.
"""
function eliminate(ode::ODE, prob = 0.99)
                                                                           
    @assert length(ode.y_vars) == 1
    minimal_poly, starting_prime = eliminate_with_love_and_support(ode, rand(2^25:2^32 - 1))
                                                                                                                
    check = min_pol -> begin
        if isone(prob)
            is_zero_mod_ode(min_pol, ode)
        else
            is_zero_mod_ode_prob(min_pol, ode, prob)
        end
    end
    
    while check(minimal_poly) == false
        @info "Incorrect minimal polynomial. Running Love & Support again (T_T)"
        starting_prime = Hecke.next_prime(starting_prime)
        minimal_poly, starting_prime = eliminate_with_love_and_support(ode, starting_prime)
    end

    return minimal_poly                                                                                    
end 


"""
    eliminate_with_love_and_support_modp(ode, x, p, rational_param, ord, possible_supp)

Computes a polynomial of the order `ord` over a finite field F_p for the `x` variable of a polynomial ODE model `ode` (without inputs and parameters) with support `possible_supp`.
If you know the rational parametrisation of the observation function, add it to rational_param

"""
function eliminate_with_love_and_support_modp(ode::ODE, p::Int, rational_param=nothing, ord::Int=minpoly_order(ode),
    possible_supp::Vector{PointVector{ZZRingElem}}=f_min_support(ode, ord); 
    info = true)
        
    @assert is_probable_prime(p) "This is not a prime number, Yulia!"  

    #setup modular enviroment
    ode_mod_p = StructuralIdentifiability.reduce_ode_mod_p(ode, p)
    R_new = ode_mod_p.poly_ring
    F = Nemo.Native.GF(p)  

    n = length(ode_mod_p.x_vars)
    m = length(ode.parameters)
    
    y_poly = first(values(ode.y_equations))
    y_poly = R_new(y_poly)

    l = length(possible_supp)
    
    high_deg = possible_supp[end][1]

    if rational_param === nothing
        if is_linear(y_poly)[1] 
            rational_param = search_rational_parametrization(y_poly, is_linear(y_poly)[2])
        else 
            rational_param = false
        end
    end

    if !(rational_param == false)
    @info "Rational parametrization case"
           possible_supp = sort_gleb_max!(possible_supp)
           dervs = lie_derivatives(y_poly, ode_mod_p, ord)
           splits = split_supp(possible_supp, high_deg)

           ker, dim, build_mat, solve_ker = solve_matrix(F, ode, n, dervs, ord, possible_supp, splits, l, rational_param; info=true) 
   
           info && @info "Matrix building took $build_mat"
           info && @info "Kernel computation took $solve_ker"
   
           result = construct_result_polynomial(ode, ker, dim, possible_supp, ord, F, info=info)
           return result, build_mat, solve_ker
    else 
        @info "General case"
                possible_supp = sort_gleb!(possible_supp)
                dervs = lie_derivatives(y_poly, ode_mod_p, ord)
                ker, dim = solve_matrix_general(F, ode, n, m, dervs, ord, possible_supp; info=true)
                result = construct_result_polynomial(ode, ker, dim, possible_supp, ord, F, info=info)
                return result, 0, 0
    end

end


"""
    eliminate_with_love_and_support(ode, x, p)

Computes the minimal polynomial for the `x` variable of a polynomial ODE model `ode` (without inputs and parameters)
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

    m = length(ode.parameters)
    degree_y_poly = total_degree(first(values(ode.y_equations)))
    # x = first(values(ode.y_equations))
    high_deg = possible_supp[end][1]
    # d = total_degree(x)

    # use_optimized = (x == ode.x_vars[1]) && (m == 0) && (high_deg > 2)
    # linear_y = (d == 1) && (m == 0) && !(x == ode.x_vars[1]) && (high_deg > 2)  
 
    linear_optimized_case = (degree_y_poly == 1) && (m == 0) && (high_deg > 2) 

    prod_of_done_primes = one(ZZ)
    prim_cnt = 0

    sol_vector = zeros(QQ, l_supp)
    crts = zeros(ZZ, l_supp)
    found_cand = falses(l_supp)
    is_stable = falses(l_supp)
    ker_t, mat_t = 0, 0
    
    is_first_prime = true
    while !all(is_stable)
        @label nxt_prm 
        

        starting_prime = Hecke.next_prime(starting_prime)
        p = ZZ(starting_prime)                                        
        prim_cnt += 1
        
        @info "Chose $prim_cnt th prime $p, $(length(findall(is_stable))) stable coefficients"
        sol_mod_p, m_t, k_t = eliminate_with_love_and_support_modp(ode, Int(p), minpoly_ord, possible_supp)  

        ker_t += k_t
        mat_t += m_t

        if is_first_prime
            filter!(exp -> !iszero(coeff(sol_mod_p, Vector{Int}(exp))), possible_supp)
            add_unit!(possible_supp) 
            l_supp = length(possible_supp)
            if linear_optimized_case
                possible_supp = sort_gleb_max!(possible_supp)
            else
                possible_supp = sort_gleb!(possible_supp)
            end
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

    println("#------------------------------------------------------------#")
    @info "Overall Matrix Building: $mat_t"
    @info "Overall Kernel Computation: $ker_t"
    println("#------------------------------------------------------------#")
    
    mons = [prod([gens(R)[k]^exp[k] for k in 1:ngens(R)]) for exp in possible_supp]
    g = sum([s * m for (s, m) in zip(sol_vector, mons)])
    return g, starting_prime
end 
                   
# -------- estimate support for f_min based on Theorem 1  -------- #

function f_min_support(ode::ODE, jacobian_rank::Int; info = true)
    n = jacobian_rank
    m = length(ode.parameters)
  
    y = first(values(ode.y_equations))

    if (y == ode.x_vars[1]) && m == 0
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

    elseif  m == 0
        d = maximum(total_degree(ode.y_equations[ode.y_vars[1]]))
        D = maximum(total_degree(ode.x_equations[x]) for x in ode.x_vars)
        D = max(D, 0)
    
        info && @info "We have d = $d and D = $D"
    
        ineq_lhs = reshape([d + (k - 1) * (D - 1) for k in 1:n+1], 1, n+1)
        ineq_rhs = [prod([d + (k - 1) * (D - 1) for k in 1:n+1])]
    
        A_final = vcat(matrix(QQ, ineq_lhs), -identity_matrix(QQ, n + 1))
        b_final = vcat(ineq_rhs, zeros(QQ, n + 1))
    else
        # Bound using Theorem 1 from https://arxiv.org/abs/2506.08824
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
        ineq_rhs_param = [sum((dp + (i - 1) * Dp) * prod(dx + (j - 1) * (Dx - 1) for j in 1:(n + 1) if j != i) for i in 1:( n + 1))]

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
    return sort_gleb_max!(collect(lattice_points(Oscar.polyhedron(A_final, b_final))))
end

