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
function eliminate(ode::ODE, x, prob = 0.99)
                                                                           
    @assert x in ode.x_vars
    minimal_poly, starting_prime = eliminate_with_love_and_support(ode, x, rand(2^25:2^32 - 1))
                                                                                                                
    check = min_pol -> begin
        if isone(prob)
            is_zero_mod_ode(min_pol, ode, x)
        else
            is_zero_mod_ode_prob(min_pol, ode, x, prob)
        end
    end
    
    while check(minimal_poly) == false
        @info "Incorrect minimal polynomial. Running Love & Support again (T_T)"
        starting_prime = Hecke.next_prime(starting_prime)
        minimal_poly, starting_prime = eliminate_with_love_and_support(ode, x, starting_prime)
    end

    return minimal_poly                                                                                    
end 


"""
    eliminate_with_love_and_support_modp(ode, x, p, ord, possible_supp)

Computes a polynomial of the order `ord` over a finite field F_p for the `x` variable of a polynomial ODE model `ode` (without inputs and parameters) with support `possible_supp`.

"""
                                                                                                                    
function eliminate_with_love_and_support_modp(ode::ODE, x, p::Int, ord::Int=minpoly_order(ode, x),
                                              possible_supp::Vector{PointVector{ZZRingElem}}=f_min_support(ode, x, ord); 
                                              info = true)
                                                                    
    @assert is_probable_prime(p) "This is not a prime number, Yulia!"

    ode_mod_p = reduce_ode_mod_p(ode, p)
    x_mod_p = switch_ring(x, ode_mod_p.poly_ring)
    n = length(ode_mod_p.x_vars)
    F = Nemo.Native.GF(p)     

    l = length(possible_supp)
    info && @info "The size of the estimates support is $(length(possible_supp))"

    tim2 = @elapsed dervs = lie_derivatives(ord, ode_mod_p, x_mod_p)
    info && @info "Computing Derivatives in: $(tim2)"

    tim2 = @elapsed ls_UL, ls_LL, ls_LR = build_matrix_multipoint(F, ode_mod_p, dervs, ord, possible_supp, info = info)
                                                                    
    info && @info "Building Linear System in: $(tim2)"

    info && @info "Linear Systems with dims $(size(ls_UL)), $(size(ls_LL)) and $(size(ls_LR))"

    ker = kernel_blocks(ls_UL, ls_LL, ls_LR)
    
    dim = size(ker)[2]
    info && @info "The dimension of the solution space is $(dim)"
    if dim > 1
        info && @info "Adding $(dim-1) rows to compensate for loss"
        strt = time()
        E = make_matrix(n, dervs, ord, possible_supp, dim - 1, l)  #Sol_space = 1 <==> No additional rows
        ker = solve_linear_combinations(E, ker)
        info && @info "Additional rows added and reduced solution space computed in $(time() - strt)"
        dim = size(ker, 2)
        info && @info "The dimension of the solution space is $(dim)"
    end


    strt = time()
    R, _ = polynomial_ring(F, [var_to_str(x), [var_to_str(x) * "^($i)" for i in 1:ord]...])
            
    mons = [prod([gens(R)[k]^exp[k] for k in 1:ngens(R)]) for exp in possible_supp]
    g = gcd([sum([s * m for (s, m) in zip(ker[:, i], mons)]) for i in 1:dim])
    
    info && @info "The resulting polynomial computes in $(time() - strt)"

    return g * (1 // Oscar.leading_coefficient(g))
end


"""
    eliminate_with_love_and_support(ode, x, p)

Computes the minimal polynomial for the `x` variable of a polynomial ODE model `ode` (without inputs and parameters)
using evaluation-interpolation approach over a finite field F_p.

"""

function eliminate_with_love_and_support(ode::ODE, x, starting_prime::Int)
    minpoly_ord = minpoly_order(ode, x) 
    possible_supp = f_min_support(ode, x, minpoly_ord)
    l_supp = length(possible_supp)
    R, _ = polynomial_ring(QQ, [var_to_str(x), [var_to_str(x) * "^($i)" for i in 1:minpoly_ord]...])
 
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
        sol_mod_p = eliminate_with_love_and_support_modp(ode, x, Int(p), minpoly_ord, possible_supp)  



        if is_first_prime
            # starting_prime = 65519
            # r = ZZ(starting_prime)                                        
            # @info "Chose $prim_cnt th prime $r, $(length(findall(is_stable))) stable coefficients"
            # sol_mod_p_2 = eliminate_with_love_and_support_modp(ode, x, Int(r), minpoly_ord, possible_supp) 
            # filter!(exp -> !iszero(coeff(sol_mod_p_2, Vector{Int}(exp))), possible_supp)
            # add_unit!(possible_supp, minpoly_ord)
            # l_supp = length(possible_supp)
            # resize!(sol_vector, l_supp)
            # resize!(crts, l_supp)
            # resize!(found_cand, l_supp)
            # resize!(is_stable, l_supp)
            # @info "updated support, new size is $(length(possible_supp))"

            filter!(exp -> !iszero(coeff(sol_mod_p, Vector{Int}(exp))), possible_supp)
            add_unit!(possible_supp, minpoly_ord)
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
                                

"""
    rand_poly(deg, vars)

Computes the polynomial of degree 'deg' in variables 'vars' 
with coefficient sampled uniformly at random from the integers in [-1000, 1000].

"""
                                
function rand_poly(deg, vars)
    result = 0
    degs = [collect(0:deg) for v in vars]

        for m in IterTools.product(degs...)
            if sum(m) <= deg
                monom = rand(1:5)
                for i in 1:length(vars)
                    monom *= vars[i]^m[i]
                end
               result += rand(1:20) * monom
            end
        end

    return result
end                                
                                
"""
    rand_ode(degs)

Computes the polynomial ODE model `ode` with the right hand side of degrees ‘degs’.
                                            
"""
                                                
function rand_ode(degs::Vector{Int}; char=0)
    n = length(degs)
    F = iszero(char) ? QQ : GF(char)
    R, vars = polynomial_ring(QQ, vcat(["x$i(t)" for i in 1:n], ["y(t)"]))
    return StructuralIdentifiability.ODE{Ptype}(
        vars[1:n],
        [vars[end]],
        Dict(vars[i] => rand_poly(degs[i], vars[1:n]) for i in 1:n),
        Dict(vars[end] => vars[1]),
        Ptype[]
    )
end


# -------- estimate support for f_min based on Theorem 1  -------- #

function f_min_support(ode::ODE, x, jacobian_rank::Int; info = true)
    n = jacobian_rank
    d1 = total_degree(ode.x_equations[x])
    @assert d1 > 0 "d1 = 0"
    D = maximum(total_degree, [eq for (v, eq) in ode.x_equations if v != x])
    D = max(D, 0)
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
        ineq_rhs = Vector{Int}(undef, n)
        for l in 0:(n-1)
            fac1 = prod(Vector{Int}([d1 + (k - 1) * (D - 1) for k in 1:l]))
            fac2 = prod(Vector{Int}([i * (d1 - 1) + l * (D - 1) + 1 for i in 1:(n - l)]))
            ineq_rhs[l+1] = fac1*fac2
        end
        A = vcat(matrix(QQ, ineq_lhs1 + ineq_lhs2), -identity_matrix(QQ, n + 1))
        b = vcat(ineq_rhs, zeros(QQ, n + 1))
    end
    return sort_gleb_max!(collect(lattice_points(Oscar.polyhedron(A, b))))
end

# -------- Functions for test of correctness -------- #
                    
function is_zero_mod_ode(pol, ode::ODE, x)
    start_time = time()                    
    n = length(ode.x_vars)
     
    dervs = lie_derivatives(n, ode, x) 
 
    res = pol(dervs...) 
                      
    @info "Checked membership deterministaically in $(time() - start_time) seconds"
    return iszero(res)
end
                        
function is_zero_mod_ode_prob(pol, ode::ODE, x, prob = 0.99) 
    start_time = time()
    n = length(ode.x_vars)
    ord = minpoly_order(ode, x) 
                                                                             
    lie_derivs = lie_derivatives(ord, ode, x) 
    D = [total_degree(d) for d in lie_derivs]                        
    deg_bnd = findmax([sum(m .* D) for m in Oscar.exponents(pol)])[1]
    
    N = Int(1 + ceil(deg_bnd / (1 - prob)))           
                                
    vec = [rand(1:N) for _ in 1:(n + 1)]
     
    evals = [deriv(vec...) for deriv in lie_derivs]                        
                       
    res = pol(evals...)
    @info "Checked membership probabilistically in $(time() - start_time) seconds"     
    return iszero(res)
end                            

# -------- Function for matrix construction for Ansatz -------- #

# This function assumes that support contains the unit vectors and is sorted by `sort_gleb_max!`
function build_matrix_multipoint(F, ode, dervs, minpoly_ord, support; info = true)
    n = length(ode.x_vars)
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1: (minpoly_ord + 1)]                #TODO: One function for 1 thing, refactor!!!!                             
                                                                                                #Make a functiont that generates the interpolation points
    lsup = length(support)  

    s1, s2, k1 = split_supp(support, 1)      #Make a new order prioritizing the first exponent while keeping sort_gleb logic in each subsupport


    A = make_matrix(n, dervs, minpoly_ord, s1, k1, k1, true)

    # v1 = kernel(A, side=:right)   #Check sol_space and adjust the additional interpolation points accordingly (extra rows for [B|C])

    # d = size(v1, 2) - 1

    BC = make_matrix(n, dervs, minpoly_ord, support, lsup - k1, lsup)  #Sol_space = 1 <==> No additional rows
    B, C = BC[:, 1:k1], BC[:, (k1 + 1):lsup]

    return A, B, C 
end

# -------- Auxiliary Functions -------- #

# return the order of the largest upper left non vanishing minor of the jacobian matrix                          
function minpoly_order(ode, x)
    n = length(ode.x_vars)  
                    
    dervs = lie_derivatives(n - 1, ode, x)
    J = jacobian_matrix(dervs)[1:n, :]
    res = rank(J)                                        
 
    return res                                                                                              
end                                            

function qq_to_mod(a::QQFieldElem, p)
    return numerator(a) * invmod(denominator(a), ZZ(p))
end


# One would suggest to hit the road...
function add_unit!(supp, jacobian_rank)
    l_supp = length(supp)
    for j in 1:(jacobian_rank + 2)         
        unit = [i == j ? one(ZZ) : zero(ZZ) for i in 1:(jacobian_rank + 1)] 
        !(unit in supp) && push!(supp, point_vector(ZZ, unit))  
    end                                                                                                                           
    l_supp < length(supp) && sort_gleb_max!(supp)
    return supp 
end                                                                                                                        

function lie_derivative(pol, ode)
    result = zero(pol)
        for v in vars(pol)
            result += derivative(pol, v) * ode.x_equations[v]
        end
   return result
end
      
        
function lie_derivatives(ord, ode, var)
    result = [var]            
    for i in 1:ord
        push!(result, lie_derivative(last(result), ode))
    end
    return result
end
            
function sort_gleb!(exp_vectors::Vector{PointVector{ZZRingElem}})
    sort!(exp_vectors, by = s -> [sum(s), s[end:-1:1]...])
end

function sort_gleb_max!(exp_vectors::Vector{PointVector{ZZRingElem}})
    sort!(exp_vectors, by = s -> (s[1] == 0 ? (0, sum(s), s[end:-1:1]...) : (1, sum(s), s[end:-1:1]...)))
end


# ————————————————— #
