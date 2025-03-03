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

    high_deg = possible_supp[end][1]

    info && @info "Highest degree of x1 in the support is $high_deg"

    ks = split_supp(possible_supp, high_deg - 1)     #Try different number of splits

    info && @info "Splitting at indices $ks for System of size $l"

    solve_ker = 0
    build_mat = 0

    for i in 1:length(ks) + 1
        if i > length(ks)
            supp = possible_supp
        else
            supp = possible_supp[1:ks[i]]
        end

        if 1 < i && i <= length(ks)         # Neither first nor last
            n_rows = ks[i] - ks[i - 1]      
        elseif 1 < i                        # Last row case
            n_rows = l - ks[i - 1]
        else                                # First row case (i == 1)
            n_rows = ks[i]
        end

        if n_rows == 0
            n_rows = 1
        end
        
        # println("Row $i is $((n_rows/l)*100)% of Total Linear System")
        strt = time()
        if i > length(ks)                   # Allows to build only each block row one by one to not overload memory.
            ls = build_matrix_multipoint(n, dervs, ord, supp, n_rows, info = info)    # Last Block row
        else
            ls = build_matrix_multipoint(n, dervs, ord, supp, n_rows, vanish_deg = Int(i), info = info)   # All other block rows
        end    
        t = time() - strt
        build_mat += t

        strt = time()
        if i == 1
            ker = kernel(ls, side=:right)     #First row
        elseif i <= length(ks)
            ker = solve_linear_combinations(ls, ker, ks[i - 1])  # All subsequent rectangular blocks
        else
            ker = solve_linear_combinations(ls, ker, ks[end])              #Last row
        end
        t = time() - strt
        solve_ker += t
        dim = size(ker)[2]
        if dim > 1
            info && @info "The dimension of the $i th solution space is $(dim)"
        end
    
    end
                                                                    
    dim = size(ker)[2]
    info && @info "The dimension of the solution space is $(dim)"
    if dim > 1
        info && @info "Adding $(dim-1) rows to compensate for loss"
        strt = time()
        E = build_matrix_multipoint(n, dervs, ord, possible_supp, dim-1, info = info)
        t = time() - strt
        build_mat += t
        info && @info "Additional rows added in $(time() - strt)"
        strt = time()
        ker = solve_linear_combinations(E, ker)
        t = time() - strt
        solve_ker += t
        info && @info "Reduced solution space computed in $t"
        dim = size(ker, 2)
        info && @info "The dimension of the new solution space is $(dim)"
    end

    info && @info "Matrix building took $build_mat"
    info && @info "Kernel computation took $solve_ker"
    strt = time()
    R, _ = polynomial_ring(F, [var_to_str(x), [var_to_str(x) * "^($i)" for i in 1:ord]...])
            
    mons = [prod([gens(R)[k]^exp[k] for k in 1:ngens(R)]) for exp in possible_supp]
    g = gcd([sum([s * m for (s, m) in zip(ker[:, i], mons)]) for i in 1:dim])
    
    info && @info "The resulting polynomial computes in $(time() - strt)"

    return g * (1 // Oscar.leading_coefficient(g)), build_mat, solve_ker
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
    ker_t, mat_t = 0, 0
    
    is_first_prime = true
    while !all(is_stable)
        @label nxt_prm 
        

        starting_prime = Hecke.next_prime(starting_prime)
        p = ZZ(starting_prime)                                        
        prim_cnt += 1
        
        @info "Chose $prim_cnt th prime $p, $(length(findall(is_stable))) stable coefficients"
        sol_mod_p, m_t, k_t = eliminate_with_love_and_support_modp(ode, x, Int(p), minpoly_ord, possible_supp)  

        ker_t += k_t
        mat_t += m_t

        if is_first_prime

            filter!(exp -> !iszero(coeff(sol_mod_p, Vector{Int}(exp))), possible_supp)
            add_unit!(possible_supp, minpoly_ord)
            l_supp = length(possible_supp)
            resize!(sol_vector, l_supp)
            resize!(crts, l_supp)
            resize!(found_cand, l_supp)
            resize!(is_stable, l_supp)
            @info "Updated Support. New Size is $(length(possible_supp))"

            is_first_prime = false
            starting_prime = Hecke.next_prime(rand(2^32:2^62))
        end
    
        sol_vector_mod_p = [coeff(sol_mod_p, Vector{Int}(exp)) for exp in possible_supp]
        for (i, a) in enumerate(sol_vector_mod_p)
            is_stable[i] && continue

            if found_cand[i]
                if Oscar.divides(denominator(sol_vector[i]), p)[1]
                    @info "Bad Prime. Restarting..."
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
               result += rand(1:2) * monom
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
function build_matrix_multipoint(n, dervs, minpoly_ord, support, n_rows; vanish_deg = false, info = true)
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1: (minpoly_ord + 1) ]                                           
    F = base_ring(parent(dervs[end]))

    support = [Vector{Int64}(p) for p in support]
    
    lsup = length(support)       

    M = Array{Any}(undef, n_rows, lsup)

    supp_to_index = Dict(s => i for (i, s) in enumerate(support))

    if vanish_deg == false
        points = generate_points_base(F, n_rows, n)
    else
        points = generate_points_dual(F, n_rows, n, vanish_deg)
    end

    # filling the columns corresponding to the derivatives
    for i in 1:n_rows
        if vanish_deg == 1 || vanish_deg == false
            M[i, 1] = F(1)
        else
            term = [F(0) for _  in 1:vanish_deg]
            term[1] = F(1)
            M[i, 1] = term
        end

        evals = evaluate_derivatives(dervs, points[i], vanish_deg)
        
        
        for j in 1:(minpoly_ord + 1)
            supp = var_to_sup(j)
            if haskey(supp_to_index, supp)
                ind = supp_to_index[supp]
                M[i, ind] = evals[j]
            end
        end
    end

    # filling the rest of the columns
    if vanish_deg != false
        for i in 1:lsup
            supp = support[i]
            supp_divisor = copy(supp)
            nonzero_ind = findfirst(x -> x > 0, supp_divisor)
            if nonzero_ind === nothing
                continue
            end
    
            supp_divisor[nonzero_ind] -= 1
            if all(x -> x == 0, supp_divisor)
                continue
            end
    
            multiplier = zeros(Int, minpoly_ord + 1)
            multiplier[nonzero_ind] += 1
            while !haskey(supp_to_index, supp_divisor)
                nonzero_ind = findfirst(x -> x > 0, supp_divisor)
                if nonzero_ind === nothing
                    break
                end
                supp_divisor[nonzero_ind] -= 1
                multiplier[nonzero_ind] += 1
            end             
            
            if !haskey(supp_to_index, supp_divisor)
                error("Unexpected Situation: Divisor not found in support. Shouldn't happen with ordering.")
            end
            
            supp_div_ind = supp_to_index[supp_divisor]
            mult_ind = get(supp_to_index, multiplier, -1)
    
            for j in 1:n_rows
                if mult_ind == -1
                    v = []
                    for k in 1:(minpoly_ord + 1)
                        supp = var_to_sup(k)
                        if haskey(supp_to_index, supp)
                            ind = supp_to_index[supp]
                            val = M[j, ind]
                            if val isa Vector{fpFieldElem}
                                push!(v, DualNumber(val, vanish_deg, F))
                            else
                                push!(v, val)
                            end
                        else
                            push!(v, Epsilon(vanish_deg, F))
                        end
                    end
                    multiplier_eval = 1
                    for k in 1:length(multiplier)
                        multiplier_eval *= v[k]^multiplier[k]
                    end
                else
                    multiplier_eval = M[j, mult_ind]
                end
               
                if multiplier_eval isa fpFieldElem
                    term = multiplier_eval * M[j, supp_div_ind]
                else
                    R, (ε,) = polynomial_ring(F, ["ε"])
                    
                    matrix_dual = if M[j, supp_div_ind] isa fpFieldElem
                        DiffMinPoly.DualNumber{fpFieldElem}(R(M[j, supp_div_ind]), vanish_deg, F)
                    elseif M[j, supp_div_ind] isa Vector{fpFieldElem}
                        DiffMinPoly.DualNumber(M[j, supp_div_ind], vanish_deg, F)
                    else
                        M[j, supp_div_ind]
                    end
                    
                    if multiplier_eval isa DiffMinPoly.DualNumber{fpFieldElem}
                        term = matrix_dual * multiplier_eval
                    else
                        multiplier_dual = if M[j, supp_div_ind] isa Vector{fpFieldElem}
                            DiffMinPoly.DualNumber(multiplier_eval, vanish_deg, F)
                        else
                            DiffMinPoly.DualNumber{fpFieldElem}(R(multiplier_eval), vanish_deg, F)
                        end
                        term = matrix_dual * multiplier_dual
                    end
                    
                    if term isa DiffMinPoly.DualNumber
                        term = DiffMinPoly.get_terms(term)
                    end
                end
                
                M[j, i] = term

            end
        end

        if !(M[1, end] isa fpFieldElem)
            for i in 1:n_rows
                for j in 1:lsup
                    if !(M[i, j] isa fpFieldElem)
                        M[i, j] = M[i, j][end]
                    else
                        continue
                    end
                end
            end
        end
        
        S = matrix_space(F, n_rows, lsup)
          
        return S(M)

    else
        for i in 1:lsup
            supp = support[i]
            supp_divisor = copy(supp)
            nonzero_ind = findfirst(x -> x > 0, supp_divisor)
            if nonzero_ind === nothing
                continue
            end
    
            supp_divisor[nonzero_ind] -= 1
            if all(x -> x == 0, supp_divisor)
                continue
            end
    
            multiplier = zeros(Int, minpoly_ord + 1)
            multiplier[nonzero_ind] += 1
            while !haskey(supp_to_index, supp_divisor)
                nonzero_ind = findfirst(x -> x > 0, supp_divisor)
                if nonzero_ind === nothing
                    break
                end
                supp_divisor[nonzero_ind] -= 1
                multiplier[nonzero_ind] += 1
            end             
            
            if !haskey(supp_to_index, supp_divisor)
                error("Unexpected Situation: Divisor not found in support. Shouldn't happen with ordering.")
            end
            
            supp_div_ind = supp_to_index[supp_divisor]
            mult_ind = get(supp_to_index, multiplier, -1)
    
            for j in 1:n_rows
                if mult_ind == -1
                    v = []
                    for k in 1:(minpoly_ord + 1)
                        supp = var_to_sup(k)
                        if haskey(supp_to_index, supp)
                            ind = supp_to_index[supp]
                            push!(v, M[j, ind])
                        else
                            println("VERY BIG PROBLEM. SHOULD NEVER HAPPEN SINCE ALL x1, x1', x1'' ... SHOULD BE IN SUPPORT.")
                            push!(v, F(0))
                        end
                    end
                    multiplier_eval = prod(v .^ multiplier)
                else
                    multiplier_eval = M[j, mult_ind]
                end
                M[j, i] = M[j, supp_div_ind] * multiplier_eval           
            end
        end

        S = matrix_space(F, n_rows, lsup)

        return S(M)  
    end  
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
    sort!(exp_vectors, by = s -> (s[1], sum(s), s[end:-1:1]...))
end


# ————————————————— #
