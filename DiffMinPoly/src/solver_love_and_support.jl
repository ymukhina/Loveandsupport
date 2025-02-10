using Oscar
using Nemo
using StructuralIdentifiability
using IterTools
import StructuralIdentifiability: reduce_ode_mod_p, var_to_str, switch_ring
using Random

const Ptype = QQMPolyRingElem


#Finding the kernel of [[A, 0],[B, C]]
function kernel_block(F, ls_UL, ls_LL, ls_LR) 
    strt = time()

    v1 = kernel(ls_UL, side=:right)
    
    #if A not full rank            # This is the chosen representation: [[A, 0], [B, C]]  a lower triangular in 4 blocks
    if rank(ls_UL) < size(ls_UL, 2)     # Start of a logic to handle linear combinations for sparse systems
                                        # Issue is that setting x1 = 0 makes many 0 appear in the matrix A
                               # Kernel(A) is then [[0,...,0],[1,0,...,0], ..., [0,...,0,1]] effectively describing all
                                      # possible vectors, hence the infinite solution triggering the error
                                                                            # "Unable to solve linear system"
        n_temp, m_temp = size(v1)
        ls_n, ls_m = size(ls_LR)
        S = matrix_space(F, ls_n + n_temp, ls_m + m_temp)
        aug = zero(S)      
        for i in 1:ls_n
            for j in 1:ls_m
                aug[i, j] = ls_LR[i, j]
            end
        end

        for i in 1:m_temp
            w_i = ls_LL * v1[:, i]
            for j in 1:n_temp
                aug[ls_n + j, ls_m + i] = w_i[j]
            end
        end     #Creating the augmented matrix [C|W] such that ker([C|W]) = [v2|lambda_coeffs]^T
                # Theoretically should satisfy Cv2 + Bv1' = 0 such that v1' = sum(lambda_i*v1_column(i))
                # Seems to be that for A not full rank, kernel of A is spanned by canonical vectors and [0,..,0]
                # This makes the set of solutions infinite, hence unsolvable
                # We can see that the v2 determined by this has linearly dependent rows
        v = kernel(aug, side=:right)
        v2, lambs = v[1:ls_m,:], v[ls_m+1:ls_m + m_temp, :]

        ker = map(1:size(lambs, 2)) do i
            v1_lamb = zero(v1[:, 1]) 
        
            for j in 1:size(v1, 2)
                v1_temp = zero(v1[:, j])
                for k in 1:size(v1, 1)    # Scale each element manually since     for j in 1:size(v1, 2)
                                                                                        #v1_lamb += lambs[j, i] * v1[:, j]
                                                                                # end
                                            # MethodError: no method matching *(::fpFieldElem, ::Vector{fpFieldElem})
                                            #Yet another naive workaround 
                    v1_temp[k] = lambs[j, i] * v1[k, j]
                end
                v1_lamb += v1_temp  
            end
        
            v2_i = v2[:, i]
            vcat(v1_lamb, v2_i)  # Recover the correct v1 in Ker(A) as a linear combination of basis vectors of ker(A)
        end


    else
        temp_LL = -(ls_LL * v1)      # Logic that used to work for dense systems
        v2 = solve(ls_LR, temp_LL, side=:right)
        ker = vcat(v1, v2)

    end

    @info "Linear System solved in $(time() - strt)"

    return ker

end


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

    # add a comment
                                                                                                                
    check = min_pol -> begin
        if isone(prob)
            is_zero_mod_ode(min_pol, ode, x)
        else
            is_zero_mod_ode_prob(min_pol, ode, x, prob)
        end
    end
    
    while check(minimal_poly) == false
        @info "Running Love & Support again. Incorrect Minimal Polynomial :("
        # error("STOP FOR TESTING")
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

    if n < 1  #Very temporary fix, but it seems to hold that Optimization1 removes too much info in non dense system of 2 variables.
        #Furthermore, the runtime improvement gained from the reduction isn't worth the loss in space and information for such systems
        #Especially considering that for small supports, e.g. here the planar case, the Original Love & Support is already very good.
        tim2 = @elapsed M = build_matrix_multipoint(F, ode_mod_p, x_mod_p, ord, possible_supp, info = info)

        info && @info "Building Linear System in: $(tim2)"

        strt = time()

        ker = kernel(M, side=:right)

        info && @info "Linear System solved in $(time() - strt)"

    else
        # This is a modified version of build_matrix_multipoint that should create only A, B and C 
        # without slicing and initializing M = [[A, 0], [B, C]]
        tim2 = @elapsed ls_UL, ls_LL, ls_LR = build_matrix_multipoint_big(F, ode_mod_p, x_mod_p, ord, possible_supp, info = info)

        info && @info "Building Linear System in: $(tim2)"

        info && @info "Linear Systems with dims $(size(ls_UL)), $(size(ls_LL)) and $(size(ls_LR))"

        ker = kernel_block(F, ls_UL, ls_LL, ls_LR)

    end


    if ndims(ker) == 1  # Check if ker is a 1D vector (n,)
        ker = reshape(ker, size(ker, 1), 1)  # Convert to (n,1)
    end

    dim = size(ker, 2)    

    info && @info "The dimension of the solution space is $(dim)"

    start_constructing_time = time()

    R, _ = polynomial_ring(F, [var_to_str(x), [var_to_str(x) * "^($i)" for i in 1:ord]...])
            
    mons = [prod([gens(R)[k]^exp[k] for k in 1:ngens(R)]) for exp in possible_supp]
    
    g = gcd([sum([s * m for (s, m) in zip(ker[:, i], mons)]) for i in 1:dim])  #This breaks for similar reason as * in kernel_block ((((((
                                                                            # I'll get to it tomorrow
    info && @info "The resulting polynomial computes in $(time() - start_constructing_time)"

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
    return sort_gleb_max!(collect(lattice_points(Oscar.polyhedron(A, b))))  # sort_gleb but with first the k1 [0, a, b] elements of the support
end                                                                         # then normal sort_gleb

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

# This function assumes that support contains the unit vectors and is sorted by `sort_gleb!`
function build_matrix_multipoint(F, ode, x, minpoly_ord, support; info = true)
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1: (minpoly_ord + 1) ]                                           
    n = length(ode.x_vars)
    @info "computing derivatives"
    dervs = lie_derivatives(minpoly_ord, ode, x)
    @info "done"

    support = [Vector{Int64}(p) for p in support]
    
    lsup = length(support)                                                    
    S = matrix_space(F, lsup, lsup)
    M = zero(S)
    supp_to_index = Dict(s => i for (i, s) in enumerate(support))

    # filling the columns corresponding to the derivatives
    for i in 1:lsup
        M[i, 1] = 1
        vec = [rand(F) for _ in 1:(n + 1) ]                                      
        evals = [derv(vec...) for derv in dervs]
        
        
        for j in 1:(minpoly_ord + 1)
            supp = var_to_sup(j)
            ind = supp_to_index[supp]
            M[i, ind] = evals[j]
        end
    end

    # filling the rest of the columns
    for i in (minpoly_ord + 3):lsup
        supp = support[i]
        supp_divisor = copy(supp)
        nonzero_ind = findfirst(x -> x > 0, supp_divisor)
        supp_divisor[nonzero_ind] -= 1                                                 
        multiplier = zeros(Int, minpoly_ord + 1)
        multiplier[nonzero_ind] += 1
        while !haskey(supp_to_index, supp_divisor)
            nonzero_ind = findfirst(x -> x > 0, supp_divisor)
            supp_divisor[nonzero_ind] -= 1
            multiplier[nonzero_ind] += 1
        end                                                    
        
        supp_div_ind = supp_to_index[supp_divisor]
        mult_ind = get(supp_to_index, multiplier, -1)
        for j in 1:lsup
            if mult_ind == -1
                multiplier_eval = prod(M[j, 2:(minpoly_ord + 2)] .^ multiplier)
            else
                multiplier_eval = M[j, mult_ind]
            end
            M[j, i] = M[j, supp_div_ind] * multiplier_eval           
        end
    end
      
  return M
end

# This function assumes that support contains the unit vectors and is sorted by `sort_gleb_max!`
function build_matrix_multipoint_big(F, ode, x, minpoly_ord, support; info = true)
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1: (minpoly_ord + 1) ]                                           
    n = length(ode.x_vars)

    tim2 = @elapsed dervs = lie_derivatives(minpoly_ord, ode, x)
    info && @info "Computing Derivatives in: $(tim2)"

    support = [Vector{Int64}(p) for p in support]
    lsup = length(support)  
    k1 = 0                     
    for i in 1:lsup
        if support[i][1] != 0
            k1 = i - 1       #last index that has [0, smthg, smthg, ...]
            break
        end
    end
    s1, s2 = support[1:k1], support[k1+1:lsup]


    S = matrix_space(F, k1, k1)
    A = zero(S)
    S = matrix_space(F, lsup - k1, k1)                 #[[A, 0], [B, C]]  a lower triangular in 4 blocks
    B = zero(S)                                     #Set up for Av1 = 0  &&  Bv1 + Cv2 = 0
    S = matrix_space(F, lsup - k1, lsup - k1)       # Issue arises when sparse sys ==> many terms in supp [0, a, b]
    C = zero(S)                                     # depend on x1, e.g. each term of x1' has x1 and each term of x1'' has x1 for example
    supp_to_index1 = Dict(s => i for (i, s) in enumerate(s1))  #first part of supp with x1^0
    supp_to_index2 = Dict(s => i for (i, s) in enumerate(s2))  #second part of supp with x1^k with k != 0


    # filling the columns corresponding to the derivatives   
    for i in 1:k1                 
        A[i, 1] = 1    
        vec = [rand(F) for _ in 1:(n + 1)]               #Similar to previous implementation
        vec[1] = F(0)                                   # with x1 = 0 

        evals = [derv(vec...) for derv in dervs]      

        for j in 2:(minpoly_ord + 1)            # since A shoudln't have terms containing x1,
            supp = var_to_sup(j)                # no need to look for the [1,0,0] in s1 
            ind = supp_to_index1[supp]   
            A[i, ind] = evals[j]          #Should fill in the [1, x1', x1'', combinations of x1' and x1'']
        end
    end

    for i in 1:(lsup - k1)                 
        B[i, 1] = 1    
        vec = [rand(F) for _ in 1:(n + 1)]         #Similar to A, but without setting x1 = 0

        evals = [derv(vec...) for derv in dervs]      

        for j in 2:(minpoly_ord + 1)          #Same logic as for A
            supp = var_to_sup(j)     
            ind = supp_to_index1[supp]   
            B[i, ind] = evals[j]        
        end
                  
        for j in 1:(lsup - k1)      #Using the same vec as for B to build C
            sup = s2[j]
            val = 1
            for s in 1:length(sup)
                val *= evals[s]^(sup[s])  #<-- should make this better eventually, very naive implementation for now
            end
            C[i, j] = val                   #Filling all of C with the same interpolation vectors as for B for consistency
        end       
    end


    # filling the rest of the columns
    for i in (minpoly_ord + 3):k1
        supp = support[i]
        supp_divisor = copy(supp)
        nonzero_ind = findfirst(x -> x > 0, supp_divisor)   #Rest of the logic is the same as before
        supp_divisor[nonzero_ind] -= 1                                                
        multiplier = zeros(Int, minpoly_ord + 1)
        multiplier[nonzero_ind] += 1
        while !(haskey(supp_to_index1, supp_divisor) || haskey(supp_to_index2, supp_divisor))
            nonzero_ind = findfirst(x -> x > 0, supp_divisor)
            supp_divisor[nonzero_ind] -= 1
            multiplier[nonzero_ind] += 1
        end                                                    

        if supp_divisor[1] == 0
            supp_div_ind = supp_to_index1[supp_divisor]
            mult_ind = get(supp_to_index1, multiplier, -1)
        else
            supp_div_ind = supp_to_index2[supp_divisor]
            mult_ind = get(supp_to_index2, multiplier, -1)
        end

        
        for j in 1:k1      #Simply treat the A and B cases. The C case was naively implemented above temporarely
            if mult_ind == -1
                multiplier_eval = prod(A[j, 2:(minpoly_ord + 2)] .^ multiplier)
            else
                multiplier_eval = A[j, mult_ind]
            end
            A[j, i] = A[j, supp_div_ind] * multiplier_eval          
        end

        for j in 1:(lsup - k1)
            if mult_ind == -1
                multiplier_eval = prod(B[j, 2:(minpoly_ord + 2)] .^ multiplier)
            else
                multiplier_eval = B[j, mult_ind]
            end
            B[j, i] = B[j, supp_div_ind] * multiplier_eval          
        end
    end

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
    sort!(exp_vectors, by = s -> (
    s[1] == 0 ? (0, s[2:end]...) : (1, sum(s), s[end:-1:1]...)
))
end

# ————————————————— #
