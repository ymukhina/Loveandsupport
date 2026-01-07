using StructuralIdentifiability

"""
    is_linear(f)

This function checks if a polinomial f linear w.r.t. one of its variables, 
returns the index of this variable or 0 otherwise.

"""

function is_linear(f)

    n = nvars(parent(f))          
    remaining = trues(n)          
    check = zeros(Int, n)

    for exps in collect(Oscar.exponents(f))
        for i in 1:n
            if remaining[i] && exps[i] > 1
                remaining[i] = false
            else 
                remaining[i] = true
                check[i] += exps[i]
            end
        end
        any(remaining) || return nothing
    end

    for i in 1:n
        if remaining[i] && check[i] > 0
            return i 
        end
    end

    return 0

end

"""
    generate_points_base(F, n_points, n_vars, set_x1 = false)

This function generates n_points random interpolation points in the field F and returns them as a vector of vectors.

"""

function generate_points_base(F, n_points, n_vars)
    vecs = Vector{Vector{fpFieldElem}}(undef, n_points)

    for i in 1:n_points
        vec = [rand(F) for _ in 1:n_vars + 1]
        vecs[i] = vec
    end
    return vecs
end

"""
    generate_points_dual(F, n_points, n_vars, set_x1 = false)

This function generates n_points random interpolation dual number points in the field F and returns them as a vector of vectors.

The set_x1 flag is used to set the first variable to a specific value, e.g. 0 or epsilon(1) for the dual number vanishing degree.

"""

function generate_points_dual(F, n_points, n_vars, set_x1 = 0)
    R, _ = polynomial_ring(F, "ε")

    if set_x1 > 1
        vecs = Vector{Vector{DualNumber{fpFieldElem}}}(undef, n_points)
    else
        vecs = Vector{Vector{fpFieldElem}}(undef, n_points)
    end

    for i in 1:n_points
        if set_x1 > 1
            vec = [Epsilon(Int(set_x1), F)]
            
            for _ in 1:n_vars
                dual_poly = R(rand(F)) + rand(F) * Epsilon(Int(set_x1), F)
                push!(vec, dual_poly)
            end
        else
            vec = [rand(F) for _ in 1:n_vars + 1]
            vec[1] = F(0)
        end
        vecs[i] = vec
    end
    return vecs

end

function solve_linear_equation(F, n_vars, coeffs)

    sum_other = zero(F)
    sol = [rand(F) for _ in 1:n_vars+1]

    i = findfirst(!iszero, coeffs)

    if i === nothing
        error("All coefficients are zero, das ist not gut")
        # @info "CATASTRIFY!!!!"
    end

        for j in 1:n_vars
            if j != i
                sol[j] = rand(F)
                sum_other += qq_to_fp(F, coeffs[j]) * sol[j]
            end
        end
        sum_other += qq_to_fp(F, coeffs[end])
        sol[i] = -sum_other / qq_to_fp(F, coeffs[i])
        

    return sol
end

function is_pole(F, r, a)
    den = denominator(r)
    return den(a) == F(0)
end

function has_denominator(r)
    try
        denominator(r)
        return true
    catch
        return false
    end
end

"""
    generate_points_rational_parametrization(F, n_points::Int, dual_deg, rp)

For the rational parametrization rp = [rp[1](t), ... , rp[n](t)] generates n_points for t 
and evaluates rp on the chosen t. Use dual_deg > 1 for the evaluation in dual numbers

"""
function generate_points_rational_parametrization(F, n_points::Int, dual_deg, rp)

    R, _ = polynomial_ring(F, "ε")
    if dual_deg > 1
        vecs = Vector{Vector{DualNumber{fpFieldElem}}}(undef, n_points)
    else
        vecs = Vector{Vector{fpFieldElem}}(undef, n_points)
    end
    
    i = 1
    while i <= n_points
        t = rand(F)
        @info t
        bad = false
        for r in rp
            if has_denominator(r) && is_pole(F, r, t)
                bad = true
                break
            end
        end
        
        if bad
            continue  
        end
        
        if dual_deg > 1
            vecs[i] = [r(t) + rand(F) * Epsilon(Int(dual_deg), F) for r in rp]
        else
            vecs[i] = [r(t) for r in rp]
        end
        
        i += 1  
    end
    
    return vecs
end


function generate_points_dual_linear(F, n_points, n_vars, dual_deg, y_poly)
    R, _ = polynomial_ring(F, "ε")

    units = unit_vectors(F, n_vars+1)


    coeffs = [Oscar.coeff(y_poly, m) for m in units]
   
    if dual_deg > 1
        vecs = Vector{Vector{DualNumber{fpFieldElem}}}(undef, n_points)

        for i in 1:n_points
            sol = solve_linear_equation(F, n_vars, coeffs) 
            vec_dual = Vector{DualNumber{fpFieldElem}}(undef, n_vars+1)
            for k in 1:n_vars+1
                if !(coeffs[k] == 0)
                    vec_dual[k] =  rand(F) * Epsilon(Int(dual_deg), F) + sol[k]
                else
                    vec_dual[k] =  DualNumber{fpFieldElem}(R(sol[k]), dual_deg)
                end
            end
            vec_dual[n_vars+1] = DualNumber{fpFieldElem}(R(F(0)), dual_deg)
            vecs[i] = vec_dual
         end
        
    else
        vecs = Vector{Vector{fpFieldElem}}(undef, n_points)
        
        for i in 1:n_points
            vecs[i] = solve_linear_equation(F, n_vars, coeffs) 
        end
    end

return vecs
end


function unit_vectors(F, n)
    units = Vector{Int}[]
    
    for i in 1:n
        vec = [fp_to_int(F(0)) for _ in 1:n]
        vec[i] = fp_to_int(F(1))
        push!(units, vec)
    end

    push!(units, [fp_to_int(F(0)) for _ in 1:n])

    return units
end

function fp_to_int(x)
    return Int(x.data)
end

function qq_to_fp(F::fpField, x::QQFieldElem)
    return F(numerator(x)) * inv(F(denominator(x)))
end

function qq_to_fp(F::fpField, x::fpFieldElem)
    return x  # Already in the right field, just return it
end


"""
    split_supp(supp, n_splits)

This function gives a list of the indices at which we can find the last element of the sections of the support that have the same highest order of x1.

For example [[0, 0, 0], [0, 1, 0], [0, 0, 1], ....., [0, ...], [1, 0, 0], [1, 1, 0], ...., and more ....] split with n_splits = 1
This will return a list containing [k1] such that support[k1] = [0, ....] and support[k1 + 1] = [1, 0, 0]

"""
# Assuming ordered by sort_gleb_max! 
function split_supp(supp, n_splits)
    ks = Int[]
    start_idx = 1
    
    for target_ord in 0:(n_splits - 1)
        next_idx = findfirst(i -> supp[i][1] > target_ord, start_idx:length(supp))
        if isnothing(next_idx)
            push!(ks, length(supp))
            break
        else
            next_idx = next_idx + start_idx - 1
            push!(ks, next_idx - 1)
            start_idx = next_idx
        end
    end
    
    return ks
end

function evaluate_polynomial(dervs, point, vanishing_deg)
   
    if (vanishing_deg < 2) || (vanishing_deg == false)
        eval = [derv(point...) for derv in dervs]
    else   
        eval = [derv(point) for derv in dervs]
    end
    

    return eval
end


function construct_result_polynomial(ode, ker, dim, possible_supp, ord, F; info=true)

    start_constructing_time = time()
    y_var = only(ode.y_vars)


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

function sort_gleb_max!(exp_vectors::Vector{PointVector{ZZRingElem}})
    sort!(exp_vectors, by = s -> (s[1], sum(s), s[end:-1:1]...))
end


# ————————————————— #



