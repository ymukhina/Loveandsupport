using StructuralIdentifiability

"""
    is_linear(f)

This function checks if a polinomial f linear w.r.t. one of its variables, 
returns the index of this variable or -1 otherwise.

"""

function is_linear(f)
    
    R = parent(f)
    n = ngens(R)

    for i in 1:n-1
        if degree(f, i) == 1
            return true, i
        end
    end

    return false, -1
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


function is_pole(F, r, a)
    den = denominator(r)
    return den(a...) == F(0)
end

function has_denominator(r)
    try
        denominator(r)
        return true
    catch
        return false
    end
end


function generate_points_rational_parametrization(F, n_points::Int, n_vars::Int, dual_deg, rp = nothing)
    if dual_deg > 1
        vecs = Vector{Vector{DualNumber{fpFieldElem}}}(undef, n_points)
    else
        vecs = Vector{Vector{fpFieldElem}}(undef, n_points)
    end

    
    R = parent(rp[1])    
    nv = ngens(R)

    i = 1
    while i <= n_points

        t = [rand(F) for _ in 1:nv]

        bad = false
        for r in rp
            if has_denominator(r) && is_pole(F, r, t)
                bad = true
                break
            end
        end
        bad && continue   

       
        if dual_deg > 1
            vecs[i] = Vector{DualNumber{fpFieldElem}}(undef, n_vars + 1)
        else
            vecs[i] = Vector{fpFieldElem}(undef, n_vars + 1)
        end
    
        for j in 1:n_vars
            if dual_deg > 1
                vecs[i][j] = rp[j](t...) + rand(F) * Epsilon(Int(dual_deg), F)
            else
                vecs[i][j] = rp[j](t...)
            end
        end

        if dual_deg > 1
            vecs[i][n_vars + 1] = rand(F) * Epsilon(Int(dual_deg), F)
        else
            vecs[i][n_vars + 1] = rand(F)
        end

        i += 1
    end

    return vecs
end

# for the rational parametrisation generates the points [sol_1 + ε * rand(F) + ε^2 * rand(F) + ... + ε^k * rand(F), ...]

function generate_truncated_dual_points(F, k::Int, n_points::Int, n_vars::Int, rp = nothing)

    vecs = Vector{Vector{DualNumber{fpFieldElem}}}(undef, n_points)
    
    R = parent(rp[1])    
    nv = ngens(R)

    i = 1
    while i <= n_points

        t = [rand(F) for _ in 1:nv]

        bad = false
        for r in rp
            if has_denominator(r) && is_pole(F, r, t)
                bad = true
                break
            end
        end
        bad && continue   

        vecs[i] = Vector{DualNumber{fpFieldElem}}(undef, n_vars + 1)    
    
        for j in 1:n_vars
            vecs[i][j] = rp[j](t...) + sum(rand(F) * Epsilon(Int(k+1), F)^i for i in 1:k+1)
        end

        vecs[i][n_vars + 1] = rand(F) * Epsilon(Int(k+1), F)

        
        i += 1
    end

    return vecs
    
end


function search_rational_parametrization(f, linear_index)
    R = parent(f)
    gens_list = gens(R)
    n = length(gens_list)
    xi = gens_list[linear_index]
    rp = Vector{Any}(undef, n-1)  

    a = derivative(f, xi)
    b = f - xi*a

    F = fraction_field(R)
    xi_rp = -F(b) / F(a)

    gens_F = [F(g) for g in gens_list]

    for i in 1:n-1
        if i == linear_index
            rp[i] = xi_rp
        else
            rp[i] = gens_F[i]
        end
    end

    return rp
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

function split_index(support, hd)

    splits = split_supp(support, hd)
    pushfirst!(splits, 0)

    delta = [splits[i] - splits[i-1] for i in 2:length(splits)]

    popfirst!(splits)
    push!(splits, splits[end] + 1)

    smart_split = [delta[i] - delta[i+1] for i in length(delta)-1:-1:1]
    pushfirst!(smart_split, delta[end])

    push!(delta, 1)
    
    return delta, splits, smart_split

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

    # @info "KERKERKER" ker

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



