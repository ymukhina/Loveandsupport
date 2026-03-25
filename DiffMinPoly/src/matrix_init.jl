using Oscar 

"""
    solve_linear_combinations(ls, sol_space, ks)

This function solves the kernel of the input linear system ls given the constraint from the previous solution space sol_space found.
Adapted to solve for the kernel of a matrix we split given the splitting index ks
"""

function solve_linear_combinations(ls, sol_space, ks=nothing)
    F = base_ring(ls)
    ls_n, ls_tot = size(ls)
    _, sol_cols = size(sol_space)
   
    if ks === nothing
        S = matrix_space(F, ls_n, sol_cols)
        aug = ls * sol_space
        
        v = kernel(aug, side=:right)
        
        if size(v, 2) > 0
            ker = sol_space * v
            return ker
        else
            return zero_matrix(F, size(sol_space, 1), 0)
        end
    end
   
    last_block_cols = ls_tot - ks
    aug_cols = last_block_cols + sol_cols
    S = matrix_space(F, ls_n, aug_cols)
    aug = zero(S)
    
    aug[:, 1:last_block_cols] = ls[:, (ks+1):end]
    aug[:, last_block_cols+1:end] = ls[:, 1:ks] * sol_space
   
    v = kernel(aug, side=:right)
   
    if size(v, 2) == 0
        return zero_matrix(F, size(sol_space, 1) + last_block_cols, 0)
    end
   
    v_last = v[1:last_block_cols, :]
    lambdas = v[(last_block_cols+1):end, :]
   
    sol_result = sol_space * lambdas
    ker = vcat(sol_result, v_last)

    return ker
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
                dual_poly = R(rand(F))
                push!(vec, DualNumber{fpFieldElem}(dual_poly, set_x1))
            end
        else
            vec = [rand(F) for _ in 1:n_vars + 1]
            vec[1] = F(0)
        end
        vecs[i] = vec
    end
    return vecs
end


function generate_points_dual_linear(F, n_points, n_vars, set_x1, y_poly)

    R, _ = polynomial_ring(F, "ε")
    poly_ring = parent(y_poly)
    total_vars = ngens(poly_ring)
    
    units = unit_vectors(F, n_vars+1)
    coeffs = [coeff(y_poly, m) for m in units]
    M = matrix(F, 1, length(coeffs), coeffs)
   
    const_term = Oscar.coeff(y_poly, [fp_to_int(F(0)) for _ in 1:n_vars+1])
    
    if !(const_term == 0)
        const_term = qq_to_fp(F, const_term)
        sol = [F(0) for _ in 1:n_vars+1]
        found = false
        for i in 1:(n_vars+1)
            if !iszero(coeffs[i])
                sol[i] = - const_term / qq_to_fp(F, coeffs[i])
                found = true
                break
            end
        end
        
        if !found && !iszero(const_term)
            error("No solution exists for the linear equation")
        end
    else
        sol = [F(0) for _ in 1:n_vars+1]
    end
        ker = kernel(M, side=:right)

    if set_x1 > 1
        #Dual number case
        vecs = Vector{Vector{DualNumber{fpFieldElem}}}(undef, n_points)
        for i in 1:n_points
            vec_base = zeros(F, total_vars)
            
            if ncols(ker) > 0
                random_coeffs = [rand(F) for _ in 1:n_vars+1]
                for j in 1:ncols(ker)
                    for k in 1:n_vars+1 
                        vec_base[k] += random_coeffs[j] * ker[k, j]
                    end
                end
            else
                vec_base = [rand(F) for _ in 1:n_vars+1]
            end

            # Convert to dual numbers by multiplying each component by epsilon
            vec_dual = Vector{DualNumber{fpFieldElem}}(undef, n_vars+1)
            for k in 1:n_vars
                vec_dual[k] =  (Epsilon(Int(set_x1), F)) * (1 - const_term) * vec_base[k]
            end
            dual_poly = R(rand(F))
            vec_dual[n_vars + 1] = DualNumber{fpFieldElem}(dual_poly, set_x1)
            vecs[i] = vec_dual
        end
    else
        # Regular field elements case
        vecs = Vector{Vector{fpFieldElem}}(undef, n_points)
        for i in 1:n_points
            vec = zeros(F, total_vars)
            
            if (ncols(ker) > 0)
                random_coeffs = [rand(F) for _ in 1:ncols(ker)]
                for j in 1:ncols(ker)
                    for k in 1:n_vars+1   
                        vec[k] += random_coeffs[j] * ker[k, j]
                    end
                end
            else
                vec = [rand(F) for _ in 1:n_vars+1]
            end
            vec = vec + sol
            vecs[i] = vec
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

   # push!(units, [fp_to_int(F(0)) for _ in 1:n])

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
    evaluate_at_point(dervs, point, set_x1 = false)

Evaluates the input interpolation point the same way as in the for loop of the original build_matrix_multipoint.
Outputs the row generated by this evaluation to set it in M.
The support is already ordered by sort_gleb_max!, i.e. ordered first by exponent of x1 in the support and then like the normal sort_gleb!
set_x1 is a flag signaling the vanishing degree of the dual number
"""
# Gleb: I think there is also an assumption that the support contains the constant term, right?
# Max: Yes, this assumption is reflected by setting row[1] = F(1) 
function evaluate_derivatives(dervs, point, set_x1 = false)
    isolate_var = dervs[1]

    if set_x1 == 1   # Epsilon(1) = F(0) for row 1
        evals = [derv(point...) for derv in dervs]
        return evals
    
    elseif set_x1 == false   # No epsilon for last row
        evals = [derv(point...) for derv in dervs]
        return evals

    else   # All other rows with different vanishing degrees
        dual_point = [p isa DualNumber ? p : point[1] for p in point]
        
        evals = [derv(dual_point) for derv in dervs]
        
        return evals
    end
end