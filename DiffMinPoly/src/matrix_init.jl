"""
    solve_linear_combinations(ls, sol_space, ks)

This function solves the kernel of the input linear system ls given the constraint from the previous solution space sol_space found.
Adapted to solve for the kernel of a matrix we split given the split indices ks

This function is used to find the kernel of all block rows in in the middle, e.g. neither the first, the last, nor the additional rows
"""
function solve_linear_combinations(ls, sol_space, ks=nothing)
    F = base_ring(ls)
    ls_n, ls_tot = size(ls)
    _, sol_cols = size(sol_space)
    
    if ks === nothing        # First row
        S = matrix_space(F, ls_n, sol_cols)
        aug = zero(S)
        
        for j in 1:sol_cols
            aug[:, j] = ls * sol_space[:, j]
        end
        
        v = kernel(aug, side=:right)
        
        ker = Matrix{eltype(sol_space)}(undef, size(sol_space, 1), size(v, 2))
        for i in 1:size(v, 2)
            ker[:, i] = sum(v[j, i] .* sol_space[:, j] for j in 1:sol_cols)
        end
        
        return matrix_space(F, size(ker, 1), size(ker, 2))(ker)
    end
    
    if isa(ks, Vector)       # Last row
        last_block_cols = ls_tot - ks[end]
        
        aug_cols = last_block_cols + sol_cols * length(ks)
        S = matrix_space(F, ls_n, aug_cols)
        aug = zero(S)
        
        aug[:, 1:last_block_cols] = ls[:, (ks[end]+1):end]
        
        col_idx = last_block_cols + 1
        
        for i in 1:length(ks)
            start_col = i > 1 ? ks[i-1] + 1 : 1
            end_col = ks[i]
            block_i = ls[:, start_col:end_col]
            
            for j in 1:sol_cols
                v_j = sol_space[start_col:end_col, j]
                aug[:, col_idx] = block_i * v_j
                col_idx += 1
            end
        end
    else                                # Middle rows
        last_block_cols = ls_tot - ks
        
        aug_cols = last_block_cols + sol_cols
        S = matrix_space(F, ls_n, aug_cols)
        aug = zero(S)
        
        aug[:, 1:last_block_cols] = ls[:, (ks+1):end]
        
        aug[:, last_block_cols+1:end] = ls[:, 1:ks] * sol_space
    end
    
    v = kernel(aug, side=:right)
    
    if size(v, 2) == 0
        return matrix_space(F, size(sol_space, 1) + last_block_cols, 0)(zeros(F, size(sol_space, 1) + last_block_cols, 0))
    end
    
    v_last = v[1:last_block_cols, :]
    lambdas = v[(last_block_cols+1):end, :]
    
    ker = Matrix{eltype(sol_space)}(undef, size(sol_space, 1) + size(v_last, 1), size(lambdas, 2))
    
    for i in 1:size(lambdas, 2)
        v1_lamb = sum(lambdas[j, i] .* sol_space[:, j] for j in 1:size(sol_space, 2))
        ker[:, i] = vcat(v1_lamb, v_last[:, i])
    end
    
    return matrix_space(F, size(ker, 1), size(ker, 2))(ker)
end


"""
    split_supp(supp, n_splits)

This function gives a list of the indices at which we can find the last element of the sections of the support that have the same highest order of x1.

For example [[0, 0, 0], [0, 1, 0], [0, 0, 1], ....., [0, ...], [1, 0, 0], [1, 1, 0], ...., and more ....] split with n_splits = 1
This will return a list containing [k1] such that support[k1] = [0, ....] and support[k1 + 1] = [1, 0, 0]

"""
 #Assuming ordered by sort_gleb_max! 
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
    if set_x1 > 1
        vecs = Vector{Vector{DiffMinPoly.DualNumber{fpFieldElem}}}(undef, n_points)
    else
        vecs = Vector{Vector{fpFieldElem}}(undef, n_points)
    end

    for i in 1:n_points
        if set_x1 > 1
            vec = [Epsilon(set_x1, F)]
            for _ in 1:n_vars
                terms = [F(0) for _ in 1:set_x1]
                terms[1] = rand(F)
                push!(vec, DualNumber(terms, set_x1, F))
            end
        else
            vec = [rand(F) for _ in 1:n_vars + 1]
            vec[1] = F(0)
        end
        vecs[i] = vec

    end
    return vecs
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
    n_vars = length(point)
    isolate_var = dervs[1]

    if set_x1 == 1   # Epsilon(1) = F(0) for row 1
        evals = [derv(point...) for derv in dervs]
        return evals
    
    elseif set_x1 == false   # No epsilon for last row
        evals = [derv(point...) for derv in dervs]
        return evals

    else   # All other rows i with epsilon(i)
        evals = [derv(point, isolate_var) for derv in dervs]
        
        e = Array{Any}(undef, length(evals))
        for i in 1:length(evals)
            a = get_terms(evals[i])
            e[i] = a
        end
        
        return e
    end
end

