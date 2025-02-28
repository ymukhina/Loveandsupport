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
    evaluate_at_point(F, dervs, point, minpoly_ord, support, var_to_sup, supp_to_index)

Evaluates the input interpolation point the same way as in the for loop of the original build_matrix_multipoint.
Outputs the row generated by this evaluation to set it in M.
The support is already ordered by sort_gleb_max!, i.e. ordered first by exponent of x1 in the support and then like the normal sort_gleb!
To save time and space, we pass var_to_sup and supp_to_index as arguments.

"""
# Gleb: I think there is also an assumption that the support contains the constant term, right?
# Max: Yes this assumption is reflected by setting row[1] = F(1) @ line 241
function evaluate_derivatives(F, dervs, point, minpoly_ord, var_to_sup, supp_to_index, set_x1 = false)
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

"""
    fill_matrix(points, dervs, minpoly_ord, support, mat_m, set_x1)

Creates and fills the block matrix M of size (length(support), mat_m).
The support is the reordered and/or reduced support containing only the relevant support vectors for the block in question.
dervs are the derivatives we pass as argument to evaluate_at_point and fill the matrix wrt their position in the support.
minpoly_ord is the order of the minimal polynomial.
set_x1 is a flag signaling which row we are handling and the vanish degree of the dual number used for that row

This function is essentially the part outside the for loop that would set up the everything necessary before evaluating and 
filling in the original build_matrix_multipoint.

Inside the following for loop, we now call a function evaluate_at_point that creates the respective row for some point and set it to be the row of M.

"""
function fill_matrix(points, dervs, minpoly_ord, support, mat_m, set_x1 = false)
    F = base_ring(parent(dervs[end]))
    support = [Vector{Int64}(p) for p in support]
    mat_n = length(points)
   
    M = Array{Any}(undef, mat_n, mat_m)

    for i in 1:mat_n
        if set_x1 == 1 || set_x1 == false
            M[i, 1] = F(1)
        else
            term = [F(0) for _  in 1:set_x1]
            term[1] = F(1)
            M[i, 1] = term
        end
    end
   
    supp_to_index = Dict(s => i for (i, s) in enumerate(support))
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1:(minpoly_ord + 1)]

    evals = evaluate_derivatives(F, dervs, points[1], minpoly_ord, var_to_sup, supp_to_index, set_x1)

    if set_x1 != false && set_x1 > 1
        evals_data = Dict{Int, Vector{Vector{elem_type(F)}}}()  # Middle block rows with dual numbers
    else
        evals_data = Dict{Int, Vector{elem_type(F)}}()  # First and Last Block Rows
    end

    evals_data[1] = evals
      
    for j in 1:(minpoly_ord + 1)
        supp = var_to_sup(j)
        if haskey(supp_to_index, supp)
            ind = supp_to_index[supp]
            M[1, ind] = evals[j]
        end
    end

    str = time()
   
    for i in 2:mat_n
        evals = evaluate_derivatives(F, dervs, points[i], minpoly_ord, var_to_sup, supp_to_index, set_x1)
        evals_data[i] = evals
        
        for j in 1:(minpoly_ord + 1)
            supp = var_to_sup(j)
            if haskey(supp_to_index, supp)
                ind = supp_to_index[supp]
                M[i, ind] = evals[j]
            end
        end
    end

    # println("Derivatives Evaluated in $(time() - str)")

    str = time()
    
    if set_x1 != false
        # Blocks with dual numbers
        for i in 1:mat_m
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
           
            while !haskey(supp_to_index, supp_divisor) && any(x -> x > 0, supp_divisor)
                nonzero_ind = findfirst(x -> x > 0, supp_divisor)
                if nonzero_ind === nothing
                    break
                end
                supp_divisor[nonzero_ind] -= 1
                multiplier[nonzero_ind] += 1
            end
           
            if !haskey(supp_to_index, supp_divisor)
                error("Unexpected situation: divisor not found in support")
            end
           
            supp_div_ind = supp_to_index[supp_divisor]
            mult_ind = get(supp_to_index, multiplier, -1)
           
            for j in 1:mat_n
                evals = evals_data[j]

                if mult_ind == -1
                    if evals isa Vector{fpFieldElem}
                        multiplier_eval = prod(evals .^ multiplier)
                    else
                        term = [F(0) for _ in 1:set_x1]
                        term[1] = F(1)
                        multiplier_eval = DiffMinPoly.DualNumber(term, set_x1, F)
                        for k in 1:length(evals)
                            term = evals[k]
                            e = DiffMinPoly.DualNumber(term, set_x1, F)
                            multiplier_eval = multiplier_eval * e ^ multiplier[k]
                        end
                    end
                else
                    multiplier_eval = M[j, mult_ind]
                end
               
                if multiplier_eval isa fpFieldElem
                    term = multiplier_eval * M[j, supp_div_ind]
                else
                    if multiplier_eval isa DiffMinPoly.DualNumber{fpFieldElem}
                        a = DualNumber(M[j, supp_div_ind], set_x1, F)
                        term = a * multiplier_eval
                    else
                        a, b = DualNumber(M[j, supp_div_ind], set_x1, F), DualNumber(multiplier_eval, set_x1, F)
                        term = a * b
                    end
                    term = get_terms(term)
                end
                M[j, i] = term
            end
        end

        if !(M[1, 1] isa fpFieldElem)
            for i in 1:mat_n
                for j in 1:mat_m
                    M[i, j] = M[i, j][end]
                end
            end
        end
        
        S = matrix_space(F, mat_n, mat_m)

        # println("Matrix filled in $(time() - str)")
        return S(M)
    else
        # No DualNumbers case
        for i in 1:mat_m
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
           
            while !haskey(supp_to_index, supp_divisor) && any(x -> x > 0, supp_divisor)
                nonzero_ind = findfirst(x -> x > 0, supp_divisor)
                if nonzero_ind === nothing
                    break
                end
                supp_divisor[nonzero_ind] -= 1
                multiplier[nonzero_ind] += 1
            end
           
            if !haskey(supp_to_index, supp_divisor)
                error("Unexpected situation: divisor not found in support")
            end
           
            supp_div_ind = supp_to_index[supp_divisor]
            mult_ind = get(supp_to_index, multiplier, -1)
           
            for j in 1:mat_n
                evals = evals_data[j]
               
                if mult_ind == -1
                    multiplier_eval = prod(evals .^ multiplier)
                else
                    multiplier_eval = M[j, mult_ind]
                end

                M[j, i] = M[j, supp_div_ind] * multiplier_eval
            end
        end
       
        S = matrix_space(F, mat_n, mat_m)
        # println("Matrix filled in $(time() - str)")

        return S(M)    
    end
end

"""
    make_matrix(n_vars, dervs, minpoly_ord, support, mat_n, mat_m, set_x1 = false)

Calls the various functions to generate the points needed for the evaluation method and creating the block matrix of size (mat_n, mat_m).
n_vars is the number of variables in the system and dervs are the Lie derivatives to evaluate at the points.
We pass the support and minpol_ord to fill_matrix to create the block matrix similarely to the original build_matrix_multipoint.
set_x1 is a flag to determine wether x1 should have a specific value, e.g. for now x1 = 0 or epsilon(1), a dual number vanishing for degree > 1.

This will be extended to allow setting x1 = epsilon(2) and build the rectangular matrices shaping the Block Lower Triangular matrix.
"""
function make_matrix(n_vars, dervs, minpoly_ord, support, mat_n, mat_m; set_x1 = false)
    F = base_ring(parent(dervs[end]))
    if set_x1 == false
        points = generate_points_base(F, mat_n, n_vars)
    else
        points = generate_points_dual(F, mat_n, n_vars, set_x1)
    end
    return fill_matrix(points, dervs, minpoly_ord, support, mat_m, set_x1)
end
