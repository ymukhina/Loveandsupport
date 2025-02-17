function kernel_blocks(ls_U, ls_L) 
    strt = time()
    v1 = kernel(ls_U, side=:right)

    ker = solve_linear_combinations(ls_L, v1, size(ls_U, 2))
    @info "Linear System solved in $(time() - strt)"

    return ker
end

# Adapt these to solve for the kernel by splitting M into D blocks for the highest order of x1 in the support

function solve_linear_combinations(ls, sol_space, k = 0)
    F = base_ring(ls)
    n_temp, m_temp = size(sol_space)
    ls_n, ls_tot = size(ls)
    
    if k > 0
        ls_L = ls[:, 1:k]
        ls_R = ls[:, (k+1):ls_tot]
        ls_m = ls_tot - k
    else
        ls_L = ls
        ls_R = zero(matrix_space(F, ls_n, 0))
        ls_m = 0
    end
    
    S = matrix_space(F, ls_n, ls_m + m_temp)
    aug = zero(S)    
    
    if k > 0
        aug[:, 1:ls_m] = ls_R
    end
    
    for j in 1:m_temp
        w_i = ls_L * sol_space[:, j]
        aug[:, ls_m + j] = w_i
    end
    
    v = kernel(aug, side=:right)
    if k > 0
        v2, lambs = v[1:ls_m,:], v[ls_m+1:ls_m + m_temp, :]
        ker = Matrix{eltype(sol_space)}(undef, size(sol_space, 1) + size(v2, 1), size(lambs, 2))
        for i in 1:size(lambs, 2)
            v1_lamb = sum(lambs[j, i] .* sol_space[:, j] for j in 1:size(sol_space, 2))
            ker[:, i] = vcat(v1_lamb, v2[:, i])
        end
    else
        ker = Matrix{eltype(sol_space)}(undef, size(sol_space, 1), size(v, 2))
        for i in 1:size(v, 2)
            ker[:, i] = sum(v[j, i] .* sol_space[:, j] for j in 1:size(sol_space, 2))
        end
    end
    
    return matrix_space(F, size(ker,1), size(ker,2))(ker)
end


 #Assuming ordered by sort_gleb_max! 
 function split_supp(supp, n_splits)
    ks = Int[]
    start_idx = 1
    
    for target_ord in 0:(n_splits)
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

function generate_points(F, n_points, n_vars, set_x1 = false)
    vecs = Vector{Vector{elem_type(F)}}(undef, n_points)
    for i in 1:n_points
        vec = [rand(F) for _ in 1:(n_vars + 1)]
        if set_x1
            vec[1] = F(0)
        end
        vecs[i] = vec
    end
    return vecs
end

function evaluate_at_point(dervs, point, minpoly_ord, support, var_to_sup, memo_eval)
    evals = [derv(point...) for derv in dervs]
    
    for j in 1:(minpoly_ord + 1)
        supp = var_to_sup(j)
        memo_eval[supp] = evals[j]
    end
    
    function decompose_multiplier(mult)
        if haskey(memo_eval, mult)
            return memo_eval[mult]
        end
        
        mult2 = zeros(Int, length(mult))
        nonzero_ind = findfirst(x -> x > 0, mult)
        
        if isnothing(nonzero_ind)
            return one(parent(evals[1]))
        end
        
        mult_copy = copy(mult)
        mult_copy[nonzero_ind] -= 1
        mult2[nonzero_ind] += 1
        
        val1 = decompose_multiplier(mult_copy)
        val2 = decompose_multiplier(mult2)
        
        result = val1 * val2
        memo_eval[mult] = result
        return result
    end
    
    row_values = Dict{Vector{Int}, elem_type(parent(evals[1]))}()
    row_values[support[1]] = one(parent(evals[1]))
    
    for supp in support[2:end]
        if haskey(memo_eval, supp)
            row_values[supp] = memo_eval[supp]
            continue
        end
        
        res = zero(parent(evals[1]))
        sup_elem = copy(supp)
        multiplier = zeros(Int, minpoly_ord + 1)
        
        while true
            nonzero_ind = findfirst(x -> x > 0, sup_elem)
            isnothing(nonzero_ind) && break
            
            sup_elem[nonzero_ind] -= 1
            multiplier[nonzero_ind] += 1
            
            if haskey(memo_eval, sup_elem)
                mult_val = decompose_multiplier(multiplier)
                res += mult_val * memo_eval[sup_elem]
                break
            end
        end
        
        memo_eval[supp] = res
        row_values[supp] = res
    end
    
    return row_values
end

function fill_matrix(points, dervs, minpoly_ord, support, mat_m)
    F = base_ring(parent(dervs[end]))
    support = [Vector{Int64}(p) for p in support]
    mat_n = length(points)
    
    S = matrix_space(F, mat_n, mat_m)
    M = zero(S)
    
    supp_to_index = Dict(s => i for (i, s) in enumerate(support))
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1:(minpoly_ord + 1)]
    
    for i in 1:mat_n
        memo_eval = Dict()
        row_values = evaluate_at_point(dervs, points[i], minpoly_ord, support, var_to_sup, memo_eval)
        
        for (supp, val) in row_values
            ind = supp_to_index[supp]
            M[i, ind] = val
        end
    end
    
    return M
end

function make_matrix(n_vars, dervs, minpoly_ord, support, mat_n, mat_m, set_x1 = false)
    F = base_ring(parent(dervs[end]))
    points = generate_points(F, mat_n, n_vars, set_x1)
    return fill_matrix(points, dervs, minpoly_ord, support, mat_m)
end