function kernel_blocks(ls_UL, ls_LL, ls_LR) 
    strt = time()

    v1 = kernel(ls_UL, side=:right)
    ker = solve_linear_combinations_blocks(ls_LL, ls_LR, v1)

    @info "Linear System solved in $(time() - strt)"

    return ker
end


function solve_linear_combinations(ls, sol_space)
    F = base_ring(ls)

    # Gleb: better variable naming
    m_temp = size(sol_space, 2)
    ls_n = size(ls, 1)
    S = matrix_space(F, ls_n, m_temp)
    aug = zero(S)     
    
    for j in 1:m_temp
        w_i = ls * sol_space[:, j]
        aug[:, j] = w_i
    end

    lambs = kernel(aug, side=:right)

    S = matrix_space(F, size(sol_space, 1), 1)
    v1_lamb = zero(S)
    ker = Matrix{eltype(sol_space)}(undef, size(sol_space, 1), size(lambs, 2))

    for i in 1:size(lambs, 2)
        v1_lamb = sum(lambs[j, i] .* sol_space[:, j] for j in 1:size(sol_space, 2))
        ker[:, i] = v1_lamb
    end

    ker = matrix_space(F, size(ker,1), size(ker,2))(ker)
    return ker
end


function solve_linear_combinations_blocks(ls_L, ls_R, sol_space)
    F = base_ring(ls_L)

    # Gleb: naming
    n_temp, m_temp = size(sol_space)
    ls_n, ls_m = size(ls_R)
    S = matrix_space(F, ls_n, ls_m + m_temp)
    aug = zero(S)     
        
    for i in 1:ls_n
        for j in 1:ls_m
            aug[i, j] = ls_R[i, j]
        end
    end

    for j in 1:m_temp
        w_i = ls_L * sol_space[:, j]
        aug[:, ls_m + j] = w_i
    end


    v = kernel(aug, side=:right)
    v2, lambs = v[1:ls_m,:], v[ls_m+1:ls_m + m_temp, :]

    S = matrix_space(F, size(sol_space, 1), 1)
    v1_lamb = zero(S)
    ker = Matrix{eltype(sol_space)}(undef, size(sol_space, 1) + size(v2, 1), size(lambs, 2))

    for i in 1:size(lambs, 2)
        v1_lamb = sum(lambs[j, i] .* sol_space[:, j] for j in 1:size(sol_space, 2))
        ker[:, i] = vcat(v1_lamb, v2[:, i])
    end

    ker = matrix_space(F, size(ker,1), size(ker,2))(ker)
    return ker
end

 #Assuming ordered by sort_gleb_max! 
 function split_supp(supp, n_splits)  #   For now support only cuts by 1 e.g. the 3 blocks A, B, C
    k1 = 0
    # Gleb: I think you can just use `findfirst` built-in function here
    for i in 1:length(supp)
        if supp[i][1] != 0
            k1 = i - 1      
            break
        end
    end
    return supp[1:k1], supp[k1+1:length(supp)], k1
end


#     WIP
# function make_interpolation_points(F, n_points, n_vars, set_x1 = false)   #Initialized n_points at random in F
#     vecs = Vector{Vector{fpFieldElem}}(undef, n_points)
#     for i in 1:n_points
#         vec = [rand(F) for _ in 1:(n_vars + 1)]
#         if set_x1
#             vec[1] = F(0)
#         end
#         vecs[i] = vec
#     end
#     return vecs
# end

function make_matrix(n_vars, dervs, minpoly_ord, support, mat_n, mat_m, set_x1 = false)     #Fix this BS pls
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1:(minpoly_ord + 1)]

    support = [Vector{Int64}(p) for p in support]

    F = base_ring(parent(dervs[end]))
    S = matrix_space(F, mat_n, mat_m)
    M = zero(S)
    supp_to_index = Dict(s => i for (i, s) in enumerate(support))
    memo_eval = Dict()

    function decompose_multiplier(mult)
        if haskey(memo_eval, mult)
            return memo_eval[mult]
        end
        
        mult2 = zeros(Int, minpoly_ord + 1)
        nonzero_ind = findfirst(x -> x > 0, mult)
        
        if isnothing(nonzero_ind)
            return F(1) 
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

    for i in 1:mat_n
        M[i, 1] = 1
        memo_eval[support[1]] = 1
        vec = [rand(F) for _ in 1:(n_vars + 1)]
        if set_x1
            vec[1] = F(0)
        end
       
        evals = [derv(vec...) for derv in dervs]
        for j in 1:(minpoly_ord + 1)
            supp = var_to_sup(j)
            memo_eval[supp] = evals[j]
            if haskey(supp_to_index, supp)
                ind = supp_to_index[supp]
                M[i, ind] = evals[j]
            end
        end

        for supp in support[2:end]
            if haskey(memo_eval, supp)
                continue
            end
            ind = supp_to_index[supp]
           
            res = 0
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
            M[i, ind] = res
        end
        empty!(memo_eval)
    end
    return M
end
