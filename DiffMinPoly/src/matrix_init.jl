function kernel_block(ls_LL, ls_LR, v1) 
    strt = time()

    F = base_ring(ls_LL)

    n_temp, m_temp = size(v1)
    ls_n, ls_m = size(ls_LR)
    S = matrix_space(F, ls_n, ls_m + m_temp)
    aug = zero(S)     
        
    for i in 1:ls_n
        for j in 1:ls_m
            aug[i, j] = ls_LR[i, j]
        end
    end

    for j in 1:m_temp
        w_i = ls_LL * v1[:, j]
        aug[:, ls_m + j] = w_i
    end


    v = kernel(aug, side=:right)
    v2, lambs = v[1:ls_m,:], v[ls_m+1:ls_m + m_temp, :]

    S = matrix_space(F, size(v1, 1), 1)
    v1_lamb = zero(S)
    ker = Matrix{eltype(v1)}(undef, size(v1, 1) + size(v2, 1), size(lambs, 2))

    for i in 1:size(lambs, 2)
        v1_lamb = sum(lambs[j, i] .* v1[:, j] for j in 1:size(v1, 2))
        ker[:, i] = vcat(v1_lamb, v2[:, i])
    end

    ker = matrix_space(F, size(ker,1), size(ker,2))(ker)

    @info "Linear System solved in $(time() - strt)"

    return ker

end

 #Assuming ordered by sort_gleb_max! 
 function split_supp(supp, n_splits)  #   For now support only cuts by 1 e.g. the 3 blocks A, B, C
    k1 = 0                     
    for i in 1:length(supp)
        if supp[i][1] != 0
            k1 = i - 1      
            break
        end
    end
    return supp[1:k1], supp[k1+1:length(supp)], k1
end


function fill_matrix(F, n_vars, dervs, minpoly_ord, support, mat_n, mat_m, set_x1 = false)
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1:(minpoly_ord + 1)]
   
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
