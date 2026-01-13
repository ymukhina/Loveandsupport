using Oscar 


function build_matrix(ode, n, dervs, minpoly_ord, support, n_rows; vanish_deg = false, info = true)

    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1: (minpoly_ord + 1) ]
    F = base_ring(parent(dervs[end]))
    # @info "parent we need", parent(dervs[end])[1]
    support = [Vector{Int64}(p) for p in support]
    lsup = length(support)       
    M = Array{Any}(undef, n_rows, lsup)

    supp_to_index = Dict(s => i for (i, s) in enumerate(support))

    if vanish_deg == false
        points = generate_points_base(F, n_rows, n)
    elseif (is_linear(dervs[1]) == -1)
         @info "Special case"
        points = generate_points_dual(F, n_rows, n, vanish_deg)
    else 
        @info dervs[1]
         @info "Proletarian case"
        points = generate_points_dual_linear(F, n_rows, n, vanish_deg, dervs[1])
    end   

 

    R, ε = polynomial_ring(F, "ε")

    # filling the columns corresponding to the derivatives
    for i in 1:n_rows
        if vanish_deg == 1 || vanish_deg == false
            M[i, 1] = F(1)
        else
            M[i, 1] = DualNumber{fpFieldElem}(R(1), vanish_deg)
        end

        evals = evaluate_polynomial(dervs, points[i], vanish_deg)
            
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
                            push!(v, val)
                        else
                            # @info "Ahoy! ", supp
                            push!(v, Epsilon(vanish_deg, F))
                        end
                    end
                    multiplier_eval = prod(v .^ multiplier)
                else
                    multiplier_eval = M[j, mult_ind]
                end
               
                M[j, i] = M[j, supp_div_ind] * multiplier_eval

            end
        end

        if M[1, end] isa DualNumber
            for i in 1:n_rows
                for j in 1:lsup
                    if !(M[i, j] isa fpFieldElem)
                        poly = M[i, j].poly
                        leading_term_exponent = degree(poly)
                        if (leading_term_exponent + 1) == vanish_deg
                            M[i, j] = AbstractAlgebra.leading_coefficient(poly)
                        elseif (leading_term_exponent + 1) < vanish_deg
                            M[i, j] = 0
                        else
                            error("Leading term exponent is greater than vanish_deg, which is impossible.")
                        end
                    else
                        # @info "Carramba", M[i, j]
                    end
                end
            end
        else
            # @info "1000 Chertey ", M[1, end], " deg ", vanish_deg
        end
        # @info "Matrix before" M
        S = matrix_space(F, n_rows, lsup)
        # @info "Matrix after" S(M)  
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


function solve_matrix(ode, n, dervs, ord, possible_supp, ks, l; info=true)
    solve_ker = 0
    build_mat = 0


    for i in 1:length(ks) + 1
        if i > length(ks)
            supp = possible_supp
        else
            supp = possible_supp[1:ks[i]]
        end

        #separate first, middle and the last blocks
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
            ls = build_matrix(ode, n, dervs, ord, supp, n_rows; info = true)
           # (n, dervs, minpoly_ord, support, n_rows, vanish_deg = false, info = true)
        else
            ls = build_matrix(ode, n, dervs, ord, supp, n_rows; vanish_deg = Int(i), info = info)
            # All other block rows
        end 

        t = time() - strt
        build_mat += t

        strt = time()
        if i == 1
            ker = kernel(ls, side=:right)
            # @info "II" size(ker)     #First row
        elseif i <= length(ks)
            ker = solve_linear_combinations(ls, ker, ks[i-1])
            # @info "AI" size(ker)  # All subsequent rectangular blocks
        else
            ker = solve_linear_combinations(ls, ker, ks[end])  
            # @info "OI" size(ker)            #Last row
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
        E = build_matrix(ode, n, dervs, ord, possible_supp, dim - 1; vanish_deg = false, info = info)
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

    return ker, dim, build_mat, solve_ker
end

#build_matrix(ode, n, dervs, minpoly_ord, support, n_rows; vanish_deg = false, info = true)
function build_matrix_general(F, ode, n, m, dervs, minpoly_ord, support; info = true)
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1:(minpoly_ord + m + 1) ]                                           

    support = [Vector{Int64}(p) for p in support]
    
    lsup = length(support)                                                    
    S = matrix_space(F, lsup, lsup)
    M = zero(S)
    supp_to_index = Dict(s => i for (i, s) in enumerate(support))

    # filling the columns corresponding to the derivatives
    for i in 1:lsup
        M[i, 1] = 1
        vec = [rand(F) for _ in 1:(n + m + 1)] 
        evals = [derv(vec...) for derv in dervs]
               
        for j in 1:(minpoly_ord + m + 1)
            supp = var_to_sup(j)
            ind = supp_to_index[supp]
            if j > m
                M[i, ind] = evals[j - m]
            else
                M[i, ind] = vec[findfirst(isequal(ode.parameters[j]), gens(parent(ode)))]
            end
        end
    end

    # filling the rest of the columns
    for i in (minpoly_ord + m + 3):lsup
        supp = support[i]
        supp_divisor = copy(supp)
        nonzero_ind = findfirst(e -> e > 0, supp_divisor)
        supp_divisor[nonzero_ind] -= 1                                                 
        multiplier = zeros(Int, minpoly_ord + m + 1)
        multiplier[nonzero_ind] += 1
        while !haskey(supp_to_index, supp_divisor)
            nonzero_ind = findfirst(e -> e > 0, supp_divisor)
            supp_divisor[nonzero_ind] -= 1
            multiplier[nonzero_ind] += 1
        end                                                    
        
        supp_div_ind = supp_to_index[supp_divisor]
        mult_ind = get(supp_to_index, multiplier, -1)
        for j in 1:lsup
            if mult_ind == -1
                multiplier_eval = prod(M[j, 2:(minpoly_ord + m + 2)] .^ multiplier)
            else
                multiplier_eval = M[j, mult_ind]
            end
            M[j, i] = M[j, supp_div_ind] * multiplier_eval           
        end
    end
    return M
end


function solve_matrix_general(F, ode, n, m, dervs, ord::Int, possible_supp; info=true)

    info && @info "The size of the estimates support is $(length(possible_supp))"

    tim2 = @elapsed ls = build_matrix_general(F, ode, n, m, dervs, ord, possible_supp; info = true)
                                                                    
    info && @info "eval method $(tim2)"

    info && @info "linear system dims $(size(ls))"
    
    system_soltime = @elapsed ker = kernel(ls, side=:right)
    info && @info "Linear system solved in $system_soltime"

    dim = size(ker)[2]
    info && @info "The dimension of the solution space is $(dim)"

    return ker, dim
end


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

