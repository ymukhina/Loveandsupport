# + Gleb: should not be here

############################################################################
# updated code for rational parametrisation
############################################################################

const MatrixElement = Union{fpFieldElem, DualNumber{fpFieldElem}}

# + Gleb: too many arguments, to discuss wrapping them into a structure
# we say adios to many arguments
function build_smart_matrix_truncated(F, odeios::ODEios, split_array::Vector{Int}, dervs, rational_param; info = true)
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1:(minpoly_ord + m + 1) ]                                           

    n = length(odeios.ode.x_vars)
    m = length(odeios.ode.parameters)
    minpoly_ord = odeios.order

    support = sort_gleb!(odeios.support)
    old_support = support
    new_support = copy(support)
    sort_gleb_max!(new_support)

    support = [Vector{Int64}(p) for p in support]

    lsup = length(support)  
    n_rows = sum(split_array) + 1

    max_k_degree = length(split_array)

    T = MatrixElement
    M = Matrix{T}(undef, n_rows, lsup)
    zero_val = F(0) 
    for i in 1:n_rows, j in 1:lsup
        M[i, j] = zero_val
    end

    supp_to_index = Dict(s => i for (i, s) in enumerate(support))

    cum_split = cumsum(split_array)

    # filling the columns corresponding to the derivatives
    for i in 1:sum(split_array)+1
        M[i, 1] = F(1)
        block = findfirst(x -> i <= x, cum_split) 

        if i == n_rows
            vec = [rand(F) for _ in 1:(n + m + 1)]
        else
            vec = generate_truncated_dual_points(F, max_k_degree - block, 1, n + m, rational_param) 
        end
     
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
        for j in 1:n_rows
            if mult_ind == -1
                multiplier_eval = prod(M[j, 2:(minpoly_ord + m + 2)] .^ multiplier)
            else
                multiplier_eval = M[j, mult_ind]
            end
            M[j, i] = M[j, supp_div_ind] * multiplier_eval           
        end
    end

    old_index = Dict(s => i for (i, s) in enumerate(old_support))
    perm = [old_index[s] for s in new_support]

    M = M[:, perm]

    return M
end

# + Gleb: what does it do? Add a docstring
"""
    submatrix_dual_matrix(M, index, size_rows, size_col)

Extract a submatrix of specified dimensions (size_rows x size_col) from a matrix of DualNumbers M, 
selecting the coefficient of the graded part at the given index.

"""

function submatrix_dual_matrix(M, index, size_rows, size_col)
    F = base_ring(M[1,2].poly)

    N = Matrix{fpFieldElem}(undef, size_rows, size_col)
    fill!(N, F(0))

    if index == -1
        last_row = M[end, 1:size_col]
        for i in 1:size_rows
            for j in 1:size_col
                N[i, j] = F(last_row[j]) 
            end
        end
    else
        for i in 1:size_rows, j in 2:size_col
            if index == 0
                N[i,1] = F(M[1,1])
            else
                N[i,1] = F(0)
            end
            N[i, j] = F(coeff(M[i, j].poly, index))  
        end
    end

    return matrix(F, N)
end


function generate_submatrix_subsequence(F, odeios::ODEios, dervs, rational_param; info = true)

    support = sort_gleb_max!(odeios.support)
    hd = support[end][1]
    splits = split_index(support, hd)
    
    ############################
    # fill the matrix with trancated epsilon ps
    #@info "old way"
    # n_rows = splits[1][1]
    # k = length(splits[2]) - 2
    # M = build_matrix_truncated(F, ode, n, m, n_rows + 1, k, dervs, minpoly_ord, support, rational_param; info = true)
    # ###########################

    ############################
    # fill the matrix with trancated epsilon ps
    @info "new way"
    M = build_smart_matrix_truncated(F, odeios::ODEios, splits[3], dervs, rational_param; info = true)
    ############################

    N = Vector{fpMatrix}(undef, length(splits[1]))

    for i in 1:(length(splits[1]) - 1)
        N[i] = submatrix_dual_matrix(M, i-1, splits[1][i], splits[2][i])
    end
    
    end_ind = length(splits[1])
    N[end_ind] = submatrix_dual_matrix(M, -1, splits[1][end_ind], splits[2][end_ind])
    

    return N

end

function mini_ker(N, sp)

    F = base_ring(N)
    n_old = nrows(sp)             
    d = ncols(sp)                 
    n_rows, n_tot = size(N)

    A = sub(N, 1:n_rows, 1:n_old)
    B = sub(N, 1:n_rows, n_old+1:n_tot)

    aug_cols = ncols(B) + d
    S = matrix_space(F, n_rows, aug_cols)
    aug = zero(S)

    aug[:, 1:ncols(B)] = B
    aug[:, ncols(B)+1:end] = A * sp

    v = kernel(aug, side = :right)

    if ncols(v) == 0
        return zero_matrix(F, n_old + ncols(B), 0)
    end

    y_part = sub(v, 1:ncols(B), 1:ncols(v))
    lambdas = sub(v, ncols(B)+1:nrows(v), 1:ncols(v))

    x_part = sp * lambdas

    return vcat(x_part, y_part)
end

function constrained_kernel(Ms...)
    sp = kernel(Ms[1], side = :right)

    for i in 2:length(Ms)  
       
        sp = mini_ker(Ms[i], sp)
    end

    return sp
end


#------- WE DO NO USE THIS FUNCTION NOW, To fill all the matrix with truncated epsilon---------------------#
function build_matrix_truncated(F, ode, n, m, n_rows, k, dervs, minpoly_ord, support, rational_param; info = true)
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1:(minpoly_ord + m + 1) ]                                           

    support = sort_gleb!(support)

    old_support = support
    new_support = copy(support)
    sort_gleb_max!(new_support)

    support = [Vector{Int64}(p) for p in support]
    
    lsup = length(support)  

    M = Array{Any}(undef, n_rows, lsup)
    for i in 1:n_rows, j in 1:lsup
        M[i, j] = 0
    end
    
    supp_to_index = Dict(s => i for (i, s) in enumerate(support))

    # filling the columns corresponding to the derivatives
    for i in 1:n_rows
        M[i, 1] = 1

        if i == n_rows
            vec = [rand(F) for _ in 1:(n + m + 1)]
        else
            vec = generate_truncated_dual_points(F, k, 1, n + m, rational_param) 
        end
     
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
        for j in 1:n_rows
            if mult_ind == -1
                multiplier_eval = prod(M[j, 2:(minpoly_ord + m + 2)] .^ multiplier)
            else
                multiplier_eval = M[j, mult_ind]
            end
            M[j, i] = M[j, supp_div_ind] * multiplier_eval           
        end
    end

    old_index = Dict(s => i for (i, s) in enumerate(old_support))
    perm = [old_index[s] for s in new_support]

    M = M[:, perm]

    return M
end

############################################################################
# code for the general case
############################################################################

#build_matrix(ode, n, dervs, minpoly_ord, support, n_rows; vanish_deg = false, info = true)
function build_matrix_general(F, ode, n, m, dervs, minpoly_ord, support; info = true)
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1:(minpoly_ord + m + 1) ]                                           

    support = [Vector{Int64}(p) for p in support]

    # @info "Support" support
    
    lsup = length(support)                                                    
    S = matrix_space(F, lsup, lsup)
    M = zero(S)
    supp_to_index = Dict(s => i for (i, s) in enumerate(support))

    # filling the columns corresponding to the derivatives
    for i in 1:lsup
        M[i, 1] = 1
        vec = [rand(F) for _ in 1:(n + m + 1)] 
        # @info "VEC" vec
        evals = [derv(vec...) for derv in dervs]
        # @info "EVALS" evals
               
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


############################################################################
# code Max
############################################################################

function build_matrix(F, ode, n, dervs, minpoly_ord, support, n_rows; vanish_deg = false, rational_param = nothing, info = true)

    # @info "Support" support

    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1: (minpoly_ord + 1) ]
    F = base_ring(parent(dervs[end]))
    # @info "parent we need", parent(dervs[end])[1]
    support = [Vector{Int64}(p) for p in support]
    lsup = length(support)       
    M = Array{Any}(undef, n_rows, lsup)

    supp_to_index = Dict(s => i for (i, s) in enumerate(support))

    if vanish_deg == false
        points = generate_points_base(F, n_rows, n)

    else (!isempty(rational_param))
        # @info "vanishing degree" vanish_deg
        points = generate_points_rational_parametrization(F, n_rows, n, vanish_deg, rational_param)
        # @info "points" points
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


function solve_matrix(F, ode, n, dervs, ord, possible_supp, ks, l, rational_param; info=true)
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
            ls = build_matrix(F, ode, n, dervs, ord, supp, n_rows, rational_param = nothing; info = true)
            # @info "Matrix1" ls
           # (n, dervs, minpoly_ord, support, n_rows, vanish_deg = false, info = true)
        else
            ls = build_matrix(F, ode, n, dervs, ord, supp, n_rows; vanish_deg = Int(i), rational_param, info = info)
            # @info "Matrix2" ls
            # All other block rows
        end 

        t = time() - strt
        build_mat += t

        strt = time()
        if i == 1
            ker = kernel(ls, side=:right)
            #  @info "II" ker    #First row
        elseif i <= length(ks)
            ker = solve_linear_combinations(ls, ker, ks[i-1])
            # @info "AI" ker  # All subsequent rectangular blocks
        else
            ker = solve_linear_combinations(ls, ker, ks[end])  
            # @info "OI" ker            #Last row
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
        E = build_matrix(F, ode, n, dervs, ord, possible_supp, dim - 1; vanish_deg = false, rational_param=nothing, info = info)
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

