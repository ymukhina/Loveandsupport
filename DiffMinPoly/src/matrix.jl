#----------------------------------------------------------------------------#
#----------------- Rational Parametrisation Case ----------------------------#
#----------------------------------------------------------------------------#

function build_smart_matrix_truncated_dual(F, odeios::ODEios, split_array::Vector{Int}, n_col::Int, rational_param; info = true)
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1:(minpoly_ord + m + 1) ]                                           

    n = length(odeios.ode.x_vars)
    m = length(odeios.ode.parameters)
    minpoly_ord = odeios.order
    dervs = odeios.dervs
   
    support = odeios.support
    support = [Vector{Int64}(p) for p in support]
    reduced_support = support[1:n_col]

    n_rows = sum(split_array)

    max_k_degree = length(split_array)


    R_eps, ε = power_series_ring(F, max_k_degree + 1, "ε")
    M = Matrix{Any}(undef, n_rows, n_col)
    zero_series = zero(R_eps)

    for i in 1:n_rows, j in 1:n_col
        M[i, j] = zero_series
    end

    supp_to_index = Dict(s => i for (i, s) in enumerate(reduced_support))


    cum_split = cumsum(split_array)

    filled_columns = falses(n_col)


    # filling the columns corresponding to the derivatives
    for i in 1:sum(split_array)
        M[i, 1] = R_eps(1)
        filled_columns[1] = true
        block = findfirst(x -> i <= x, cum_split) 

        vec_ps = generate_truncated_dual_points_new(F, max_k_degree - block, 1, n + m, rational_param)[1] 
        # @info "vector to evaluate" vec_ps
        evals = [derv(vec_ps[1:n+1]...) for derv in dervs]
        # @info "evals" evals
               
        for j in 1:(minpoly_ord + m + 1)
            supp = var_to_sup(j)
            if haskey(supp_to_index, supp)  
                ind = supp_to_index[supp]
                if j > m
                    M[i, ind] = evals[j - m]
                else
                    M[i, ind] = vec[findfirst(isequal(ode.parameters[j]), gens(parent(ode)))]
                end
            filled_columns[ind] = true
            end        
        end
    end


    # filling the rest of the columns
    remaining_indices = findall(.!filled_columns)

    for i in remaining_indices 
        supp = reduced_support[i]
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
    # @info "our matrix" M
    return M
end

function build_smart_matrix_constant_part(F, odeios::ODEios, n_rows::Int; info = true)
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1:(minpoly_ord + m + 1) ]                                           

    n = length(odeios.ode.x_vars)
    m = length(odeios.ode.parameters)
    minpoly_ord = odeios.order
    dervs = odeios.dervs
   
    support = odeios.support
    support = [Vector{Int64}(p) for p in support]
    n_col = length(support) 
   
    M = zeros(F, n_rows, n_col)

    supp_to_index = Dict(s => i for (i, s) in enumerate(support))

    filled_columns = falses(n_col)


    # filling the columns corresponding to the derivatives
    for i in 1:n_rows
        M[i, 1] = F(1)
        filled_columns[1] = true
        
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
            filled_columns[ind] = true
        end
    end


    # filling the rest of the columns
    remaining_indices = findall(.!filled_columns)

    for i in remaining_indices 
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

    return M
end


function generate_submatrix_subsequence_smart(F, odeios::ODEios, rational_param; info = true)

    support = odeios.support
    hd = support[end][1]
    #splits = split_index(support, hd)
    splits = smarter_split_index(support, hd)
    # @info splits
    @info "SO so smart, masters of knowledge" 
    ############################
    # fill the matrix with trancated epsilon ps
    @info "very new way with power series"
    M = build_smart_matrix_truncated_dual(F, odeios::ODEios, splits[3], splits[2][end - 1], rational_param; info = true)
    ############################

    N = Vector{fpMatrix}(undef, length(splits[1]))

    for i in 1:(length(splits[1]) - 1)
        N[i] = submatrix_dual_matrix(M, i-1, splits[1][i], splits[2][i])
    end
    
    end_ind = length(splits[1])
    constant_matrix = build_smart_matrix_constant_part(F, odeios, splits[1][end_ind]; info = true)
    N[end_ind] = matrix(F, constant_matrix)
    # @info "Size last" size(N[end_ind])
    return N

end

#-------------------Good but we can do better as it seems-----------------------#
function build_smart_matrix_truncated(F, odeios::ODEios, split_array::Vector{Int}, rational_param; info = true)
    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1:(minpoly_ord + m + 1) ]                                           

    n = length(odeios.ode.x_vars)
    m = length(odeios.ode.parameters)
    minpoly_ord = odeios.order
    dervs = odeios.dervs
   
    support = odeios.support
    support = [Vector{Int64}(p) for p in support]

    lsup = length(support)  
    n_rows = sum(split_array) + 1

    max_k_degree = length(split_array)


    R_eps, ε = power_series_ring(F, max_k_degree + 1, "ε")
    M = Matrix{Any}(undef, n_rows, lsup)
    zero_series = zero(R_eps)

    for i in 1:n_rows, j in 1:lsup
        M[i, j] = zero_series
    end

    supp_to_index = Dict(s => i for (i, s) in enumerate(support))

    cum_split = cumsum(split_array)
    filled_columns = falses(lsup)


    # filling the columns corresponding to the derivatives
    for i in 1:sum(split_array)+1
        M[i, 1] = R_eps(1)
        filled_columns[1] = true
        block = findfirst(x -> i <= x, cum_split) 

        if i == n_rows
            vec = [rand(F) for _ in 1:(n + m + 1)]
            vec_ps = [R_eps(v) for v in vec]
        else
            vec_ps = generate_truncated_dual_points_new(F, max_k_degree - block, 1, n + m, rational_param)[1] 
        end
     
        evals = [derv(vec_ps[1:n+1]...) for derv in dervs]
               
        for j in 1:(minpoly_ord + m + 1)
            supp = var_to_sup(j)
            ind = supp_to_index[supp]
            if j > m
                M[i, ind] = evals[j - m]
            else
                M[i, ind] = vec[findfirst(isequal(ode.parameters[j]), gens(parent(ode)))]
            end
            filled_columns[ind] = true
        end
    end


    # filling the rest of the columns
    remaining_indices = findall(.!filled_columns)

    for i in remaining_indices 
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

    return M
end


"""
    submatrix_dual_matrix(M, index, size_rows, size_col)

Extract a submatrix of specified dimensions (size_rows x size_col) from a matrix of DualNumbers M, 
selecting the coefficient of the graded part at the given index.

"""
function submatrix_dual_matrix(M, index, size_rows, size_col)
    
    F = base_ring(M[1,1])

    N = Matrix{fpFieldElem}(undef, size_rows, size_col)
    fill!(N, F(0))

    if index == -1
        last_row = M[end, 1:size_col]
        for i in 1:size_rows
            for j in 1:size_col
                N[i, j] = coeff(last_row[j], 0) 
            end
        end
    else
        if index == 0
            for i in 1:size_rows
                N[i, 1] = one(F)
            end
        end
        for i in 1:size_rows, j in 2:size_col
            N[i, j] = coeff(M[i, j], index)  
        end
    end

    return matrix(F, N)
end


function generate_submatrix_subsequence(F, odeios::ODEios, rational_param; info = true)

    support = odeios.support
    hd = support[end][1]
    splits = split_index(support, hd)
    # splits = smart_split_index(support, hd)
    ############################
    # fill the matrix with trancated epsilon ps
    @info "very new way with power series"
    M = build_smart_matrix_truncated(F, odeios::ODEios, splits[3], rational_param; info = true)
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

#---------------------------------------------------------------------------------------#
#--------------------------------- General Case ----------------------------------------#
#---------------------------------------------------------------------------------------#

function build_matrix(F, odeios::ODEios; info = true)

    ode = odeios.ode
    support = odeios.support
    minpoly_ord = odeios.order
    dervs = odeios.dervs
    
    n = length(ode.x_vars)
    m = length(ode.parameters)

    var_to_sup = var_ind -> [(k == var_ind) ? 1 : 0 for k in 1:(minpoly_ord + m + 1) ]                                           

    support = [Vector{Int64}(p) for p in support]

    # @info "Support" support
    
    lsup = length(support)                                                    
    S = matrix_space(F, lsup, lsup)
    M = zero(S)
    supp_to_index = Dict(s => i for (i, s) in enumerate(support))

    filled_columns = falses(lsup)

    # filling the columns corresponding to the derivatives
    for i in 1:lsup
        M[i, 1] = 1
        filled_columns[1] = true

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
            filled_columns[ind] = true
        end
    end
    
    # filling the rest of the columns
    remaining_indices = findall(.!filled_columns)
    for i in remaining_indices
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

function solve_matrix(F, odeios::ODEios; info=true)

    # info && @info "The size of the estimates support is $(length(odeios.support))"

    tim2 = @elapsed ls = build_matrix(F, odeios; info = true)
                                                                    
    # info && @info "eval method $(tim2)"

    # info && @info "linear system dims $(size(ls))"
    
    system_soltime = @elapsed ker = kernel(ls, side=:right)
    # info && @info "Linear system solved in $system_soltime"

    dim = size(ker)[2]
    # info && @info "The dimension of the solution space is $(dim)"

    return ker, dim, tim2, system_soltime
end

