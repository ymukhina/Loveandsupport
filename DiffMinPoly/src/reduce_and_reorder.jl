using Oscar
using Nemo

function find_new_order(M::fpMatrix, pivot::Int)   #O(l_supp)
    non_zero_cols, zero_cols = Int[], Int[]

    for i in 1:size(M, 2) 
        if !iszero(M[pivot, i])             
            push!(non_zero_cols, i)         
        else
            push!(zero_cols, i)
        end

    end
    n_ord = vcat(non_zero_cols, zero_cols)
    return n_ord
end

#In place reordering of the LinSys to save space:   Time: O(n^2)  Space: O(n)  with n the size of LinSys
function reorder_ls(M::fpMatrix, new_order::Vector{Int}, field::T) where T
    n, m = size(M, 1), length(new_order)
    S = matrix_space(field, n, m)                  
    ls = zero(S)
    for i in 1:m
        ls[:, i] = M[:, new_order[i]]
    end
    M = nothing
    return ls
end


function reorder_supp(supp::Vector{PointVector{T}}, new_order::Vector{Int}) where T

    @assert length(supp) == length(new_order)
    reordered_supp = Vector{PointVector{T}}(undef, length(new_order))
    
    for i in 1:length(new_order)
        reordered_supp[i] = supp[new_order[i]]
    end
    supp = nothing
    return reordered_supp
end