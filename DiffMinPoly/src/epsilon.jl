using Oscar
using Nemo

struct DualNumber{T}
    terms::Vector{T}
    vanishing_degree::Integer
    field::Any

    function DualNumber(terms::Vector{T}, vanishing_degree::Integer, field=nothing) where T
        if length(terms) != vanishing_degree
            zero_elem = isnothing(field) ? zero(T) : field(0)
            new_terms = Vector{T}(undef, vanishing_degree)
            for i in 1:vanishing_degree
                new_terms[i] = i <= length(terms) ? terms[i] : zero_elem
            end
            terms = new_terms
        end
        new{T}(terms, vanishing_degree, field)
    end
end

function Epsilon(vanishing_degree::Integer, field=nothing)
    if vanishing_degree < 1
        throw(ArgumentError("Vanishing degree must be positive"))
    end
    zero_elem = isnothing(field) ? 0 : field(0)
    one_elem = isnothing(field) ? 1 : field(1)
    terms = [zero_elem for _ in 1:vanishing_degree]
    if vanishing_degree != 1
        terms[2] = one_elem
    end
    return DualNumber(terms, vanishing_degree, field)
end

function convert(x::fpMPolyRingElem, F::fpField)
    if x == 0
        return F(0)
    end

    exps = collect(exponent_vectors(x))
    if isempty(exps) || (length(exps) == 1 && exps[1] == zeros(length(exps[1])))
        return coeff(x, 0)
    else
        throw(ArgumentError("Cannot convert non-constant polynomial to field element"))
    end
end

Base.:*(a::DualNumber{T}, b::DualNumber{T}) where T = begin
    # println("(",a,")*(", b,")")
    degree = min(a.vanishing_degree, b.vanishing_degree)
    zero_elem = isnothing(a.field) ? zero(T) : a.field(0)
    result_terms = [zero_elem for _ in 1:degree]
    
    for i in 1:degree
        for j in 1:i
            k = i - j + 1
            result_terms[i] += a.terms[j] * b.terms[k]
        end
    end
    
    return DualNumber(result_terms, degree, a.field)
end


Base.:*(a::DualNumber{T}, b::fpMPolyRingElem) where T = begin
    # println(a,"*", b)
    if !isnothing(a.field)
        b_converted = convert(b, a.field)
        new_terms = [term * b_converted for term in a.terms]
        return DualNumber(new_terms, a.vanishing_degree, a.field)
    end
    throw(ArgumentError("Cannot multiply DualNumber without field with polynomial"))
end

Base.:*(a::DualNumber{T}, b::Any) where T = begin
    # println(a,"*", b)
    new_terms = [term * b for term in a.terms]
    return DualNumber(new_terms, a.vanishing_degree, a.field)
end

Base.:*(b::fpMPolyRingElem, a::DualNumber{T}) where T = begin
    a * b
end

Base.:*(b::Any, a::DualNumber{T}) where T = begin
    a * b
end


Base.:+(a::DualNumber{T}, b::DualNumber{T}) where T = begin
    # println(a,"+", b)
    degree = min(a.vanishing_degree, b.vanishing_degree)
    result_terms = [a.terms[i] + b.terms[i] for i in 1:degree]
    return DualNumber(result_terms, degree, a.field)
end

Base.:+(a::DualNumber{T}, b::fpMPolyRingElem) where T = begin
    # println(a,"+", b)
    if !isnothing(a.field)
        b_converted = convert(b, a.field)
        new_terms = copy(a.terms)
        new_terms[1] += b_converted
        return DualNumber(new_terms, a.vanishing_degree, a.field)
    end
    throw(ArgumentError("Cannot add DualNumber without field with polynomial"))
end

Base.:+(a::DualNumber{T}, b::Any) where T = begin
    # println(a, "+", b)
    a.terms[1] = a.terms[1] + b
    return DualNumber(a.terms, a.vanishing_degree, a.field)
end

Base.:+(b::fpMPolyRingElem, a::DualNumber{T}) where T = begin
    a + b
end

Base.:+(b::Any, a::DualNumber{T}) where T = begin
    a + b
end


Base.:-(a::DualNumber{T}, b::DualNumber{T}) where T = begin
    # println(a, "-", b)
    degree = min(a.vanishing_degree, b.vanishing_degree)
    result_terms = [a.terms[i] - b.terms[i] for i in 1:degree]
    return DualNumber(result_terms, degree, a.field)
end

Base.:-(b::fpMPolyRingElem, a::DualNumber{T}) where T = begin
    - a + b
end


Base.:-(a::DualNumber{T}, b::fpMPolyRingElem) where T = begin
    # println(a,"-", b)
    if !isnothing(a.field)
        b_converted = convert(b, a.field)
        a.terms[1] = a.terms[1] - b_converted
        return DualNumber(a.terms, a.vanishing_degree, a.field)
    end
    throw(ArgumentError("Cannot add DualNumber without field with polynomial"))
end


Base.:-(a::DualNumber{T}, b::Any) where T = begin
    # println(a, "-", b)
    a.terms[1] = a.terms[1] - b
    DualNumber(a.terms, a.vanishing_degree, a.field)
end

Base.:-(b::Any, a::DualNumber{T}) where T = begin
    # println(b,"-", a)
    a.terms[1] = b - a.terms[1]
    for i in 2:a.vanishing_degree
        a.terms[i] = - a.terms[i]
    end
    DualNumber(a.terms, a.vanishing_degree, a.field)
end

Base.:-(a::DualNumber{T}) where T = begin
    # println("-", a)
    for i in 1:a.vanishing_degree
        a.terms[i] = -a.terms[i]
    end
    DualNumber(a.terms, a.vanishing_degree, a.field)
end

Base.:^(a::DualNumber{T}, n::Integer) where T = begin
    # println(a, "^", n)
    if n == 0
        zero_elem = isnothing(a.field) ? zero(T) : a.field(0)
        one_elem = isnothing(a.field) ? one(T) : a.field(1)
        return DualNumber([one_elem; [zero_elem for _ in 2:a.vanishing_degree]], a.vanishing_degree, a.field)
    elseif n == 1
        return a
    end

    result = a
    for _ in 2:n
        result = result * a
    end
    return result
end

Base.show(io::IO, d::DualNumber) = begin
    terms_str = String[]
    for (i, term) in enumerate(d.terms)
        if term != 0
            if i == 1
                push!(terms_str, string(term))
            elseif term == 1
                push!(terms_str, "ε" * (i > 2 ? "^$(i-1)" : ""))
            else
                push!(terms_str, string(term, "ε" * (i > 2 ? "^$(i-1)" : "")))
            end
        end
    end
    
    if isempty(terms_str)
        print(io, "0")
    else
        print(io, join(terms_str, " + "))
    end
end

# @inline function set_in_mat!(M::fpMatrix, a::DiffMinPoly.DualNumber{fpFieldElem}, i::Int, j::Int, n::Int)
#     if i <= n && j <=n
#         (base_ring(M) != parent(a)) && error("Parent objects must coincide")
#         set_indexraw!(M, a, i, j)
#     else
#         error("Index ($i, $j) out of range for Matrix of size ($n x $n)")
#     end
# end

# @inline fucn


function (f::fpMPolyRingElem)(args::Vector{DiffMinPoly.DualNumber{fpFieldElem}})
    if !isa(args[1], DualNumber)
        throw(ArgumentError("First argument must be a DualNumber"))
    end
    
    result = zero(parent(f))
    
    exps = collect(exponent_vectors(f))
    coeffs = [coeff(f, i) for i in 1:length(exps)]
    
    for (coeff, exp) in zip(coeffs, exps)
        if exp[1] >= args[1].vanishing_degree   #Need Dual in first var
            continue
        end
        
        term = coeff
        for (i, p) in enumerate(exp)
            if i == 1
                if p > 0
                    term *= args[1]^p
                end
            elseif i <= length(args)
                term *= args[i]^p
            end
        end
        result += term
    end
    
    return result
end

# p = 3300697969
# F = Nemo.Native.GF(p)  
# R, (x1, x2, x3) = polynomial_ring(F, ["x1", "x2", "x3"])

# f1 = x1

# f2 = 2*x1 + 10*x2 + 6*x3 + 8

# f3 = 10*x1^2 + 50*x1*x2 + 80*x1*x3 + 94*x1 + 80*x2^2 + 60*x2*x3 + 58*x2 + 30*x3^2 + 74*x3 + 128

# f4 = 50*x1^3 + 410*x1^2*x2 + 460*x1^2*x3 + 798*x1^2 + 1200*x1*x2^2 + 1880*x1*x2*x3 + 2190*x1*x2 + 630*x1*x3^2 + 1814*x1*x3 + 1726*x1 + 1280*x2^3 + 1440*x2^2*x3 + 1464*x2^2 + 840*x2*x3^2 + 2668*x2*x3 + 3398*x2 + 180*x3^3 + 1074*x3^2 + 2362*x3 + 1480

# a, b = rand(F), rand(F)
# point = [ε, a, b]
# println(point)
# r1 = f1(point, F)
# r2 = f2(point, F)
# r3 = f3(point, F)
# r4 = f4(point, F)
# println(r1)
# println("##")
# println(r2)
# println(10*a + 6*b + 8)
# println(r3)
# println(50*a + 80*b + 94)
# println(80*a^2 + 60*a*b + 58*a + 30*b^2 + 74*b + 128)
# println(r4)
# println(1200a^2 + 1800*a*b + 2190*a + 630*b^2 + 1814*b + 1726)
# println(1280*a^3 + 1440*a^2*b + 1464*a^2 + 840*a*b^2 + 2668*a*b + 3398*a + 180*b^3 + 1074*b^2 + 2362*a + 1480)

