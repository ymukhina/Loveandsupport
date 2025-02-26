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

function get_terms(d::DualNumber{T}) where T
    return d.terms
end

function convert(x::fpMPolyRingElem, F::fpField)::fpFieldElem
    if x == 0
        return F(0)
    end

    exps = collect(exponent_vectors(x))

    if isempty(exps) || (length(exps) == 1 && exps[1] == zeros(length(exps[1])))
        a = constant_coefficient(x)
        return a
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
            if j <= length(a.terms) && (i-j+1) <= length(b.terms)
                result_terms[i] += a.terms[j] * b.terms[i-j+1]
            end
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



function (f::fpMPolyRingElem)(args::Vector{DiffMinPoly.DualNumber{fpFieldElem}})
    if !isa(args[1], DualNumber)
        throw(ArgumentError("First argument must be a DualNumber"))
    end
    
    result = zero(parent(f))
    
    exps = collect(exponent_vectors(f))
    coeffs = [coeff(f, i) for i in 1:length(exps)]
    
    for (coeff, exp) in zip(coeffs, exps)
        if exp[1] >= args[1].vanishing_degree  
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

