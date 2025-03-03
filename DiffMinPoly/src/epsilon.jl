using Oscar

struct DualNumber{T}
    poly::fpMPolyRingElem
    vanishing_degree::Integer
    # Gleb: I think that the ground field can be deduced from `base_ring` of the polynomial
    field::Any

    function DualNumber{T}(poly::fpMPolyRingElem, vanishing_degree::Integer, field=nothing) where T
        if isnothing(field)
            throw(ArgumentError("Field must be provided for polynomial-based DualNumber"))
        end
        new{T}(poly, vanishing_degree, field)
    end

    function DualNumber(terms::Vector{T}, vanishing_degree::Integer, field=nothing) where T
        if isnothing(field)
            throw(ArgumentError("Field must be provided for polynomial-based DualNumber"))
        end
        
        R, (ε,) = polynomial_ring(field, ["ε"])
        
        poly = R(0)
        for (i, coef) in enumerate(terms)
            if i <= vanishing_degree && coef != field(0)
                poly += coef * ε^(i-1)
            end
        end
        
        new{T}(poly, vanishing_degree, field)
    end
end


# ------------------------------ Dual Number Utils ------------------------------

function Epsilon(vanishing_degree::Integer, field=nothing)
    if vanishing_degree < 1
        throw(ArgumentError("Vanishing degree must be positive"))
    end
    if isnothing(field)
        throw(ArgumentError("Field must be provided for polynomial-based DualNumber"))
    end
    
    R, (ε,) = polynomial_ring(field, ["ε"])
    
    poly = ε
    
    return DualNumber{fpFieldElem}(poly, vanishing_degree, field)
end

function get_terms(d::DualNumber{T}) where T
    # Gleb: I think you can do this with a single list comprehension or using `coefficinents` function
    terms = Vector{T}(undef, d.vanishing_degree)
    
    for i in 0:(d.vanishing_degree-1)
        c = coeff(d.poly, [i])
        terms[i+1] = c
    end
    
    return terms
end

# Gleb: there is a special function for truncation https://nemocas.github.io/AbstractAlgebra.jl/stable/polynomial/#Base.truncate-Tuple{PolyRingElem,%20Int64}
function truncate_poly(poly::fpMPolyRingElem, vanishing_degree::Integer)
    R = parent(poly)
    result = R(0)

    exps = collect(exponent_vectors(poly))
    coeffs = [coeff(poly, i) for i in 1:length(exps)]
    
    for (c, e) in zip(coeffs, exps)
        if e[1] < vanishing_degree
            result += c * gen(R, 1)^e[1]
        end
    end
    
    return result
end

function convert_to_field(x::fpMPolyRingElem, F::fpField)::fpFieldElem
    if x == 0
        return F(0)
    end

    exps = collect(exponent_vectors(x))

    if isempty(exps) || (length(exps) == 1 && all(e -> e == 0, exps[1]))
        a = constant_coefficient(x)
        return a
    else
        throw(ArgumentError("Cannot convert non-constant polynomial to field element"))
    end
end

# ------------------------------ Dual Number Math ------------------------------
Base.:*(a::DualNumber{T}, b::DualNumber{T}) where T = begin
    degree = min(a.vanishing_degree, b.vanishing_degree)
    result_poly = truncate_poly(a.poly * b.poly, degree)
    return DualNumber{T}(result_poly, degree, a.field)
end

Base.:+(a::DualNumber{T}, b::DualNumber{T}) where T = begin
    degree = min(a.vanishing_degree, b.vanishing_degree)
    result_poly = truncate_poly(a.poly + b.poly, degree)
    return DualNumber{T}(result_poly, degree, a.field)
end

Base.:-(a::DualNumber{T}, b::DualNumber{T}) where T = begin
    degree = min(a.vanishing_degree, b.vanishing_degree)
    result_poly = truncate_poly(a.poly - b.poly, degree)
    return DualNumber{T}(result_poly, degree, a.field)
end

Base.:-(a::DualNumber{T}) where T = begin
    result_poly = -a.poly
    return DualNumber{T}(result_poly, a.vanishing_degree, a.field)
end

Base.:^(a::DualNumber{T}, n::Integer) where T = begin
    if n == 0
        R = parent(a.poly)
        return DualNumber{T}(R(1), a.vanishing_degree, a.field)
    elseif n == 1
        return a
    end
    
    result_poly = truncate_poly(a.poly^n, a.vanishing_degree)
    return DualNumber{T}(result_poly, a.vanishing_degree, a.field)
end

Base.:*(a::DualNumber{T}, b::fpMPolyRingElem) where T = begin
    if !isnothing(a.field)
        b_converted = convert_to_field(b, a.field)
        result_poly = truncate_poly(a.poly * b_converted, a.vanishing_degree)
        return DualNumber{T}(result_poly, a.vanishing_degree, a.field)
    end
    throw(ArgumentError("Cannot multiply DualNumber without field with polynomial"))
end

Base.:*(b::fpMPolyRingElem, a::DualNumber{T}) where T = begin
    return a * b
end

Base.:*(a::DualNumber{T}, b::Any) where T = begin
    result_poly = truncate_poly(a.poly * b, a.vanishing_degree)
    return DualNumber{T}(result_poly, a.vanishing_degree, a.field)
end

Base.:*(b::Any, a::DualNumber{T}) where T = begin
    return a * b
end

Base.:+(a::DualNumber{T}, b::fpMPolyRingElem) where T = begin
    if !isnothing(a.field)
        b_converted = convert_to_field(b, a.field)
        R = parent(a.poly)
        result_poly = truncate_poly(a.poly + R(b_converted), a.vanishing_degree)
        return DualNumber{T}(result_poly, a.vanishing_degree, a.field)
    end
    throw(ArgumentError("Cannot add DualNumber without field with polynomial"))
end

Base.:+(b::fpMPolyRingElem, a::DualNumber{T}) where T = begin
    return a + b
end

Base.:+(a::DualNumber{T}, b::Any) where T = begin
    R = parent(a.poly)
    result_poly = truncate_poly(a.poly + R(b), a.vanishing_degree)
    return DualNumber{T}(result_poly, a.vanishing_degree, a.field)
end

Base.:+(b::Any, a::DualNumber{T}) where T = begin
    return a + b
end

Base.:-(a::DualNumber{T}, b::fpMPolyRingElem) where T = begin
    if !isnothing(a.field)
        b_converted = convert_to_field(b, a.field)
        R = parent(a.poly)
        result_poly = truncate_poly(a.poly - R(b_converted), a.vanishing_degree)
        return DualNumber{T}(result_poly, a.vanishing_degree, a.field)
    end
    throw(ArgumentError("Cannot subtract DualNumber without field with polynomial"))
end

Base.:-(b::fpMPolyRingElem, a::DualNumber{T}) where T = begin
    return -a + b
end

Base.:-(a::DualNumber{T}, b::Any) where T = begin
    R = parent(a.poly)
    result_poly = truncate_poly(a.poly - R(b), a.vanishing_degree)
    return DualNumber{T}(result_poly, a.vanishing_degree, a.field)
end

Base.:-(b::Any, a::DualNumber{T}) where T = begin
    R = parent(a.poly)
    result_poly = truncate_poly(R(b) - a.poly, a.vanishing_degree)
    return DualNumber{T}(result_poly, a.vanishing_degree, a.field)
end

Base.show(io::IO, d::DualNumber) = begin
    terms = get_terms(d)
    terms_str = String[]
    
    for (i, term) in enumerate(terms)
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

function evaluate_poly_with_dual(f::fpMPolyRingElem, args::Vector{DualNumber{fpFieldElem}}, isolate_var)
    field = args[1].field
    vanishing_degree = args[1].vanishing_degree
    
    R, (ε,) = polynomial_ring(field, ["ε"])
    result_poly = R(0)
    
    exps = collect(exponent_vectors(f))
    coeffs = [coeff(f, i) for i in 1:length(exps)]

    for (coeff, exp) in zip(coeffs, exps)
        if exp[1] >= vanishing_degree
            continue
        end
        
        term_poly = R(coeff)
        
        for (i, p) in enumerate(exp)
            if p > 0
                if i == 1 && !isnothing(args[1])
                    dual_pow = args[1]^p
                    term_poly *= dual_pow.poly
                elseif i <= length(args) && !isnothing(args[i])
                    dual_pow = args[i]^p
                    term_poly *= dual_pow.poly
                end
            end
        end
        
        term_poly = truncate_poly(term_poly, vanishing_degree)
        result_poly += term_poly
    end
    
    result_poly = truncate_poly(result_poly, vanishing_degree)
    return DualNumber{fpFieldElem}(result_poly, vanishing_degree, field)
end


function (f::fpMPolyRingElem)(args::Vector{DiffMinPoly.DualNumber{fpFieldElem}}, isolate_var)
    return DiffMinPoly.evaluate_poly_with_dual(f, args, isolate_var)
end
