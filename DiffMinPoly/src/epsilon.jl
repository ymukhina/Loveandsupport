using Oscar

struct DualNumber{T}
    poly::fpPolyRingElem
    vanishing_degree::Integer

    function DualNumber{T}(poly::fpPolyRingElem, vanishing_degree::Integer) where T
        new{T}(poly, vanishing_degree)
    end

    function DualNumber{T}(poly::fpFieldElem, vanishing_degree::Integer) where T
        F = parent(poly)
        R, _ = polynomial_ring(F, "ε")
        new{T}(R(poly), vanishing_degree)
    end
end

# ------------------------------ Dual Number Utils ------------------------------

function Epsilon(vanishing_degree::Integer, field::fpField)
    vanishing_degree < 1 && throw(ArgumentError("Vanishing degree must be positive"))
    R, ε = polynomial_ring(field, "ε")
    return DualNumber{fpFieldElem}(R(ε), vanishing_degree)
end

function truncate_poly(poly::fpPolyRingElem, vanishing_degree::Integer)
    return truncate(poly, vanishing_degree)
end

function convert_to_field(x::fpPolyRingElem)::fpFieldElem
    if iszero(x)
        return base_ring(x)(0)
    end
    if degree(x) == 0
        return coeff(x, 0)
    else
        throw(ArgumentError("Cannot convert non-constant polynomial to field element"))
    end
end

# ------------------------------ Dual Number Math ------------------------------

Base.:*(a::DualNumber{T}, b::DualNumber{T}) where T = begin
    degree = min(a.vanishing_degree, b.vanishing_degree)
    result_poly = truncate_poly(a.poly * b.poly, degree)
    return DualNumber{T}(result_poly, degree)
end

Base.:+(a::DualNumber{T}, b::DualNumber{T}) where T = begin
    degree = min(a.vanishing_degree, b.vanishing_degree)
    result_poly = truncate_poly(a.poly + b.poly, degree)
    return DualNumber{T}(result_poly, degree)
end

Base.:-(a::DualNumber{T}, b::DualNumber{T}) where T = begin
    degree = min(a.vanishing_degree, b.vanishing_degree)
    result_poly = truncate_poly(a.poly - b.poly, degree)
    return DualNumber{T}(result_poly, degree)
end

Base.:-(a::DualNumber{T}) where T = begin
    result_poly = -a.poly
    return DualNumber{T}(result_poly, a.vanishing_degree)
end

Base.:^(a::DualNumber{T}, n::Integer) where T = begin
    if n == 0
        R = parent(a.poly)
        return DualNumber{T}(R(1), a.vanishing_degree)
    elseif n == 1
        return a
    end
    result_poly = truncate_poly(a.poly^n, a.vanishing_degree)
    return DualNumber{T}(result_poly, a.vanishing_degree)
end

Base.:*(a::DualNumber{T}, b::fpFieldElem) where T = begin
    result_poly = truncate_poly(a.poly * b, a.vanishing_degree)
    return DualNumber{T}(result_poly, a.vanishing_degree)
end

Base.:*(b::fpFieldElem, a::DualNumber{T}) where T = a * b

Base.:*(a::DualNumber{T}, b::Any) where T = begin
    result_poly = truncate_poly(a.poly * b, a.vanishing_degree)
    return DualNumber{T}(result_poly, a.vanishing_degree)
end

Base.:*(b::Any, a::DualNumber{T}) where T = a * b

Base.:+(a::DualNumber{T}, b::fpFieldElem) where T = begin
    R = parent(a.poly)
    result_poly = truncate_poly(a.poly + R(b), a.vanishing_degree)
    return DualNumber{T}(result_poly, a.vanishing_degree)
end

Base.:+(b::fpFieldElem, a::DualNumber{T}) where T = a + b

Base.:+(a::DualNumber{T}, b::Any) where T = begin
    R = parent(a.poly)
    result_poly = truncate_poly(a.poly + R(b), a.vanishing_degree)
    return DualNumber{T}(result_poly, a.vanishing_degree)
end

Base.:+(b::Any, a::DualNumber{T}) where T = a + b

Base.:-(a::DualNumber{T}, b::fpFieldElem) where T = begin
    R = parent(a.poly)
    result_poly = truncate_poly(a.poly - R(b), a.vanishing_degree)
    return DualNumber{T}(result_poly, a.vanishing_degree)
end

Base.:-(b::fpFieldElem, a::DualNumber{T}) where T = -a + b

Base.:-(a::DualNumber{T}, b::Any) where T = begin
    R = parent(a.poly)
    result_poly = truncate_poly(a.poly - R(b), a.vanishing_degree)
    return DualNumber{T}(result_poly, a.vanishing_degree)
end

Base.:-(b::Any, a::DualNumber{T}) where T = begin
    R = parent(a.poly)
    result_poly = truncate_poly(R(b) - a.poly, a.vanishing_degree)
    return DualNumber{T}(result_poly, a.vanishing_degree)
end

function get_terms(d::DualNumber{T}) where T                # USED ONLY FOR PRETTY-PRINTING PURPOSES, POINTLESS TO IMPROVE
    [coeff(d.poly, i) for i in 0:(d.vanishing_degree-1)]
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

function evaluate_poly_with_dual(f::fpMPolyRingElem, args::Vector{DualNumber{fpFieldElem}})
    field = base_ring(args[1].poly)
    vanishing_degree = args[1].vanishing_degree

    R, _ = polynomial_ring(field, "ε")
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
    return DualNumber{fpFieldElem}(result_poly, vanishing_degree)
end

function (f::fpMPolyRingElem)(args::Vector{DualNumber{fpFieldElem}})
    return evaluate_poly_with_dual(f, args)
end