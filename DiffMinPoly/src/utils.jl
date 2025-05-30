"""
    subtotal_degree(poly, vars)

Computes the total degree of a polynomial wrt a subset of variables
"""

function subtotal_degree(poly, vars)
    mask = [(v in vars) for v in gens(parent(poly))]
    return maximum(map(e -> sum(e[mask]), exponent_vectors(poly)))
end

# -------- Generating random dense models -------- #

function all_monomials(deg, vars)
    result = []
    degs = [collect(0:deg) for v in vars]

    for m in IterTools.product(degs...)
        if reduce(+, m, init = 0) <= deg
            push!(result, reduce(* , [v^e for (v, e) in zip(vars, m)], init = 1))
        end
    end

    return result
end

"""
    rand_poly(deg, vars)

Computes the polynomial of multi-degree 'deg' in groups of variables 'vars' 
with coefficient sampled uniformly at random from the integers in [-10, 10].

"""
                                
function rand_poly(deg, vars)
    result = 0

    for m in IterTools.product([all_monomials(d, v) for (d, v) in zip(deg, vars)]...)
        c = rand(-10:10)
        c = (c == 0) ? 1 : c
        monom = c * prod(m)
        result += monom
    end

    return result
end                                
                                
"""
    rand_ode_x(degs)

Computes the polynomial ODE model `ode` with the right hand side of degrees ‘degs’ and output equal to x1.
                                            
"""
                                                
function rand_ode_x(degs::Vector{Int}; char=0)
    n = length(degs)
    F = iszero(char) ? QQ : GF(char)
    R, vars = polynomial_ring(QQ, vcat(["x$i(t)" for i in 1:n], ["y(t)"]))
    return StructuralIdentifiability.ODE{Ptype}(
        vars[1:n],
        [vars[end]],
        Dict(vars[i] => rand_poly([degs[i]], [vars[1:n]]) for i in 1:n),
        Dict(vars[end] => vars[1]),
        Ptype[]
    )
end

"""
    rand_ode(degs)

Computes the polynomial ODE model `ode` with the right hand side of (bi-)degrees ‘degs’ such that y(t) has the (bi-)degree degs[end].
                                            
"""

function rand_ode(degs::Vector{Int}; char=0)
    return rand_ode([(d, 0) for d in degs])
end

function rand_ode(degs::Vector{Tuple{Int, Int}}; char=0, num_params=0)
    n = length(degs)
    F = iszero(char) ? QQ : GF(char)
    R, vars = polynomial_ring(QQ, vcat(["x$i(t)" for i in 1:n-1], ["a$i" for i in 1:num_params], ["y(t)"]))
    return StructuralIdentifiability.ODE{Ptype}(
        vars[1:n-1],
        [vars[end]],
        Dict(vars[i] => rand_poly(degs[i], [vars[1:n-1], vars[(n:n + num_params - 1)]]) for i in 1:n-1),
        Dict(vars[end] => rand_poly(degs[n], [vars[1:n-1], vars[n:n + num_params - 1]])),
        Ptype[]
    )
end


