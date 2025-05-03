# -------- Generating random dense models -------- #

"""
    rand_poly(deg, vars)

Computes the polynomial of degree 'deg' in variables 'vars' 
with coefficient sampled uniformly at random from the integers in [-1000, 1000].

"""
                                
function rand_poly(deg, vars)
    result = 0
    degs = [collect(0:deg) for v in vars]

        for m in IterTools.product(degs...)
            if sum(m) <= deg
                monom = rand(1:5)
                for i in 1:length(vars)
                    monom *= vars[i]^m[i]
                end
               result += rand(1:20) * monom
            end
        end

    return result
end                                
                                
"""
    rand_ode(degs)

Computes the polynomial ODE model `ode` with the right hand side of degrees ‘degs’.
                                            
"""
                                                
function rand_ode(degs::Vector{Int}; char=0)
    n = length(degs)
    F = iszero(char) ? QQ : GF(char)
    R, vars = polynomial_ring(QQ, vcat(["x$i(t)" for i in 1:n], ["y(t)"]))
    return StructuralIdentifiability.ODE{Ptype}(
        vars[1:n],
        [vars[end]],
        Dict(vars[i] => rand_poly(degs[i], vars[1:n]) for i in 1:n),
        Dict(vars[end] => vars[1]),
        Ptype[]
    )
end

"""
    rand_ode_y(degs)

Computes the polynomial ODE model `ode` with the right hand side of degrees ‘degs’ such that y(t) has the degree degs[end].
                                            
"""

function rand_ode_y(degs::Vector{Int}; char=0)
    n = length(degs)
    F = iszero(char) ? QQ : GF(char)
    R, vars = polynomial_ring(QQ, vcat(["x$i(t)" for i in 1:n-1], ["y(t)"]))
    return StructuralIdentifiability.ODE{Ptype}(
        vars[1:n-1],
        [vars[end]],
        Dict(vars[i] => rand_poly(degs[i], vars[1:n-1]) for i in 1:n-1),
        Dict(vars[end] => rand_poly(degs[n], vars[1:n-1])),
        Ptype[]
    )
end


