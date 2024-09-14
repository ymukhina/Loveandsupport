# This script produces the table for Section 4.2 of the paper aiming
# at experimental exploration of an alternative approach via tropical implicitization

using TropicalImplicitization, Oscar, IterTools

#Randomize general polynomial of degree d over the ring R
function rand_poly(deg, vars)
    result = 0
    degs = [collect(0:deg) for v in vars]

        for m in IterTools.product(degs...)
            if sum(m) <= deg
                monom = 1 
                for i in 1:length(vars)
                    monom *= vars[i]^m[i]
                end
               result += rand(1:200) * monom
            end
        end

    return result
end

function diff_elimination_polytope(polys)
    vs = gens(parent(polys[1]))
    poly_map = [vs[1], polys[1]]
    for i in 1:(length(polys) - 1)
        newp = 0
        for (i, v) in enumerate(vs)
            newp += derivative(poly_map[end], v) * polys[i]
        end
        push!(poly_map, newp)
    end
    newton_pols = [newton_polytope(f) for f in poly_map]
    cone_list, weight_list = get_tropical_cycle(newton_pols)
    Delta = get_polytope_from_cycle(cone_list, weight_list)
    return Delta
end

cases = []

R, (x1, x2) = PolynomialRing(QQ, ["x1", "x2"])
push!(
    cases,
    Dict(
        :name => "Special system [2, 1]",
        :polys => [x2^2 + x1 * x2 + x1^2 + 1, x2],
    )
)
push!(
    cases,
    Dict(
        :name => "Generic system [2, 1]",
        :polys => [rand_poly(2, [x1, x2]), rand_poly(1, [x1, x2])],
    )
)

R, (x1, x2, x3) = PolynomialRing(QQ, ["x1", "x2", "x3"])

push!(
    cases,
    Dict(
        :name => "Generic system [2, 1, 1]",
        :polys => [rand_poly(2, [x1, x2, x3]), rand_poly(1, [x1, x2, x3]), rand_poly(1, [x1, x2, x3])],
    )
)

push!(
    cases,
    Dict(
        :name => "Generic system [3, 1, 1]",
        :polys => [rand_poly(3, [x1, x2, x3]), rand_poly(1, [x1, x2, x3]), rand_poly(1, [x1, x2, x3])],
    )
)

for c in cases
    polys = c[:polys]
    P =  diff_elimination_polytope(polys)
    vs = vertices(P) 
    newton_polytope = convex_hull(vs)
    np_size = length(lattice_points(newton_polytope))
    np_equations = facets(newton_polytope)    
    @info "Newton polytope for $(c[:name]) is $np_equations, size $np_size"
end
