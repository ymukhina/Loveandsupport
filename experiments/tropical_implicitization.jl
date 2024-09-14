# This script produces the table for Section 4.2 of the paper aiming
# at experimental exploration of an alternative approach via tropical implicitization

using TropicalImplicitization, Oscar, IterTools

#Randomize general polynomial of degree d over the ring R
function rand_poly(deg, vars)
    result = 0
    degs = [collect(0:deg) for v in vars]

        for m in IterTools.product(degs...)
            if sum(m) <= deg
                monom = rand(1:5)
                for i in 1:length(vars)
                    monom *= vars[i]^m[i]
                end
               result += rand(1:200) * monom
            end
        end

    return result
end

#Find Newton polytope for polnomials f1, f2
function diff_elimination_polytope_1(f1, f2)
           f3 = derivative(f1,gens(parent(f1))[1]) * f1 + derivative(f1,gens(parent(f1))[2]) * f2
           t_list = gens(parent(f1))       
           newton_pols = [newton_polytope(f) for f in [t_list[1], f1, f3]]
           cone_list, weight_list = get_tropical_cycle(newton_pols)
           Delta = get_polytope_from_cycle(cone_list, weight_list)
           return Delta
end

#Find Newton polytope for polnomials f1, f2, f3
function diff_elimination_polytope_2(f1, f2, f3)
                  f4 = derivative(f1,gens(parent(f1))[1]) * f1 + derivative(f1,gens(parent(f1))[2]) * f2 +  derivative(f1,gens(parent(f1))[3]) * f3
                  f5 = derivative(f4,gens(parent(f1))[1]) * f1 + derivative(f4,gens(parent(f1))[2]) * f2 +  derivative(f4,gens(parent(f1))[3]) * f3
           
                  t_list = gens(parent(f1))        
                  newton_pols = [newton_polytope(f) for f in [t_list[1], f1, f4, f5]]
                  cone_list, weight_list = get_tropical_cycle(newton_pols)
                  Delta = get_polytope_from_cycle(cone_list, weight_list)
                  return Delta
end

function diff_elimination_polytope(polys)
    if length(polys) == 2
        return diff_elimination_polytope_1(polys...)
    elseif length(polys) == 3
        return diff_elimination_polytope_2(polys...)
    end
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

for c in cases
    polys = c[:polys]
    P =  diff_elimination_polytope(polys)
    vs = vertices(P) 
    newton_polytope = convex_hull(vs)
    np_size = length(lattice_points(newton_polytope))
    np_equations = facets(newton_polytope)    
    @info "Newton polytope for $(c[:name]) is $np_equations, size $np_size"
end
