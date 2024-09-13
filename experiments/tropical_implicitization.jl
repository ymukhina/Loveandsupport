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

#cases
r1 = x2^2 + x1*x2 + x1^2 + 1
r2 = x2 

g1 = rand_poly(2, [x1, x2])
g2 = rand_poly(1, [x1, x2])

b1 = rand_poly(2, [x1, x2, x3])
b2 = rand_poly(1, [x1, x2, x3])
b3 = rand_poly(1, [x1, x2, x3])



R =  diff_elimination_polytope_1(r1, r2)
vertices_1 = vertices(R) 
newton_polytope_1 = convex_hull(vertices_1)
np_size_1 = length(lattice_points(newton_polytope_1))
np_equations_1 = facets(newton_polytope_1)    

G =  diff_elimination_polytope_1(g1, g2)
vertices_2 = vertices(G) 
newton_polytope_2 = convex_hull(vertices_2)
np_size_2 = length(lattice_points(newton_polytope_2))
np_equations_2 = facets(newton_polytope_2)

B =  diff_elimination_polytope_2(b1, b2, b3)
vertices_3 = vertices(B)
newton_polytope_3 = convex_hull(vertices_3)
np_size_3 = length(lattice_points(newton_polytope_3))
np_equations_3 = facets(newton_polytope_3)


@info "Newton polytope for the spesial case of the system [2,1] $np_equations_1, size $np_size_1"    
@info "Newton polytope for the system [2,1] $np_equations_2, size $np_size_2" 
@info "Newton polytope for the system [2,1,1] $np_equations_3, size $np_size_3" 