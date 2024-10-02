using DiffMinPoly
using Oscar
using Nemo
using StructuralIdentifiability
using IterTools
import GLPK
import StructuralIdentifiability: _reduce_mod_p, reduce_ode_mod_p, power_series_solution, ps_diff 

using DynamicPolynomials
using MixedSubdivisions

const Ptype = QQMPolyRingElem

# -------- Test function -------- #

# test ansatz against IO enemy equation 
function check_ansatz_modp(ode::ODE, p::Int)

    x = first(sort(ode.x_vars, rev = true))
    ord = DiffMinPoly.minpoly_order(ode, x) 
    possible_supp = DiffMinPoly.f_min_support(ode, x, ord)
    
    @info "Solving with love and support!"
    tim = @elapsed io_tocheck = eliminate_with_love_and_support_modp(ode, x, p, ord, possible_supp; info = false)
    @info "time: $(tim) seconds"

    @info "Solving without love and support :("
    tim = @elapsed io_correct = first(values(find_ioequations(ode)))
    io_correct = _reduce_mod_p(io_correct, p)
    io_correct *= Oscar.leading_coefficient(io_correct)^(-1)
    @info "time: $(tim) seconds"

    R = parent(io_tocheck)
    S = parent(io_correct)

    # the variables of io_correct
    n = ngens(R)
    println(gens(R))
    println(gens(S))
    # the worst thing in the known universe
    svnames = (string).(S.S)
    y_index = findfirst(vname -> vname == "y(t)_0", svnames)
    ys = gens(S)[y_index:y_index+(n-1)]
    @info "IO variables $(ys)"

    phi = hom(R, S, ys)

    quot, rem = divrem(phi(io_tocheck), io_correct)
    return iszero(rem) && iszero(total_degree(quot))
end

# Gleb: there are no tests for the main function!
# Yulia: Changed to eliminate below
function check_ansatz(ode::ODE)
        
    x = first(sort(ode.x_vars, rev = true))    

    @info "Solving with love and support!"
    tim = @elapsed io_tocheck = DiffMinPoly.eliminate(ode, x)
    @info "time: $(tim) seconds"
   

    @info "Solving without love and support :("
    tim = @elapsed io_correct = first(values(find_ioequations(ode)))
    io_correct *= (Oscar.leading_coefficient(io_correct)^(-1))
    @info "time: $(tim) seconds"
   

    R = parent(io_tocheck)
    S = parent(io_correct)
        
    n = ngens(R) 
  
    svnames = (string).(S.S)
    y_index = findfirst(vname -> vname == "y(t)_0", svnames)
    ys = gens(S)[y_index:y_index+(n-1)]
    @info "IO variables $(ys)"    

    phi = hom(R, S, ys)
    quot, rem = divrem(phi(io_tocheck), io_correct)
    return iszero(rem) && iszero(total_degree(quot))
end
        
# To circumvent a StackOverflow in Oscar
function my_convex_hull(points)
    pts = permutedims(hcat([Oscar.homogenize(v, 1) for v in points]...))
    ptype = Oscar._scalar_type_to_polymake(QQFieldElem)
    lin = zero_matrix(QQ, 0, size(pts, 2))
    return Oscar.Polyhedron{QQFieldElem}(Polymake.polytope.Polytope{ptype}(; VERTICES = pts, LINEALITY_SPACE=lin))
end
   
"""
test_table_2()

Computes the Newton polytope of the minimal polynomial for the systems in the Table 4.2
            
Code for the Section 4.2 of the article. 
        
"""            
            
function test_table_2()   
   
    cases = [
       rand_ode([2,1])
       rand_ode([2,1,1])         
    ]
                
    push!(
    cases,
    @ODEmodel(
        x1'(t) = x1(t)^2 + x2(t)^2,
        x2'(t) = x2(t) + 1,
        y(t) = x1(t)
    )
)            
                
    
 
     for c in cases
        ode = c
        x = first(sort(ode.x_vars, rev = true))
        jac = DiffMinPoly.minpoly_order(ode, x)                
                        
        supp = DiffMinPoly.f_min_support(ode, x, jac)
                    
        l = my_convex_hull(supp) 
        F = facets(l)            
        @info "Newton polytope for the system $c is"
        @info "$F"            
     end          
end
            
            

"""
test_table_3()

Computes the Mixed Fiber Polytopes for the systems (2,1,1), (3,1,1) and 
computes the size for its support.

Code for the Section 4.3 of the article.            
        
"""               
            
function test_table_3()
        
        cases = [      
            #[2, 1, 1],
            #[3, 1, 1],
            [3, 2, 2],        
]         
                
        R, (x3, x2, x1, x11, x111, x1111) = polynomial_ring(GF(65521), ["x3","x2", "x1", "x11", "x111", "x1111"])

        for vars in cases        
                
        g1 = rand_poly(vars[1], [x1, x2, x3])
        g2 = rand_poly(vars[2], [x1, x2, x3])
        g3 = rand_poly(vars[3], [x1, x2, x3])
                       

        polys = [x11 - g1, x111 - (derivative(g1, x1) * x11 + derivative(g1, x2) * g2 + derivative(g1, x3) * g3)]
        push!(polys, x1111 - (derivative(polys[2], x1) * x11 + derivative(polys[2], x2) * g2 + derivative(polys[2], x3) * g3 + derivative(polys[2], x11) * x111))

        polys_new = [sum([rand(1:1000)*m for m in Oscar.monomials(p)]) for p in polys]    

        S = parent(first(polys))
        gb = groebner_basis_f4(Oscar.ideal(S, polys_new), eliminate = 2, complete_reduction=true, info_level=2)
        m = collect(exponent_vectors(gb[1]))
        l = convexhull(m...)
        L = removevredundancy(l, GLPK.Optimizer)
        size = length(lattice_points(my_convex_hull(m)))
                    
        @info "Bound for mixed Fiber polytope for the case $vars is $L"
        @info "Number of points inside Mixed Fiber Polytope | $size " 
          
        end                                   
end            
                    
# Counting the integer points inside guessed polytopes and the ones from the theorem

function bound_difference( d1::Int, d2::Int, p::Int) # counting for systems [d1, d2, 1:p]
    for i in 1:p
        @info "System | Num of the int points inside the polytope | Num of monomials in the min pol"
        ode = rand_ode([d1,d2,i])
        x = first(sort(ode.x_vars, rev = true))                
        s = size(DiffMinPoly.f_min_support(ode, x, minpoly_order(ode, x)))[1]
               
        l = length(eliminate_with_love_and_support_modp(ode, x, Int(rand_bits_prime(ZZ, 32)), info = true))            
        @info "[$d1,$d2,$i] | $s | $l"
     end      
end
