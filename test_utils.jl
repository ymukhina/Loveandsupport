using Oscar
using Nemo
using StructuralIdentifiability
using IterTools
import StructuralIdentifiability: _reduce_mod_p, reduce_ode_mod_p, power_series_solution, ps_diff 

const Ptype = QQMPolyRingElem

# -------- Test function -------- #

# test ansatz against IO enemy equation 
function check_ansatz_modp(ode::ODE, p::Int)

    @info "Solving with love and support!"
    tim = @elapsed io_tocheck = eliminate_with_love_and_support_modp(ode, p)
    println(io_tocheck)
    @info "time: $(tim) seconds"

    @info "Solving without love and support :("
    tim = @elapsed io_correct = first(values(find_ioequations(ode)))
    io_correct = _reduce_mod_p(io_correct, p)
    io_correct *= Oscar.leading_coefficient(io_correct)^(-1)
    println(io_correct)
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
    
function check_ansatz(ode::ODE)

    @info "Solving with love and support!"
    tim = @elapsed io_tocheck = eliminate_with_love_and_support(ode, rand(1:2^30))[1]
            println(io_tocheck)
    @info "time: $(tim) seconds"
   

    @info "Solving without love and support :("
    tim = @elapsed io_correct = first(values(find_ioequations(ode)))
    io_correct *= (Oscar.leading_coefficient(io_correct)^(-1))
        println(io_correct)
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
    # all(c -> c in Oscar.coefficients(io_tocheck), Oscar.coefficients(io_correct))
    #Oscar.coefficients(io_tocheck) .- Oscar.coefficients(io_correct)
    return iszero(rem) && iszero(total_degree(quot))
end
                
                    
                    
                    
# Counting the integer points inside guessed polytopes and the ones from the theorem

function bound_difference(d1::Int, d2::Int, p::Int) # counting for systems [d1, d2, 1:p]
    for i in 1:p
        @info "System | Num of the int points inside the polytope | Num of monomials in the min pol"
        ode = rand_ode([d1,d2,i])
        s = size(f_min_support(ode, jacobian_check(ode)))[1]
        l = length(eliminate_with_love_and_support_modp(ode, Int(rand_bits_prime(ZZ, 32)), info = true))            
        @info "[$d1,$d2,$i] | $s | $l"
     end      
end



# ------------------------------- #


# -------- Helper functions -------- #

#Randomize general polynomial
function rand_poly(deg, vars)
    result = 0
    degs = [collect(0:deg) for v in vars]

    for m in IterTools.product(degs...)
        if sum(m) <= deg
            monom = rand(1:5)
            for i in 1:length(vars)
                monom *= vars[i]^m[i]
            end
            result += rand(1:1000) * monom
        end
    end

    return result
end

function rand_ode(degs::Vector{Int})
    n = length(degs)
    R, vars = polynomial_ring(QQ, vcat(["x$i(t)" for i in 1:n], ["y(t)"]))
    return StructuralIdentifiability.ODE{Ptype}(
        Dict(vars[i] => rand_poly(degs[i], vars[1:n]) for i in 1:n),
        Dict(vars[end] => vars[1]),
        Ptype[]
    )
end
 
