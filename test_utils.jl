using Oscar
using Nemo
using StructuralIdentifiability
using IterTools
import StructuralIdentifiability: _reduce_mod_p, reduce_ode_mod_p, power_series_solution, ps_diff 

const Ptype = QQMPolyRingElem

# -------- Test function -------- #
# test interpolation_ansatz against IO enemy equation 
function check_interpol_ansatz_modp(ode::ODE, p::Int)

    @info "Solving with love and support!"
    tim = @elapsed io_tocheck = interpolation_with_love_and_support_modp(ode, p)
    @info "time: $(tim) seconds"

    @info "Solving without love and support :("
    tim = @elapsed io_correct = first(values(find_ioequations(ode)))
    io_correct = _reduce_mod_p(io_correct, p)
    @info "time: $(tim) seconds"

    R = parent(io_tocheck)
    S = parent(io_correct)

    # the variables of io_correct
    n = ngens(R)
    ys = gens(S)[2 * (n - 1) + 2 : 3 * (n - 1) + 2]
    @info "IO variables $(ys)"

    phi = hom(R, S, ys)

    quot, rem = divrem(phi(io_tocheck), io_correct)
    return iszero(rem) && iszero(total_degree(quot))
end
    



# test ansatz against IO enemy equation 
function check_ansatz_modp(ode::ODE, p::Int)

    @info "Solving with love and support!"
    tim = @elapsed io_tocheck = eliminate_with_love_and_support_modp(ode, p)
    @info "time: $(tim) seconds"

    @info "Solving without love and support :("
    tim = @elapsed io_correct = first(values(find_ioequations(ode)))
    io_correct = _reduce_mod_p(io_correct, p)
    @info "time: $(tim) seconds"

    R = parent(io_tocheck)
    S = parent(io_correct)

    # the variables of io_correct
    n = ngens(R)
    ys = gens(S)[2 * (n - 1) + 2 : 3 * (n - 1) + 2]
    @info "IO variables $(ys)"

    phi = hom(R, S, ys)

    quot, rem = divrem(phi(io_tocheck), io_correct)
    return iszero(rem) && iszero(total_degree(quot))
end
    
function check_ansatz(ode::ODE)

    @info "Solving with love and support!"
    tim = @elapsed io_tocheck = eliminate_with_love_and_support(ode)
    @info "time: $(tim) seconds"
    println(Oscar.terms(io_tocheck))

    @info "Solving without love and support :("
    tim = @elapsed io_correct = first(values(find_ioequations(ode)))
    io_correct *= (Nemo.leading_coefficient(io_correct)^(-1))
    @info "time: $(tim) seconds"
    println(Oscar.terms(io_correct))

    R = parent(io_tocheck)
    S = parent(io_correct)

    # the variables of io_correct
    n = ngens(R) 
    ys = gens(S)[2 * (n - 1) + 2 : 3 * (n - 1) + 2]
    @info "IO variables $(ys)"

    phi = hom(R, S, ys)
    quot, rem = divrem(phi(io_tocheck), io_correct)
    # all(c -> c in Oscar.coefficients(io_tocheck), Oscar.coefficients(io_correct))
    Oscar.coefficients(io_tocheck) .- Oscar.coefficients(io_correct)
end
                
                    
                    
                    
# Counting the integer points inside guessed polytopes and the ones from the theorem

function bound_difference(d1::Int, d2::Int, p::Int) # counting for systems [d1, d2, 1:p]
    for i in 1:p
        @info "System | Num of the int points inside the polytope | Num of monomials in the min pol"
        ode = rand_ode([d1,d2,i])
        s = size(f_min_support(ode))[1]
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
 
