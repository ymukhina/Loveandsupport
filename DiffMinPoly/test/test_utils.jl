using DiffMinPoly
using Oscar
using Nemo
using StructuralIdentifiability
import StructuralIdentifiability: _reduce_mod_p, reduce_ode_mod_p, power_series_solution, ps_diff 

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

    @info "Solving with Love & Support!"
    tim = @elapsed io_tocheck = DiffMinPoly.eliminate(ode, x)
    @info "Time: $(tim) seconds"
   

    @info "Solving without Love & Support :("
    tim = @elapsed io_correct = first(values(find_ioequations(ode)))
    io_correct *= (Oscar.leading_coefficient(io_correct)^(-1))
    @info "Time: $(tim) seconds"
   

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
        
