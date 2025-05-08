using Oscar
using Nemo
using StructuralIdentifiability
import StructuralIdentifiability: _reduce_mod_p, reduce_ode_mod_p, power_series_solution, ps_diff, switch_ring

using DiffMinPoly
using DiffMinPoly: f_min_support

const Ptype = QQMPolyRingElem

# -------- Test function -------- #

# test ansatz against IO enemy equation 
function check_ansatz_modp(ode::ODE, p::Int)

    ord = minpoly_order(ode) 
    possible_supp = f_min_support(ode, ord)
    
    @info "Solving with love and support!"
    tim = @elapsed io_tocheck = eliminate_with_love_and_support_modp(ode, p, ord, possible_supp; info = false)
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
    # the worst thing in the known universe
    svnames = (string).(S.S)
    y_index = findfirst(isequal("y(t)_0"), svnames)
    ys = gens(S)[y_index:y_index+(n-1)]
    params = [switch_ring(p, S) for p in ode.parameters]

    phi = hom(R, S, vcat(params, ys))

    quot, rem = divrem(phi(io_tocheck), io_correct)
    return iszero(rem) && iszero(total_degree(quot)) && !iszero(quot)
end


function check_ansatz(ode::ODE)

    @info "Solving with love and support!"
    tim = @elapsed io_tocheck = DiffMinPoly.eliminate(ode)
    @info "time: $(tim) seconds"
   

    @info "Solving without love and support :("
    tim = @elapsed io_correct = first(values(find_ioequations(ode)))
    io_correct *= (Oscar.leading_coefficient(io_correct)^(-1))
    @info "time: $(tim) seconds"
   

    R = parent(io_tocheck)
    S = parent(io_correct)
        
    n = ngens(R) - length(ode.parameters) 
  
    svnames = (string).(S.S)
    y_index = findfirst(vname -> vname == "y(t)_0", svnames)
    ys = gens(S)[y_index:y_index+(n-1)]
    params = [switch_ring(p, S) for p in ode.parameters]

    phi = hom(R, S, vcat(params, ys))
    quot, rem = divrem(phi(io_tocheck), io_correct)
    return iszero(rem) && iszero(total_degree(quot)) && !iszero(quot)
end

