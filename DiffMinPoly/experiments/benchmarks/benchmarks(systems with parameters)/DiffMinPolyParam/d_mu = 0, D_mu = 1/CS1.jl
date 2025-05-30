using StructuralIdentifiability

include("../../../../../src/solver_love_and_support.jl")
include("../../../../../src/utils.jl")

ode_prec = rand_ode([2,1])
eliminate(ode_prec)

ode = @ODEmodel(
               x1'(t) = x1(t) * (0.1 + 0.9 * x1(t) + a1 * x2(t)),
               x2'(t) = x2(t) * (0.5 - a1 * x1(t) + 0.2 * x2(t)),
               y(t) =   0.9 * x1(t)^2 + 0.2 * x2(t)^2
           )

@elapsed h = eliminate(ode)