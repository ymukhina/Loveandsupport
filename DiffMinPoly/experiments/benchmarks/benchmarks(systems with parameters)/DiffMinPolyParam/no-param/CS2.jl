using StructuralIdentifiability
include("../../../../../src/solver_love_and_support.jl")
include("../../../../../src/utils.jl")

ode_prec = rand_ode([2,1])
eliminate(ode_prec)

ode = @ODEmodel(
               x1'(t) = x1(t) * (0.3 + 0.9 * x1(t) + 0.1 * x2(t)),
               x2'(t) = x2(t) * (0.4 + 0.8 * x1(t) + 0.3 * x2(t)),
               y(t) =   x1(t)^3 + x2(t)^3
           )


@elapsed h = eliminate(ode)