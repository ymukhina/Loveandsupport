using StructuralIdentifiability
include("../../../../../src/solver_love_and_support.jl")
include("../../../../../src/utils.jl")

ode_prec = rand_ode([2,1])
eliminate(ode_prec)

ode = @ODEmodel(
               x1'(t) =14*x1(t)*a1 - 45*x1(t) + 48*x2(t)*a1 - 11*x2(t) - 4*a1 - 76,
               x2'(t) = 71*x1(t)*a1 - 58*x1(t) + 54*x2(t)*a1 - 65*x2(t) + 67*a1 + 85,
               y(t) =   10*x1(t)^3 - 35*x1(t)^2*x2(t) - 64*x1(t)^2 + 3*x1(t)*x2(t)^2 + 67*x1(t)*x2(t) + 43*x1(t) + 3*x2(t)^3 + 86*x2(t)^2 - 3*x2(t) + 62
           )

@elapsed h = eliminate(ode)