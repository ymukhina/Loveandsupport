using StructuralIdentifiability
include("../../../../../src/solver_love_and_support.jl")
include("../../../../../src/utils.jl")

ode_prec = rand_ode([2,1])
eliminate(ode_prec)

ode = @ODEmodel(
               x1'(t) = -3*x1(t)^2 - 84*x1(t)*x2(t) - 17*x1(t) + 6*x2(t)^2 + 60*x2(t) - 60,
               x2'(t) =  100*x1(t)^2 - 96*x1(t)*x2(t) + 61*x1(t) - 8*x2(t)^2 + 78*x2(t) - 80,
               y(t) =  -87*x1(t)^3 + 8*x1(t)^2*x2(t) - 53*x1(t)^2 - 69*x1(t)*x2(t)^2 - 72*x1(t)*x2(t) - 85*x1(t) - 71*x2(t)^3 - 53*x2(t)^2 - 21*x2(t) - 79
           )

@elapsed h = eliminate(ode)
