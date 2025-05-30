using StructuralIdentifiability

include("../../../../../src/solver_love_and_support.jl")
include("../../../../../src/utils.jl")

ode_prec = rand_ode([2,1])
first(values(find_ioequations(ode_prec)))

ode = @ODEmodel(
               x1'(t) =65*x1(t)^3 - 16*x1(t)^2*x2(t) - 18*x1(t)^2 - 62*x1(t)*x2(t)^2 - 77*x1(t)*x2(t) - 94*x1(t) + 65*x2(t)^3 + 55*x2(t)^2 - 49*x2(t) + 6,
               x2'(t) =  52*x1(t)^3 - 3*x1(t)^2*x2(t) - 57*x1(t)^2 + 96*x1(t)*x2(t)^2 + 69*x1(t)*x2(t) - 63*x1(t) - 35*x2(t)^3 + 47*x2(t)^2 + 84*x2(t) + 10,
               y(t) =  -70*x1(t)^3 + x1(t)^2*x2(t) + 2*x1(t)^2 + 39*x1(t)*x2(t)^2 + 84*x1(t)*x2(t) + 15*x1(t) - 85*x2(t)^3 + 63*x2(t)^2 - 19*x2(t) + 51
           )

tim = @elapsed io_correct = first(values(find_ioequations(ode)))   