using StructuralIdentifiability
include("../../../../../src/solver_love_and_support.jl")
include("../../../../../src/utils.jl")

ode_prec = rand_ode([2,1])
first(values(find_ioequations(ode_prec)))

ode = @ODEmodel(
               x1'(t) = -8*x1(t)*a1 + 85*x1(t)*a2 + 96*x1(t) + 86*x2(t)*a1 - x2(t)*a2 + 47*x2(t) - 18*a1 - 50*a2 - 49,
               x2'(t) = 44*x1(t)*a1 - 23*x1(t)*a2 + 74*x1(t) - 22*x2(t)*a1 - 57*x2(t)*a2 - 13*x2(t) - 92*a1 - 48*a2 + 92,
               y(t) =   5*x1(t)^2*a1 + 95*x1(t)^2*a2 - 21*x1(t)^2 + 75*x1(t)*x2(t)*a1 + 49*x1(t)*x2(t)*a2 - 44*x1(t)*x2(t) + 50*x1(t)*a1 + 64*x1(t)*a2 + 23*x1(t) + 80*x2(t)^2*a1 + 52*x2(t)^2*a2 - 71*x2(t)^2 - 18*x2(t)*a1 - 4*x2(t)*a2 - 14*x2(t) - 19*a1 - 30*a2 + 43
           )

tim = @elapsed io_correct = first(values(find_ioequations(ode)))