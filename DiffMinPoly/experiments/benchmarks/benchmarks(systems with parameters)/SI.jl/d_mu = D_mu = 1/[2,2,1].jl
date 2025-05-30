using StructuralIdentifiability
include("../../../../../src/solver_love_and_support.jl")
include("../../../../../src/utils.jl")

ode_prec = rand_ode([2,1])
first(values(find_ioequations(ode_prec)))

ode = @ODEmodel(
               x1'(t) = -51*x1(t)^2*a1 - 96*x1(t)^2 + 44*x1(t)*x2(t)*a1 - 68*x1(t)*x2(t) - 43*x1(t)*a1 + 32*x1(t) - 80*x2(t)^2*a1 - 72*x2(t)^2 - 49*x2(t)*a1 - 7*x2(t) - 39*a1 - 24,
               x2'(t) = 5*x1(t)^2*a1 - 31*x1(t)^2 - 5*x1(t)*x2(t)*a1 - 8*x1(t)*x2(t) - 57*x1(t)*a1 - 33*x1(t) - 13*x2(t)^2*a1 - 53*x2(t)^2 - 73*x2(t)*a1 + 28*x2(t) + 87*a1 + 88,
               y(t) =   -91*x1(t)^2*a1 + 96*x1(t)^2 - 59*x1(t)*x2(t)*a1 - 81*x1(t)*x2(t) - 91*x1(t)*a1 - 22*x1(t) + 84*x2(t)^2*a1 - 30*x2(t)^2 - 3*x2(t)*a1 - 81*x2(t) - 59*a1 - 73
           )

tim = @elapsed io_correct = first(values(find_ioequations(ode)))