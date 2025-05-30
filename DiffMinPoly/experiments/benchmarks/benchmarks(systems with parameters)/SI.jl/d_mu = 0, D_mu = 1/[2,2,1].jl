using StructuralIdentifiability
include("../../../../../src/solver_love_and_support.jl")
include("../../../../../src/utils.jl")

ode_prec = rand_ode([2,1])
first(values(find_ioequations(ode_prec)))

ode = @ODEmodel(
               x1'(t) =62*x1(t)^2*a1 + 22*x1(t)^2 - 47*x1(t)*x2(t)*a1 + 81*x1(t)*x2(t) - 64*x1(t)*a1 + 33*x1(t) - 77*x2(t)^2*a1 - 56*x2(t)^2 - 17*x2(t)*a1 - 59*x2(t) - 22*a1 - 71,
               x2'(t) = 47*x1(t)^2*a1 - 90*x1(t)^2 - 12*x1(t)*x2(t)*a1 - 92*x1(t)*x2(t) - 17*x1(t)*a1 + 72*x1(t) - 30*x2(t)^2*a1 - 22*x2(t)^2 + 88*x2(t)*a1 + 13*x2(t) - 21*a1 + 85,
               y(t) =   8*x1(t)^2 - 49*x1(t)*x2(t) - 97*x1(t) + 25*x2(t)^2 - 34*x2(t) - 55
           )

tim = @elapsed io_correct = first(values(find_ioequations(ode)))    