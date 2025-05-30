using StructuralIdentifiability
include("../../../../../src/solver_love_and_support.jl")
include("../../../../../src/utils.jl")

ode_prec = rand_ode([2,1])
first(values(find_ioequations(ode_prec)))

ode = @ODEmodel(
               x1'(t) = -12*x1(t)*a1 + 43*x1(t) + 12*x2(t)*a1 + 40*x2(t) - 92*a1 + 49,
               x2'(t) = -22*x1(t)*a1 + 88*x1(t) - 86*x2(t)*a1 + 91*x2(t) + 55*a1 + 85,
               y(t) =   -41*x1(t)^3*a1 - 76*x1(t)^3 + 27*x1(t)^2*x2(t)*a1 + 21*x1(t)^2*x2(t) - 92*x1(t)^2*a1 + 38*x1(t)^2 - 42*x1(t)*x2(t)^2*a1 - 94*x1(t)*x2(t)^2 - 3*x1(t)*x2(t)*a1 + 77*x1(t)*x2(t) + 20*x1(t)*a1 + 46*x1(t) + 40*x2(t)^3*a1 - 91*x2(t)^3 + 38*x2(t)^2*a1 + 12*x2(t)^2 - 73*x2(t)*a1 - 89*x2(t) - 24*a1 - 8
           )

tim = @elapsed io_correct = first(values(find_ioequations(ode)))