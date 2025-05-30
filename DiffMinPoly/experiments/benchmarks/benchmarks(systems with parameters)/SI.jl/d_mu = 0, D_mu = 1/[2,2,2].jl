using StructuralIdentifiability
include("../../../../../src/solver_love_and_support.jl")
include("../../../../../src/utils.jl")

ode_prec = rand_ode([2,1])
first(values(find_ioequations(ode_prec)))

ode = @ODEmodel(
               x1'(t) = -75*x1(t)^2*a1 + 60*x1(t)^2*a2 + 93*x1(t)^2 - 26*x1(t)*x2(t)*a1 + 57*x1(t)*x2(t)*a2 - 84*x1(t)*x2(t) + 69*x1(t)*a1 - 60*x1(t)*a2 + 46*x1(t) + 84*x2(t)^2*a1 + 74*x2(t)^2*a2 + 47*x2(t)^2 + 46*x2(t)*a1 + 100*x2(t)*a2 - 44*x2(t) + 95*a1 + 93*a2 - 95,
               x2'(t) = 50*x1(t)^2*a1 + 62*x1(t)^2*a2 + 59*x1(t)^2 + 51*x1(t)*x2(t)*a1 - 73*x1(t)*x2(t)*a2 - 51*x1(t)*x2(t) + 63*x1(t)*a1 + 83*x1(t)*a2 + 63*x1(t) + 61*x2(t)^2*a1 - 73*x2(t)^2*a2 + 39*x2(t)^2 - 43*x2(t)*a1 + 30*x2(t)*a2 + 30*x2(t) - 47*a1 - 56*a2 - 54,
               y(t) =   8*x1(t)^2 - 29*x1(t)*x2(t) + 45*x1(t) + 26*x2(t)^2 + 36*x2(t) + 42
           )

tim = @elapsed io_correct = first(values(find_ioequations(ode)))    
