include("Love and support.jl")

#a = 18, b = 71
cata = @ODEmodel(
    x1'(t) = (2 + 18 - 10*(x1(t)^2 + x2(t)^2))*x1(t) + x2(t)^2 + 2*x2(t) + x3(t)^2,
    x2'(t) = - x3(t)^3 - (1+x2(t))*(x2(t)^2 + 2*x2(t) + x3(t)^2) - 4*x1(t) + 18*x2(t), 
    x3'(t) = (1 + x2(t))*x3(t)^2 + x1(t)^2 - 71,
    y(t) = x1(t)
)

russional_reconstruction(cata)

#check_ansatz(rand_ode([3, 3, 3]), 2^31 - 1)
