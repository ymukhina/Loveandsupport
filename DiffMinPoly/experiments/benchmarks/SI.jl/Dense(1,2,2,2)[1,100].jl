using StructuralIdentifiability, DiffMinPoly

ode = @ODEmodel(
    x1'(t) = 60*x1(t) + 65*x2(t) + 18*x3(t) + 33*x4(t) + 24,
    x2'(t) = 18*x1(t)^2 + 12*x1(t)*x2(t) + 44*x1(t)*x3(t) + 30*x1(t)*x4(t) + 24*x1(t) + 20*x2(t)^2 + 24*x2(t)*x3(t) + 17*x2(t)*x4(t) + 80*x2(t) + 14*x3(t)^2 + 4*x3(t)*x4(t) + 9*x3(t) + 34*x4(t)^2 + 18*x4(t) + 16,
    x3'(t) = 54*x1(t)^2 + 5*x1(t)*x2(t) + 5*x1(t)*x3(t) + 35*x1(t)*x4(t) + 70*x1(t) + 40*x2(t)^2 + 20*x2(t)*x3(t) + 36*x2(t)*x4(t) + 14*x2(t) + 39*x3(t)^2 + 45*x3(t)*x4(t) + 45*x3(t) + 20*x4(t)^2 + 35*x4(t) + 76,
    x4'(t) = 85*x1(t)^2 + 24*x1(t)*x2(t) + 60*x1(t)*x3(t) + 48*x1(t)*x4(t) + 2*x1(t) + 36*x2(t)^2 + 10*x2(t)*x3(t) + 32*x2(t)*x4(t) + 25*x2(t) + 50*x3(t)^2 + 68*x3(t)*x4(t) + 10*x3(t) + 44*x4(t)^2 + 20*x4(t) + 15,
    y(t) = x1(t)
)

ode_prec = rand_ode([2,1])
first(values(find_ioequations(ode_prec)))

tim = @elapsed io_correct = first(values(find_ioequations(ode)))
exit()


