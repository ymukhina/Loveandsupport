using StructuralIdentifiability, DiffMinPoly

ode = @ODEmodel(
    x1'(t) = 12*x1(t)^2 + 14*x1(t)*x2(t) + 39*x1(t)*x3(t) + 21*x1(t) + 24*x2(t)^2 + 76*x2(t)*x3(t) + 32*x2(t) + 14*x3(t)^2 + 34*x3(t) + 56,
    x2'(t) = 48*x1(t)^3 + 54*x1(t)^2*x2(t) + 9*x1(t)^2*x3(t) + 22*x1(t)^2 + 37*x1(t)*x2(t)^2 + 50*x1(t)*x2(t)*x3(t) + 76*x1(t)*x2(t) + 32*x1(t)*x3(t)^2 + 2*x1(t)*x3(t) + 49*x1(t) + 33*x2(t)^3 + 18*x2(t)^2*x3(t) + 7*x2(t)^2 + 36*x2(t)*x3(t)^2 + 13*x2(t)*x3(t) + 36*x2(t) + 12*x3(t)^3 + 31*x3(t)^2 + 38*x3(t) + 16,
    x3'(t) = 13*x1(t)^3 + 2*x1(t)^2*x2(t) + 14*x1(t)^2*x3(t) + 27*x1(t)^2 + 34*x1(t)*x2(t)^2 + 42*x1(t)*x2(t)*x3(t) + 41*x1(t)*x2(t) + 22*x1(t)*x3(t)^2 + 14*x1(t)*x3(t) + 54*x1(t) + 38*x2(t)^3 + 38*x2(t)^2*x3(t) + 45*x2(t)^2 + 49*x2(t)*x3(t)^2 + 82*x2(t)*x3(t) + 42*x2(t) + 74*x3(t)^3 + 34*x3(t)^2 + 27*x3(t) + 36,
    y(t) = x1(t)
    )

ode_prec = rand_ode([2,1])
first(values(find_ioequations(ode_prec)))

tim = @elapsed io_correct = first(values(find_ioequations(ode)))
exit()

