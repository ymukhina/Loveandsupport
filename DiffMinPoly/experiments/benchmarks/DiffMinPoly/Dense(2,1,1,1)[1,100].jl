using StructuralIdentifiability, DiffMinPoly

ode = @ODEmodel(
    x1'(t) = 30*x1(t)^2 + 60*x1(t)*x2(t) + 9*x1(t)*x3(t) + x1(t)*x4(t) + 6*x1(t) + 2*x2(t)^2 + 4*x2(t)*x3(t) + 55*x2(t)*x4(t) + 51*x2(t) + 64*x3(t)^2 + 48*x3(t)*x4(t) + 4*x3(t) + 35*x4(t)^2 + 18*x4(t) + 16,
    x2'(t) = 32*x1(t) + 95*x2(t) + 60*x3(t) + 25*x4(t) + 56,
    x3'(t) = 33*x1(t) + 28*x2(t) + 20*x3(t) + 60*x4(t) + 36,
    x4'(t) = 64*x1(t) + 17*x2(t) + 20*x3(t) + 3*x4(t) + 45,
    y(t) = x1(t)
)

ode_prec = rand_ode([2,1])
x = first(sort(ode_prec.x_vars, rev = true))
eliminate(ode_prec, x, 0.99)


x = first(sort(ode.x_vars, rev = true))
@elapsed h = eliminate(ode, x, 0.99)
exit()
  