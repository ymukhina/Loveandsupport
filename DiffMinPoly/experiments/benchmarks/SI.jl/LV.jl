using StructuralIdentifiability, DiffMinPoly

ode = @ODEmodel(
    x1'(t) = x1(t) * (1 - 0.1 * x1(t) - 0.3 * x2(t) - 0.5 * x3(t)) + 0.4 * x2(t)^2 + 0.7 * x3(t)^2,
    x2'(t) = x2(t) * (1 - 0.2 * x1(t) - 0.8 * x2(t) - 0.5 * x3(t) + 0.6 * x3(t)^3),
    x3'(t) = x3(t) * (1 - 0.9 * x1(t) - 0.7 * x2(t) - 0.2 * x3(t) + 0.2 * x3(t)^3),
    y(t) = x1(t)
)

ode_prec = rand_ode([2,1])
first(values(find_ioequations(ode_prec)))

tim = @elapsed io_correct = first(values(find_ioequations(ode)))
exit()
