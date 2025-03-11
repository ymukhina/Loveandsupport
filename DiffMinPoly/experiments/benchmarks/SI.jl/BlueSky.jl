using StructuralIdentifiability, DiffMinPoly

ode = @ODEmodel(
               x1'(t) = (2 + 0.456 - 10 * (x1(t)^2 + x2(t)^2) ) * x1(t) + x2(t)^2 + 2 * x2(t) + x3(t)^2,
               x2'(t) = - x3(t)^3 - (1 + x2(t)) * (x2(t)^2 + 2 * x2(t) + x3(t)^2) - 4 * x1(t) + 0.456 * x2(t),
               x3'(t) = (1 + x2(t)) * x3(t)^2 + x1(t)^2 - 0.0357,
               y(t) = x1(t)
           )
           
ode_prec = rand_ode([2,1])
first(values(find_ioequations(ode_prec)))

tim = @elapsed io_correct = first(values(find_ioequations(ode)))
exit()
