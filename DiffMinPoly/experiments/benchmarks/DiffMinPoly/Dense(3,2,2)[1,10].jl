using StructuralIdentifiability, DiffMinPoly

ode = @ODEmodel(x1'(t) = 4*x1(t)^3 + 10*x1(t)^2*x2(t) + 5*x1(t)^2*x3(t) + 4*x1(t)^2 + x1(t)*x2(t)^2 + 5*x1(t)*x2(t)*x3(t) + 2*x1(t)*x2(t) + 3*x1(t)*x3(t)^2 + x1(t)*x3(t) + 8*x1(t) + 5*x2(t)^3 + 2*x2(t)^2*x3(t) + 2*x2(t)^2 + 2*x2(t)*x3(t)^2 + 10*x2(t)*x3(t) + x2(t) + 5*x3(t)^3 + 4*x3(t)^2 + 5*x3(t) + 8,
       x2'(t) = 2*x1(t)^2 + 10*x1(t)*x2(t) + 8*x1(t)*x3(t) + 10*x1(t) + x2(t)^2 + 6*x2(t)*x3(t) + 4*x2(t) + 2*x3(t)^2 + 4*x3(t) + 10,
       x3'(t) = 8*x1(t)^2 + 3*x1(t)*x2(t) + 3*x1(t)*x3(t) + 4*x1(t) + 10*x2(t)^2 + 2*x2(t)*x3(t) + 2*x2(t) + 5*x3(t)^2 + 4*x3(t) + 6,
       y(t) = x1(t)
       )
       
ode_prec = rand_ode([2,1])
x = first(sort(ode_prec.x_vars, rev = true))
eliminate(ode_prec, x, 0.99)


x = first(sort(ode.x_vars, rev = true))
@elapsed h = eliminate(ode, x, 0.99)
exit()