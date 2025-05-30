import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, y, a = function ('x1, x2, y, a')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2], y, a])

sys = [ 
 Derivative(x2(t), t) - (71*x1(t)*a(t) - 58*x1(t) + 54*x2(t)*a(t) - 65*x2(t) + 67*a(t) + 85),
 Derivative(x1(t), t) - (14*x1(t)*a(t) - 45*x1(t) + 48*x2(t)*a(t) - 11*x2(t) - 4*a(t) - 76),
 Derivative(a(t), t) - 0,   
 y(t) - (10*x1(t)**3 - 35*x1(t)**2*x2(t) - 64*x1(t)**2 + 3*x1(t)*x2(t)**2 + 67*x1(t)*x2(t) + 43*x1(t) + 3*x2(t)**3 + 86*x2(t)**2 - 3*x2(t) + 62)   
]


start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time)) 