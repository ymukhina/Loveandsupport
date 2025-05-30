import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, y = function ('x1, x2, y')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2], y])

sys = [ 
 Derivative(x2(t), t) - (15*x1(t)**4 + 75*x1(t)**3*x2(t) - 36*x1(t)**3 - 40*x1(t)**2*x2(t)**2 + 17*x1(t)**2*x2(t) + 89*x1(t)**2 - 47*x1(t)*x2(t)**3 + 73*x1(t)*x2(t)**2 + 91*x1(t)*x2(t) - 79*x1(t) + x2(t)**4 - 85*x2(t)**3 - 47*x2(t)**2 - 73*x2(t) - 68),
 Derivative(x1(t), t) - (52*x1(t)**4 - 43*x1(t)**3*x2(t) + 45*x1(t)**3 - 5*x1(t)**2*x2(t)**2 + x1(t)**2*x2(t) + 100*x1(t)**2 - 43*x1(t)*x2(t)**3 + 29*x1(t)*x2(t)**2 - 95*x1(t)*x2(t) + 24*x1(t) + 48*x2(t)**4 - 36*x2(t)**3 - 33*x2(t)**2 + 90*x2(t) - 58),
y(t) - (78*x1(t)**2 - 33*x1(t)*x2(t) + 21*x1(t) + 93*x2(t)**2 + 51*x2(t) + 24)
]


start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time)) 