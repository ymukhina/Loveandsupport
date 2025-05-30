import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, y, a = function ('x1, x2, y, a')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2], y, a])

sys = [ 
 10 * Derivative(x2(t), t) - (x2(t) * (5 - 10 * a(t) * x1(t) + 2 * x2(t))),
 10 * Derivative(x1(t), t) - (x1(t) * (1 + 9 * x1(t) + 10 * a(t) * x2(t))),
 Derivative(a(t), t) - 0,   
 10 * y(t) - (9 * x1(t)**2 + 2 * x2(t)**2)   
]


start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time)) 