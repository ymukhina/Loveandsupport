import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, y = function ('x1, x2, y')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2], y])

sys = [ 
 10 * Derivative(x2(t), t) - (x2(t) * (4 + 8 * x1(t) + 3 * x2(t))),
 10 * Derivative(x1(t), t) - (x1(t) * (3 + 9 * x1(t) + 1 * x2(t))),
y(t) - (x1(t)**3 + x2(t)**3)
]


start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time)) 