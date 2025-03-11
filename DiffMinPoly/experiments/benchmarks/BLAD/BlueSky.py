import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, x3, y = function ('x1, x2, x3, y')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2, x3], y])

sys = [
 Derivative(x3(t), t) - ((1 + x2(t)) * x3(t)**2 + x1(t)**2 - 357/10000), 
 Derivative(x2(t), t) - (- x3(t)**3 - (1 + x2(t)) * (x2(t)**2 + 2 * x2(t) + x3(t)**2) - 1 * x1(t) + 456 * x2(t)/1000),
 Derivative(x1(t), t) - ((2 + 456/1000 - 10 * (x1(t)**2 + x2(t)**2) ) * x1(t) + x2(t)**2 + 2 * x2(t) + x3(t)**2),
y(t) - (x1(t))
]


start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time))