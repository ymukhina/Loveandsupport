import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, x3, y = function ('x1, x2, x3, y')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2, x3], y])

sys = [ 
 Derivative(x3(t), t) - (x3(t) * (1 - 9*x1(t)/10 - 7*x2(t)/10 - 2*x3(t)/10 +  2*(x3(t)**3)/10 )), 
 Derivative(x2(t), t) - (x2(t) * (1 -  2*x1(t)/10 - 8*x2(t)/10 -  5*x3(t)/10 +  6 * (x3(t)**3)/10)),
 Derivative(x1(t), t) - (x1(t) * (1 -  x1(t)/10 - 3 * x2(t)/10 - 5 * x3(t)/10) + 4 * (x2(t)**2)/10 + 7 * (x3(t)**2)/10),
y(t) - (x1(t))
]

start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time))