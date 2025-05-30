import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, y = function ('x1, x2, y')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2], y])

sys = [ 
 Derivative(x2(t), t) - (-43*x1(t)**2 - 94*x1(t)*x2(t) + 66*x1(t) - 32*x2(t)**2 + x2(t) - 76),
 Derivative(x1(t), t) - (9*x1(t)**2 + 32*x1(t)*x2(t) + 93*x1(t) + 20*x2(t)**2 - 89*x2(t) - 63),
y(t) - (7*x1(t)**4 + 43*x1(t)**3*x2(t) + 90*x1(t)**3 + x1(t)**2*x2(t)**2 + 86*x1(t)**2*x2(t) - 21*x1(t)**2 - 55*x1(t)*x2(t)**3 + 32*x1(t)*x2(t)**2 + 76*x1(t)*x2(t) + 69*x1(t) - 41*x2(t)**4 + 61*x2(t)**3 + 23*x2(t)**2 - 68*x2(t) - 22)
]


start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time)) 