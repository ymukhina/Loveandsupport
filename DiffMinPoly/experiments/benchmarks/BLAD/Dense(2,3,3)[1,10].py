import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, x3, y = function ('x1, x2, x3, y')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2, x3], y])

sys = [ 
 Derivative(x3(t), t) - (4*x1(t)**3 + 6*x1(t)**2*x2(t) + 2*x1(t)**2*x3(t) + 2*x1(t)**2 + 2*x1(t)*x2(t)**2 + 2*x1(t)*x2(t)*x3(t) + 4*x1(t)*x2(t) + 4*x1(t)*x3(t)**2 + 3*x1(t)*x3(t) + 4*x1(t) + 2*x2(t)**3 + x2(t)**2*x3(t) + 5*x2(t)**2 + 2*x2(t)*x3(t)**2 + 3*x2(t)*x3(t) + 2*x2(t) + 6*x3(t)**3 + x3(t)**2 + 4*x3(t) + 5), 
 Derivative(x2(t), t) - (6*x1(t)**3 + 3*x1(t)**2*x2(t) + 6*x1(t)**2*x3(t) + 10*x1(t)**2 + 6*x1(t)*x2(t)**2 + 4*x1(t)*x2(t)*x3(t) + 8*x1(t)*x2(t) + 4*x1(t)*x3(t)**2 + 3*x1(t)*x3(t) + 6*x1(t) + x2(t)**3 + 6*x2(t)**2*x3(t) + 2*x2(t)**2 + 10*x2(t)*x3(t)**2 + 3*x2(t)*x3(t) + 4*x2(t) + 4*x3(t)**3 + 4*x3(t)**2 + 4*x3(t) + 4),
 Derivative(x1(t), t) - (8*x1(t)**2 + 2*x1(t)*x2(t) + 3*x1(t)*x3(t) + 6*x1(t) + 8*x2(t)**2 + 4*x2(t)*x3(t) + 2*x2(t) + 3*x3(t)**2 + 8*x3(t) + 1),
y(t) - (x1(t))
]

start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time))