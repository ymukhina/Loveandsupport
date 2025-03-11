import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, x3, x4, y = function ('x1, x2, x3, x4, y')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2, x3, x4], y])

sys = [
 Derivative(x4(t), t) - (10*x1(t) + 4*x2(t) + 8*x3(t) + 2*x4(t) + 8),  
 Derivative(x3(t), t) - (4*x1(t) + x2(t) + 4*x3(t) + 4*x4(t) + 2), 
 Derivative(x2(t), t) - (6*x1(t) + x2(t) + 4*x3(t) + 5*x4(t) + 5),
 Derivative(x1(t), t) - (4*x1(t)**2 + 6*x1(t)*x2(t) + 4*x1(t)*x3(t) + 2*x1(t)*x4(t) + 6*x1(t) + 2*x2(t)**2 + 4*x2(t)*x3(t) + 6*x2(t)*x4(t) + 6*x2(t) + 2*x3(t)**2 + 8*x3(t)*x4(t) + 3*x3(t) + 2*x4(t)**2 + 10*x4(t) + 4),
y(t) - (x1(t))
]

start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time))