import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, y = function ('x1, x2, y')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2], y])

sys = [ 
 Derivative(x2(t), t) - (100*x1(t)**2 - 96*x1(t)*x2(t) + 61*x1(t) - 8*x2(t)**2 + 78*x2(t) - 80),
 Derivative(x1(t), t) - (-3*x1(t)**2 - 84*x1(t)*x2(t) - 17*x1(t) + 6*x2(t)**2 + 60*x2(t) - 60),
y(t) - (-87*x1(t)**3 + 8*x1(t)**2*x2(t) - 53*x1(t)**2 - 69*x1(t)*x2(t)**2 - 72*x1(t)*x2(t) - 85*x1(t) - 71*x2(t)**3 - 53*x2(t)**2 - 21*x2(t) - 79)
]

start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time)) 