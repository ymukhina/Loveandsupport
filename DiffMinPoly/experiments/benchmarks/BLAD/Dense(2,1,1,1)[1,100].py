import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, x3, x4, y = function ('x1, x2, x3, x4, y')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2, x3, x4], y])

sys = [
 Derivative(x4(t), t) - (64*x1(t) + 17*x2(t) + 20*x3(t) + 3*x4(t) + 45),  
 Derivative(x3(t), t) - (33*x1(t) + 28*x2(t) + 20*x3(t) + 60*x4(t) + 36), 
 Derivative(x2(t), t) - (32*x1(t) + 95*x2(t) + 60*x3(t) + 25*x4(t) + 56),
 Derivative(x1(t), t) - (30*x1(t)**2 + 60*x1(t)*x2(t) + 9*x1(t)*x3(t) + x1(t)*x4(t) + 6*x1(t) + 2*x2(t)**2 + 4*x2(t)*x3(t) + 55*x2(t)*x4(t) + 51*x2(t) + 64*x3(t)**2 + 48*x3(t)*x4(t) + 4*x3(t) + 35*x4(t)**2 + 18*x4(t) + 16),
y(t) - (x1(t))
]

start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time))