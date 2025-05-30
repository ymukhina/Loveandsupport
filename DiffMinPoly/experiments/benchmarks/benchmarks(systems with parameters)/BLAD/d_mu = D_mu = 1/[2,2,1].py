import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, y, a = function ('x1, x2, y, a')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2], y, a])

sys = [ 
 Derivative(x2(t), t) - (5*x1(t)**2*a(t) - 31*x1(t)**2 - 5*x1(t)*x2(t)*a(t) - 8*x1(t)*x2(t) - 57*x1(t)*a(t) - 33*x1(t) - 13*x2(t)**2*a(t) - 53*x2(t)**2 - 73*x2(t)*a(t) + 28*x2(t) + 87*a(t) + 88),
 Derivative(x1(t), t) - (-51*x1(t)**2*a(t) - 96*x1(t)**2 + 44*x1(t)*x2(t)*a(t) - 68*x1(t)*x2(t) - 43*x1(t)*a(t) + 32*x1(t) - 80*x2(t)**2*a(t) - 72*x2(t)**2 - 49*x2(t)*a(t) - 7*x2(t) - 39*a(t) - 24),
 Derivative(a(t), t) - 0,   
 y(t) - ( -91*x1(t)**2*a(t) + 96*x1(t)**2 - 59*x1(t)*x2(t)*a(t) - 81*x1(t)*x2(t) - 91*x1(t)*a(t) - 22*x1(t) + 84*x2(t)**2*a(t) - 30*x2(t)**2 - 3*x2(t)*a(t) - 81*x2(t) - 59*a(t) - 73)   
]


start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time)) 