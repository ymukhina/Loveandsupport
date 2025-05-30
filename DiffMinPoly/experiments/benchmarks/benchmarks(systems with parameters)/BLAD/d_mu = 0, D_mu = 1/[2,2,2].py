import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, y, a, b = function ('x1, x2, y, a, b')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2], y, a, b])

sys = [ 
 Derivative(x2(t), t) - (50*x1(t)**2*a(t) + 62*x1(t)**2*b(t) + 59*x1(t)**2 + 51*x1(t)*x2(t)*a(t) - 73*x1(t)*x2(t)*b(t) - 51*x1(t)*x2(t) + 63*x1(t)*a(t) + 83*x1(t)*b(t) + 63*x1(t) + 61*x2(t)**2*a(t) - 73*x2(t)**2*b(t) + 39*x2(t)**2 - 43*x2(t)*a(t) + 30*x2(t)*b(t) + 30*x2(t) - 47*a(t) - 56*b(t) - 54),
 Derivative(x1(t), t) - (-75*x1(t)**2*a(t) + 60*x1(t)**2*b(t) + 93*x1(t)**2 - 26*x1(t)*x2(t)*a(t) + 57*x1(t)*x2(t)*b(t) - 84*x1(t)*x2(t) + 69*x1(t)*a(t) - 60*x1(t)*b(t) + 46*x1(t) + 84*x2(t)**2*a(t) + 74*x2(t)**2*b(t) + 47*x2(t)**2 + 46*x2(t)*a(t) + 100*x2(t)*b(t) - 44*x2(t) + 95*a(t) + 93*b(t) - 95),
 Derivative(a(t), t) - 0, 
 Derivative(b(t), t) - 0,   
 y(t) - (8*x1(t)**2 - 29*x1(t)*x2(t) + 45*x1(t) + 26*x2(t)**2 + 36*x2(t) + 42)   
]


start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time)) 