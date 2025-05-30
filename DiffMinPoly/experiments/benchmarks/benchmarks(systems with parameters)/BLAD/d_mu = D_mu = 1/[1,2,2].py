import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, y, a, b = function ('x1, x2, y, a, b')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2], y, a, b])

sys = [ 
 Derivative(x2(t), t) - (44*x1(t)*a(t) - 23*x1(t)*b(t) + 74*x1(t) - 22*x2(t)*a(t) - 57*x2(t)*b(t) - 13*x2(t) - 92*a(t) - 48*b(t) + 92),
 Derivative(x1(t), t) - (-8*x1(t)*a(t) + 85*x1(t)*b(t) + 96*x1(t) + 86*x2(t)*a(t) - x2(t)*b(t) + 47*x2(t) - 18*a(t) - 50*b(t) - 49),
 Derivative(a(t), t) - 0, 
 Derivative(b(t), t) - 0,
 y(t) - (5*x1(t)**2*a(t) + 95*x1(t)**2*b(t) - 21*x1(t)**2 + 75*x1(t)*x2(t)*a(t) + 49*x1(t)*x2(t)*b(t) - 44*x1(t)*x2(t) + 50*x1(t)*a(t) + 64*x1(t)*b(t) + 23*x1(t) + 80*x2(t)**2*a(t) + 52*x2(t)**2*b(t) - 71*x2(t)**2 - 18*x2(t)*a(t) - 4*x2(t)*b(t) - 14*x2(t) - 19*a(t) - 30*b(t) + 43)   
]


start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time)) 