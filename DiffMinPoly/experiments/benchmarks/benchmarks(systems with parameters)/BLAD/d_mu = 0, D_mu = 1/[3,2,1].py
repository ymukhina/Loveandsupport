import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, y, a = function ('x1, x2, y, a')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2], y, a])

sys = [ 
 Derivative(x2(t), t) - (-67*x1(t)**3*a(t) + 64*x1(t)**3 - 54*x1(t)**2*x2(t)*a(t) + 50*x1(t)**2*x2(t) + 61*x1(t)**2*a(t) + 56*x1(t)**2 + 56*x1(t)*x2(t)**2*a(t) - 35*x1(t)*x2(t)**2 - 62*x1(t)*x2(t)*a(t) - 34*x1(t)*x2(t) + 42*x1(t)*a(t) + 49*x1(t) + 98*x2(t)**3*a(t) + 78*x2(t)**3 - 65*x2(t)**2*a(t) - 96*x2(t)**2 - 41*x2(t)*a(t) - 89*x2(t) - 19*a(t) + 63),
 Derivative(x1(t), t) - (90*x1(t)**3*a(t) + 7*x1(t)**3 - 18*x1(t)**2*x2(t)*a(t) - 97*x1(t)**2*x2(t) + 9*x1(t)**2*a(t) - 68*x1(t)**2 + 95*x1(t)*x2(t)**2*a(t) + 65*x1(t)*x2(t)**2 - 36*x1(t)*x2(t)*a(t) + 91*x1(t)*x2(t) + 70*x1(t)*a(t) + 3*x1(t) + 4*x2(t)**3*a(t) + 47*x2(t)**3 + 82*x2(t)**2*a(t) - 6*x2(t)**2 - 61*x2(t)*a(t) - 32*x2(t) - 8*a(t) + 83),
 Derivative(a(t), t) - 0,   
 y(t) - (-60*x1(t)**2 + 6*x1(t)*x2(t) - 60*x1(t) + 32*x2(t)**2 + 42*x2(t) - 73)   
]


start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time)) 