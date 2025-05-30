import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, y, a = function ('x1, x2, y, a')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2], y, a])

sys = [ 
 Derivative(x2(t), t) - (47*x1(t)**2*a(t) - 90*x1(t)**2 - 12*x1(t)*x2(t)*a(t) - 92*x1(t)*x2(t) - 17*x1(t)*a(t) + 72*x1(t) - 30*x2(t)**2*a(t) - 22*x2(t)**2 + 88*x2(t)*a(t) + 13*x2(t) - 21*a(t) + 85),
 Derivative(x1(t), t) - (62*x1(t)**2*a(t) + 22*x1(t)**2 - 47*x1(t)*x2(t)*a(t) + 81*x1(t)*x2(t) - 64*x1(t)*a(t) + 33*x1(t) - 77*x2(t)**2*a(t) - 56*x2(t)**2 - 17*x2(t)*a(t) - 59*x2(t) - 22*a(t) - 71),
 Derivative(a(t), t) - 0,   
 y(t) - (8*x1(t)**2 - 49*x1(t)*x2(t) - 97*x1(t) + 25*x2(t)**2 - 34*x2(t) - 55)   
]


start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time)) 