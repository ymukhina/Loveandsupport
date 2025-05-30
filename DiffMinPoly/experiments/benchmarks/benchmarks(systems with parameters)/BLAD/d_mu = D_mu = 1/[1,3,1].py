import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, y, a = function ('x1, x2, y, a')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2], y, a])

sys = [ 
 Derivative(x2(t), t) - (-22*x1(t)*a(t) + 88*x1(t) - 86*x2(t)*a(t) + 91*x2(t) + 55*a(t) + 85),
 Derivative(x1(t), t) - (-12*x1(t)*a(t) + 43*x1(t) + 12*x2(t)*a(t) + 40*x2(t) - 92*a(t) + 49),
 Derivative(a(t), t) - 0,   
 y(t) - (-41*x1(t)**3*a(t) - 76*x1(t)**3 + 27*x1(t)**2*x2(t)*a(t) + 21*x1(t)**2*x2(t) - 92*x1(t)**2*a(t) + 38*x1(t)**2 - 42*x1(t)*x2(t)**2*a(t) - 94*x1(t)*x2(t)**2 - 3*x1(t)*x2(t)*a(t) + 77*x1(t)*x2(t) + 20*x1(t)*a(t) + 46*x1(t) + 40*x2(t)**3*a(t) - 91*x2(t)**3 + 38*x2(t)**2*a(t) + 12*x2(t)**2 - 73*x2(t)*a(t) - 89*x2(t) - 24*a(t) - 8)   
]


start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time)) 