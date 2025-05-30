import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, y = function ('x1, x2, y')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2], y])

sys = [ 
 Derivative(x2(t), t) - (60*x1(t)**3 + 34*x1(t)**2*x2(t) - 74*x1(t)**2 - 27*x1(t)*x2(t)**2 + 66*x1(t)*x2(t) + 67*x1(t) + 87*x2(t)**3 - 86*x2(t)**2 - 12*x2(t) - 1),
 Derivative(x1(t), t) - (-61*x1(t)**3 + 75*x1(t)**2*x2(t) - 25*x1(t)**2 - 38*x1(t)*x2(t)**2 + 84*x1(t)*x2(t) - 66*x1(t) + 25*x2(t)**3 + 35*x2(t)**2 + 55*x2(t) + 52),
y(t) - (24*x1(t)**4 + 88*x1(t)**3*x2(t) - 64*x1(t)**3 + 24*x1(t)**2*x2(t)**2 - 15*x1(t)**2*x2(t) + 99*x1(t)**2 - 7*x1(t)*x2(t)**3 - 48*x1(t)*x2(t)**2 - 33*x1(t)*x2(t) - 32*x1(t) - 42*x2(t)**4 - 67*x2(t)**3 + 82*x2(t)**2 + 97*x2(t) + 92)
]


start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time)) 