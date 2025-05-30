import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, y = function ('x1, x2, y')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2], y])

sys = [ 
 Derivative(x2(t), t) - (52*x1(t)**3 - 3*x1(t)**2*x2(t) - 57*x1(t)**2 + 96*x1(t)*x2(t)**2 + 69*x1(t)*x2(t) - 63*x1(t) - 35*x2(t)**3 + 47*x2(t)**2 + 84*x2(t) + 10),
 Derivative(x1(t), t) - (65*x1(t)**3 - 16*x1(t)**2*x2(t) - 18*x1(t)**2 - 62*x1(t)*x2(t)**2 - 77*x1(t)*x2(t) - 94*x1(t) + 65*x2(t)**3 + 55*x2(t)**2 - 49*x2(t) + 6),
y(t) - (-70*x1(t)**3 + x1(t)**2*x2(t) + 2*x1(t)**2 + 39*x1(t)*x2(t)**2 + 84*x1(t)*x2(t) + 15*x1(t) - 85*x2(t)**3 + 63*x2(t)**2 - 19*x2(t) + 51)
]


start_time = time.time()
ideal = R.RosenfeldGroebner(sys)
print("--- %s seconds ---" % (time.time() - start_time)) 