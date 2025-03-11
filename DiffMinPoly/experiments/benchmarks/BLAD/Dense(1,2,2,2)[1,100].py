import time
from sympy import *
from DifferentialAlgebra import *
init_printing ()

t = var ('t')
x1, x2, x3, x4, y = function ('x1, x2, x3, x4, y')

R = DifferentialRing(derivations = [t], blocks = [[x1, x2, x3, x4], y])

sys = [
 Derivative(x4(t), t) - (85*x1(t)**2 + 24*x1(t)*x2(t) + 60*x1(t)*x3(t) + 48*x1(t)*x4(t) + 2*x1(t) + 36*x2(t)**2 + 10*x2(t)*x3(t) + 32*x2(t)*x4(t) + 25*x2(t) + 50*x3(t)**2 + 68*x3(t)*x4(t) + 10*x3(t) + 44*x4(t)**2 + 20*x4(t) + 15),  
 Derivative(x3(t), t) - (54*x1(t)**2 + 5*x1(t)*x2(t) + 5*x1(t)*x3(t) + 35*x1(t)*x4(t) + 70*x1(t) + 40*x2(t)**2 + 20*x2(t)*x3(t) + 36*x2(t)*x4(t) + 14*x2(t) + 39*x3(t)**2 + 45*x3(t)*x4(t) + 45*x3(t) + 20*x4(t)**2 + 35*x4(t) + 76), 
 Derivative(x2(t), t) - (18*x1(t)**2 + 12*x1(t)*x2(t) + 44*x1(t)*x3(t) + 30*x1(t)*x4(t) + 24*x1(t) + 20*x2(t)**2 + 24*x2(t)*x3(t) + 17*x2(t)*x4(t) + 80*x2(t) + 14*x3(t)**2 + 4*x3(t)*x4(t) + 9*x3(t) + 34*x4(t)**2 + 18*x4(t) + 16),
 Derivative(x1(t), t) - (60*x1(t) + 65*x2(t) + 18*x3(t) + 33*x4(t) + 24),
y(t) - (x1(t))
]


f = open("./dense_1222_1100_py.txt", "w+")  