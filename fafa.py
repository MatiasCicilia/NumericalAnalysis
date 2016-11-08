import numpy
import math

def runge_kutta(f, iterations, h, x, y):
    value = runge_kutta_formula(f,h,x,y)
    for i in range(1, iterations):
        value = runge_kutta_formula(f, h, value[0], value[1])
    return value

def runge_kutta_formula(f, h, x, y):
    x_new = x + h
    k1 = h * f(x, y)
    k2 = h * f(x + (h / 2), y + (k1 / 2))
    k3 = h * f(x + (h / 2), y + (k2 / 2))
    k4 = h * f(x + h, y + k3)
    m = (k1 + 2*k2 + 2*k3 + k4)/6
    y_new = y + m
    return [x_new,y_new]

def romberg(f, a, b, p,error):
    matrix = numpy.zeros((p, p))
    for k in range(0, p):

        matrix[k, 0] = trapezoidal_rule(f, a, b, 2 ** k)

        for j in range(0, k):
            matrix[k, j + 1] = (4 ** (j + 1) * matrix[k, j] - matrix[k - 1, j]) / (4 ** (j + 1) - 1)

        if((abs(matrix[k,k]-matrix[k,k-1]) < error)) :
           return k,matrix

    return p,matrix

def trapezoidal_rule(f, a, b, n):
    h = (b - a) / n
    x = a

    sum = f(a)
    for k in range(1, n):
        x = x + h
        sum += 2 * f(x)

    return (sum + f(b)) * h * 0.5

"""
MAIN.................................................................................
"""

f = lambda x: math.sin(x)/x
g = lambda x: math.log(x, math.e)
h = lambda x,y: -2.0*x*y

#-------------------------------TRAPECIOS-------------------------------
for i in range(1,25):
     print "Resultado de trapecios para f (%d) = %f " % (i,trapezoidal_rule(f,2.0,3.0,i))

for i in range(1,25):
     print "Resultado de trapecios para g (%d) = %f " % (i, trapezoidal_rule(g, 1.0, 3.0, i))

#-------------------------------ROMBERG-------------------------------
result = romberg(f, 2.0, 3.0, 15, 10**-6)
last = result[0]
print "Result romberg %d , %f" % (result[0],result[1][last-1][last-1])

result = romberg(g, 1.0, 3.0, 15,10**-6)
last = result[0]
print "Result romberg %d , %f" % (result[0],result[1][last-1][last-1])

#-------------------------------RUNGE KUTTA-------------------------------
for i in range(1, 11):
    result = runge_kutta(h, i, 0.1, 0, 1)
    print "X,Y(%d) = %f,%f" % (i, result[0],result[1])