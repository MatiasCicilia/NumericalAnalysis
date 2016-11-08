import numpy
import math

def fun(x, y):
    return -2*(x*y)


def runge_kutta(h, x0, y0):
    x = x0 + h
    k1 = h * fun(x0, y0)
    k2 = h * fun(x0 + (h/2), y0 + (k1/2))
    k3 = h * fun(x0 + (h/2), y0 + (k2/2))
    k4 = h * fun(x0 + h, y0 + k3)
    m = (k1 + 2*k2 + 2*k3 + k4)/6
    y = y0 + m
    return [x, y]


def runge_kutta_iter(iterations, h, x0, y0):
    result = runge_kutta(h, x0, y0)
    for i in range(1, iterations):
        result = runge_kutta(h, result[0], result[1])
    return result


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

for i in range(1,25):
     print "Resultado de trapecios para f (%d) = %f " % (i,trapezoidal_rule(f,2.0,3.0,i))

for i in range(1,25):
     print "Resultado de trapecios para g (%d) = %f " % (i, trapezoidal_rule(g, 1.0, 3.0, i))

result = romberg(f, 2.0, 3.0, 15, 10**-6)
last = result[0]
print "Result romberg %d , %f" % (result[0],result[1][last-1][last-1])

result = romberg(g, 1.0, 3.0, 15,10**-6)
last = result[0]
print "Result romberg %d , %f" % (result[0],result[1][last-1][last-1])

h = 0.1
y0 = 1
x0 = 0
for i in range(1, 11):
    result = runge_kutta_iter(i,h,x0,y0)
    print "X(%d) = %f" % (i, result[0])
    print "Y(%d) = %f" % (i, result[1])
