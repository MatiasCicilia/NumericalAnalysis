import numpy

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

def romberg( f, a, b, n, cota):
    """
    INPUT:
        f       - function to integrate,
        [a, b]  - the interval of integration,
        n       - number of levels of recursion

    OUTPUT:
        numpy float array - Romberg integration array; most accurate
                            answer should be at bottom of right-most column,
    """

    matrix = numpy.array( [[0] * (n+1)] * (n+1), float )
    h = b - a
    matrix[0,0] = 0.5 * h * ( f(a) + f(b) )

    powerOf2 = 1
    for i in range( 1, n + 1 ):
        h = 0.5 * h

        sum = 0.0
        powerOf2 = 2 * powerOf2
        for k in range( 1, powerOf2, 2 ):
            sum = sum + f( a + k * h )

        matrix[i,0] = 0.5 * matrix[i-1,0] + sum * h

        powerOf4 = 1
        for j in xrange( 1, i + 1 ):
            powerOf4 = 4 * powerOf4
            matrix[i,j] = matrix[i,j-1] + ( matrix[i,j-1] - matrix[i-1,j-1] ) / ( powerOf4 - 1 )

    return matrix


def trapezoidal_rule(f, a, b, n):
    h = (b - a) / n
    s = f(a) + f(b)
    for i in xrange(1, n):
        s += 2 * f(a + i * h)
    return s * h / 2


"""
MAIN.................................................................................
"""
h = 0.1
y0 = 1
x0 = 0
for i in range(1, 11):
    result = runge_kutta_iter(i,h,x0,y0)
    print "X(%d) = %f" % (i, result[0])
    print "Y(%d) = %f" % (i, result[1])
