import numpy as np
import matplotlib.pyplot as plt

#Точное решение: y = c*e^x -x/3 -1/3
x0 = 0.5
y0= 0.6

def getC():
    return (y0 + x0/3 +1/3)/np.exp(x0)

def getExactY(x,c):
    return c*np.exp(x) -x/3 -1/3

def f(x, y):
    return x/3 +y

def calcErrorLastPoint(y, b, p):
    m = len(y)-1
    err = abs(getExactY(b, getC())-y[m])/((2**p)-1)
    return err

def calculateError(y1, y2, p):
    n = len(y1)-1
    m = len(y2)-1
    err = abs(y1[n]-y2[m])/((2**p)-1)
    return err

def eulerMethod(a, b, y0, n):
    h = (b-a)/n
    x = np.linspace(a, b, n+1)
    y = np.zeros(n+1)
    y[0] = y0
    for i in range(n):
        y[i+1] = y[i] + h*f(x[i], y[i])
    return y

def runge_kutta_method(a, b, y0, n):
    h = (b-a)/n
    x = np.linspace(a, b, n+1)
    y = np.zeros(n+1)
    y[0] = y0
    for i in range(n):
        k1 = h*f(x[i], y[i])
        k2 = h*f(x[i] + 1/2*h, y[i] + 1/2*k1)
        k3 = h*f(x[i] + 3/4*h, y[i] + 3/4*k2)
        y[i+1] = y[i] + (2*k1 + 3*k2 + 4*k3)/9
    return y

def requiredError(a, b, y0, p, e, method):
    n = 1
    prevErr = -1
    h=0
    
    while True:
        y1 = method(a, b, y0, n)
        y2 = method(a, b, y0, 2*n)
        h = (b-a)/2*n
        err = calculateError(y1, y2, p)
        if err <= e:
            return y2
        if ((prevErr - err) == 0):
            print('процесс решения прекращен, т.к. с уменьшением шага погрешность не уменьшается')
            return y2
        if (h == 0):
            print('процесс решения прекращен, т.к. значение шага стало недопустимо малым')
            return y2
        prevErr = err
        n *= 2



a = 0.5
b = 4
n = 16
e = 0.0000000001

#y = requiredError(a, b, y0, 1, e, eulerMethod)
y2 = requiredError(a, b, y0, 3, e, runge_kutta_method)

error = calcErrorLastPoint(y2, b, 3)

print(error)










    