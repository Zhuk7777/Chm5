import numpy as np
import matplotlib.pyplot as plt

def f1(x,u1,u2):
    return 3*u1 + 2*u2 + 3*np.exp(2*x)
    
def f2(x,u1,u2):
    return u1 + 2*u2 +np.exp(2*x)

def calculate_error(u1_n, u1_2n, p):
    err = abs(u1_n - u1_2n)/((2**p)-1)
    return err

def runge_kutta_method(a, b, u10, u20, n):
    
        x = np.linspace(a, b, n+1)
        u1 = np.zeros(n+1)
        u2 = np.zeros(n+1)
        u1[0], u2[0] = u10, u20
        h = (b - a) / n
        for i in range(n):
            k11, k12 = h*f1(x[i], u1[i], u2[i]), h*f2(x[i], u1[i], u2[i])
            k21, k22 = h*f1(x[i] + h/4, u1[i] + 1/4 * k11, u2[i] + 1/4 * k12), h*f2(x[i] + h/4, u1[i] + 1/4 * k11, u2[i] + 1/4 * k12)
            k31, k32 = h*f1(x[i] + h/2, u1[i] + 1/2 * k21, u2[i] + 1/2 * k22), h*f2(x[i] + h/2, u1[i] + 1/2 * k21, u2[i] + 1/2 * k22)
            k41, k42 = h*f1(x[i] + h, u1[i] + k11 - 2*k21 + 2*k31, u2[i] + k12 - 2*k22 + 2*k32), h*f2(x[i] + h, u1[i] + k11 - 2*k21 + 2*k31, u2[i] + k12 - 2*k22 + 2*k32)
            u1[i+1] = u1[i] + 1/6 * (k11 + 4*k31 + k41)
            u2[i+1] = u2[i] + 1/6 * (k12 + 4*k32 + k42)
        return u1[n], u2[n]
    
def requiredError(a, b, u10, u20, p, e):
    n = 1
    prevErr = -1
    h=0
    
    while True:
        u1_n, u2_n = runge_kutta_method(a, b, u10, u20, n)
        u1_2n, u2_2n = runge_kutta_method(a, b, u10, u20, 2*n)
        h = (b-a)/2*n
        err1 = calculate_error(u1_n, u1_2n, p)
        err2 = calculate_error(u2_n, u2_2n, p)
        err = max(err1, err2)
        if err <= e:
            return err
        if ((prevErr - err) == 0):
            print('процесс решения прекращен, т.к. с уменьшением шага погрешность не уменьшается')
            return err
        if (h == 0):
            print('процесс решения прекращен, т.к. значение шага стало недопустимо малым')
            return err
        prevErr = err
        n *= 2
    
    
a = 0
b = 2
u10 = 0
u20 = -2
e = 0.00001

print('погрешность: ',requiredError(a, b, u10, u20, 4, e))
