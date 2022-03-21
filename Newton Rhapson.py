import numpy as np

f = lambda x: x**2 - 3*x + 2
fprime = lambda x: 2*x - 3

tol = 1e-7
x0 = 0.0
i = 1
n = 20

while True:
    if abs(f(x0)) < tol:
        print("Root is ", x0)
        break
    if fprime(x0) == 0:
        print("Mathematical error \n")
        break
    else:
        x1 = x0 - f(x0)/fprime(x0)
        print("For iteration ", i, " value of root is ", x1, "\n")
        x0 = x1
    i = i+1
    if i > n:
        print("Not convergent \n")
        break
