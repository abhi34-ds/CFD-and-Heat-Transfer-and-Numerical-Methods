import numpy as np
import matplotlib.pyplot as plt
import math

x = math.atan(1)

y, y_prev = 0.0, 0.0
flag = 1
n = 1

while flag == 1:
    y = y + (2 * ((-1) ** (n + 1)) * math.sin(n * x)) / n

    print("n ", n, " y ", y)
    if (abs((y - y_prev)) / y) * 100 > 5e-5:
        n = n + 1
        y_prev = y
    else:
        flag = 0


print("Value of f(x) at x = pi/4 accurate up to 5 decimal places is ", y)
print("Value of f(x) = y = x ", y)
print("Value of f(x) ", x)


