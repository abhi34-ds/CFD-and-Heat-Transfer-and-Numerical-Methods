import numpy as np
import matplotlib.pyplot as plt
import math

alpha = 0.01
p = 4.0 * math.atan(1)

X = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
T = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 50.0]

for i in range(11):
    for j in range(10):
        n = 0
        t, t_prev = 0.0, 0.0
        flag = 1

        while flag == 1:
            t = t + ((800 / (p * p)) * ((-1.0) ** n) * math.sin((2 * n + 1) * p * X[i]) * math.exp(
                (-1) * (2 * n + 1) * (2 * n + 1) * p * p * alpha * T[j])) / ((2 * n + 1) * (2 * n + 1))
            n = n + 1

            if t != 0:
                if abs(t - t_prev) * 100.0 / t > 5e-5:
                    t_prev = t
                else:
                    flag = 0
            else:
                break


        print("T(", X[i], ",", T[j], ") = ", t)