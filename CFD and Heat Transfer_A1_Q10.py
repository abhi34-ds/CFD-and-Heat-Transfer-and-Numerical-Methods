import numpy as np
import matplotlib.pyplot as plt
import math

p = 4.0 * math.atan(1)

X = [0.0, 1.25, 2.5, 3.75, 5.0]
Y = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0]

for i in range(5):
    for j in range(15):
        n = 0
        t, t_prev = 0.0, 0.0
        flag = 1

        while flag == 1:
            t = t + (400.0/p)*(math.sin(((2*n+1)*p*X[i])/10.0)*math.sinh(((2*n+1)*p*Y[j])/10.0))/((2*n+1)*math.sinh(((2*n+1)*p*15.0)/10.0))

            n = n + 1

            if t != 0:
                if abs(t - t_prev) * 100.0 / t > 5e-5:
                    t_prev = t
                else:
                    flag = 0
            else:
                break


        print("T(", X[i], ",", Y[j], ") = ", t)