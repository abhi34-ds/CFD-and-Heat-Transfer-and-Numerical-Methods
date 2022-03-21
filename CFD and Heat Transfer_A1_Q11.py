import numpy as np
import matplotlib.pyplot as plt



s = 2.0 + ( 16.0 / 9.0)
i, j = 3, 4
m, n = 8, 9
Y1 = []
Y2 = []
X = []

t = 16.0 / 9.0
flag = 1
x = 2
while flag == 1:
    t = t * ((m + i)/(n + j))

    if t > 1e-5:
        x = x + 1
        s = s + t
        m = m + i
        n = n + j
        Y1.append(t)
        Y2.append(s)
        X.append(x)

    else:
        flag = 0

plt.plot(X, Y1, color = 'r', label = 'Term')
plt.plot(X, Y2, color = 'g', label = 'Sum')
plt.xticks(X)
plt.xlabel("Number of terms")
plt.legend()
plt.show()