import numpy as np
import matplotlib.pyplot as plt
import math

# Part a
x_guess = 0.5
x_i = 0
y_i = 0.0
step = 0
flag = 1
X1 = []
Y1 = []

while flag ==1 and step <50:
    x_i = 1.0/math.sin(x_guess)

    if abs((x_i - x_guess)/ x_i) > 1e-5:
        step = step + 1
        x_guess = x_i
        X1.append(step)
        Y1.append(x_i)
    else:
        flag = 0


# Part b

x_guess = 0.5
x_i = 0
y_i = 0.0
step = 0
flag = 1
X2 = []
Y2 = []

while flag == 1 and step < 50:
    x_i = 1.0/math.tan(x_guess)

    if abs((x_i - x_guess)/ x_i) > 1e-5:
        step = step + 1
        x_guess = x_i
        X2.append(step)
        Y2.append(x_i)
    else:
        flag = 0



# Plotting

figure, axis = plt.subplots(1,2)

axis[0].plot(X1, Y1)
axis[0].set_title(" x = 1/sin(x)")

axis[1].plot(X2, Y2)
axis[1].set_title(" x = 1/tan(x)")

plt.show()
