import numpy as np
import matplotlib.pyplot as plt


y_guess = 0.0
x_i = 0.0
y_i = 0.0
y_t = 0.0
X1 = []
Y1 = []
flag = 1
step = 0

# Part a
while flag == 1 and step < 50:
    x_i = y_guess
    y_i = (3-x_i)*0.5

    if abs((y_i - y_t) / y_i) > 1e-9:
        step = step + 1
        y_guess = y_i
        y_t = y_i
        X1.append(step)
        Y1.append(y_i)
    else:
        flag = 0

# Part b

X2 = []
Y2 = []
y_guess = 0.0
y_t = 0.0
flag = 1

x_i = 0.0
y_i = 0.0
step = 0

while flag == 1 and step < 50:
    x_i = 3 - 2*y_guess
    y_i = x_i

    if abs((y_i - y_t) / y_i) > 1e-9:
        step = step + 1
        y_guess = y_i
        y_t = y_i
        X2.append(step)
        Y2.append(y_i)
    else:
        flag = 0

figure, axis = plt.subplots(1, 2)

axis[0].plot(X1, Y1)
axis[0].set_title("Part A")
axis[0].set_xlabel('x-axis')
axis[0].set_ylabel('y-axis')

axis[1].plot(X2, Y2)
axis[1].set_title("Part B")
axis[1].set_xlabel('x-axis')
axis[1].set_ylabel('y-axis')

print(Y1, "\n")
print(Y2)

plt.show()


