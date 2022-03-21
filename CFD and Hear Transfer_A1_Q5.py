import numpy as np
import matplotlib.pyplot as plt
import math

# Part A
print("Part A exp(1) \n")

f, sum, sum_prev = 1, 1.0, 1.0
x_i = 1

flag = 1
step = 1
X = []
# Y = []
# Z = []

while flag == 1:
    f = f * (x_i/step)
    sum = sum + f
    if abs((sum - sum_prev)/sum)*100 > 5e-5:
        step = step + 1
        X.append(sum)
        # Y.append(sum-sum_prev)
        # Z.append((sum-sum_prev)/sum)
        sum_prev = sum
    else:
        flag = 0


print(X)
# print(Y)
# print(Z)

print("Value using in - built function accurate upto 5 ddecimal places ", '%.5f'%math.exp(1))
print("Value using Taylor Series", X[-1])

print("Value using Taylor Series", '%.5f'%X[-1])
print("Value using Taylor Series", "{:.5}".format(X[-1]))

print("\n")

# Part B
print("Part B exp(-1) \n")

f, sum, sum_prev = 1, 0.0, 0.0
x_i = -1

flag = 1
step = 2
X = []
# Y = []
# Z = []
f = x_i

while flag == 1:
    f = f * (x_i/step)
    sum = sum + f
    if abs((sum - sum_prev)/sum)*100 > 5e-5:
        step = step + 1
        X.append(sum)
        # Y.append(sum-sum_prev)
        # Z.append((sum-sum_prev)/sum)
        sum_prev = sum
    else:
        flag = 0


print(X)
# print(Y)
# print(Z)

print("Value using in - built function accurate up to 5 decimal places ", '%.5f'%math.exp(-1))
print("Value using Taylor Series", X[-1])

print("Value using Taylor Series", '%.5f'%X[-1])
print("Value using Taylor Series", "{:.5}".format(X[-1]))




