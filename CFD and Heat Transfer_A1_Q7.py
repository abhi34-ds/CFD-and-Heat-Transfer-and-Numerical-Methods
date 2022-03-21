import numpy as np
import matplotlib.pyplot as plt

m = 2
p = 10

m_i = 10 ** (m - 1)

dummy = 0.0
flag = 1
X = []

while m_i < 10 ** m:

    m_d = m_i
    sum = 0
    while m_d != 0:
        sum = sum + m_d % 10
        m_d = m_d // 10


    if sum % p == 0:
        X.append(m_i)

    m_i = m_i + 1

print("List of integers satisfying the condition ", X)

# Part b

n = len(X)
y1, y2, y3 = 0, 0, 0

for i in range(n):
    for j in range(n):
        y1 = y1 + X[i]*X[j]

        if i != j:
            y2 = y2 + X[i]*X[j]
        if i < j:
            y3 = y3 + X[i]*X[j]

print("y1", " ", "y2", " ", "y3")
print(y1, " ", y2, " ", y3, " ")


