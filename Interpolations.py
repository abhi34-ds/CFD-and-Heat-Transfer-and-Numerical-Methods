import numpy as np

n = 6

X = np.asarray([250, 300, 350, 400, 450, 500]).astype("float64")
Y = np.asarray([1.003, 1.005, 1.008, 1.013, 1.020, 1.029]).astype("float64")

op = int(input("Enter order of polynomial : "))
xp = int(input("Enter the point to be interpolated : "))

yp = 0.0

# Lagrange Interpolation
for i in range(n):
    p = 1

    for j in range(n):
        if i != j:
            p = p * (xp - X[j])/(X[i]-X[j])

    yp = yp + p * Y[i];
print("Interpolated value as per Lagrange Interpolation is : ", yp)

# Newton Interpolation

s = (xp-X[0])/(X[1]-X[0])
y = 0.0
if op == n-1:
    print("Data Insufficient ")
else:

    if op == 1:
        y = Y[0] + s * (Y[1]-Y[0])
    if op == 2:
        y = Y[0] + s * (Y[1] - Y[0]) + ((s*s - s)/2.0) * (Y[2] - 2 * Y[1] + Y[0])
    if op == 3:
        y = Y[0] + s * (Y[1] - Y[0]) + ((s * s - s) / 2.0) * (Y[2] - 2 * Y[1] + Y[0]) + ((s * s * s-3 * s * s + 2 * s) / 6.0) * (Y[3] - 3 * Y[2] + 3 * Y[1] - Y[0])

    print("interpolated value by Newton Interpolation Rule is :", y)