import numpy as np
import sys

n = int(input("Enter number size of array : "))

aug = np.zeros((n, n+1))

X = np.zeros(n)

print("Enter the aug matrix coefficients : ")

for i in range (n):
    for j in range(n+1):
        aug[i][j]=float(input("a["+str(i)+"]["+str(j)+"]= "))

# Gauss Elimination

for i in range(n):
    if aug[i][i]==0.0:
       sys.exit("Mathematical error")

    for j in range(i+1,n):
        r = aug[j][i]/aug[i][i]

        for k in range(n+1):
            aug[j][k] = aug[j][k] - r*aug[i][k]

# Back Substitution

X[n-1] = aug[n-1][n]/aug[n-1][n-1]

for i in range(n-2, -1, -1):
    X[i] = aug[i][n]
    for j in range (i+1, n):
        X[i] = X[i] - aug[i][j]*X[j]
    X[i] = X[i]/aug[i][i]

# Display Solution

print("\n Solution is : ")

for i in range (n):
    print('X%d = %0.2f' % (i, X[i]), end='\t')
