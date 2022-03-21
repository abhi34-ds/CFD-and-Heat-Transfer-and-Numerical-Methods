import numpy as np
import sys

def my_LU_Decomposition(A, B):
    n = B.size

    X = np.zeros(n)
    L = np.array([[0.0]*n]*n)
    U = np.array([[0.0]*n]*n)

    # Get L and U Matrix
    for i in range(n):
        for j in range(n):
            if j < i:
                L[j][i] = 0.0
            else:
                L[j][i] = A[j][i]
                for k in range(i):
                    L[j][i] = L[j][i] -L[j][k] * U[k][i]

    for j in range(n):
        if j < i:
            U[i][j] = 0.0
        elif j == i:
            U[i][j] = 1.0
        else:
            U[i][j] = A[i][j]/L[i][i]
            for k in range(i):
                U[i][j] = U[i][j] - ((L[i][k]*U[k][j])/L[i][i])


    for i in range(n):
        for j in range(n):
            print(L[i][j]," ")
        print("\n")


    # Forward on L and b Matrix to get d

    d = np.zeros(n)

    for i in range(n-1):
        if L[i][i] == 0.0:
            print("Division by zeros")
        for j in range(i+1, n):
            r = L[j][i]/L[i][i]
            for k in range(n):
                L[j][k] = L[j][k] - r*L[i][k]
            B[j] = B[j] - r * B[i]

    # Back Substitution
    d[n-1] = B[n-1]/L[n-1][n-1]

    sum = 0.0

    for i in range(n-1, -1, -1):
        sum = B[i]

        for j in range(i+1, n):
            sum = sum - L[i][j]*d[j]

        d[i] = sum/L[i][i]

    # Forward on U and d matrix to get X

    X = np.zeros(n)

    for i in range(n - 1):
        if U[i][i] == 0.0:
            print("Division by zeros")
        for j in range(i + 1, n):
            r = U[j][i] / U[i][i]
            for k in range(n ):
                U[j][k] = U[j][k] - r * U[i][k]
            d[j] = d[j] - r * d[i]

    # Back Substitution
    X[n - 1] = d[n - 1] / U[n - 1][n - 1]

    s = 0.0

    for i in range(n - 1, -1, -1):
        s = d[i]

        for j in range(i + 1, n):
            s = s - U[i][j] * X[j]

        X[i] = s / U[i][i]

    return X


def main():
    a = np.array([[2, -1, -2], [-4, 6, 3], [-4, -2, 8]]).astype('float64')
    b = np.array([2, 4, 6]).astype('float64')
    X_sol = my_LU_Decomposition(a, b)

    print("Final Solution is: \n")

    for i in range(X_sol.size):
        print('X%d = %0.2f' % (i, X_sol[i]), end='\t')

if __name__ == "__main__":
    main()
