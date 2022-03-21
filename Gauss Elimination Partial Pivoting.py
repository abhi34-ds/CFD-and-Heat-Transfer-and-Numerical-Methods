import numpy as np
import sys

def my_GaussElimination_PartialPivoting(Aug):
    n, m = np.shape(Aug)
    X = np.zeros(n)

    # Forward Elimination

    for k in range(n):
        i_max = k
        v_max = Aug[i_max][k]

        for i in range(k+1, n):
            if abs(Aug[i][k]) > v_max:
                v_max = Aug[i][k]
                i_max = i

        if i_max != k:
            for m in range (n+1):
                temp = Aug[k][m]
                Aug[k][m] = Aug[i_max][m]
                Aug[i_max][m] = temp

        for i in range(k+1, n):
            f = Aug[i][k]/Aug[k][k]

            for j in range(k+1, n+1):
                Aug[i][j] = Aug[i][j] - f*Aug[k][j]
            Aug[i][k] = 0.0

    # Back Substitution

    s = 0.0

    for i in range(n-1, -1, -1):
        X[i] = Aug[i][n]

        for j in range(i+1, n):
            X[i] = X[i] - Aug[i][j]*X[j]

        X[i] = X[i]/Aug[i][i]

    # Print Solution

    for i in range(n):
        print('X%d = %0.2f' % (i, X[i]), end='\t')

    print("\n")

    return X

def main():
    aug = np.array([[1, 2, 4], [1, -1, 1]]).astype('float64')
    X_sol = my_GaussElimination_PartialPivoting(aug)

    print("Final Solution is: \n")

    for i in range(X_sol.size):
        print('X%d = %0.2f' % (i, X_sol[i]), end='\t')

if __name__ == "__main__":
    main()
