import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def Input():
    l, d, n = 0.25, 0, 5
    s = l / (2.0 * (n - 1))
    return s, n

def gridGeneration(s, n):
    Xw, Xe, X, dxw, dxe, dx = ([] for i in range(6))

    Xw.append(0.0)

    for i in range(1, n):
        if i == 1:
            Xw.append(Xw[i - 1] + s)
        else:
            Xw.append(Xw[i - 1] + 2 * s)

    Xe.append(s)

    for i in range(1, n):
        if i == n - 1:
            Xe.append(Xe[i - 1] + s)
        else:
            Xe.append(Xe[i - 1] + 2 * s)

    X.append(0.0)

    for i in range(1, n):
         X.append(X[i - 1] + 2 * s)

    dxw.append(0.0)
    for i in range(1, n):
        dxw.append(2 * s)

    for i in range(0, n - 1):
        dxe.append(2 * s)
    dxe.append(0.0)

    for i in range(n):
        if i == 0 or i == n - 1:
            dx.append(s)
        else:
            dx.append(2 * s)

    print("Xw array is ", Xw)
    print("Xe array is ", Xe)
    print("X array is ", X)
    print("dxw array is ", dxw)
    print("dxe array is ", dxe)
    print("dx array is ", dx)

    return Xw, Xe, X, dxw, dxe, dx


import math


def area(X, Xw, n):
    Af = []
    A = []
    Y = []
    Yf = []
    P = []
    Pf = []
    m, c = 1.0, math.sqrt(0.0000000061) / 2.0
    for i in range(n):
        Y.append(m * X[i] + c)
        Yf.append(m * Xw[i] + c)

    Af.append(0.0)
    for i in range(n):
        A.append(Y[i] * Y[i])
        Af.append(Yf[i] * Yf[i])
    Af.append(0.0)
    Pf.append(0.0)
    for i in range(n):
        P.append(8 * Y[i])
        Pf.append(8 * Yf[i])
    Pf.append(0.0)

    print("Area of nodes : ", A)
    print("Area of faces : ", Af)
    print("Perimeter at nodes : ", P)
    print("Perimeter at faces : ", Pf)


def coefficient_matrix(dxw, dxe, dx, Xw, Xe, X, n):
    psi, q, h_inf, t_inf, mf, density, Cp, dt = 1.0, 0.0, 64.0, 20.0, 0.0, 0.0, 0.0, 0.01

    Ap, Af, kf, peri, AE, AW, AQ, AC, CAE, CAEO, CAW, CAWO, CAQ, CAQO, CAC, CACO, AT, ATO, CAP, CAPO = ([] for i in range(20))

    for i in range(n):
        Af.append(0.0000000061)
        Ap.append(0.0000000061)
        peri.append(0.0003125)
        kf.append(50.0)
    Af.append(0.0000000061)
    kf.append(50.0)

    AQ.insert(0, Ap[0] * q * dx[0])
    AC.insert(0, peri[0] * h_inf * dx[0])
    CAE.insert(0, 0.0)
    CAEO.insert(0, 0.0)
    CAW.insert(0, 0.0)
    CAWO.insert(0, 0.0)
    CAQ.insert(0, 0.0)
    CAQO.insert(0, (1.0 - psi) * AQ[0])
    CAC.insert(0, 0.0)
    CACO.insert(0, (1.0 - psi) * AC[0])
    AT.insert(0, 0.0)
    ATO.insert(0, 0.0)
    AE.insert(0, (Af[0] * kf[0]) / dxe[0])
    AW.insert(0, 0.0)

    for i in range(1, n-1):
        AE.insert(i, (Af[i+1] * kf[i+1]) / dxe[i])
        AW.insert(i, (Af[i] * kf[i]) / dxw[i])
        AQ.insert(i, Ap[i] * q * dx[i])
        AC.insert(i, peri[i] * h_inf * dx[i])
        CAE.insert(i, psi * AE[i])
        CAEO.insert(i, (1.0 - psi) * AE[i])
        CAW.insert(i, psi * AW[i])
        CAWO.insert(i, (1.0 - psi) * AW[i])
        CAQ.insert(i, psi * AQ[i])
        CAQO.insert(i, (1.0 - psi) * AQ[i])
        CAC.insert(i, psi * AC[i])
        CACO.insert(i, (1.0 - psi) * AC[i])
        AT.insert(i, 0.0)  #(density * Ap[i] * Cp * dx[i]) / dt)
        ATO.insert(i, 0.0)

    AE.insert(n-1, 0.0)
    AW.insert(n-1, (Af[n-1] * kf[n-1]) / dxw[n-1])
    AQ.insert(n-1, Ap[n-1] * q * dx[n-1])
    AC.insert(n-1, peri[n-1] * h_inf * dx[n-1])
    CAE.insert(n-1, 0.0)
    CAEO.insert(n-1, 0.0)
    CAW.insert(n-1, 0.0)
    CAWO.insert(n-1, 0.0)
    CAQ.insert(n-1, 0.0)
    CAQO.insert(n-1, (1.0 - psi) * AQ[n-1])
    CAC.insert(n-1, 0.0)
    CACO.insert(n-1, (1.0 - psi) * AC[n-1])
    AT.insert(n-1, 0.0)
    ATO.insert(n-1, 0.0)
    CAP.insert(0, 1.0)
    CAPO.insert(0, 0.0)

    for i in range(1, n-1):
        CAP.insert(i, AT[i]+CAW[i]+CAE[i]+CAC[i])
        CAPO.insert(i, ATO[i]-(CAWO[i]+CAEO[i]+CACO[i]))

    CAP.insert(n-1, 1.0)
    CAPO.insert(n-1, 0.0)

    return CAPO, CAWO, CAEO, CAQO, CAQ, CAC, CACO, CAP, CAW, CAE, t_inf

def Source(CAPO, CAWO, CAEO, CAQO, CAQ, CAC, CACO, n, T_inf):

    S = []
    T = []
    S.insert(0, 120.0)
    for i in range(1,n-1):
        S.insert(i, CAC[i]*T_inf)
    S.insert(n-1, 30.0)
    return S

def TDMA(CAE, CAW, CAP, S, n):
    A = []
    B = []
    T = []
    A.insert(0, 0.0)
    B.insert(0, S[0])

    for i in range(1, n):
        A.insert(i, CAE[i]/(CAP[i] - (CAW[i]*A[i-1])))
        B.insert(i, (S[i]+(CAW[i]*B[i-1]))/(CAP[i]-(CAW[i]*A[i-1])))

    for i in range(n):
        T.append(0)

    T[n-1] = B[n-1]

    for i in range(n-2, -1, -1):
        T[i] = A[i]*T[i+1] + B[i]

    print("Ai = ", A)
    print("Bi = ", B)

    return T

def FileGeneration(CAP, CAW, CAE, CAC, CAQ, S, T, Xw):
    df = pd.DataFrame()
    df['Xw'] = Xw
    df['CAP'] = CAP
    df['CAW'] = CAW
    df['CAE'] = CAE
    df['CAC'] = CAC
    df['CAQ'] = CAQ
    df['S'] = S
    df['T'] = T
    df.to_excel('result2.xlsx', index=True)

def plottingfunc(T, n, Xw):
    ypoints = np.array(T)
    xpoints = np.array(Xw)
    plt.plot(xpoints, ypoints, linestyle="solid")
    plt.show()


def main():
    s, n = Input()
    Xw, Xe , X, dxw, dxe, dx = gridGeneration(s, n)

    CAPO, CAWO, CAEO, CAQO, CAQ, CAC, CACO, CAP, CAW, CAE, t_inf = coefficient_matrix(dxw, dxe, dx, Xw, Xe, X, n)

    S = Source(CAPO, CAWO, CAEO, CAQO, CAQ, CAC, CACO, n, t_inf)
    print("Source Vector is", S)
    T = TDMA(CAE, CAW, CAP, S, n)
    print("Temperature distribution at node points ", T)
    FileGeneration(CAP, CAW, CAE, CAC, CAQ, S, T, Xw)
    plottingfunc(T, n, Xw)
    area(X, Xw, n)

if __name__ == "__main__":
    main()

