import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

def Input():
    l, d, n = 0.1, 0, 11
    s = l / (2.0 * (n - 1))
    tr, tl, newman_r, newman_l, c_r, c_l, ins_r, ins_l = 1, 0, 0, 1, 0, 0, 0, 0
    return s, n, tr, tl, newman_r, newman_l, c_r, c_l, ins_r, ins_l

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
    print("****** For n =", n, "nodes ******")
    print("Xw array is ", Xw)
    print("Xe array is ", Xe)
    print("X array is ", X)
    print("dxw array is ", dxw)
    print("dxe array is ", dxe)
    print("dx array is ", dx)

    return Xw, Xe, X, dxw, dxe, dx


def coefficient_matrix(dxw, dxe, dx, Xw, Xe, X, n, tr, tl, newman_r, newman_l, c_r, c_l, ins_r, ins_l):
    psi, mf, density, Cp, dt = 1.0, 0.0, 2000, 1000.0, 0.01
    ql = -1e4
    Ap, Af, kf, p_node, AE, AW, AQ, AC, CAE, CAEO, CAW, CAWO, CAQ, CAQO, CAC, CACO, AT, ATO, CAP, CAPO = ([] for i in range(20))
    R = 0.002
    TR = 20.0

    q = []
    S = []
    h_inf = []
    t_inf = []
    for i in range(n):
        Af.append(1.0)
        Ap.append(1.0)
        p_node.append(1.0)
        h_inf.append(0.0)
        t_inf.append(0.0)
        q.append(4e5)
        kf.append(10.0)
    Af.append(1.0)
    kf.append(10.0)

    AT.insert(0, density * Ap[0] * Cp * dx[0] * mf / dt)
    ATO.insert(0, density * Ap[0] * Cp * dx[0] * mf / dt)

    # Boundary Condition at Left Face
    # Temperature BC


    # Insulated Left

    # Heat FLux BC
    if newman_l == 1:
        CAE.insert(0, psi * Af[1] * kf[1] / dxe[0])
        CAEO.insert(0, (1 - psi) * Af[1] * kf[1] / dxe[0])
        CAW.insert(0, 0)
        CAWO.insert(0, 0)
        CAQ.insert(0, psi * Ap[0] * q[0] * dx[0])
        CAQO.insert(0, (1 - psi) * Ap[0] * q[0] * dx[0])
        CAC.insert(0, psi * p_node[0] * dx[0] * h_inf[0])
        CACO.insert(0, (1 - psi) * p_node[0] * dx[0] * h_inf[0])
        CAP.insert(0, AT[0] + CAW[0] + CAE[0] + CAC[0])
        CAPO.insert(0, ATO[0] - CAWO[0] - CAEO[0] - CACO[0])
        S.insert(0,  CAQ[0])
        S[0] = S[0] + psi * Af[0] * ql + (1 - psi) * Af[0] * ql

    # Convection BC


    # Intermediate Control Volume
    for i in range(1, n - 1):
        CAE.insert(i, psi * Af[i + 1] * kf[i + 1] / dxe[i])
        CAEO.insert(i, (1 - psi) * Af[i + 1] * kf[i + 1] / dxe[i])
        CAW.insert(i, psi * Af[i] * kf[i] / dxw[i])
        CAWO.insert(i, (1 - psi) * Af[i] * kf[i] / dxw[i])
        CAQ.insert(i, psi * Ap[i] * q[i] * dx[i])
        CAQO.insert(i, (1 - psi) * Ap[i] * q[i] * dx[i])
        CAC.insert(i, psi * p_node[i] * dx[i] * h_inf[i])
        CACO.insert(i, (1 - psi) * p_node[i] * dx[i] * h_inf[i])
        AT.insert(i, density * Ap[i] * Cp * dx[i] * mf / dt)
        ATO.insert(i, density * Ap[i] * Cp * dx[i] * mf / dt)
        CAP.insert(i, AT[i] + CAW[i] + CAE[i] + CAC[i])
        CAPO.insert(i, ATO[i] - CAWO[i] - CAEO[i] - CACO[i])
        S.insert(i,  CAQ[i])
    # Intermediate Control Volume Ends

    AT.insert(n - 1, density * Ap[n - 1] * Cp * dx[n - 1] * mf / dt)
    ATO.insert(n - 1, density * Ap[n - 1] * Cp * dx[n - 1] * mf / dt)

    # Right Boundary Control Volume

    # Temporal BC
    if tr == 1:
        CAE.insert(n - 1, 0.0)
        CAEO.insert(n - 1, 0.0)
        CAW.insert(n - 1, 0.0)
        CAWO.insert(n - 1, 0.0)
        CAQ.insert(n - 1, psi * Ap[n-1] * q[n-1] * dx[n-1])
        CAQO.insert(n - 1, (1.0 - psi) * Ap[n-1] * q[n-1] * dx[n-1])
        CAC.insert(n - 1, 0.0)
        CACO.insert(n - 1, 0.0)
        CAP.insert(n - 1, 1.0)
        CAPO.insert(n - 1, 0.0)
        S.insert(n - 1, TR)



    # Heat Flux BC


    # Convection Heat Flux BC


    return CAPO, CAWO, CAEO, CAQO, CAQ, CAC, CACO, CAP, CAW, CAE, S

def TDMA(CAE, CAW, CAP, S, n):
    A = []
    B = []
    T = []
    A.insert(0, CAE[0]/CAP[0])
    B.insert(0, S[0]/CAP[0])

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



def FileGeneration(Xw, Xe , X, dxw, dxe, dx, CAP, CAW, CAE, CAC, CAQ, S, T):
    df = pd.DataFrame()
    df['Xw'] = Xw
    df['Xe'] = Xe
    df['X'] = X
    df['dxw'] = dxw
    df['dxe'] = dxe
    df['dx'] = dx
    df['CAP'] = CAP
    df['CAW'] = CAW
    df['CAE'] = CAE
    df['CAC'] = CAC
    df['CAQ'] = CAQ
    df['S'] = S
    df['T'] = T
    df.to_excel('Result MSEA .xlsx', index=True)

def plottingfunc(T, n):
    ypoints = np.array(T)

    plt.plot(ypoints, linestyle="solid")
    plt.show()


def main():
    s, n, tr, tl, newman_r, newman_l, c_r, c_l, ins_r, ins_l = Input()
    Xw, Xe , X, dxw, dxe, dx = gridGeneration(s, n)

    CAPO, CAWO, CAEO, CAQO, CAQ, CAC, CACO, CAP, CAW, CAE, S = coefficient_matrix(dxw, dxe, dx, Xw, Xe, X, n, tr, tl, newman_r, newman_l, c_r, c_l, ins_r, ins_l)


    print("Source Vector is", S)
    T = TDMA(CAE, CAW, CAP, S, n)
    print("Temperature distribution at node points taking ", n, "points", T)
    FileGeneration(Xw, Xe , X, dxw, dxe, dx, CAP, CAW, CAE, CAC, CAQ, S, T)
    plottingfunc(T, n)

    print("QR at x = b = 0.1m is :", (10.0*(T[n-2]-20.0))/(dx[n-1]))

if __name__ == "__main__":
    main()

