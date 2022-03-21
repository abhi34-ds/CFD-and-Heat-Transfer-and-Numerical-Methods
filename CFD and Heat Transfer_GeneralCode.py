import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

def Input():
    d = 0
    L = [0.04]
    M = [11]
    density, Cp, TL, TR, dt = 1000.0, 10.0, 60.0, 80.0, 1
    Af, Ap, kf, kp, p_node, t_inf, h_inf = ([] for i in range(7))
    ql, qr = 0.0, 0.0
    q = []
    R = 0.002

    for i in range(M[0]):
        Af.append(math.pi * R * R)
        Ap.append(math.pi * R * R)
        p_node.append(math.pi * 2.0 * R)
        kf.append(5.0)
        q.append(0.25e6)
    Af.append(math.pi * R * R)
    kf.append(5.0)

    h_inf.append(12.5)
    t_inf.append(20.0)
    for i in range(1, M[0]-1):
        h_inf.append(12.5)
        t_inf.append(20.0)
    h_inf.append(12.5)
    t_inf.append(20.0)

    tr, tl, newman_r, newman_l, c_r, c_l = 1, 1, 0, 0, 0, 0

    return L, M, d, Af, Ap, kf, kp, p_node, t_inf, h_inf, density, Cp, q, ql, qr, TL, TR, dt, tr, tl, newman_r, newman_l, c_r, c_l

def gridGeneration(L, M, d):
    Xw, Xe, X, dxw, dxe, dx = ([] for i in range(6))
    s = []

    if d == 0:
        n = M[0]
        Xw.append(0.0)

        s.append(L[0] / (2.0 * (n - 1)))

        for i in range(1, n):
            if i == 1:
                Xw.append(Xw[i - 1] + s[0])
            else:
                Xw.append(Xw[i - 1] + 2 * s[0])

        Xe.append(s[0])

        for i in range(1, n):
            if i == n - 1:
                Xe.append(Xe[i - 1] + s[0])
            else:
                Xe.append(Xe[i - 1] + 2 * s[0])

        X.append(0.0)

        for i in range(1, n):
            X.append(X[i - 1] + 2 * s[0])

        dxw.append(0.0)
        for i in range(1, n):
            dxw.append(2 * s[0])

        for i in range(0, n - 1):
            dxe.append(2 * s[0])
        dxe.append(0.0)

        for i in range(n):
            if i == 0 or i == n - 1:
                dx.append(s[0])
            else:
                dx.append(2 * s[0])

    if d == 1:

        # Step Calculation
        # First Region
        s.insert(0, L[0] / (2 * M[0] - 1))

        # Last Region
        s.append((L[-1] - s[-1])/(2 * (M[-1] - 1)))

        k = 1
        Xw.append(0.0)
        for i in range(d + 1):
            for j in range(M[i]):
                if i == 0 and j == 0:
                    Xw.append(s[i])
                    k = k + 1
                elif i != 0 and j == 0:
                    Xw.append(Xw[k - 1] + s[i - 1] + s[i])
                    k = k + 1
                elif i == d and j == M[i] - 1:
                    break
                else:
                    Xw.append(Xw[k - 1] + 2 * s[i])
                    k = k + 1

        k = 0
        for i in range(d + 1):
            for j in range(M[i]):
                if i == 0 and j == 0:
                    Xe.append(s[i])
                    k = k + 1
                elif i != 0 and j == 0:
                    Xe.append(Xe[k - 1] + s[i - 1] + s[i])
                    k = k + 1
                elif i == d and j == M[i] - 1:
                    Xe.append(0.0)
                    k = k + 1
                else:
                    Xe.append(Xe[k - 1] + 2 * s[i])
                    k = k + 1

        X.append(0.0)
        k = 1
        for i in range(d + 1):
            for j in range(M[i]):
                if i == d and j == M[i] - 1:
                    break
                else:
                    X.append(X[k - 1] + 2 * s[i])
                    k = k + 1

        for i in range(d + 1):
            for j in range(M[i]):
                if i == 0 and j == 0:
                    dxw.append(0.0)
                elif i != 0 and j == 0:
                    dxw.append(s[i - 1] + s[i])
                else:
                    dxw.append(2 * s[i])

        for i in range(d + 1):
            for j in range(M[i]):
                if i != 0 and j == 0:
                    dxe.append(s[i - 1] + s[i])
                elif i == d and j == M[i] - 1:
                    dxe.append(0.0)
                else:
                    dxe.append(2 * s[i])

        for i in range(d + 1):
            for j in range(M[i]):
                if (i == 0 and j == 0) or (j == M[i] - 1 and i == d):
                    dx.append(s[i])
                elif (i != 0 or i != d) and j == 0:
                    dx.append(s[i - 1] + s[i])
                else:
                    dx.append(2 * s[i])

    elif d >= 2:
        # Step Calculation
        # First Region
        s.insert(0, L[0] / (2 * M[0] - 1))

        # Intermediate Regions
        j = 1
        for i in range(d - 1):
            s.append((L[j] - s[j - 1]) / (2 * M[j] - 1))
            j = j + 1

        # Last Region
        s.append((L[-1] - s[-2]) / (2 * (M[-1] - 1)))
        # Step Calculation ends

        # Grid Generation Begins

        k = 1
        Xw.append(0.0)
        for i in range(d + 1):
            for j in range(M[i]):
                if i == 0 and j == 0:
                    Xw.append(s[i])
                    k = k + 1
                elif i != 0 and j == 0:
                    Xw.append(Xw[k - 1] + s[i - 1] + s[i])
                    k = k + 1
                elif i == d and j == M[i] - 1:
                    break
                else:
                    Xw.append(Xw[k - 1] + 2 * s[i])
                    k = k + 1

        k = 0
        for i in range(d + 1):
            for j in range(M[i]):
                if i == 0 and j == 0:
                    Xe.append(s[i])
                    k = k + 1
                elif i != 0 and j == 0:
                    Xe.append(Xe[k - 1] + s[i - 1] + s[i])
                    k = k + 1
                elif i == d and j == M[i] - 1:
                    Xe.append(0.0)
                    k = k + 1
                else:
                    Xe.append(Xe[k - 1] + 2 * s[i])
                    k = k + 1

        X.append(0.0)
        k = 1
        for i in range(d + 1):
            for j in range(M[i]):
                if i == d and j == M[i] - 1:
                    break
                else:
                    X.append(X[k - 1] + 2 * s[i])
                    k = k + 1

        for i in range(d + 1):
            for j in range(M[i]):
                if i == 0 and j == 0:
                    dxw.append(0.0)
                elif i != 0 and j == 0:
                    dxw.append(s[i - 1] + s[i])
                else:
                    dxw.append(2 * s[i])

        for i in range(d + 1):
            for j in range(M[i]):
                if i != 0 and j == 0:
                    dxe.append(s[i - 1] + s[i])
                elif i == d and j == M[i] - 1:
                    dxe.append(0.0)
                else:
                    dxe.append(2 * s[i])

        for i in range(d + 1):
            for j in range(M[i]):
                if (i == 0 and j == 0) or (j == M[i] - 1 and i == d):
                    dx.append(s[i])
                elif (i != 0 or i != d) and j == 0:
                    dx.append(s[i - 1] + s[i])
                else:
                    dx.append(2 * s[i])

    print("Step array is ", s)
    print("Xw array is ", Xw)
    print("Xe array is ", Xe)
    print("X array is ", X)
    print("dxw array is ", dxw)
    print("dxe array is ", dxe)
    print("dx array is ", dx)

    return Xw, Xe, X, dxw, dxe, dx, s, M


def coefficient_matrix(dxw, dxe, dx, Xw, Xe, X, n, Af, Ap, kf, kp, p_node, t_inf, h_inf, density, Cp, q, ql, qr, TL, TR, dt, T_old, tr, tl, newman_r, newman_l, c_r, c_l):
    psi, mf = 1.0, 1.0
    AT, ATO, AE, AW, AQ, AC = ([] for i in range(6))
    CAE, CAEO, CAW, CAWO, CAQ, CAQO, CAC, CACO, CAP, CAPO = ([] for i in range(10))
    S = []
    AT.insert(0, density*Ap[0]*Cp*dx[0]*mf/dt)
    ATO.insert(0, density*Ap[0]*Cp*dx[0]*mf/dt)

    # Boundary Condition at Left Face
    # Temperature BC
    if tl == 1:
        CAE.insert(0, 0.0)
        CAEO.insert(0, 0.0)
        CAW.insert(0, 0.0)
        CAWO.insert(0, 0.0)
        CAQ.insert(0, 0.0)
        CAQO.insert(0, 0.0)
        CAC.insert(0, 0.0)
        CACO.insert(0, 0.0)
        CAP.insert(0, 1.0)
        CAPO.insert(0, 0.0)
        S.insert(0, TL)


    # Heat FLux BC
    if newman_l == 1:
        CAE.insert(0, psi * Af[1] * kf[1] / dxe[0])
        CAEO.insert(0, (1 - psi) * Af[1] * kf[1] / dxe[0])
        CAW.insert(0, 0)
        CAWO.insert(0, 0)
        CAQ.insert(0, psi * Ap[0] * ql * dx[0])
        CAQO.insert(0, (1 - psi) * Ap[0] * ql * dx[0])
        CAC.insert(0, psi * p_node[0] * dx[0] * h_inf[0])
        CACO.insert(0, (1 - psi) * p_node[0] * dx[0] * h_inf[0])
        CAP.insert(0, AT[0] + CAW[0] + CAE[0] + CAC[0])
        CAPO.insert(0, ATO[0] - CAWO[0] - CAEO[0] - CACO[0])
        S.insert(0, CAPO[0] * T_old[0] + CAEO[0] * T_old[1] + CAQ[0] + CAQO[0] + CAC[0] * t_inf + CACO[0] * t_inf)
        S[0] = S[0] + psi * Af[0] * ql + (1 - psi) * Af[0] * ql

    #Convection BC
    if c_l == 1:
        CAE.insert(0, psi * Af[1] * kf[1] / dxe[0])
        CAEO.insert(0, (1 - psi) * Af[1] * kf[1] / dxe[0])
        CAW.insert(0, 0)
        CAWO.insert(0, 0)
        CAQ.insert(0, psi * Ap[0] * ql * dx[0])
        CAQO.insert(0, (1 - psi) * Ap[0] * ql * dx[0])
        CAC.insert(0, psi * p_node[0] * dx[0] * h_inf[0])
        CACO.insert(0, (1 - psi) * p_node[0] * dx[0] * h_inf[0])
        CAP.insert(0, AT[0] + CAW[0] + CAE[0] + CAC[0])
        CAP[0] = CAP[0] + psi * Af[0] * h_inf[0]
        CAPO.insert(0, ATO[0] - CAWO[0] - CAEO[0] - CACO[0])
        CAPO[0] = CAPO[0] - (1 - psi) * Af[0] * h_inf[0]
        S.insert(0, CAPO[0] * T_old[0] + CAEO[0] * T_old[1] + CAQ[0] + CAQO[0] + CAC[0] * t_inf[0] + CACO[0] * t_inf[0])
        S[0] = S[0] + psi * Af[0] * h_inf[0] * t_inf[0] + (1 - psi) * Af[0] * h_inf[0] * t_inf[0]

    # Intermediate Control Volume
    for i in range(1, n-1):
        CAE.insert(i, psi*Af[i+1]*kf[i+1]/dxe[i])
        CAEO.insert(i, (1-psi)*Af[i+1]*kf[i+1]/dxe[i])
        CAW.insert(i, psi*Af[i]*kf[i]/dxw[i])
        CAWO.insert(i, (1 - psi)*Af[i]*kf[i]/dxw[i])
        CAQ.insert(i, psi*Ap[i]*q[i]*dx[i])
        CAQO.insert(i, (1-psi)*Ap[i]*q[i]*dx[i])
        CAC.insert(i, psi*p_node[i]*dx[i]*h_inf[i])
        CACO.insert(i, (1-psi)*p_node[i]*dx[i]*h_inf[i])
        AT.insert(i, density*Ap[i]*Cp*dx[i]*mf/dt)
        ATO.insert(i, density*Ap[i]*Cp*dx[i]*mf/dt)
        CAP.insert(i, AT[i] + CAW[i] + CAE[i] + CAC[i])
        CAPO.insert(i, ATO[i] - CAWO[i] - CAEO[i] - CACO[i])
        S.insert(i, CAPO[i]*T_old[i] + CAWO[i]*T_old[i-1] + CAEO[i]*T_old[i+1] + CAQ[i] + CAQO[i] + CAC[i]*t_inf[i] + CACO[i]*t_inf[i])
    # Intermediate Control Volume Ends

    AT.insert(n-1, density * Ap[n-1] * Cp * dx[n-1] * mf / dt)
    ATO.insert(n-1, density * Ap[n-1] * Cp * dx[n-1] * mf / dt)

    # Right Boundary Control Volume

    # Temporal BC
    if tr == 1:
        CAE.insert(n - 1, 0.0)
        CAEO.insert(n - 1, 0.0)
        CAW.insert(n - 1, 0.0)
        CAWO.insert(n - 1, 0.0)
        CAQ.insert(n - 1, 0.0)
        CAQO.insert(n - 1, 0.0)
        CAC.insert(n - 1, 0.0)
        CACO.insert(n - 1, 0.0)
        CAP.insert(n - 1, 1.0)
        CAPO.insert(n - 1, 0.0)
        S.insert(n-1, TR)

    # Heat Flux BC
    if newman_r == 1:
        CAW.insert(n - 1, psi * Af[n - 1] * kf[n - 1] / dxw[n - 1])
        CAWO.insert(n - 1, (1 - psi) * Af[n - 1] * kf[n - 1] / dxw[n - 1])
        CAE.insert(n - 1, 0)
        CAEO.insert(n - 1, 0)
        CAQ.insert(n - 1, psi * Ap[n - 1] * qr * dx[n - 1])
        CAQO.insert(n - 1, (1 - psi) * Ap[n - 1] * qr * dx[n - 1])
        CAC.insert(n - 1, psi * p_node[n - 1] * dx[n - 1] * h_inf[n-1])
        CACO.insert(n - 1, (1 - psi) * p_node[n - 1] * dx[n - 1] * h_inf[n-1])
        CAP.insert(n - 1, AT[n - 1] + CAW[n - 1] + CAE[n - 1] + CAC[n - 1])
        CAPO.insert(n - 1, ATO[n - 1] - CAWO[n - 1] - CAEO[n - 1] - CACO[n - 1])
        S.insert(n - 1, CAPO[n - 1] * T_old[n - 1] + CAWO[n - 1] * T_old[n - 1] + CAQ[n - 1] + CAQO[n - 1] + CAC[
            n - 1] * t_inf[n-1] + CACO[n - 1] * t_inf[n-1])
        S[n - 1] = S[n - 1] - psi * Af[n] * qr + (1 - psi) * Af[n] * qr

    # Convection Heat Flux BC
    if c_r == 1:
        CAW.insert(n - 1, psi * Af[n - 1] * kf[n - 1] / dxw[n - 1])
        CAWO.insert(n - 1, (1 - psi) * Af[n - 1] * kf[n - 1] / dxw[n - 1])
        CAE.insert(n - 1, 0)
        CAEO.insert(n - 1, 0)
        CAQ.insert(n - 1, psi * Ap[n - 1] * qr * dx[n - 1])
        CAQO.insert(n - 1, (1 - psi) * Ap[n - 1] * qr * dx[n - 1])
        CAC.insert(n - 1, psi * p_node[n - 1] * dx[n - 1] * h_inf[n-1])
        CACO.insert(n - 1, (1 - psi) * p_node[n - 1] * dx[n - 1] * h_inf[n-1])
        CAP.insert(n - 1, AT[n - 1] + CAW[n - 1] + CAE[n - 1] + CAC[n - 1])
        CAP[n - 1] = CAP[n - 1] + psi * Af[n] * h_inf[n-1]
        CAPO.insert(n - 1, ATO[n - 1] - CAWO[n - 1] - CAEO[n - 1] - CACO[n - 1])
        CAPO[n - 1] = CAPO[n - 1] - (1 - psi) * Af[n] * h_inf[n-1]
        S.insert(n - 1, CAPO[n - 1] * T_old[n - 1] + CAWO[n - 1] * T_old[n - 1] + CAQ[n - 1] + CAQO[n - 1] + CAC[
            n - 1] * t_inf[n-1] + CACO[n - 1] * t_inf[n-1])
        S[n - 1] = S[n - 1] + psi * Af[n] * h_inf[n-1] * t_inf[n-1] + (1 - psi) * Af[n] * h_inf[n-1] * t_inf[n-1]

    return CAPO, CAWO, CAEO, CAQO, CAQ, CAC, CACO, CAP, CAW, CAE, S


# Matrix Inversion Function
def TDMA(CAE, CAW, CAP, S, n):
    A = []
    B = []
    T = []
    A.insert(0, CAE[0] / CAP[0])
    B.insert(0, S[0] / CAP[0])

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


def TransientSolver(Xw, Xe, X, dxw, dxe, dx, s, n, l, Af, Ap, kf, kp, p_node, t_inf, h_inf, density, Cp, q, ql, qr, TL, TR, dt, tr, tl, newman_r, newman_l, c_r, c_l):
    time, N, dt, T_inf, T_inf_old = 0.0, 0, 1.0, 20.0, 20.0
    temp_l = TL
    temp_r = TR
    temp_old = []
    temp = []
    s = []
    T = []
    x = []
    time_list = []
    cap, caw, cae, cac, caq = ([] for i in range(5))
    temp_old.insert(0, temp_l)
    temp.insert(0, temp_l)
    for i in range(1, n-1):
        temp_old.insert(i, temp_l+(temp_r - temp_l)*X[i]/l)
        temp.insert(i, temp_l+(temp_r - temp_l)*X[i]/l)
    temp_old.insert(n, temp_r)
    temp.insert(n, temp_r)

    print("Temp_Old[i] = ", temp_old)

    while N <= 10:
        N = N+1
        time = time+dt

        print("n = ", N, "Time = ", time)
        print("Temp_old[i] = ", temp_old)

        CAPO, CAWO, CAEO, CAQO, CAQ, CAC, CACO, CAP, CAW, CAE, S = coefficient_matrix(dxw, dxe, dx, Xw, Xe, X, n, Af, Ap, kf, kp, p_node, t_inf, h_inf, density, Cp, q, ql, qr, TL, TR, dt, temp_old, tr, tl, newman_r, newman_l, c_r, c_l)
        cap, caw, cae, cac, caq = CAP, CAW, CAE, CAC, CAQ
        temp = TDMA(CAE, CAW, CAP, S, n)
        print("Temp[i] = ", temp)
        print("N is : ", N)
        T.append(temp)
        x.append(X)
        time_list.append(time)
        temp_old = temp
        s = S

    plt.xlabel("X-axis")
    plt.ylabel("Y - axis")
    for i in range(N):
        plt.plot(x[i], T[i], label="Time " + str(time_list[i]))

    plt.legend()
    plt.show()

    return cap, caw, cae, cac, caq, s, temp_old

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
    df.to_excel('ResultQ1_candd.xlsx', index=True)

def plottingfunc(T, n):
    ypoints = np.array(T)

    plt.plot(ypoints, linestyle="solid")
    plt.show()


def main():
    L, M, d,  Af, Ap, kf, kp, p_node, t_inf, h_inf, density, Cp, q, ql, qr, TL, TR, dt, tr, tl, newman_r, newman_l, c_r, c_l = Input()
    Xw, Xe, X, dxw, dxe, dx, s, M = gridGeneration(L, M, d)
    n = M[0]
    s = s[0]
    ll = L[0]

    #CAPO, CAWO, CAEO, CAQO, CAQ, CAC, CACO, CAP, CAW, CAE, S = coefficient_matrix(dxw, dxe, dx, Xw, Xe, X, n)
    CAP, CAW, CAE, CAC, CAQ, S, T = TransientSolver(Xw, Xe, X, dxw, dxe, dx, s, n, ll, Af, Ap, kf, kp, p_node, t_inf, h_inf, density, Cp, q, ql, qr, TL, TR, dt, tr, tl, newman_r, newman_l, c_r, c_l)
    #print("Final Temp", T)

    # print("Source Vector is", S)
    # T = TDMA(CAE, CAW, CAP, S, n)
    # print("Temperature distribution at node points taking ", n, "points", T)

    FileGeneration(Xw, Xe, X, dxw, dxe, dx, CAP, CAW, CAE, CAC, CAQ, S, T)
    #plottingfunc(T, n)

if __name__ == "__main__":
    main()

