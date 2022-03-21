import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import random


def Input():
    d = 2
    L = [9, 11, 10]
    M = [5, 3, 2]
    A = [10, 10, 10]
    # if problem is linear
    Q = [0, 0, 0]
    P = [1, 1, 1]
    density, Cp, TL, TR, dt = 2000.0, 750.0, 150, 75, 1.0
    t_inf, h_inf, ql, qr = [20] * sum(M), [10] * sum(M), 0, 0

    tr, tl, newman_r, newman_l, c_r, c_l, ins_r, ins_l = 0, 0, 0, 1, 1, 0, 0, 0

    return L, M, d, A, Q, P, density, Cp, TL, TR, dt, t_inf, h_inf, ql, qr, tr, tl, newman_r, newman_l, c_r, c_l, ins_r, ins_l


def gridGeneration(L, M, d, A, T_prev_face, T_prev_node, Q, P):
    Xw, Xe, X, dxw, dxe, dx = ([] for i in range(6))
    s = []
    c = sum(M)
    print("d is", d)
    Ap = []
    Af = []
    kf = []
    kp = []
    q = []
    AKeff = []
    p_node = []

    for i in range(3):
        for j in range(M[i]):
            p_node.append(P[i])

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

        for j in range(M[0]+1):
            if j != M[0]:
                Ap.append(A[0])
                q.append(Q[0])
                kp.append(0.05*(1+0.008*T_prev_node[j]))
            Af.append(A[0])
            kf.append(0.05*(1+0.008*T_prev_face[j]))

    if d == 1:

        # Step Calculation
        # First Region
        s.insert(0, L[0] / (2 * M[0] - 1))

        # Last Region
        s.append((L[-1] - s[-1]) / (2 * (M[-1] - 1)))

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

        pn = 0  # Area , K, Q for non linear
        pf = 0
        for i in range(2):
            for j in range(M[i]):
                Ap.append(A[i])
                q.append(Q[i])
                if i == 0:
                    kp.append(0.05*(1+0.008*T_prev_node[pn]))
                    pn = pn + 1
                    Af.append(A[i])
                    kf.append(0.05*(1+0.008*T_prev_face[pf]))
                    pf = pf + 1
                if i == 1:
                    kp.append(0.04*(1+0.0075*T_prev_node[pn]))
                    pn = pn + 1
                    Af.append(A[i])
                    kf.append(0.04*(1+0.0075*T_prev_face[pf]))
                    pf = pf + 1
        AKeff.append((2 * Ap[4] * Ap[5] * kp[4] * kp[5]) / (Ap[4] * kp[4] + Ap[5] * kp[5]))

        print("Ap", Ap)
        print("q", q)
        print("Af", Af)
        print('kf', kf)

    if d >= 2:
        # Step Calculation
        # First Region
        s.insert(0, L[0] / (2 * M[0] - 1))

        # Intermediate Regions
        j = 1
        for i in range(d - 1):
            s.append((L[j] - s[j - 1]) / (2 * M[j] - 1))
            j = j + 1

        # Last Region
        s.append((L[-1] - s[-1]) / (2 * (M[-1] - 1)))
        # Step Calculation ends

        # Grid Generation Begins

        k = 1
        Xw.append(0.0)
        for i in range(d + 1):  # no of discontinuity
            for j in range(M[i]):  # in a particular continuity
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
                    dxw.append(2 * s[i - 1])
                else:
                    dxw.append(2 * s[i])

        for i in range(d + 1):
            for j in range(M[i]):

                if i == d and j == M[i] - 1:
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

        pn = 0  # Area , K, Q for non linear
        pf = 0
        for i in range(3):      # change it as per no of discontinuity
            for j in range(M[i]):
                Ap.append(A[i])
                q.append(Q[i])
                if i == 0:
                    kp.append(0.05*(1+0.008*T_prev_node[pn]))
                    pn = pn + 1
                    Af.append(A[i])
                    kf.append(0.05*(1+0.008*T_prev_face[pf]))
                    pf = pf + 1
                if i == 1:
                    kp.append(0.04*(1+0.0075*T_prev_node[pn]))
                    pn = pn + 1
                    Af.append(A[i])
                    kf.append(0.04*(1+0.0075*T_prev_face[pf]))
                    pf = pf + 1
                if i == 2:
                    kp.append(0.04*(1+0.0075*T_prev_node[pn]))   # as per continuity
                    Af.append(A[i])
                    kf.append(0.04*(1+0.0075*T_prev_face[pf]))   # as per continuity
                    pf = pf + 1
        AKeff.append((2 * Ap[4] * Ap[5] * kp[4] * kp[5]) / (Ap[4] * kp[4] + Ap[5] * kp[5]))
        AKeff.append((2 * Ap[7] * Ap[8] * kp[7] * kp[8]) / (Ap[7] * kp[7] + Ap[8] * kp[8]))
        # AKeff.append((2 * Ap[9] * Ap[10] * kp[7] * kp[8]) / (Ap[7] * kp[7] + Ap[8] * kp[8]))

        print("Ap", Ap)
        print("q", q)
        print("Af", Af)
        print('kf', kf)


    Akeff = []
    t = 0
    z = 0
    for i in range(c + 1):
        if i == 0 or i == c:
            Akeff.append(Af[z] * kf[z])
            z = z + 1
        elif i == 5 or i == 8: # as per continuity
            Akeff.append(AKeff[t])
            t = t + 1
        else:
            Akeff.append(Af[z] * kf[z])
            z = z + 1

    print("Akeff", Akeff)

    return Xw, Xe, X, dxw, dxe, dx, s, M, Af, kf, kp, Ap, Akeff, q, p_node


def coefficient_matrix(Xw, Xe, X, dxw, dxe, dx, s, M, Af, kf, Ap, kp, Akeff, q, p_node, t_inf, h_inf, density, Cp, ql,
                       qr, TL, TR, dt, tr, tl, newman_r, newman_l, c_r, c_l, ins_r, ins_l):
    psi, mf = 1.0, 0.0
    n = sum(M)
    AT, ATO, AE, AW, AQ, AC = ([] for i in range(6))
    CAE, CAEO, CAW, CAWO, CAQ, CAQO, CAC, CACO, CAP, CAPO = ([] for i in range(10))
    S = []
    AT.insert(0, density * Ap[0] * Cp * dx[0] * mf / dt)
    ATO.insert(0, density * Ap[0] * Cp * dx[0] * mf / dt)

    # Boundary Condition at Left Face
    # Temperature BC
    if tl == 1:
        CAE.insert(0, 0.0)
        CAEO.insert(0, 0.0)
        CAW.insert(0, 0.0)
        CAWO.insert(0, 0.0)
        CAQ.insert(0, psi * Ap[0] * q[0] * dx[0])
        CAQO.insert(0, (1 - psi) * Ap[0] * q[0] * dx[0])
        CAC.insert(0, psi * p_node[0] * dx[0] * h_inf[0])
        CACO.insert(0, (1 - psi) * p_node[0] * dx[0] * h_inf[0])
        CAP.insert(0, AT[0] + CAW[0] + CAE[0] + CAC[0])
        CAPO.insert(0, ATO[0] - CAWO[0] - CAEO[0] - CACO[0])
        S.insert(0, TL)

    # Insulated Left
    if ins_l == 1:
        ql = 0.0
        CAE.insert(0, psi * Akeff[1] / dxe[0])
        CAEO.insert(0, (1 - psi) * Akeff[1] / dxe[0])
        CAW.insert(0, 0)
        CAWO.insert(0, 0)
        CAQ.insert(0, psi * Ap[0] * q[0] * dx[0])
        CAQO.insert(0, (1 - psi) * Ap[0] * q[0] * dx[0])
        CAC.insert(0, psi * p_node[0] * dx[0] * h_inf[0])
        CACO.insert(0, (1 - psi) * p_node[0] * dx[0] * h_inf[0])
        CAP.insert(0, AT[0] + CAW[0] + CAE[0] + CAC[0])
        CAPO.insert(0, ATO[0] - CAWO[0] - CAEO[0] - CACO[0])
        S.insert(0, CAQ[0] + CAQO[0] + CAC[0] * t_inf[0] + CACO[0] * t_inf[0])
        S[0] = S[0] + psi * Af[0] * ql + (1 - psi) * Af[0] * ql

    # Heat FLux BC
    if newman_l == 1:
        CAE.insert(0, psi * Akeff[1] / dxe[0])
        CAEO.insert(0, (1 - psi) * Akeff[1] / dxe[0])
        CAW.insert(0, 0)
        CAWO.insert(0, 0)
        CAQ.insert(0, psi * Ap[0] * q[0] * dx[0])
        CAQO.insert(0, (1 - psi) * Ap[0] * q[0] * dx[0])
        CAC.insert(0, psi * p_node[0] * dx[0] * h_inf[0])
        CACO.insert(0, (1 - psi) * p_node[0] * dx[0] * h_inf[0])
        CAP.insert(0, AT[0] + CAW[0] + CAE[0] + CAC[0])
        CAPO.insert(0, ATO[0] - CAWO[0] - CAEO[0] - CACO[0])
        S.insert(0, CAQ[0] + CAQO[0] + CAC[0] * t_inf[0] + CACO[0] * t_inf[0])
        S[0] = S[0] + psi * Af[0] * ql + (1 - psi) * Af[0] * ql

    # Convection BC
    if c_l == 1:
        CAE.insert(0, psi * Akeff[1] / dxe[0])
        CAEO.insert(0, (1 - psi) * Akeff[1] / dxe[0])
        CAW.insert(0, 0)
        CAWO.insert(0, 0)
        CAQ.insert(0, psi * Ap[0] * q[0] * dx[0])
        CAQO.insert(0, (1 - psi) * Ap[0] * q[0] * dx[0])
        CAC.insert(0, psi * p_node[0] * dx[0] * h_inf[0])
        CACO.insert(0, (1 - psi) * p_node[0] * dx[0] * h_inf[0])
        CAP.insert(0, AT[0] + CAW[0] + CAE[0] + CAC[0])
        CAP[0] = CAP[0] + psi * Af[0] * h_inf[0]
        CAPO.insert(0, ATO[0] - CAWO[0] - CAEO[0] - CACO[0])
        CAPO[0] = CAPO[0] - (1 - psi) * Af[0] * h_inf[0]
        S.insert(0, CAQ[0] + CAQO[0] + CAC[0] * t_inf[0] + CACO[0] * t_inf[0])
        S[0] = S[0] + psi * Af[0] * h_inf[0] * t_inf[0] + (1 - psi) * Af[0] * h_inf[0] * t_inf[0]

    # Intermediate Control Volume
    for i in range(1, n - 1):
        CAE.insert(i, psi * Akeff[i + 1] / dxe[i])
        CAEO.insert(i, (1 - psi) * Akeff[i + 1] / dxe[i])
        CAW.insert(i, psi * Akeff[i] / dxw[i])
        CAWO.insert(i, (1 - psi) * Akeff[i] / dxw[i])
        CAQ.insert(i, psi * Ap[i] * q[i] * dx[i])
        CAQO.insert(i, (1 - psi) * Ap[i] * q[i] * dx[i])
        CAC.insert(i, psi * p_node[i] * dx[i] * h_inf[i])
        CACO.insert(i, (1 - psi) * p_node[i] * dx[i] * h_inf[i])
        AT.insert(i, density * Ap[i] * Cp * dx[i] * mf / dt)
        ATO.insert(i, density * Ap[i] * Cp * dx[i] * mf / dt)
        CAP.insert(i, AT[i] + CAW[i] + CAE[i] + CAC[i])
        CAPO.insert(i, ATO[i] - CAWO[i] - CAEO[i] - CACO[i])
        S.insert(i, CAQ[i] + CAQO[i] + CAC[i] *
                 t_inf[i] + CACO[i] * t_inf[i])
    # Intermediate Control Volume Ends

    AT.insert(n - 1, density * Ap[n - 1] * Cp * dx[n - 1] * mf / dt)
    ATO.insert(n - 1, density * Ap[n - 1] * Cp * dx[n - 1] * mf / dt)

    # Right Boundary Control Volume

    # Temporal BC
    if tr == 1:
        CAE.insert(n - 1, 0.0)
        CAEO.insert(n - 1, 0.0)
        CAW.insert(n - 1, psi * Akeff[n - 1] / dxw[n - 1])
        CAWO.insert(n - 1, (1 - psi) * Akeff[n - 1] / dxw[n - 1])
        CAQ.insert(n - 1, psi * Ap[n - 1] * q[n - 1] * dx[n - 1])
        CAQO.insert(n - 1, (1.0 - psi) * Ap[n - 1] * q[n - 1] * dx[n - 1])
        CAC.insert(n - 1, psi * p_node[n - 1] * dx[n - 1] * h_inf[n - 1])
        CACO.insert(n - 1, (1 - psi) * p_node[n - 1] * dx[n - 1] * h_inf[n - 1])
        CAP.insert(n - 1, AT[n - 1] + CAW[n - 1] + CAE[n - 1] + CAC[n - 1])
        CAPO.insert(n - 1, ATO[n - 1] - CAWO[n - 1] - CAEO[n - 1] - CACO[n - 1])
        S.insert(n - 1, TR)

    if ins_r == 1:
        qr = 0.0
        CAW.insert(n - 1, psi * Akeff[n - 1] / dxw[n - 1])
        CAWO.insert(n - 1, (1 - psi) * Akeff[n - 1] / dxw[n - 1])
        CAE.insert(n - 1, 0)
        CAEO.insert(n - 1, 0)
        CAQ.insert(n - 1, psi * Ap[n - 1] * q[n - 1] * dx[n - 1])
        CAQO.insert(n - 1, (1 - psi) * Ap[n - 1] * q[n - 1] * dx[n - 1])
        CAC.insert(n - 1, psi * p_node[n - 1] * dx[n - 1] * h_inf[n - 1])
        CACO.insert(n - 1, (1 - psi) * p_node[n - 1] * dx[n - 1] * h_inf[n - 1])
        CAP.insert(n - 1, AT[n - 1] + CAW[n - 1] + CAE[n - 1] + CAC[n - 1])
        CAPO.insert(n - 1, ATO[n - 1] - CAWO[n - 1] - CAEO[n - 1] - CACO[n - 1])
        S.insert(n - 1, CAQ[n - 1] + CAQO[n - 1] + CAC[
            n - 1] * t_inf[n - 1] + CACO[n - 1] * t_inf[n - 1])
        S[n - 1] = S[n - 1] - psi * Af[n] * qr + (1 - psi) * Af[n] * qr

    # Heat Flux BC
    if newman_r == 1:
        CAW.insert(n - 1, psi * Akeff[n - 1] / dxw[n - 1])
        CAWO.insert(n - 1, (1 - psi) * Akeff[n - 1] / dxw[n - 1])
        CAE.insert(n - 1, 0)
        CAEO.insert(n - 1, 0)
        CAQ.insert(n - 1, psi * Ap[n - 1] * q[n - 1] * dx[n - 1])
        CAQO.insert(n - 1, (1 - psi) * Ap[n - 1] * q[n - 1] * dx[n - 1])
        CAC.insert(n - 1, psi * p_node[n - 1] * dx[n - 1] * h_inf[n - 1])
        CACO.insert(n - 1, (1 - psi) * p_node[n - 1] * dx[n - 1] * h_inf[n - 1])
        CAP.insert(n - 1, AT[n - 1] + CAW[n - 1] + CAE[n - 1] + CAC[n - 1])
        CAPO.insert(n - 1, ATO[n - 1] - CAWO[n - 1] - CAEO[n - 1] - CACO[n - 1])
        S.insert(n - 1, CAQ[n - 1] + CAQO[n - 1] + CAC[
            n - 1] * t_inf[n - 1] + CACO[n - 1] * t_inf[n - 1])
        S[n - 1] = S[n - 1] - psi * Af[n] * qr + (1 - psi) * Af[n] * qr

    # Convection Heat Flux BC
    if c_r == 1:
        CAW.insert(n - 1, psi * Akeff[n - 1] / dxw[n - 1])
        CAWO.insert(n - 1, (1 - psi) * Akeff[n - 1] / dxw[n - 1])
        CAE.insert(n - 1, 0)
        CAEO.insert(n - 1, 0)
        CAQ.insert(n - 1, psi * Ap[n - 1] * q[n - 1] * dx[n - 1])
        CAQO.insert(n - 1, (1 - psi) * Ap[n - 1] * q[n - 1] * dx[n - 1])
        CAC.insert(n - 1, psi * p_node[n - 1] * dx[n - 1] * h_inf[n - 1])
        CACO.insert(n - 1, (1 - psi) * p_node[n - 1] * dx[n - 1] * h_inf[n - 1])
        CAP.insert(n - 1, AT[n - 1] + CAW[n - 1] + CAE[n - 1] + CAC[n - 1])
        CAP[n - 1] = CAP[n - 1] + psi * Af[n-1] * h_inf[n - 1] #CAP[n - 1] = CAP[n - 1] + psi * Af[n] * h_inf[n - 1]
        CAPO.insert(n - 1, ATO[n - 1] - CAWO[n - 1] - CAEO[n - 1] - CACO[n - 1])
        CAPO[n - 1] = CAPO[n - 1] - (1 - psi) * Af[n-1] * h_inf[n - 1] #CAPO[n - 1] = CAPO[n - 1] - (1 - psi) * Af[n] * h_inf[n - 1]
        S.insert(n - 1,  CAQ[n - 1] + CAQO[n - 1] + CAC[
            n - 1] * t_inf[n - 1] + CACO[n - 1] * t_inf[n - 1])
        S[n - 1] = S[n - 1] + psi * Af[n-1] * h_inf[n - 1] * t_inf[n - 1] + (1 - psi) * Af[n-1] * h_inf[n - 1] * t_inf[
            n - 1]
        #S[n - 1] = S[n - 1] + psi * Af[n] * h_inf[n - 1] * t_inf[n - 1] + (1 - psi) * Af[n] * h_inf[n - 1] * t_inf[n - 1]

    return CAPO, CAWO, CAEO, CAQO, CAQ, CAC, CACO, CAP, CAW, CAE, S, AT, ATO


# Matrix Inversion Function
def TDMA(CAE, CAW, CAP, S, n):
    A = []
    B = []
    T = []
    A.insert(0, CAE[0] / CAP[0])
    B.insert(0, S[0] / CAP[0])

    for i in range(1, n):
        A.insert(i, CAE[i] / (CAP[i] - (CAW[i] * A[i - 1])))
        B.insert(i, (S[i] + (CAW[i] * B[i - 1])) / (CAP[i] - (CAW[i] * A[i - 1])))

    for i in range(n):
        T.append(0)

    T[n - 1] = B[n - 1]

    for i in range(n - 2, -1, -1):
        T[i] = A[i] * T[i + 1] + B[i]

    print("Ai = ", A)
    print("Bi = ", B)
    return T

def FileGeneration(Xw, Xe, X, dxw, dxe, dx, CAP, CAW, CAE, CAC, CAQ, S, T, AT, ATO, CAPO):
    df = pd.DataFrame()
    df['Xw'] = Xw
    df['Xe'] = Xe
    df['X'] = X
    df['dxw'] = dxw
    df['dxe'] = dxe
    df['dx'] = dx
    df['AT'] = AT
    df['ATO'] = ATO
    df['CAP'] = CAP
    df['CAPO'] = CAPO
    df['CAW'] = CAW
    df['CAE'] = CAE
    df['CAC'] = CAC
    df['CAQ'] = CAQ
    df['S'] = S
    df['T'] = T
    df.to_excel('CFD and Heat Transfer_1D Non Linear and Discontinuous and Steady.xlsx', index=True)


def plottingfunc(T, n):
    ypoints = np.array(T)

    plt.plot(ypoints, linestyle="solid")
    plt.show()


def main():
    L, M, d, A, Q, P, density, Cp, TL, TR, dt, t_inf, h_inf, ql, qr, tr, tl, newman_r, newman_l, c_r, c_l, ins_r, ins_l = Input()
    si = sum(M)
    y = 0
    #T_prev_face = random.sample(range(TL,TR), sum(M+1))
    T_prev_face = list(np.random.randint(low = TR, high = TL, size = si + 1))
    #T_prev = random.sample(range(TL,TR), sum(M))
    T_prev = list(np.random.randint(low  = TR, high = TL, size = si))

    print("Initial T at face randomly generated", T_prev_face)
    print("Initial T at points randomly generated", T_prev)

    n = sum(M)
    print("Number of points/ nodes", n)

    while True:
        Xw, Xe, X, dxw, dxe, dx, s, M, Af, kf, kp, Ap, Akeff, q, p_node = gridGeneration(L, M, d, A, T_prev_face, T_prev, Q, P)

        CAPO, CAWO, CAEO, CAQO, CAQ, CAC, CACO, CAP, CAW, CAE, S, AT, ATO = coefficient_matrix(Xw, Xe, X, dxw, dxe, dx, s, M, Af, kf, Ap, kp, Akeff, q, p_node, t_inf, h_inf, density, Cp, ql,
                             qr, TL, TR, dt, tr, tl, newman_r, newman_l, c_r, c_l, ins_r, ins_l)

        T = TDMA(CAE, CAW, CAP, S, n)

        e = []
        for i in range(n):
            e.insert(i, (T[i] - T_prev[i]) / T_prev[i])

        print("Error matrix at ", y+1, " iteration is", e)

        if max(e) < 1e-5:
            print("Final Temperature Profile ", T)
            print("Successful after ", y, " times")
            plottingfunc(T, n)
            break
        else:
            w = 0.5
            y = y + 1
            T_prev = T
            for i in range(n+1):
                if i == 0:
                    T_prev_face[i] = T[i]
                elif i == n:
                    T_prev_face[i] = T[n-1]
                else:
                    T_prev_face[i] = ((T[i-1]+T[i])/2.0)

            print("************************************************************************************************************************************************")

    # print("Final Temp", T)

        # print("Source Vector is", S)
        #
        # print("Temperature distribution at node points taking ", n, "points", T)
    print("Non Linear Iteration is run", y, " times")
    FileGeneration(Xw, Xe, X, dxw, dxe, dx, CAP, CAW, CAE, CAC, CAQ, S, T, AT, ATO, CAPO)



if __name__ == "__main__":
    main()
