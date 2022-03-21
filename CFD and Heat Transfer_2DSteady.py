import numpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

def Input():
    d = 0
    lx, ly = [10], [10]
    nx, ny = [3], [3]

    return lx, ly, nx, ny, d


def gridGeneration(lx, ly, nx, ny, d):
    Xw, Xe, X, dxw, dxe, dx = ([] for i in range(6))
    Yn, Ys, Y, dyn, dys, dy = ([] for i in range(6))

    s = []

    if d == 0:
        sx = lx[0] / (2.0 * (nx[0] - 1))
        sy = ly[0] / (2.0 * (ny[0] - 1))

        Xw.append(0.0)

        for i in range(1, nx[0]):
            if i == 1:
                Xw.append(Xw[i - 1] + sx)
            else:
                Xw.append(Xw[i - 1] + 2 * sx)

        Xe.append(sx)

        for i in range(1, nx[0]):
            if i == nx[0] - 1:
                Xe.append(Xe[i - 1] + sx)
            else:
                Xe.append(Xe[i - 1] + 2 * sx)

        X.append(0.0)

        for i in range(1, nx[0]):
            X.append(X[i - 1] + 2 * sx)

        dxw.append(0.0)
        for i in range(1, nx[0]):
            dxw.append(2 * sx)

        for i in range(0, nx[0] - 1):
            dxe.append(2 * sx)
        dxe.append(0.0)

        for i in range(nx[0]):
            if i == 0 or i == nx[0] - 1:
                dx.append(sx)
            else:
                dx.append(2 * sx)

        Ys.append(0.0)

        for i in range(1, ny[0]):
            if i == 1:
                Ys.append(Ys[i - 1] + sy)
            else:
                Ys.append(Ys[i - 1] + 2 * sy)

        Yn.append(sy)

        for i in range(1, ny[0]):
            if i == ny[0] - 1:
                Yn.append(Yn[i - 1] + sy)
            else:
                Yn.append(Yn[i - 1] + 2 * sy)

        Y.append(0.0)

        for i in range(1, ny[0]):
            Y.append(Y[i - 1] + 2 * sy)

        dys.append(0.0)

        for i in range(1, ny[0]):
            dys.append(2 * sy)

        for i in range(0, ny[0] - 1):
            dyn.append(2 * sy)
        dyn.append(0.0)

        for i in range(ny[0]):
            if i == 0 or i == ny[0] - 1:
                dy.append(sy)
            else:
                dy.append(2 * sy)

        print("****** For nx =", nx[0], "nodes ******")
        print("Xw array is ", Xw)
        print("Xe array is ", Xe)
        print("X array is ", X)
        print("dxw array is ", dxw)
        print("dxe array is ", dxe)
        print("dx array is ", dx)
        print("\n")
        print("****** For ny =", ny[0], "nodes ******")
        print("Ys array is ", Ys)
        print("Yn array is ", Yn)
        print("Y array is ", Y)
        print("dys array is ", dys)
        print("dyn array is ", dyn)
        print("dy array is ", dy)
        print("\n")

    return Xw, Xe, X, dxw, dxe, dx, Ys, Yn, Y, dys, dyn, dy


def coefficient_matrix(Xw, Xe, X, dxw, dxe, dx, Ys, Yn, Y, dys, dyn, dy, lx, ly, nx, ny, d, T_old, TW, TE, TN, TS, qw, qe, qn, qs):
    psi, q, mf, density, Cp, dt = 1.0, 0.0, 0.0, 100.0, 10.0, 1.0

    AT, ATO, AS, AN, AW, AE, AQ, CAS, CAN, CAW, CAE, CAQ, CAP, CAPO, S = ([] for i in range(15))
    kcf, k = [], []
    AT_ls, ATO_ls, AS_ls, AN_ls, AW_ls, AE_ls, AQ_ls, CAS_ls, CAN_ls, CAW_ls, CAE_ls, CAQ_ls, CAP_ls, CAPO_ls = ([] for i in range(14))
    S_ls, S = [], []

    for i in range(ny[0]+1):
        for j in range(nx[0]+1):
            k.append(1.0)
        kcf.append(k)
    tw, te, ts, tn = 0, 1, 1, 1
    newman_w, newman_e, newman_s, newman_n = 1, 0, 0, 0
    c_w, c_e, c_s, c_n = 0, 0, 0, 0

    h_inf_w, h_inf_e, h_inf_s, h_inf_n = 0, 0, 0, 0
    t_inf_w, t_inf_e, t_inf_s, t_inf_n = 0, 0, 0, 0

    # CVs
    for j in range(ny[0]): # rows
        AT, ATO, AS, AN, AW, AE, AQ, CAS, CAN, CAW, CAE, CAQ, CAP, CAPO, S = ([] for i in range(15))
        for i in range(nx[0]): # columns
            if i == 0 and j == 0:  # Bottom Left corner
                if tw == 1 and ts == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, (TW+TS)/2.0)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_w == 1 and newman_s == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, 0)
                    CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + qw * dy[j] * 1 + qs * dx[i] * 1
                if c_w == 1 and c_s == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, 0)
                    CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + dy[j] * 1 * h_inf_w + dx[i] * 1 * h_inf_s
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + dy[j] * 1 * h_inf_w * t_inf_w + dx[i] * 1 * h_inf_s * t_inf_s
                if tw == 1 and newman_s == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TW)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if tw == 1 and c_s == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TW)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_w == 1 and ts == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TS)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if c_w == 1 and ts == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TS)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_w == 1 and c_s == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, 0)
                    CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + h_inf_s * dx[i] * 1
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + qw * dy[j] * 1 + h_inf_s * dx[i] * 1 * t_inf_s
                if c_w == 1 and newman_s == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, 0)
                    CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + dy[j] * 1 * h_inf_w
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + qs * dx[i] * 1 + h_inf_w * t_inf_w * dy[j] * 1
            if i == 0 and j == ny[0] - 1:  # Top Left corner
                if tw == 1 and tn == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, (TW+TN)/2.0)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_w == 1 and newman_n == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                    CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + qw * dy[j] * 1 + qs * dx[i] * 1
                if c_w == 1 and c_n == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                    CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + dy[j] * 1 * h_inf_w + dx[i] * 1 * h_inf_n
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + dy[j] * 1 * h_inf_w * t_inf_w + dx[i] * 1 * h_inf_n * t_inf_n
                if tw == 1 and newman_n == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TW)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])

                if tw == 1 and c_n == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TW)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_w == 1 and tn == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TN)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if c_w == 1 and tn == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TN)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_w == 1 and c_n == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                    CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + dx[i] * 1 * h_inf_n
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + h_inf_n * t_inf_n * 1 * dx[i] + qw * 1 * dy[j]
                if c_w == 1 and newman_n == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                    CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + dy[j] * 1 * h_inf_w
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + qn * 1 * dx[i] + h_inf_w * t_inf_w * 1 * dy[j]
            if i == nx[0] - 1 and j == 0:  # Bottom Right Corner
                if te == 1 and ts == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, (TE+TS)/2.0)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_e == 1 and newman_s == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + qe * dy[j] * 1 + qs * dx[i] * 1
                if c_e == 1 and c_s == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + h_inf_e * 1 * dy[j] + h_inf_s * 1 * dx[i]
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + h_inf_e * t_inf_e * dy[j] * 1 + h_inf_s * t_inf_s * dx[i] * 1
                if te == 1 and newman_s == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TE)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if te == 1 and c_s == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TE)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_e == 1 and ts == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TS)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if c_e == 1 and ts == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TS)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_e == 1 and c_s == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + h_inf_s * dx[i] * 1
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] + S[i] + h_inf_s * t_inf_s * dx[i] * 1 + qe * 1 * dy[j]
                if c_e == 1 and newman_s == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + h_inf_e * dy[j] * 1
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + h_inf_e * t_inf_e * dy[j] * 1 + qs * dx[i] * 1
            if i == nx[0] - 1 and j == ny[0] - 1:  # Top Right Corner
                if te == 1 and tn == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, (TE+TN)/2.0)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_e == 1 and newman_n == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                    CAE.insert(i, 0)
                    CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + qe * dy[j] * 1 + qn * dx[i] * 1
                if c_e == 1 and c_n == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                    CAE.insert(i, 0)
                    CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + h_inf_n * dx[i] * 1 + h_inf_e * dy[j] * 1
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + h_inf_n * t_inf_n * dx[i] * 1 + h_inf_e * t_inf_e * dy[j] * 1
                if te == 1 and newman_n == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TE)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if te == 1 and c_n == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TE)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_e == 1 and tn == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TN)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if c_e == 1 and tn == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TN)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_e == 1 and c_n == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                    CAE.insert(i, 0)
                    CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + h_inf_n * dx[i] * 1
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + h_inf_n * dx[i] * t_inf_n * 1 + qe * dy[j] * 1
                if c_e == 1 and newman_n == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                    CAE.insert(i, 0)
                    CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP [i] + h_inf_e * dy[j] * 1
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + h_inf_e * dy[j] * 1 * t_inf_e + qn * dx[i] * 1
            if i == 0 and j != 0 and j != ny[0] - 1:  # Left Boundary
                if tw == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TW)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_w == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                    CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + qw * dy[j] * 1
                if c_w == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                    CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + h_inf_w * 1 * dy[j]
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + h_inf_w * t_inf_w * 1 * dy[j]
            if j == 0 and i != 0 and i != nx[0] - 1:  # Bottom Boundary
                if ts == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TS)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_s == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, 0)
                    CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                    CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + qs * 1 * dx[i]
                if c_s == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, 0)
                    CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                    CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + h_inf_s * dx[i] * 1
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + h_inf_s * t_inf_s * dx[i] * 1
            if i == nx[0] - 1 and j != 0 and j != ny[0] - 1:  # Right Boundary
                if te == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TE)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_e == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                    CAE.insert(i, 0)
                    CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + qe * 1 * dy[j]
                if c_e == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                    CAE.insert(i, 0)
                    CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + h_inf_e * dy[j] * 1
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + h_inf_e * dy[j] * t_inf_e * 1
            if j == ny[0] - 1 and i != 0 and i != nx[0] - 1:  # Top Boundary
                if tn == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, 0)
                    CAE.insert(i, 0)
                    CAW.insert(i, 0)
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, 1)
                    CAPO.insert(i, ATO[i])
                    S.insert(i, TN)  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                if newman_n == 1:
                    CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                    CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                    CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                    CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + qn * 1 * dx[i]
                if c_n == 1:
                    CAN.insert(i, 0)
                    CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                    CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                    CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                    CAQ.insert(i, psi * q * dx[i] * dy[j])
                    AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                    CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                    CAP[i] = CAP[i] + h_inf_n * dx[i] * 1
                    CAPO.insert(i, ATO[i])
                    S.insert(i, CAQ[i])  # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
                    S[i] = S[i] + h_inf_n * t_inf_n * 1 * dx[i]
            if i != 0 and i != nx[0] - 1 and j != 0 and j != ny[0] - 1:
                CAN.insert(i, (psi * kcf[i][j + 1] * dx[i]) / dyn[j])
                CAS.insert(i, (psi * kcf[i][j] * dx[i]) / dys[j])
                CAE.insert(i, (psi * kcf[i + 1][j] * dy[j]) / dxe[i])
                CAW.insert(i, (psi * kcf[i][j] * dy[j]) / dxw[i])
                CAQ.insert(i, psi * q * dx[i] * dy[j])
                AT.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                ATO.insert(i, (density * Cp * mf * dx[i] * dy[j]) / dt)
                CAP.insert(i, AT[i] + CAE[i] + CAW[i] + CAN[i] + CAS[i])
                CAPO.insert(i, ATO[i])
                S.insert(i, CAQ[i]) # S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
        print("CAN", CAN)
        CAN_ls.append( CAN)
        print("CAN_ls", CAN_ls)
        print("CAS", CAS)
        CAS_ls.append( CAS)
        print("CAS_ls", CAS_ls)
        print("CAE", CAE)
        CAE_ls.append( CAE)
        print("CAE_ls", CAE_ls)
        print("CAW", CAW)
        CAW_ls.append( CAW)
        print("CAW_ls", CAW_ls)
        CAQ_ls.append( CAQ)
        AT_ls.append( AT)
        ATO_ls.append( ATO)
        CAP_ls.append( CAP)
        CAPO_ls.append( CAPO)
        S_ls.append( S)

        #CAN_ls.insert(j, CAN)
        #CAS_ls.insert(j, CAS)
        #CAE_ls.insert(j, CAE)
        #CAW_ls.insert(j, CAW)
        #CAQ_ls.insert(j, CAQ)
        #AT_ls.insert(j, AT)
        #ATO_ls.insert(j, ATO)
        #CAP_ls.insert(j, CAP)
        #CAPO_ls.insert(j, CAPO)
        #S_ls.insert(j, S)


    return CAN_ls, CAS_ls, CAE_ls, CAW_ls, CAQ_ls, AT_ls, ATO_ls, CAP_ls, CAPO_ls, S_ls


def ADIscheme (CAN_ls, CAS_ls, CAE_ls, CAW_ls, CAQ_ls, AT_ls, ATO_ls, CAP_ls, CAPO_ls, S_ls, nx, ny, T_old):
    T_old1 = T_old
    r_m = 1
    while r_m > 1e-5:
        T_ls = []
        sx = []
        for j in range(ny[0]):
            sxls = []
            for i in range(nx[0]):
                if j == 0:
                    sxls.insert(i, S_ls[j][i] + CAN_ls[j][i] * T_old1[j][i+1])
                else:
                    sxls.insert(i, S_ls[j][i] + CAN_ls[j][i] * T_old1[j][i+1] + CAS_ls[j][i] * T_old1[j][i-1])
            sx.append(sxls)

        print('***********************SX************************')
        print("Shape of SX array is", np.array(sx).shape)
        for j in range(ny[0] - 1, -1, -1):
            for i in range(nx[0]):
                print(sx[j][i], end=' ')
            print("\n")

        # j sweep

        for j in range(ny[0]):
            A, B, T = [], [], []
            A.insert(0, CAE_ls[0][j] / CAP_ls[0][j])
            B.insert(0, sx[0][j] / CAP_ls[0][j])

            for i in range(1, nx[0]):
                A.insert(i, CAE_ls[i][j] / (CAP_ls[i][j] - (CAW_ls[i][j] * A[i - 1])))
                B.insert(i, (sx[i][j] + (CAW_ls[i][j] * B[i - 1])) / (
                            CAP_ls[i][j] - (CAW_ls[i][j] * A[i - 1])))

            for i in range(nx[0]):
                T.append(0)

            T[nx[0] - 1] = B[nx[0] - 1]

            for i in range(nx[0] - 2, -1, -1):
                T[i] = A[i] * T[i + 1] + B[i]

            T_ls.append(T)
        print('***********************T_ls after J sweep ************************')
        print("Shape of SX array is", np.array(sx).shape)
        for j in range(ny[0] - 1, -1, -1):
            for i in range(nx[0]):
                print(T_ls[j][i], end=' ')
            print("\n")


        sy = []
        for j in range(ny[0]):
            syls = []
            for i in range(nx[0]):
                if i == 0:
                    syls.insert(i, S_ls[i][j] + CAE_ls[i][j] * T_old1[i+1][j])
                else:
                    syls.insert(i, S_ls[i][j] + CAE_ls[i][j] * T_old1[i+1][j] + CAW_ls[i][j] * T_old1[i-1][j])
            sy.append(syls)
        print('***********************SY************************')
        print("Shape of SY array is", np.array(sy).shape)
        for j in range(ny[0] - 1, -1, -1):
            for i in range(nx[0]):
                print(sy[j][i], end=' ')
            print("\n")


        # i sweep

        T_ls = []
        for i in range(nx[0]):
            C, D, Ty = [], [], []
            C.insert(0, CAN_ls[i][0] / CAP_ls[i][0])
            D.insert(0, sy[i][0] / CAP_ls[i][0])

            for j in range(1, ny[0]):
                C.insert(j, CAN_ls[i][j] / (CAP_ls[i][j] - (CAS_ls[i][j] * C[j - 1])))
                D.insert(j, (sy[i][j] + (CAS_ls[i][j] * D[j - 1])) / (
                        CAP_ls[i][j] - (CAS_ls[i][j] * C[j - 1])))

            for j in range(ny[0]):
                Ty.append(0)

            Ty[ny[0] - 1] = D[ny[0] - 1]

            for j in range(ny[0] - 2, -1, -1):
                Ty[j] = C[j] * Ty[j + 1] + D[j]

            T_ls.append(Ty)

        print('***********************T_ls after i sweep************************')
        print("Shape of T_ls array is", np.array(T_ls).shape)
        for j in range(ny[0] - 1, -1, -1):
            for i in range(nx[0]):
                print(T_ls[j][i], end=' ')
            print("\n")

        T_ls = numpy.transpose(T_ls)
        print("T_ls after i sweep after transpose ")
        print(np.matrix(T_ls))

        r, r_ls = [], []

        for j in range(ny[0]):
            r = []
            for i in range(nx[0]):
                if i == 0 and j == 0:  # Bottom Left corner
                    r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAN_ls[i][j]*T_ls[i][j+1] - S_ls[i][j])
                if i == 0 and j == ny[0] - 1:  # Top Left corner
                    r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])
                if i == nx[0] - 1 and j == 0:  # Bottom Right Corner
                    r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAN_ls[i][j]*T_ls[i][j+1] - S_ls[i][j])
                if i == nx[0] - 1 and j == ny[0] - 1:  # Top Right Corner
                    r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])
                if i == 0 and j != 0 and j != ny[0] - 1:  # Left Boundary
                    r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAN_ls[i][j]*T_ls[i][j+1] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])
                    #r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAN_ls[i][j]*T_ls[i][j+1] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])
                if j == 0 and i != 0 and i != nx[0] - 1:  # Bottom Boundary
                    r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAN_ls[i][j]*T_ls[i][j+1] - S_ls[i][j])
                    #r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAN_ls[i][j]*T_ls[i][j+1] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])
                if i == nx[0] - 1 and j != 0 and j != ny[0] - 1:  # Right Boundary
                    #r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAN_ls[i][j]*T_ls[i][j+1] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])
                    r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAN_ls[i][j]*T_ls[i][j+1] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])
                if j == ny[0] - 1 and i != 0 and i != nx[0] - 1:  # Top Boundary
                    r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])
                    #r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAN_ls[i][j]*T_ls[i][j+1] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])
                if i != 0 and i != nx[0] - 1 and j != 0 and j != ny[0] - 1: # Intermediate CVs
                    r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAN_ls[i][j]*T_ls[i][j+1] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])

            r_ls.insert(j, r)
        print("********************* r matrix ************************")
        print(np.matrix(r_ls))
        print("\n")

        a = np.array(r_ls)
        indices = np.where(a == a.max())
        print(a[indices])
        r_m = a[indices]
        r_m = 0
         # just to get out of loop after running for first time - will change in transient code

    return (T_ls)

def FileGeneration(Xw, Xe, X, dxw, dxe, dx, Ys, Yn, Y, dys, dyn, dy, CAN_ls, CAS_ls, CAE_ls, CAW_ls, CAQ_ls, AT_ls, ATO_ls, CAP_ls, CAPO_ls, S_ls):
    df = pd.DataFrame()
    df['Xw'] = Xw
    df['Xe'] = Xe
    df['X'] = X
    df['dxw'] = dxw
    df['dxe'] = dxe
    df['dx'] = dx

    df1 = pd.DataFrame()
    df1['Ys'] = Ys
    df1['Yn'] = Yn
    df1['Y'] = Y
    df1['dys'] = dys
    df1['dyn'] = dyn
    df1['dy'] = dy

    #df2 = pd.DataFrame(CAN_ls, columns =["CAN_ls"]*21)



    with pd.ExcelWriter('CFD and Heat Transfer_2D Steady.xlsx') as writer:
        df.to_excel(writer, sheet_name="X grid")
        df1.to_excel(writer, sheet_name="Y grid")
        #df2.to_excel(writer, sheet_name="Coefficient")

def main():
    lx, ly, nx, ny, d = Input()
    Xw, Xe, X, dxw, dxe, dx, Ys, Yn, Y, dys, dyn, dy = gridGeneration(lx, ly, nx, ny, d)
    TW, TE, TN, TS, qw, qe, qn, qs = 0.0, 200.0, 100.0, 200.0, 10.0, 0.0, 0.0, 0.0

    T_old = []
    for j in range(ny[0]+1):
        T_oldls = []
        for i in range(nx[0]+1):
            T_oldls.insert(i, 200.0)
        T_old.append(T_oldls)

    CAN_ls, CAS_ls, CAE_ls, CAW_ls, CAQ_ls, AT_ls, ATO_ls, CAP_ls, CAPO_ls, S_ls = coefficient_matrix(Xw, Xe, X, dxw, dxe, dx, Ys, Yn, Y, dys, dyn, dy, lx, ly, nx, ny, d, T_old, TW, TE, TN, TS, qw, qe, qn, qs)
    print('***********************CAN_ls ************************')
    print("Shape of CAN_ls array is", np.array(CAN_ls).shape)
    for j in range(ny[0]-1, -1, -1):
        for i in range(nx[0]):
            print(CAN_ls[j][i], end = ' ')
        print("\n")

    print('***********************CAS_ls ************************')
    print("Shape of CAS_ls array is", np.array(CAS_ls).shape)
    for j in range(ny[0]-1, -1, -1):
        for i in range(nx[0]):
            print(CAS_ls[j][i], end = ' ')
        print("\n")

    print('***********************CAW_ls ************************')
    print("Shape of CAW_ls array is", np.array(CAW_ls).shape)
    for j in range(ny[0]-1, -1, -1):
        for i in range(nx[0]):
            print(CAW_ls[j][i], end = ' ')
        print("\n")

    print('***********************CAE_ls ************************')
    print("Shape of CAE_ls array is", np.array(CAE_ls).shape)
    for j in range(ny[0]-1, -1, -1):
        for i in range(nx[0]):
            print(CAE_ls[j][i], end = ' ')
        print("\n")

    print('***********************CAP_ls ************************')
    print("Shape of CAP_ls array is", np.array(CAP_ls).shape)
    for j in range(ny[0]-1, -1, -1):
        for i in range(nx[0]):
            print(CAP_ls[j][i], end = ' ')
        print("\n")


    print('***********************CAQ_ls ************************')
    print("Shape of CAP_ls array is", np.array(CAQ_ls).shape)
    for j in range(ny[0]-1, -1, -1):
        for i in range(nx[0]):
            print(CAQ_ls[j][i], end = ' ')
        print("\n")

    print('***********************S_ls ************************')
    print("Shape of S_ls array is", np.array(S_ls).shape)
    for j in range(ny[0]-1, -1, -1):
        for i in range(nx[0]):
            print(S_ls[j][i], end = ' ')
        print("\n")

    print("***********************Coefficient Printing Ends***************************")
    print("\n")


    T = ADIscheme(CAN_ls, CAS_ls, CAE_ls, CAW_ls, CAQ_ls, AT_ls, ATO_ls, CAP_ls, CAPO_ls, S_ls, nx, ny, T_old)
    FileGeneration(Xw, Xe, X, dxw, dxe, dx, Ys, Yn, Y, dys, dyn, dy, CAN_ls, CAS_ls, CAE_ls, CAW_ls, CAQ_ls, AT_ls, ATO_ls, CAP_ls, CAPO_ls, S_ls)

    print("\n")
    print("***********************Final Temperature Distribution***************************")
    for j in range(ny[0]-1, -1, -1):
        for i in range(nx[0]):
            print(T[i][j], end = " ")
        print("\n")



if __name__ == "__main__":
    main()