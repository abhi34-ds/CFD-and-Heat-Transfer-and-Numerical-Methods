import numpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

def Input():
    d = 0
    lx, ly = [0.04], [0.08]
    nx, ny = [11], [21]

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
        print("****** For ny =", ny[0], "nodes ******")
        print("Ys array is ", Ys)
        print("Yn array is ", Yn)
        print("Y array is ", Y)
        print("dys array is ", dys)
        print("dyn array is ", dyn)
        print("dy array is ", dy)

    return Xw, Xe, X, dxw, dxe, dx, Ys, Yn, Y, dys, dyn, dy


def coefficient_matrix(Xw, Xe, X, dxw, dxe, dx, Ys, Yn, Y, dys, dyn, dy, lx, ly, nx, ny, d, T_old, TW, TE, TN, TS, qw, qe, qn, qs):
    psi, q, mf, density, Cp, dt = 1.0, 0.25e6, 1.0, 2000, 0.75, 0.1

    AT, ATO, AS, AN, AW, AE, AQ, CAS, CAN, CAW, CAE, CAQ, CAP, CAPO, S = ([] for i in range(15))
    kcf, k= [], []
    AT_ls, ATO_ls, AS_ls, AN_ls, AW_ls, AE_ls, AQ_ls, CAS_ls, CAN_ls, CAW_ls, CAE_ls, CAQ_ls, CAP_ls, CAPO_ls = ([] for i in range(14))
    S_ls, S = [], []

    for i in range(ny[0]):
        for j in range(nx[0]):
            k.append(20.0)
        kcf.append(k)
    tw, te, ts, tn = 1, 1, 1, 1
    newman_w, newman_e, newman_s, newman_n = 0, 0, 0, 0
    c_w, c_e, c_s, c_n = 0, 0, 0, 0

    h_inf_w, h_inf_e, h_inf_s, h_inf_n = 1, 1, 1, 1
    t_inf_w, t_inf_e, t_inf_s, t_inf_n = 1, 1, 1, 1

    # CVs
    for j in range(ny[0]):
        for i in range(nx[0]):
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                    S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])
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
                S.insert(i, CAPO[i] * T_old[i][j] + CAQ[i])

        CAN_ls.insert(j, CAN)
        CAS_ls.insert(j, CAS)
        CAE_ls.insert(j, CAE)
        CAW_ls.insert(j, CAW)
        CAQ_ls.insert(j, CAQ)
        AT_ls.insert(j, AT)
        ATO_ls.insert(j, ATO)
        CAP_ls.insert(j, CAP)
        CAPO_ls.insert(j, CAPO)
        S_ls.insert(j, S)


    return CAN_ls, CAS_ls, CAE_ls, CAW_ls, CAQ_ls, AT_ls, ATO_ls, CAP_ls, CAPO_ls, S_ls


def ADIscheme (CAN_ls, CAS_ls, CAE_ls, CAW_ls, CAQ_ls, AT_ls, ATO_ls, CAP_ls, CAPO_ls, S_ls, nx, ny):

    r_m = 1
    while r_m > 1e-5:
        T, T_ls = [], []

        # j sweep
        A, B = [], []
        for j in range(ny[0]):
            A.insert(0, CAE_ls[0][j] / CAP_ls[0][j])
            B.insert(0, S_ls[0][j] / CAP_ls[0][j])

            for i in range(1, nx[0]):
                A.insert(i, CAE_ls[i][j] / (CAP_ls[i][j] - (CAW_ls[i][j] * A[i - 1])))
                B.insert(i, (S_ls[i][j] + CAN_ls[i][j] + CAS_ls[i][j] + (CAW_ls[i][j] * B[i - 1])) / (
                            CAP_ls[i][j] - (CAW_ls[i][j] * A[i - 1])))

            for i in range(nx[0]):
                T.append(0)

            T[nx[0] - 1] = B[nx[0] - 1]

            for i in range(nx[0] - 2, -1, -1):
                T[i] = A[i] * T[i + 1] + B[i]

            T_ls.append(T)

        print("T_ls after j sweep ", T_ls)

        # i sweep
        Ty = []
        C, D = [], []
        for i in range(nx[0]):
            C.insert(0, CAN_ls[i][0] / CAP_ls[i][0])
            D.insert(0, S_ls[i][0] / CAP_ls[i][0])

            for j in range(1, ny[0]):
                C.insert(j, CAN_ls[i][j] / (CAP_ls[i][j] - (CAS_ls[i][j] * C[i - 1])))
                D.insert(j, (S_ls[i][j] + CAE_ls[i][j] + CAW_ls[i][j] + (CAS_ls[i][j] * D[i - 1])) / (
                        CAP_ls[i][j] - (CAS_ls[i][j] * C[i - 1])))

            for j in range(ny[0]):
                Ty.append(0)

            Ty[ny[0] - 1] = D[ny[0] - 1]

            for j in range(ny[0] - 2, -1, -1):
                Ty[j] = C[j] * Ty[j + 1] + D[j]

            T_ls.insert(i, Ty)

        print("T_ls after i sweep ", T_ls)
        T_ls = numpy.transpose(T_ls)
        print("T_ls after i sweep after transpose ", T_ls)

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
                    r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAN_ls[i][j]*T_ls[i][j+1] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])
                if j == 0 and i != 0 and i != nx[0] - 1:  # Bottom Boundary
                    r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAN_ls[i][j]*T_ls[i][j+1] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])
                if i == nx[0] - 1 and j != 0 and j != ny[0] - 1:  # Right Boundary
                    r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAN_ls[i][j]*T_ls[i][j+1] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])
                if j == ny[0] - 1 and i != 0 and i != nx[0] - 1:  # Top Boundary
                    r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAN_ls[i][j]*T_ls[i][j+1] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])
                if i != 0 and i != nx[0] - 1 and j != 0 and j != ny[0] - 1: # Intermediate CVs
                    r.insert(i, CAP_ls[i][j]*T_ls[i][j] - CAE_ls[i][j]*T_ls[i+1][j] - CAW_ls[i][j]*T_ls[i-1][j] - CAN_ls[i][j]*T_ls[i][j+1] - CAS_ls[i][j]*T_ls[i][j-1] - S_ls[i][j])

            r_ls.insert(j, r)

        a = np.array(r_ls)
        indices = np.where(a == a.max())
        print(a[indices])
        r_m = a[indices]
        r_m = 0.0


    return (T_ls)

def TransientSolver_2D(Xw, Xe, X, dxw, dxe, dx, Ys, Yn, Y, dys, dyn, dy, lx, ly, nx, ny, d):
    T_prev = np.random.randint(1, 5, size=(nx[0], ny[0]))
    time, dt, N = 0.0, 1.0, 0
    x, y = [], []
    time_list, temp_list = [], []
    TW, TE, TN, TS, qw, qe, qn, qs = 10, 10, 10, 10, 20, 20, 20, 20
    CAN_ls1, CAS_ls1, CAE_ls1, CAW_ls1, CAQ_ls1, AT_ls1, ATO_ls1, CAP_ls1, CAPO_ls1, S_ls1 = ([] for i in range(10))

    time_list.append(time)
    temp_list.append(T_prev)
    x.append(X)
    y.append(Y)

    while N <= 100:
        N = N + 1
        time = time + dt

        print("n = ", N, "Time = ", time)
        print("T_prev = ", T_prev)

        CAN_ls, CAS_ls, CAE_ls, CAW_ls, CAQ_ls, AT_ls, ATO_ls, CAP_ls, CAPO_ls, S_ls = coefficient_matrix(Xw, Xe, X, dxw, dxe, dx, Ys, Yn, Y, dys, dyn, dy, lx, ly, nx, ny, d, T_prev, TW, TE, TN, TS, qw, qe, qn, qs)
        T_curr = ADIscheme(CAN_ls, CAS_ls, CAE_ls, CAW_ls, CAQ_ls, AT_ls, ATO_ls, CAP_ls, CAPO_ls, S_ls, nx, ny)
        CAN_ls1, CAS_ls1, CAE_ls1, CAW_ls1, CAQ_ls1, AT_ls1, ATO_ls1, CAP_ls1, CAPO_ls1, S_ls1 = CAN_ls, CAS_ls, CAE_ls, CAW_ls, CAQ_ls, AT_ls, ATO_ls, CAP_ls, CAPO_ls, S_ls
        print("T = ", T_curr)

        if time == 30.0:
            time_list.append(time)
            temp_list.append(T_prev)
            x.append(X)
            y.append(Y)

        T_prev = T_curr

    return CAN_ls1, CAS_ls1, CAE_ls1, CAW_ls1, CAQ_ls1, AT_ls1, ATO_ls1, CAP_ls1, CAPO_ls1, S_ls1

def main():
    lx, ly, nx, ny, d = Input()
    Xw, Xe, X, dxw, dxe, dx, Ys, Yn, Y, dys, dyn, dy = gridGeneration(lx, ly, nx, ny, d)
    CAN_ls1, CAS_ls1, CAE_ls1, CAW_ls1, CAQ_ls1, AT_ls1, ATO_ls1, CAP_ls1, CAPO_ls1, S_ls1 = TransientSolver_2D(Xw, Xe, X, dxw, dxe, dx, Ys, Yn, Y, dys, dyn, dy, lx, ly, nx, ny, d)



if __name__ == "__main__":
    main()