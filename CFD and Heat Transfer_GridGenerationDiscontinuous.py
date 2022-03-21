def Input():
    d = 2
    L = [9, 11, 10]
    M = [5, 3, 3]
    A = [10, 5, 70]
    K = [20, 20, 20]
    return L, M, d, A, K

def gridGeneration(L, M, d, A, K):
    Xw, Xe, X, dxw, dxe, dx = ([] for i in range(6))
    s = []

    Ap = []
    Af = []
    kf = []
    kp = []
    AKeff = []

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

        for i in range(d+1):
            for j in range(M[i]):
                Ap.append(A[i])
                kp.append(K[i])
            if i != 0 and i != d:
                for j in range(M[i]-1):
                    Af.append(A[i])
                    kf.append(K[i])
            else:
                for j in range(M[i]):
                    Af.append(A[i])
                    kf.append(K[i])
            if i != d:
                Af.append((2*(A[i]*A[i+1]))/(A[i]+A[i+1]))
                kf.append((2*K[i]*K[i+1])/(K[i]+K[i+1]))
                AKeff.append((2*A[i]*A[i+1]*K[i]*K[i+1])/(A[i]*K[i]+A[i+1]*K[i+1]))


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

        for i in range(d+1):
            for j in range(M[i]):
                Ap.append(A[i])
                kp.append(K[i])
            if i != 0 and i != d:
                for j in range(M[i]-1):
                    Af.append(A[i])
                    kf.append(K[i])
            else:
                for j in range(M[i]):
                    Af.append(A[i])
                    kf.append(K[i])
            if i != d:
                Af.append((2*(A[i]*A[i+1]))/(A[i]+A[i+1]))
                kf.append((2*K[i]*K[i+1])/(K[i]+K[i+1]))
                AKeff.append((2*A[i]*A[i+1]*K[i]*K[i+1])/(A[i]*K[i]+A[i+1]*K[i+1]))

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
        s.append((L[-1] - s[-1]) / (2 * (M[-1] - 1)))
        # Step Calculation ends

        # Grid Generation Begins

        k = 1
        Xw.append(0.0)
        for i in range(d + 1): # no of discontinuity
            for j in range(M[i]): # in a particular continuity
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

        for i in range(d+1):
            for j in range(M[i]):
                Ap.append(A[i])
                kp.append(K[i])
            if i != 0 and i != d:
                for j in range(M[i]-1):
                    Af.append(A[i])
                    kf.append(K[i])
            else:
                for j in range(M[i]):
                    Af.append(A[i])
                    kf.append(K[i])
            if i != d:
                Af.append((2*(A[i]*A[i+1]))/(A[i]+A[i+1]))
                kf.append((2*K[i]*K[i+1])/(K[i]+K[i+1]))
                AKeff.append((2*A[i]*A[i+1]*K[i]*K[i+1])/(A[i]*K[i]+A[i+1]*K[i+1]))

    Akeff = []
    t = 0
    c = sum(M)
    for i in range(c+1):
        if i == 0 or i == c:
            Akeff.append(Af[i] * kf[i])
        elif i == 5 or i == 8:
            Akeff.append(AKeff[t])
            t = t + 1
        else:
            Akeff.append(Af[i] * kf[i])



    print("Step array is ", s)
    print("Xw array is ", Xw)
    print("Xe array is ", Xe)
    print("X array is ", X)
    print("dxw array is ", dxw)
    print("dxe array is ", dxe)
    print("dx array is ", dx)
    print("Length of X is ", len(X))
    print("Length of dxw is ", len(dxw))
    print("Length of dxe is ", len(dxe))
    print("Area at face", Af)
    print("Area at point", Ap)
    print("K at face", kf)
    print("K at point", kp)
    print("Kf length", len(kf))
    print(len(Af))
    print("AK eff", AKeff)
    print("Ak eff", Akeff)


L, M, d, A, K = Input()
print(d)
gridGeneration(L, M, d, A, K)