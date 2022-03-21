import math
def area(X, Xw, n):
    Af = []
    A = []
    Y = []
    Yf = []
    P = []
    Pf = []
    m , c = 1.0, math.sqrt(0.0000000061)/2.0
    for i in range(n):
        Y.append(m*X[i]+c)
        Yf.append(m*Xw[i]+c)

    Af.append(0.0)
    for i in range(n):
        A.append(Y[i]*Y[i])
        Af.append(Yf[i]*Yf[i])
    Af.append(0.0)
    Pf.append(0.0)
    for i in range(n):
        P.append(8*Y[i])
        Pf.append(8*Yf[i])
    Pf.append(0.0)

    print("Area of nodes : ", A)
    print("Area of faces : ", Af)
    print("Perimeter at nodes : ", P)
    print("Perimeter at faces : ", Pf)
