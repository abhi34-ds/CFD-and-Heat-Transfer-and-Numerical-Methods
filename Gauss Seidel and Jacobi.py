import numpy as np

def f1(x, y, z):
    return (17-y+2*z)/20.0
def f2(x, y, z):
    return (-18-3*x+z)/20.0
def f3(x, y, z):
    return(25-2*x+3*y)/20.0

def GaussSiedel(x0, y0, z0, e):

    x1 = f1(x0, y0, z0)
    y1 = f2(x1, y0, z0)
    z1 = f3(x1, y1, z0)

    e1 = abs(x0-x1)
    e2 = abs(y0-y1)
    e3 = abs(z0-z1)

    while e1>e and e2>e and e3>e:
        x0 = x1
        y0 = y1
        z0 = z1

        x1 = f1(x0, y0, z0)
        y1 = f2(x1, y0, z0)
        z1 = f3(x1, y1, z0)

        e1 = abs(x0 - x1)
        e2 = abs(y0 - y1)
        e3 = abs(z0 - z1)
    print("Gauss Siedel \n")
    print("Solution: x = ", x1, " y = ", y1, " z = ", z1, "\n")


def GaussJacobi(x0, y0, z0, e):
    x1 = f1(x0, y0, z0)
    y1 = f2(x0, y0, z0)
    z1 = f3(x0, y0, z0)

    e1 = abs(x0 - x1)
    e2 = abs(y0 - y1)
    e3 = abs(z0 - z1)

    while e1 > e and e2 > e and e3 > e:
        x0 = x1
        y0 = y1
        z0 = z1

        x1 = f1(x0, y0, z0)
        y1 = f2(x0, y0, z0)
        z1 = f3(x0, y0, z0)

        e1 = abs(x0 - x1)
        e2 = abs(y0 - y1)
        e3 = abs(z0 - z1)
    print("Gauss Jacobi \n")
    print("Solution: x = ", x1, " y = ", y1, " z = ", z1, "\n")


def main ():
    X0 = 0.0
    Y0 = 0.0
    Z0 = 0.0

    err = float(input("Enter tolerable error"))
    GaussSiedel(X0, Y0, Z0, err)
    GaussJacobi(X0, Y0, Z0, err)

if __name__ == "__main__":
    main()

