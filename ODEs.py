import numpy as np

def func(x, y):
    return ((( x * x )-( y * y ))/( 2 * x * y))

def Euler(x0, y0, h):
    yn = 0.0
    print("Euler Method \n")

    while x0 <= 10:
        yn = y0 + h * func(x0, y0)
        y0 = yn
        print(" y( ", x0, ")", yn,"\n")
        x0 = x0 + h

def RK4(x0, y0, h):
    I = 0.0
    yn = 0.0
    print("RK4 method \n")

    for i in range(10):
         k1 = h * func(x0, y0)
         k2 = h * func((x0 + h / 2), (y0 + k1 / 2))
         k3 = h * func((x0 + h / 2), (y0 + k2 / 2))
         k4 = h * func((x0 + h), (y0 + k3))
         k = (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
         yn = y0 + k
         print("y( ", x0, ") = ", yn, "\n")
         x0 = x0 + h
         y0 = yn


def main():
    h = int(input("Enter value of h : "))
    X0 = 1.0
    Y0 = 1.0

    Euler(X0, Y0, h)
    RK4(X0, Y0, h)

if __name__ == "__main__":
    main()