import numpy as np
import math

def func(X, A, B):
    y=((A+B)+(B-A)*X)/2.0
    return ((math.log(1+y))/(1+(y*y)))

def main():
    m = int(input("Enter the number of segments : "))
    x = int(input("Enter the number of points : "))

    u = 1.0
    l = 0.0
    a = 0.0

    n = (u-l)/float(m)
    b = a + n
    I = 0.0

    for i in range(m):
        if x == 1:
            w1 = 2.0
            s1 = 0.0
            I = I + ((b - a) * (w1 * func(s1, a, b))) / 2.0

        if x == 2:
            w1 = 1.0
            w2 = 1.0
            s1 = 0.577
            s2 = -0.577
            I = I + ((b - a) * (w1 * func(s1, a, b) + w2 * func(s2, a, b))) / 2.0
        if x == 3:
            w1 = 0.55
            w2 = 0.88
            w3 = 0.55
            s1 = 0.774
            s2 = 0.0
            s3 = -0.774
            I = I + ((b - a) * (w1 * func(s1, a, b) + w2 * func(s2, a, b) + w3 * func(s3, a, b))) / 2.0

        if x == 4:
            w1 = 0.347855
            w2 = 0.652145
            w3 = 0.652145
            w4 = 0.347855
            s1 = 0.8611363
            s2 = 0.3399810
            s3 = -0.3399810
            s4 = -0.8611363
            I = I + ((b - a) * (w1 * func(s1, a, b) + w2 * func(s2, a, b) + w3 * func(s3, a, b) + w4 * func(s4, a, b))) / 2.0

        a = b;
        b = a + n;

    print("Integration using GQ Rule is ", I)

if __name__ == "__main__":
    main()