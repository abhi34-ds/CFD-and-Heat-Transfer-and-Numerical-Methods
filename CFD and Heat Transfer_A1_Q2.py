import numpy as np
import matplotlib.pyplot as plt

def f(w):
    x_guess = 1.0
    x_i = 0.0
    x_t = 0.0
    step = 0
    X = []
    Y = []
    flag = 1

    while flag == 1 and step < 50:
        x_i = (w * (4.0 / x_guess) + (1.0-w) * x_guess)

        if abs((x_i - x_t) / x_i) > 1e-5:
            step = step + 1
            x_guess = x_i
            x_t = x_i
            X.append(step)
            Y.append(x_i)
        else:
            flag = 0
    return X, Y

def main():
    X1, Y1 = f(1.0)
    X2, Y2 = f(0.5)
    X3, Y3 = f(0.75)
    X4, Y4 = f(0.25)

    figure, axis = plt.subplots(2, 2)

    axis[0, 0].plot(X1, Y1)
    axis[0, 0].set_title("w=1")
    axis[0, 1].plot(X2, Y2)
    axis[0, 1].set_title("w=0.5")
    axis[1, 0].plot(X3, Y3)
    axis[1, 0].set_title("w=0.75")
    axis[1, 1].plot(X4, Y4)
    axis[1, 1].set_title("w=0.25")

    plt.show()

if __name__ == "__main__":
    main()