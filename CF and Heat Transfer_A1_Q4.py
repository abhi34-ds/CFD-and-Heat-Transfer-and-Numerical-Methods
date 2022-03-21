import numpy as np
import matplotlib.pyplot as plt

eps = 1.0
X = []
while 1+eps != 1:
    eps = eps / 2.0
    X.append(eps)

print("My machine epsilon is -> ", X[-1])
