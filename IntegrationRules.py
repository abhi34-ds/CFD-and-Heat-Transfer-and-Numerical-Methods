import numpy as np

n = 11
h = int(input("Enter value of h "))

f = np.array([-1, 2, 23, 86, 215, 434, 767, 1238, 1871, 2690, 3719]).astype('float64')

# Trapezoidal Rule
s_t = f[0]
for i in range(1, n-1):
    s_t += f[i]
s_t = s_t + f[n-1]

print("Integral using Trapezoidal rule is : ", (s_t*h)/2.0)

# Simpsons 1/3 rule

s_s3 = 0.0
for i in range(1, n-2, 2):
    s_s3 = s_s3 + 4*f[i] + 2*f[i+1]

s_s3 = s_s3 + 4*f[n-2] + f[n-1]
print("Integral using Simpsons 1/3 rule is : ", (h*s_s3)/3.0)

# Simpsons 3/8 rule

s_s8 = 0.0
for i in range(1, n-1):
    if i % 3 == 0:
        s_s8 = s_s8 + 2*f[i]
    else:
        s_s8 = s_s8 + 3*f[i]

s_s8 = s_s8 + f[0] + f[n-1]

print("Integral using Simpsons 3/8 rule is : ", (s_s8*h*3)/8.0)

