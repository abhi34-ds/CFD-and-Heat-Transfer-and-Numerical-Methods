import math

# Part A
print("Part A")

y1 = 0.0
yp1 = 0.0
flag = 1
n = 1

while flag == 1:

    y1 = y1 + (1.0/(n*n))
    if ((y1-yp1)/y1)*100 > 5e-2:
        n = n+1
        yp1 = y1
    else:
        flag = 0
print("Highest value of n for which accuracy up to 2 place will be achieved is ", n)
print("Value of y1 is ", y1)
# Part B
print("Part B")

y2 = 0.0
yp2 = 0.0
flag = 1
n = 1

while flag == 1:

    y2 = y2 + (1.0/(n*n))
    if abs((y2-yp2)/y2)*100 > 5e-5:
        n = n+1
        yp2 = y2
    else:
        flag = 0

print("Highest value of n for which accuracy up to 5 place will be achieved is ", n)
print("Value of y2 is ", y2)

print("\n")

# Part C

print("Part C")
print("|y1-pi^2/6| is ", abs(y1-((math.pi*math.pi)/6.0)))
print("|y2-pi^2/6| is ", abs(y2-((math.pi*math.pi)/6.0)))


