import matplotlib.pyplot as plt
from math import exp, sqrt
data_file = open("result.txt")
y = []
x = []
lines = 0

for line in open("result.txt"):
    lines += 1
    x.append(lines)
for line in data_file:
	line = line.replace('\n', '')
	line = line.replace(',', '')
	y.append(float(line))

plt.plot(x,y)
plt.xlabel("T")
plt.ylabel("function(x)")
plt.show()
