import math
import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 3:
        print 'Need to specify the file containing the function and the output file.'
        exit(1)

with open(sys.argv[1]) as f:
        lines = f.readlines()

x = []
y = []
y1 = []

for i in xrange(len(lines)):
        s = lines[i].split()
        x.append(float(s[0]))
        y.append(float(s[1]))
        xVal = float(s[0])
        expected = np.exp(- xVal * xVal / 10) / np.sqrt(10 * math.pi)
        y1.append(expected)

plt.plot(x, y, 'b', linewidth = 2)
plt.plot(x, y1, 'g--', linewidth = 2)
plt.xlabel('$x$', fontsize = 20)
#plt.ylabel('$y$', fontsize = 15)

plt.axis([-8, 8, 0, 0.2])
plt.gcf().set_size_inches(4.0, 3.0)
plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig(sys.argv[2], format = 'eps', dpi = 1000)

