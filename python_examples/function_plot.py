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

for i in xrange(len(lines)):
        s = lines[i].split()
        x.append(float(s[0]))
        y.append(float(s[1]))

plt.plot(x, y)
plt.xlabel('$x$', fontsize = 15)
#plt.ylabel('$y$', fontsize = 15)

plt.savefig(sys.argv[2], format = 'eps', dpi = 1000)

