import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys
import scipy.interpolate

if len(sys.argv) < 4:
        print 'Need to specify the 2D distribution file, a file with contour levels, and the output file.'
        exit(1)

with open(sys.argv[1]) as f:
        lines = f.readlines()

levels = []

with open(sys.argv[2]) as f1:
        lines1 = f1.readlines()

for l1 in lines1:
        levels.append(float(l1.rstrip()))

x = []
y = []

s = lines[0].split()
for ss in s:
        x.append(float(ss))

s = lines[1].split()
for ss in s:
        y.append(float(ss))

z = np.zeros((len(x), len(y)))

if len(lines) - 2 != len(x):
        print 'Disaster'
        exit(1)

for i in xrange(2, len(lines)):
        s = lines[i].split()
        if len(s) != len(y):
                print 'Disaster'
                exit(1)
        for j in xrange(len(s)):
                z[i - 2, j] = float(s[j])

X, Y = np.meshgrid(x, y)

plt.figure()

#plot the distribution itself
#plt.pcolor(X, Y, z)
#plt.colorbar()

#plot the contours
plt.contour(X, Y, z, levels)

#put labels on the contours
#plt.clabel(CS, inline=1, fontsize=10)

plt.xlabel('$x$', fontsize = 15)
plt.ylabel('$y$', fontsize = 15)

#plt.show()

plt.savefig(sys.argv[3], format = 'eps', dpi = 1000)
