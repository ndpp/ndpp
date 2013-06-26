#!/usr/bin/env python2

import ndpp_data as ndpp
import numpy as np
import scipy.special as ss

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

h1 = ndpp.NDPP_lib('./1001.70c.g2','binary')

bins = 21
mu = np.linspace(-1.0, 1.0, bins)

one_group = np.zeros((h1.NE_scatt, h1.scatt_order))
data = np.zeros((h1.NE_scatt, bins))

for iE in xrange(h1.NE_scatt):
    if h1.scatter[iE].gmin == h1.scatter[iE].gmax:
        one_group[iE][:] += h1.scatter[iE].outgoing[h1.scatter[iE].gmin][:]        
    else:
        for g in xrange(h1.scatter[iE].gmin, h1.scatter[iE].gmax + 1):
            one_group[iE][:] += h1.scatter[iE].outgoing[g][:]
    # Convert to a functional expansion
    for l in xrange(h1.scatt_order):
        data[iE][:] = data[iE][:] + (float(l) + 0.5) * ss.eval_legendre(l, mu[:]) * one_group[iE][l]
    print data[iE][:]

logEin = np.log(h1.Ein_scatt[:])



X, Y = np.meshgrid(h1.Ein_scatt, mu)
Z = data.reshape(np.shape(X))
print np.shape(X), np.shape(Y), np.shape(Z)
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, antialiased = True)
plt.show()