#!/usr/bin/env python2

import ndpp_data as ndpp
import numpy as np
import scipy.special as ss

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

#~ h1 = ndpp.NDPP_lib('./1001.70c.g1','binary')
#~ h1 = ndpp.NDPP_lib('./1002.84c.g1','binary')
#~ h1 = ndpp.NDPP_lib('./8016.70c.g1','binary')
h1 = ndpp.NDPP_lib('./92235.70c.g1','binary')

bins = 201
mu = np.linspace(-1.0, 1.0, bins)

one_group = h1.condense_outgoing([0])
data = h1.expand(one_group, bins)


lethargy = np.log(h1.Ein_scatt[h1.NE_scatt - 1] / h1.Ein_scatt[:])
print lethargy
from mayavi import mlab
s = mlab.surf(mu, lethargy[::-1], data.transpose(), representation = 'wireframe')
mlab.axes(ranges = [-1, 1, 0.0, lethargy[0], -1, 3], xlabel= 'mu' ,ylabel = 'Energy', 
          zlabel='Distribution')
mlab.show()
