#!/usr/bin/env python2

import ndpp_data as ndpp
import numpy as np
from mayavi import mlab

h1 = ndpp.NDPP_lib('./1001.70c.g60','binary')
#~ h1 = ndpp.NDPP_lib('./1001.70c.g3','binary')
#~ h1 = ndpp.NDPP_lib('./1002.84c.g1','binary')
#~ h1 = ndpp.NDPP_lib('./8016.70c.g1','binary')
#~ h1 = ndpp.NDPP_lib('./8016.70c.g3','binary')
#~ h1 = ndpp.NDPP_lib('./92235.70c.g3','binary')
#~ h1 = ndpp.NDPP_lib('./92235.70c.g1','binary')
#~ h1 = ndpp.NDPP_lib('./94240.70c.g3','binary')

bins = 21
print h1.test_scatt_positivity(num_mu_pts = bins)
one_group = h1.condense_outgoing_scatt([0, 1, 2])
glow = h1.condense_outgoing_scatt([0])
gmid = h1.condense_outgoing_scatt([1])
ghigh = h1.condense_outgoing_scatt([2])
(data_1g, mu) = h1.expand_scatt(one_group, bins, order = 6)
lethargy = np.log(h1.Ein_scatt[h1.NE_scatt - 1] / h1.Ein_scatt[:])

x,y=np.meshgrid(mu,lethargy)

s = mlab.mesh(x,y, data_1g, colormap = 'Spectral')
mlab.axes(ranges = [-1, 1, 0.0, lethargy[0], -1, 3], xlabel= 'mu' ,ylabel = 'Energy', 
          zlabel='Distribution')
mlab.outline()
mlab.colorbar()
mlab.show()
mlab.clf()

(data_g0, mu) = h1.expand_scatt(glow, bins, order = 11)
(data_g1, mu) = h1.expand_scatt(gmid, bins, order = 11)
(data_g2, mu) = h1.expand_scatt(ghigh, bins, order = 11)

s = mlab.mesh(x,y, data_g2, colormap = 'Spectral')
mlab.axes(ranges = [-1, 1, 0.0, lethargy[0], -1, 3], xlabel= 'mu' ,ylabel = 'Energy', 
          zlabel='Distribution')
mlab.outline()
mlab.colorbar()
mlab.show()

                    

