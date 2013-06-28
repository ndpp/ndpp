#!/usr/bin/env python2

import ndpp_data as ndpp
import numpy as np
from xml.dom import minidom
import sys

mu_bins = 21

# Read command line arguments
if len(sys.argv) < 2:
    sys.exit("Usage: validate_library.py <dir of ndpp_lib.xml>")
working_dir = sys.argv[1] + '/'
ndpp_lib_xml = working_dir + 'ndpp_lib.xml'

library = ndpp.lib_xml(ndpp_lib_xml)

neg_grps = [[] for i in xrange(library.n_nuclides)]
neg_grps_flag = [True for i in xrange(library.n_nuclides)]
min_value = [None for i in xrange(library.n_nuclides)]

for (i, nuc) in enumerate(library.nuc_info):
    nuc_data = ndpp.NDPP_lib(nuc.path, library.filetype)
    (neg_grps_flag[i], neg_grps_temp, min_value[i]) =  \
        nuc_data.test_scatt_positivity(num_mu_pts = mu_bins)
    neg_grps[i].append(neg_grps_temp)
    print('Testing ' + nuc.alias + ' Library:')
    if neg_grps_flag[i]:
        print('\tLibrary maintains positivity.')
    else:
        print('\tLibrary becomes negative, with minimum value of ' + 
            str(min_value[i]))
        print('\tOffending incoming energies and outgoing groups:')
        for iE in xrange(len(neg_grps_temp)):
            print('\t\tIncoming Energy = ' + str(neg_grps_temp[iE][0]) + 
                  '\tOutgoing Group = ' + str((neg_grps_temp[iE][1]) + 1))
        
