#!/usr/bin/env python2

import ndpp_data as ndpp
import numpy as np
from xml.dom import minidom
import sys


# Read command line arguments
if len(sys.argv) < 2:
    sys.exit("Usage: validate_library.py <dir of ndpp_lib.xml>")
working_dir = sys.argv[1] + '/'
ndpp_lib_xml = working_dir + 'ndpp_lib.xml'

# Read xsdata and create XML document object
xmldoc = minidom.parse(ndpp_lib_xml)

# Get type of library
filetype = xmldoc.getElementsByTagName('filetype')[0].toxml().replace('<filetype>','').replace('</filetype>','').strip()

# Make list of nuclides
nuc_list = xmldoc.getElementsByTagName('ndpp_table')

for nuc in nuc_list:
    nuc_file = working_dir + nuc.attributes['path'].value
    nuc_data = ndpp.NDPP_lib(nuc_file,filetype)
    print nuc.attributes['alias'].value + ': ', nuc_data.test_scatt_positivity(num_mu_pts = 2001)
