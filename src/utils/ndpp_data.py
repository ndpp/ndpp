#!/usr/bin/env python2

import struct
import numpy as np
import scipy.special as ss
from xml.dom import minidom

SCATT_TYPE_LEGENDRE = 0
SCATT_TYPE_TABULAR  = 1

class ndpp_table(object):
    def __init__(self, nuc):
        self.alias = nuc.attributes['alias'].value
        self.awr = nuc.attributes['awr'].value
        self.location = nuc.attributes['location'].value
        self.name = nuc.attributes['name'].value
        self.path = nuc.attributes['path'].value
        self.kT = nuc.attributes['temperature'].value
        self.zaid = nuc.attributes['zaid'].value

class lib_xml(object):
    def __init__(self, lib_file):
        # Read xsdata and create XML document object
        xmldoc = minidom.parse(lib_file)

        # Get type of library
        self.filetype = self.get_attrib(xmldoc, "filetype")

        # Get directory of library
        self.directory = self.get_attrib(xmldoc, "directory")

        # Get number of entries
        self.n_nuclides = int(self.get_attrib(xmldoc, "entries"))

        # Get flag if chi is present
        self.chi_present = bool(self.get_attrib(xmldoc, "chi_present"))

        # Get scattering type
        self.scatt_type = int(self.get_attrib(xmldoc, "scatt_type"))

        # Get scattering order
        self.scatt_order = int(self.get_attrib(xmldoc, "scatt_order"))

        # Get printing tolerance
        self.print_tol = float(self.get_attrib(xmldoc, "print_tol"))

        # Get thinning tolerance
        self.thin_tol = float(self.get_attrib(xmldoc, "thin_tol"))

        # Get mu_bins
        self.mu_bins = int(self.get_attrib(xmldoc, "mu_bins"))

        # Get energy_bins
        self.energy_bins = []
        e_bins_tmp = ((self.get_attrib(xmldoc, "energy_bins")).split())
        for i in xrange(len(e_bins_tmp)):
            self.energy_bins.append(float(e_bins_tmp[i]))
        self.energy_bins = np.asarray(self.energy_bins, dtype=np.double)

        # Make list of nuclides XML info
        nuc_list = xmldoc.getElementsByTagName('ndpp_table')

        # Create lists for zaid, filenames, alias, ...
        self.nuc_info = []
        for nuc in nuc_list:
            self.nuc_info.append(ndpp_table(nuc))

    def get_attrib(self, xmldoc, tag):
        start = '<' + tag + '>'
        start.strip()
        end = '</' + tag + '>'
        end.strip()

        return xmldoc.getElementsByTagName(tag.strip())[0].toxml().replace( \
            start,'').replace(end,'').strip()

class scatt_data(object):
    def __init__(self, gmin, gmax, order):
        # Set gmin and gmax, but subtract by one to move to python indexing
        self.gmin = gmin - 1
        self.gmax = gmax - 1
        # Initialize the storage of the outgoing data to order
        self.outgoing = np.zeros(shape = (gmax - gmin + 1, order))

class NDPP_lib(object):
    def __init__(self, filename, filetype):

        self.filetype = filetype.lower()

        if self.filetype == 'hdf5':
            import h5py
            self._f = h5py.File(filename, 'r')
            self._hdf5 = True
        elif self.filetype == 'binary':
            self._f = open(filename, 'rb')
            self._hdf5 = False
        elif self.filetype == 'ascii':
            self._f = open(filename, 'r')
            self._hdf5 = False

        # Set flags for what data  was read
        self._metadata = False
        self._data = False

        # Initialize the metadata
        self.name = ""
        self.kT = 0.0
        self.NG = 0
        self.E_bins = []
        self.nuscatter_present = False
        self.chi_present = False
        self.thin_tol = 0.0
        self.print_tol = 0.0
        self.mu_bins = 0

        # Initialize arrays for Energies, group structures, scatter and chi
        # data.
        self.NE_scatt = 0
        self.Ein_scatt = []
        self.gmin = []
        self.gmax = []
        self.scatter = []
        self.nuscatter = []
        self.NE_chi = 0
        self.Ein_chi = []
        self.chi = []

        # Read all metadata
        self._read_metadata()

        # Read scattering data
        self._read_scatt()

        # Read chi
        if self.chi_present:
            self._read_chi()

    def _read_metadata(self):
        if self._metadata == True:
            return

        # Read nuclide name
        self.name = self._get_string(10, path='name')

        # Read kT value
        self.kT = self._get_double(path='kT')[0]

        # Read Energy Groups
        self.NG = self._get_int(path='NG')[0]

        # Read Group Structure
        self.E_bins = self._get_double(n = self.NG + 1, path='E_bins')

        # Read scattering type
        self.scatt_type = self._get_int(path='scatt_type')[0]

        # Read scatering order
        self.scatt_order = self._get_int(path='scatt_order')[0]

        if self.scatt_type == SCATT_TYPE_LEGENDRE:
            self.scatt_order += 1

        # Get flag for if nuscatter is present
        self.nuscatter_present = bool(self._get_int(path='nuscatter')[0])

        # Get flag for if chi is present
        self.chi_present = bool(self._get_int(path='chi_present')[0])

        # Get thinning tolerance
        self.thin_tol = self._get_double(path='thin_tol')[0]

        # Get mu_bins
        self.mu_bins = self._get_int(path='mu_bins')[0]

    def _read_scatt(self):
        # Get NE_scatt
        self.NE_scatt = self._get_int(path='scatt/NE_scatt')[0]

        # Get Ein
        self.Ein_scatt = np.asarray(self._get_double(n=self.NE_scatt,
                                                     path='scatt/Ein_scatt'),
                                    dtype = np.double)

        # Get scattering data for each Ein
        base = 'scatt/scatt_data/'

        for iE in xrange(self.NE_scatt):
            iE_base = base + str(iE)+ '/'
            gmin = self._get_int(path=iE_base+'gmin')[0]
            gmax = self._get_int(path=iE_base+'gmax')[0]
            self.scatter.append(scatt_data(gmin, gmax, self.scatt_order))
            if (gmin > 0):
                for g in xrange(gmax - gmin + 1):
                    # this base name needs some fixing for Hdf5
                    iE_g_base = iE_base + str(g + gmin + 1)
                    self.scatter[iE].outgoing[g][:] = np.asarray( \
                        self._get_double(n=self.scatt_order, path=iE_g_base))

        if self.nuscatter_present:
            for iE in xrange(self.NE_scatt):
                iE_base = base + str(iE)+ '/'
                gmin = self._get_int(path=iE_base+'nu_gmin')[0]
                gmax = self._get_int(path=iE_base+'nu_gmax')[0]
                self.nuscatter.append(scatt_data(gmin, gmax, self.scatt_order))
                if (gmin > 0):
                    for g in xrange(gmax - gmin + 1):
                        # this base name needs some fixing for Hdf5
                        iE_g_base = iE_base + str(g + gmin + 1)
                        self.nuscatter[iE].outgoing[g][:] = np.asarray( \
                            self._get_double(n=self.scatt_order, path=iE_g_base))

    def _read_chi(self):
        # NOT YET IMPLEMENTED!!
        pass

    def condense_outgoing_scatt(self, nu=False, groups=None):
        if groups is None:
            groups = range(self.NG)
        condensed = np.zeros((self.NE_scatt, self.scatt_order))

        if nu:
            if self.nuscatter_present:
                scatter = self.nuscatter
            else:
                raise ValueError("Error: Nu-Scatter data requested, but none present!")
        else:
            scatter = self.scatter

        for iE in xrange(self.NE_scatt):
            gmin = scatter[iE].gmin
            gmax = scatter[iE].gmax

            if gmin == gmax:
                if gmin in groups:
                    condensed[iE][:] += scatter[iE].outgoing[gmin- gmin][:]
            else:
                for g in xrange(gmin, gmax + 1):
                    if g in groups:
                        condensed[iE][:] += scatter[iE].outgoing[g-gmin][:]
        return condensed

    def expand_scatt(self, outgoing, num_mu_pts = 201, order = None):
        # Outgoing is the set of legendre moments vs incoming energies
        # It is only for one group (or an already condensed set of groups)
        # num_mu_pts is the number of pts to use on the mu variable to set up
        # the functional data.
        # (Perhaps this could be moved to a sympy function in the future instead
        # of discrete pts)

        expanded = np.zeros((self.NE_scatt, num_mu_pts))
        mu = np.linspace(-1.0, 1.0, num_mu_pts)

        # Set maxL equal to whatever is smaller, order, or self.scatt_order
        # Note that this will only be used if scatt_type == SCATT_TYPE_LEGENDRE
        if order != None:
            maxL = min(order, self.scatt_order)
        else:
            maxL = self.scatt_order

        if self.scatt_type == SCATT_TYPE_LEGENDRE:
            for iE in xrange(self.NE_scatt):
                for l in xrange(maxL):
                    expanded[iE][:] = expanded[iE][:] + \
                        (float(l) + 0.5) * ss.eval_legendre(l, mu[:]) * \
                        outgoing[iE][l]

        elif self.scatt_type == SCATT_TYPE_TABULAR:
            pass
        return (expanded, mu)

    def test_scatt_positivity(self, nu=False, num_mu_pts = 201, order = None):
        # This function will pass through each Ein and outgoing group and
        # will return a flag if the data set is negative and will also return
        # a list of negative iE and groups

        if nu:
            if self.nuscatter_present:
                scatter = self.nuscatter[iE]
            else:
                raise ValueError("Error: Nu-Scatter data requested, but none present!")
        else:
            scatter = self.scatter[iE]

        # Initialize the return values
        positivity = True
        negativity_list = []
        min_value = 1.0E50

        # Set maxL equal to whatever is smaller, order, or self.scatt_order
        # Note that this will only be used if scatt_type == SCATT_TYPE_LEGENDRE
        if order != None:
            maxL = min(order, self.scatt_order)
        else:
            maxL = self.scatt_order

        mu = np.linspace(-1.0, 1.0, num_mu_pts)

        if self.scatt_type == SCATT_TYPE_LEGENDRE:
          for iE in xrange(self.NE_scatt):
              gmin = scatter[iE].gmin
              gmax = scatter[iE].gmax

              for g in xrange(gmin, gmax + 1):
                  expanded = np.zeros(num_mu_pts)
                  for l in xrange(maxL):
                    expanded[:] = expanded[:] + \
                        (float(l) + 0.5) * ss.eval_legendre(l, mu[:]) * \
                        scatter[iE].outgoing[g - gmin][l]

                  minval = min(expanded)
                  if minval < min_value:
                      min_value = minval
                  if minval < 0.0:
                      positivity = False
                      negativity_list.append((iE, g + gmin))
        elif self.scatt_type == SCATT_TYPE_TABULAR:
            pass

        return (positivity, negativity_list, min_value)

    def _get_data(self, n, typeCode, size):
        return list(struct.unpack('={0}{1}'.format(n,typeCode),
                                  self._f.read(n*size)))

    def _get_int(self, n=1, path=None):
        if self._hdf5:
            return [int(v) for v in self._f[path].value]
        else:
            return [int(v) for v in self._get_data(n, 'i', 4)]

    def _get_long(self, n=1, path=None):
        if self._hdf5:
            return [long(v) for v in self._f[path].value]
        else:
            return [long(v) for v in self._get_data(n, 'q', 8)]

    def _get_float(self, n=1, path=None):
        if self._hdf5:
            return [float(v) for v in self._f[path].value]
        else:
            return [float(v) for v in self._get_data(n, 'f', 4)]

    def _get_double(self, n=1, path=None):
        if self._hdf5:
            return [float(v) for v in self._f[path].value]
        else:
            return [float(v) for v in self._get_data(n, 'd', 8)]

    def _get_string(self, n=1, path=None):
        if self._hdf5:
            return str(self._f[path].value)
        else:
            return str(self._get_data(n, 's', 1)[0])
