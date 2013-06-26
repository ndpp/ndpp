#!/usr/bin/env python2

import struct
import numpy as np

SCATT_TYPE_LEGENDRE = 0
SCATT_TYPE_TABULAR  = 1

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
                                               path='scatt/Ein_scatt'))
        
        # Get scattering data for each Ein  
        base = 'scatt/scatt_data/'
                                             
        for iE in xrange(self.NE_scatt):
            iE_base = base + str(iE)+ '/'
            gmin = self._get_int(path=iE_base+'gmin')[0]
            gmax = self._get_int(path=iE_base+'gmax')[0]
            self.scatter.append(scatt_data(gmin, gmax, self.scatt_order))
            for g in xrange(gmax - gmin + 1):
                iE_g_base = iE_base + str(g + gmin + 1)
                self.scatter[iE].outgoing[g][:] = np.asarray( \
                    self._get_double(n=self.scatt_order, path=iE_g_base))
                    
    def _read_chi(self):
        # NOT YET IMPLEMENTED!!        
        pass
        
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