"""Contains the object which is the main workhorse for scattering data
"""

class ScattData(object):
    """
    """

    @property
    def ne(self):
        return self._ne

    @ne.setter(self, ne):
        self._ne = ne

    @property
    def e_grid(self):
        return self._e_grid

    @e_grid.setter
    def e_grid(self,e_grid):
        self._e_grid = e_grid






    def __init___(self):
        self._is_init =