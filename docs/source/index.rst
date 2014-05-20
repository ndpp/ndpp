===================================================
The NDPP Continuous Energy Data Pre-Processing Code
===================================================

NDPP [Nuclear Data PreProcessor] is a continuous-energy neutron data
preprocessing code which converts continuous-energy ACE data to
information that is more readily tallied by a Monte Carlo solver for the
generation of multi-group cross-sections. Specifically, it integrates the
outgoing energy and angular distributions of a reaction to produce a format
(Legendre moments or tabular data) and energy group structure defined by the
user.  NDPP utilizes code written for the OpenMC Monte Carlo neutron
transport code.

NDPP is developed by the Adam Nelson and contributors.  For
more information, feel free to contact the lead developer, `Adam Nelson`_.

The development of OpenMC is led by the `Computational Reactor Physics Group`_
at the `Massachusetts Institute of Technology`_.

.. _University of Michigan: http://www-ners.engin.umich.edu
.. _Adam Nelson: mailto:nelsonag@umich.edu
.. _Computational Reactor Physics Group: http://crpg.mit.edu
.. _Massachusetts Institute of Technology: http://web.mit.edu

.. toctree::
    :maxdepth: 1

    quickinstall
    usersguide/index
    methods/index
    releasenotes/index
    publications
    license
