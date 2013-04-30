=============================================
The NDPP Monte Carlo Data Pre-Processing Code
=============================================

NDPP [Nuclear Data PreProcessor] is a Monte Carlo neutron transport data 
preprocessing code which converts continuous-energy ACE data in to 
information that is more readily tallied by a Monte Carlo for the generation
of multi-group cross-sections. Specifically, it integrates the outgoing energy
and angular distributions of a reaction over the outgoing energies and angles
in to a format and energy group structure defined by the user.  NDPP uses 
subroutines written for the OpenMC Monte Carlo particle transport code.

NDPP is developed by the `University of Michigan`_.  For 
more information, feel free to contact `Adam Nelson`_.

The development of OpenMC is led by the `Computational Reactor Physics Group`_
at the `Massachusetts Institute of Technology`_. 

.. _University of Michigan: http://www-ners.engin.umich.edu
.. _Adam Nelson: mailto:nelsonag@umich.edu
.. _Computational Reactor Physics Group: http://crpg.mit.edu
.. _Massachusetts Institute of Technology: http://web.mit.edu
.. _Paul Romano: mailto:paul.k.romano@gmail.com

--------
Contents
--------

.. toctree::
    :maxdepth: 1

    quickinstall
    releasenotes/index
    methods/index
    usersguide/index
    publications
    license
