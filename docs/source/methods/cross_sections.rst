.. _methods_cross_sections:

============================
Cross Section Representation
============================

The data governing the interaction of neutrons with various nuclei are
represented using the ACE format which is used by OpenMC_, MCNP_, and Serpent_. 
ACE-format data can be generated from `ENDF data`_ with the NJOY_ nuclear 
data processing system. The use of a standard cross section format allows for a
direct comparison of OpenMC with other codes since the same cross section
libraries can be used.

The ACE format contains continuous-energy cross sections for the following types
of reactions: elastic scattering, fission (or first-chance fission,
second-chance fission, etc.), inelastic scattering, :math:`(n,xn)`,
:math:`(n,\gamma)`, and various other absorption reactions. For those reactions
with one or more neutrons in the exit channel, secondary angle and energy
distributions may be provided. In addition, fissionable nuclides have total,
prompt, and/or delayed :math:`\nu` as a function of energy and neutron precursor
distributions. Many nuclides also have probability tables to be used for
accurate treatment of self-shielding in the unresolved resonance range. For
bound scatterers, separate tables with :math:`S(\alpha,\beta,T)` scattering law
data can be used.


.. _OpenMC: http://mit-crpg.github.io/openmc/
.. _MCNP: http://mcnp.lanl.gov
.. _Serpent: http://montecarlo.vtt.fi
.. _NJOY: http://t2.lanl.gov/codes.shtml
.. _ENDF data: http://www.nndc.bnl.gov/endf
