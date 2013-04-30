==========================================
Nuclear Data PreProcessor (NDPP) Code
==========================================

The NDPP program is a code which extracts nuclear data from ACE format 
data files (used by the Monte Carlo neutron transport codes MCNP_, Serpent_, and 
OpenMC_, among others), and pre-processes the data to a new data library
which has outgoing particle angle and energy probability distributions from both scattering and 
fission reactions in a form that is more readily tallied by the Monte Carlo 
codes with a significantly increased convergence rate.  The project started 
under the Fission Systems and Radiation Transport group at the University of
Michigan.

------------
Installation
------------

Detailed installation instructions can be found in the User's Guide.

---------------
Troubleshooting
---------------

If you run into problems compiling, installing, or running NDPP, first check
the Troubleshooting section in the User's Guide. If you are not able to find
a solution to your problem there, please contact the `developer`_.

--------------
Reporting Bugs
--------------

NDPP is hosted on BitBucket and all bugs are reported and tracked through the
Issues_ feature on BitBucket. However, BitBucket Issues should not be used for 
common troubleshooting purposes. If you are having trouble installing the code 
or getting your model to run properly, you should first send a message to the
`developer`_. If it turns out your issue really is a bug in the
code, an issue will then be created on BitBucket. If you want to request that a
feature be added to the code, you may create an Issue on BitBucket.

-------
License
-------

NDPP is distributed under the MIT/X license.


.. _MCNP: http://mcnp.lanl.gov
.. _Serpent: http://montecarlo.vtt.fi
.. _OpenMC: http://mit-crpg.github.io/openmc/index.html
.. _developer: mailto:nelsonag@umich.edu
