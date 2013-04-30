.. _notes_0.5.0:

==============================
Release Notes for NDPP 0.1.0
==============================

-------------------
System Requirements
-------------------

There are no special requirements for running the NDPP code other than
those required by the OpenMC code. As of this release, NDPP and OpenMC have 
been tested on a variety of Linux distributions, Mac OS X, and Microsoft 
Windows 7. Memory requirements will vary depending on the size of
the problem at hand (mostly on the number of nuclides in the problem).

------------
New Features
------------

- All user input options that formerly accepted "off" or "on" should now be
  "false" or "true" (the proper XML schema datatype).

---------
Bug Fixes
---------

- 737b90_: Coincident surfaces from separate universes / particle traveling
  tangent to a surface.


.. _737b90: https://github.com/mit-crpg/openmc/commit/737b90
