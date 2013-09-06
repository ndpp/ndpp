.. _quickinstall:

===================
Quick Install Guide
===================

This quick install guide outlines the basic steps needed to install NDPP on
your computer. For more detailed instructions on configuring and installing
NDPP, see :ref:`usersguide_install` in the User's Manual.

-------------
Prerequisites
-------------

In order to compile and run NDPP, a number of prerequisite software packages
and libraries may be needed. These include:

- A Fortran compiler such as gfortran_ (version 4.6 or newer)
- git_ version control software for obtaining source code (*optional*)
- HDF5_ Library for portable binary output format (*optional* and *not yet implemented*)

.. _gfortran: http://gcc.gnu.org/wiki/GFortran
.. _git: http://git-scm.com
.. _HDF5: http://www.hdfgroup.org/HDF5/

--------------------
Obtaining the Source
--------------------

All NDPP source code is hosted on BitBucket_. The latest version of the NDPP
source code can be downloaded from GitHub with the following command::

    git clone git@bitbucket.org:nelsonag/ndpp.git

If you don't have git installed, you can download the source as a zip or tar
file directly from BitBucket_.

.. _BitBucket: https://github.com/nelsonag/ndpp

-------------------------------
Compiling on Linux and Mac OS X
-------------------------------

To compile NDPP on Linux or Max OS X, run the following commands from within
the root directory of the source code:

.. code-block:: sh

    cd src
    make

This will build an executable named ``ndpp``.

--------------------
Compiling on Windows
--------------------

To compile NDPP on a Windows operating system, you will need to first install
Cygwin_, a Linux-like environment for Windows. When configuring Cygwin, make
sure you install both the gfortran_ compiler (version 4.6 or newer) as well as 
git. Once you have obtained the source code, run the following commands from 
within the source code root directory:

.. code-block:: sh

    cd src
    make

This will build an executable named ``ndpp``.

.. _Cygwin: http://www.cygwin.com/
