.. _usersguide_install:

==============================
Installation and Configuration
==============================

--------------------
Building from Source
--------------------

Prerequisites
-------------

.. admonition:: Required

    * A Fortran compiler such as gfortran_

      In order to compile NDPP, you will need to have a Fortran compiler
      installed on your machine. Since a number of Fortran 2003/2008 features
      are used in the code, it is recommended that you use the latest version of
      whatever compiler you choose. For gfortran_, it is necessary to use
      version 4.6.0 or above.

      If you are using Debian or a Debian derivative such as Ubuntu, you can
      install the gfortran compiler using the following command::

          sudo apt-get install gfortran

.. admonition:: Optional

    * An MPI implementation for distributed-memory parallel runs
    * This feature is not yet implemented*

      To compile with support for parallel runs on a distributed-memory
      architecture, you will need to have a valid implementation of MPI
      installed on your machine. The code has been tested and is known to work
      with the latest versions of both OpenMPI_ and MPICH_. Note that if using
      OpenMPI, make sure that --with-mpi-f90-size is not set to medium or large
      since this may prevent MPI calls from completing successfully in
      NDPP. OpenMPI and/or MPICH can be installed on Debian derivatives
      with::

          sudo apt-get install mpich2 libmpich2-dev
          sudo apt-get install openmpi1.6-bin libopenmpi1.6-dev

    * HDF5_ Library for portable binary output format

      To compile with support for HDF5_ output (highly recommended), you will
      need to have HDF5 installed on your computer. The installed version will
      need to have been compiled with the same compiler you intend to compile
      NDPP with.

    * git_ version control software for obtaining source code

.. _gfortran: http://gcc.gnu.org/wiki/GFortran
.. _OpenMPI: http://www.open-mpi.org
.. _MPICH: http://www.mpich.org
.. _HDF5: http://www.hdfgroup.org/HDF5/

Obtaining the Source
--------------------

All NDPP source code is hosted on GitHub_. You can download the source code
directly from GitHub or, if you have the git_ version control software installed
on your computer, you can use git to obtain the source code. The latter method
has the benefit that it is easy to receive updates directly from the GitHub
repository. With git installed and setup, the following command will download
the full source code from the GitHub repository::

    git clone git@ndpp.org:ndpp/ndpp.git

By default, the cloned repository will be set to the development branch. To
switch to the source of the latest stable release, run the following commands::

    cd ndpp/src
    git checkout master

.. _GitHub: https://github.com/ndpp/ndpp
.. _git: http://git-scm.com

Build Configuration
-------------------

All configuration for NDPP is done within the Makefile located in
``src/Makefile``. In the Makefile, you will see that there are a number of User
Options which can be changed. It is recommended that you do not change anything
else in the Makefile unless you are experienced with compiling and building
software using Makefiles. The following parameters can be set from the User
Options sections in the Makefile:

COMPILER
  This variable tells the Makefile which compiler to use. Valid options are
  gnu, and intel The default is gnu (gfortran).

DEBUG
  Enables debugging when compiling. The flags added are dependent on which
  compiler is used.

PROFILE
  Enables profiling using the GNU profiler, gprof.

OPTIMIZE
  Enables high-optimization using compiler-dependent flags. For gfortran and
  Intel Fortran, this compiles with -O3.

MPI
  Enables parallel runs using the Message Passing Interface. The MPI_DIR
  variable should be set to the base directory of the MPI implementation.

OPENMP
  Enables shared-memory parallelism using the OpenMP API. The Fortran compiler
  must support OpenMP.

HDF5
  Enables HDF5 output in addition to normal screen and text file output. The
  HDF5_DIR variable should be set to the base directory of the HDF5
  installation.
  * This feature is not yet implemented*

It is also possible to change these options from the command line itself. For
example, if you want to compile with DEBUG turned on without actually change the
Makefile, you can enter the following from a terminal::

    make DEBUG=yes

Compiling on Linux and Mac OS X
-------------------------------

To compile NDPP on Linux or Max OS X, run the following commands from within
the root directory of the source code:

.. code-block:: sh

    cd src
    make
    sudo make install

This will build an executable named ``ndpp`` and install it (by default in
/usr/local/bin).

Compiling on Windows
--------------------

Using Cygwin
++++++++++++

One option for compiling NDPP on a Windows operating system is to use Cygwin_,
a Linux-like environment for Windows. You will need to first `install
Cygwin`_. When you are asked to select packages, make sure the following are
selected:

* Devel: gcc4-core
* Devel: gcc4-fortran
* Devel: make

If you plan on obtaining the source code directly using git, select the
following packages:

* Devel: git
* Devel: git-completion (Optional)
* Devel: gitk (Optional)

In order to use the Python scripts provided with NDPP, you will also need to
install Python. This can be done within Cygwin or directly in Windows. To
install within Cygwin, select the following packages:

* Python: python (Version > 2.7 recommended)

Once you have obtained the source code, run the following commands from within
the source code root directory:

.. code-block:: sh

    cd src
    make

This will build an executable named ``ndpp``.

.. _Cygwin: http://cygwin.com/
.. _install Cygwin: http://cygwin.com/setup.exe

Using MinGW
+++++++++++

An alternate option for installing NDPP on Windows is using MinGW_, which
stands for Minimalist GNU for Windows. An executable for installing the MinGW
distribution is available on SourceForge_. When installing MinGW, make sure the
following components are selected:

* MinGW Compiler Suite: Fortran Compiler
* MSYS Basic System

Once MinGW is installed, copy the NDPP source distribution to your MinGW home
directory (usually C:\\MinGW\\msys\\1.0\\home\\YourUsername). Once you have
the source code in place, run the following commands from within the MinGW shell
in the root directory of the NDPP distribution:

.. code-block:: sh

    cd src
    make

This will build an executable named ``ndpp``.

.. _MinGW: http://www.mingw.org
.. _SourceForge: http://sourceforge.net/projects/mingw

---------------------------
Cross Section Configuration
---------------------------


In order to run a simulation with NDPP, you will need cross section data for
each nuclide in your problem. Since NDPP uses ACE format cross sections, you
can use nuclear data that was processed with NJOY_, such as that distributed
with MCNP_ or Serpent_. Several sources provide free processed ACE data as
described below. The TALYS-based evaluated nuclear data library, TENDL_, is also
openly available in ACE format.

In the following discussion, note that the ``cross_sections.xml`` file can be
the same as is used by OpenMC_.  Therefore, these steps can be skipped if
OpenMC_ is currently installed on the system.


Using ENDF/B-VII.1 Cross Sections from NNDC
-------------------------------------------

The NNDC_ provides ACE data from the ENDF/B-VII.1 neutron and thermal scattering
sublibraries at four temperatures processed using NJOY_. To use this data with
NDPP, a script is provided with NDPP that will automatically download,
extract, and set up a configuration file:

.. code-block:: sh

    cd ndpp/data
    python get_nndc_data.py

At this point, you should set the :envvar:`CROSS_SECTIONS` environment variable
to the absolute path of the file ``ndpp/data/nndc/cross_sections.xml``.

Using JEFF Cross Sections from OECD/NEA
---------------------------------------

The NEA_ provides processed ACE data from the JEFF_ nuclear library upon
request. A DVD of the data can be requested here_. To use this data with NDPP,
the following steps must be taken:

1. Copy and unzip the data on the DVD to a directory on your computer.
2. In the root directory, a file named ``xsdir``, or some variant thereof,
   should be present. This file contains a listing of all the cross sections and
   is used by MCNP. This file should be converted to a ``cross_sections.xml``
   file for use with NDPP. A Python script is provided in the NDPP
   distribution for this purpose:

   .. code-block:: sh

       ndpp/src/utils/convert_xsdir.py xsdir31 cross_sections.xml

3. In the converted ``cross_sections.xml`` file, change the contents of the
   <directory> element to the absolute path of the directory containing the
   actual ACE files.
4. Additionally, you may need to change any occurrences of upper-case "ACE"
   within the ``cross_sections.xml`` file to lower-case.
5. Either set the :ref:`cross_sections` in a settings.xml file or the
   :envvar:`CROSS_SECTIONS` environment variable to the absolute path of the
   ``cross_sections.xml`` file.

Using Cross Sections from MCNP
------------------------------

To use cross sections distributed with MCNP, change the <directory> element in
the ``cross_sections.xml`` file in the root directory of the NDPP distribution
to the location of the MCNP cross sections. Then, either set the
:ref:`cross_sections` in a settings.xml file or the :envvar:`CROSS_SECTIONS`
environment variable to the absolute path of the ``cross_sections.xml`` file.

Using Cross Sections from Serpent
---------------------------------

To use cross sections distributed with Serpent, change the <directory> element
in the ``cross_sections_serpent.xml`` file in the root directory of the NDPP
distribution to the location of the Serpent cross sections. Then, either set the
:ref:`cross_sections` in a settings.xml file or the :envvar:`CROSS_SECTIONS`
environment variable to the absolute path of the ``cross_sections_serpent.xml``
file.

.. _NJOY: http://t2.lanl.gov/nis/codes.shtml
.. _NNDC: http://www.nndc.bnl.gov/endf/b7.1/acefiles.html
.. _NEA: http://www.oecd-nea.org
.. _JEFF: http://www.oecd-nea.org/dbdata/jeff/
.. _here: http://www.oecd-nea.org/dbdata/pubs/jeff312-cd.html
.. _MCNP: http://mcnp.lanl.gov
.. _Serpent: http://montecarlo.vtt.fi
.. _TENDL: ftp://ftp.nrg.eu/pub/www/talys/tendl2012/tendl2012.html
.. _OpenMC: https://github.com/mit-crpg/openmc

------------
Running NDPP
------------

Once an NDPP input file has been created (see :ref:`usersguide_input`), you can
either run the ``ndpp`` executable directly from the directory containing your
XML input files, or you can specify as a command-line argument the directory containing
the XML input files. For example, if the path of your NDPP executable is
``/home/username/ndpp/src/ndpp`` and your XML input files are in the
directory ``/path/to/somemodel``, one way to run the simulation would be:

.. code-block:: sh

    cd /path/to/somemodel
    ndpp

Alternatively, you could run from any directory:

.. code-block:: sh

    ndpp /path/to/someplace

Note that in the latter case, any output files will be placed in the present
working directory which may be different from ``/path/to/somemodel``.