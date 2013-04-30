.. _usersguide_troubleshoot:

====================
Troubleshooting NDPP
====================

-------------------------
Problems with Compilation
-------------------------

If you are experiencing problems trying to compile OpenMC and NDPP, first check if the
error you are receiving is among the following options.

Fatal Error: File 'xml_data_settings_t.mod' opened at (1) is not a GFORTRAN module file
***************************************************************************************

When OpenMC and NDPP compile, the first thing it needs to do is compile source in the
xml-fortran subdirectory. If you compiled everything with a compiler other than
gfortran, performed a :program:`make clean`, and then tried to :program:`make`
with gfortran, the xml-fortran modules would have been compiled with a different
compiler. To fix this, try clearing out all modules and object files with
:program:`make distclean` and then recompiling.

gfortran: unrecognized option '-cpp'
************************************

You are probably using a version of the gfortran compiler that is too
old. Download and install the latest version of gfortran_.

f951: error: unrecognized command line option "-fbacktrace"
***********************************************************

You are probably using a version of the gfortran compiler that is too
old. Download and install the latest version of gfortran_.


make[1]: ifort: Command not found
*********************************

You tried compiling with the Intel Fortran compiler and it was not found on your
:envvar:`PATH`. If you have the Intel compiler installed, make sure the shell
can locate it (this can be tested with :program:`which ifort`).

make[1]: pgf90: Command not found
*********************************

You tried compiling with the PGI Fortran compiler and it was not found on your
:envvar:`PATH`. If you have the PGI compiler installed, make sure the shell can
locate it (this can be tested with :program:`which pgf90`).

-------------------------
Problems with Simulations
-------------------------

Segmentation Fault
******************

A segmentation fault occurs when the program tries to access a variable in
memory that was outside the memory allocated for the program. The best way to
debug a segmentation fault is to re-compile OpenMC with debug options turned
on. First go to your ``openmc/src`` directory where OpenMC was compiled and type
the following commands:

.. code-block:: sh

    make distclean
    make DEBUG=yes

Now when you re-run your problem, it should report exactly where the program
failed. If after reading the debug output, you are still unsure why the program
failed, send an email to the OpenMC `developers
<mailto:paul.k.romano@gmail.com>`_.

ERROR: No cross_sections.xml file was specified in settings.xml or in the CROSS_SECTIONS environment variable.
**************************************************************************************************************

OpenMC and NDPP need to know where to find cross section data for each
nuclide. Information on what data is available and in what files is summarized
in a cross_sections.xml file. You need to tell OpenMC and NDPP where to find the
cross_sections.xml file either with the :ref:`cross_sections` in settings.xml or
with the :envvar:`CROSS_SECTIONS` environment variable. It is recommended to add
a line in your ``.profile`` or ``.bash_profile`` setting the
:envvar:`CROSS_SECTIONS` environment variable.

.. _gfortran: http://gcc.gnu.org/wiki/GFortran

