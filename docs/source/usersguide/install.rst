.. _usersguide_install:

==============================
Installation and Configuration
==============================

-------------------
Build Configuration
-------------------

The NDPP source package includes the source and build system necessary
to build the NDPP pre-processor.  NDPP can be compiled by entering the following
from a terminal in the ``ndpp/src/`` directory::

    make

All other details of compiling are the same for NDPP as for OpenMC and the reader
should reference to the OpenMC_ manual for further details.
    
--------------
Running OpenMC
--------------

Once an NDPP input file has been created (see :ref:`usersguide_input`), you can 
either run the ``ndpp`` executable directly from the directory containing your 
XML input files, or you can specify as a command-line argument the directory containing
the XML input files. For example, if the path of your NDPP executable is
``/home/username/ndpp/src/ndpp`` and your XML input files are in the
directory ``/home/username/somemodel/``, one way to run the simulation would be:

.. code-block:: sh

    cd /home/username/somemodel
    ndpp

Alternatively, you could run from any directory:

.. code-block:: sh

    ndpp /home/username/somemodel

Note that in the latter case, any output files will be placed in the present
working directory which may be different from ``/home/username/somemodel``.

.. _OpenMC: http://mit-crpg.github.io/openmc/
