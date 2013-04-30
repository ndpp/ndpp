.. _usersguide_install:

==============================
Installation and Configuration
==============================

-------------------
Build Configuration
-------------------

The OpenMC source package includes the source and build system necessary
to build the NDPP pre-processor.  NDPP can be compiled when OpenMC is compiled
by setting the NDPP flag to yes in the OpenMC Makefile or you can enter the following
from a terminal in the ``openmc/src/`` directory::

    make NDPP=yes

All other details of compiling are the same for NDPP as for OpenMC and the reader
should reference to the OpenMC manual for further details.
    
--------------
Running OpenMC
--------------

Once an NDPP input file has been created (see :ref:`usersguide_input`), you can 
either run the ``ndpp`` executable directly from the directory containing your 
XML input files, or you can specify as a command-line argument the directory containing
the XML input files. For example, if the path of your NDPP executable is
``/home/username/openmc/ndpp/src/ndpp`` and your XML input files are in the
directory ``/home/username/somemodel/``, one way to run the simulation would be:

.. code-block:: sh

    cd /home/username/somemodel
    ndpp

Alternatively, you could run from any directory:

.. code-block:: sh

    ndpp /home/username/somemodel

Note that in the latter case, any output files will be placed in the present
working directory which may be different from ``/home/username/somemodel``.

-----------------------------------------------------
Configuring Input Validation with GNU Emacs nXML mode
-----------------------------------------------------

The `GNU Emacs`_ text editor has a built-in mode that extends functionality for
editing XML files. One of the features in nXML mode is the ability to perform
real-time `validation`_ of XML files against a `RELAX NG`_ schema. The OpenMC
source contains RELAX NG schemas for each type of user input file. In order for
nXML mode to know about these schemas, you need to tell emacs where to find a
"locating files" description. Adding the following lines to your ``~/.emacs``
file will enable real-time validation of XML input files:

.. code-block:: common-lisp

    (require 'rng-loc)
    (add-to-list 'rng-schema-locating-files "~/openmc/schemas.xml")

Make sure to replace the last string on the second line with the path to the
schemas.xml file in your own OpenMC source directory.

.. _GNU Emacs: http://www.gnu.org/software/emacs/
.. _validation: http://en.wikipedia.org/wiki/XML_validation
.. _RELAX NG: http://relaxng.org/
