.. _usersguide_input:

=======================
Writing XML Input Files
=======================

The input files for OpenMC and NDPP are structured in a set of XML_ files. XML,
which stands for eXtensible Markup Language, is a simple format that allows data
to be exchanged efficiently between different programs and interfaces.  See the
OpenMC manual for further information on XML files.

-----------------
Overview of Files
-----------------

To provide the input parameters for NDPP, there is one required input file, and another
which is optional. The first, required, input file is ndpp.xml which 
provides NDPP with the options the user wishes NDPP to use when processing a
set of nuclides.  The set of nuclides to evaluate is provided in the second, optional,
file, cross_sections.xml. This file is of the exact same format as the 
cross_sections.xml file used for OpenMC.  If an explicit cross_sections.xml 
file is not provided (and referenced in ``ndpp.xml), then the default 
cross_sections.xml referenced by the `CROSS_SECTIONS` environment variable
will be used.  This effectively means that if you do not want to process every nuclide, 
temperature, and evaluation provided in your default cross_sections.xml file, a 
separate file must be provided with the explicit nuclides requested.

--------------------------------------
NDPP Options Specification -- ndpp.xml
--------------------------------------

All calculational parameters and desired output options are specified in the
ndpp.xml file.

``<scatt_type>`` Element
----------------------------------

The ``<scatt_type>`` element has no attributes and has an accepted
value of "histogram" or "legendre". If set to "histogram", the incoming 
energy-to-group scattering distributions will be output in histogram format
with the number of bins defined by the ``<scatt_order>`` element.  If this
element is set to "legendre", the Legendre Moments of the angular distribution
will be output with the scattering order defined in the ``<scatt_order>`` element.

  *Default*: legendre

``<scatt_order>`` Element
----------------------------------

The ``<scatt_order>`` element has no attributes and contains a single integer. 
As discussed in the ``<scatt_type>`` section, the ``<scatt_order>`` element 
provides either the number of histogram bins or the scattering order to NDPP, 
depending on the value of ``<scatt_type>``.

  *Default*: 5
  
``<mu_bins>`` Element
----------------------------------

The ``<mu_bins>`` element has no attributes and contains a single integer.  This
value represents the number of angular points used when converting the ACE-
format angular distributions to an intermediate tabular representation.  
Increasing the number of points increases the output accuracy, however, it
comes at a cost of an increase in NDPP memory usage and computational cost.
This parameter has no impact on the runtime memory or computational costs of
the target Monte Carlo code.

  *Default*: 2001
  
.. _cross_sections:

``<cross_sections>`` Element
----------------------------

The ``<cross_sections>`` element has no attributes and simply indicates the path
to an XML cross section listing file (usually named cross_sections.xml). If this
element is absent from the settings.xml file, the :envvar:`CROSS_SECTIONS`
environment variable will be used to find the path to the XML cross section
listing.

``<energy_bins>`` Element
-------------------------

The ``<energy_bins>`` element provides the energy group structure to NDPP.
``<energy_bins>`` simply contains a monotonically increasing list of 
bounding energies for a number of groups. For example, if this element is specified as
``<energy_bins> 0.0 1.0 20.0 </energy_bins>``, then two energy groups
will be created, one with energies between 0 and 1 MeV and the other with
energies between 1 and 20 MeV.

``<integrate_chi>`` Element
-----------------------

The ``<integrate_chi>`` element has no attributes and has an accepted value of
"true" or "false". If set to "true", all fissionable nuclides will have their
fission neutron spectrum (chi) integrated over the provided energy group structure
and writen to the output files.  If "false", then the chi integration will not
be performed.

  *Default*: true

``<thinning_tol>`` Element
------------------

The ``<thinning_tol>`` element has no attributes and accepts a single
floating-point number.  This element is used to set the percent tolerance for 
thinning the energy grid of the calculated data (chi and the scattering 
distributions). The larger this value is the smaller the memory footprint is of the
resultant data, but the larger the inaccuracy introduced is.

  *Default*: 0.2%
  
``<print_tol>`` Element
------------------

The ``<print_tol>`` element has no attributes and accepts a single
floating-point number.  This element is used to set the minimum value for 
group-to-group transfer that will be printed.  Increasing this value 
decreases the output file size but can reduce accuracy when the libraries are
used in a Monte Carlo code.

  *Default*: 1.0E-8
  
``<output_format>`` Element
--------------------

The ``<output_format>`` element determines what type of output library/libraries
should be written to disk during the run. This element has no attributes and is
simply a string.  Valid options are "ascii", "binary", "hdf5", and "none".  
If "ascii" is specified, an output library will be written for each entry in the 
cross_sections.xml file which contains the requested data in human-readable 
ASCII code. If "binary" is specified, the same will be written, but in a 
machine-readable binary format.  If "hdf5" is specified, a single binary HDF5 
library will be created which contains the data for all the cross_sections.xml
file entries. Finally, if "none" is specified, then no library will be written.

  *Default*: "ascii"
