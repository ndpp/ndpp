.. _usersguide_input:

=======================
Writing XML Input Files
=======================

The input files for NDPP are structured in a set of XML_ files. XML,
which stands for eXtensible Markup Language, is a simple format that allows data
to be exchanged efficiently between different programs and interfaces.  

-----------------
Overview of Files
-----------------

To provide the input parameters for NDPP, there is one required input file, 
`ndpp.xml`, and another, `cross_sections.xml`, which is optional. 
The required `ndpp.xml` file provides NDPP with the options the user wishes to 
use when processing a set of nuclides.  The set of nuclides to evaluate is 
provided in the optional file, `cross_sections.xml`. This file is of the exact 
same format as the cross_sections.xml file used for OpenMC_.  If an explicit 
`cross_sections.xml` file is not provided (and referenced in `ndpp.xml`), then 
the default `cross_sections.xml` referenced by the `CROSS_SECTIONS` environment variable
will be used.  In effect, this means that the user's entire cross-section library
(in the default `cross_sections.xml` file) will be processed, unless a reduced
set is provided with a local `cross_sections.xml.'

--------------------------------------
NDPP Options Specification -- ndpp.xml
--------------------------------------

All calculational parameters and desired output options are specified in the
ndpp.xml file.  In the following discussion, if a parameter has a default value
then its presence is optional; if it is omitted, the default will be used.

``<scatt_type>`` Element
------------------------

The ``<scatt_type>`` element has no attributes and has an accepted
value of "tabular" or "legendre". If set to "tabular", the incoming 
energy-to-group scattering distributions will be output in a tabular format
with the number of bins defined by the ``<scatt_order>`` element.  If this
element is set to "legendre", the Legendre moments of the angular distribution
will be output with the scattering order defined in the ``<scatt_order>`` element.

  *Default*: legendre
  *NOTE*: The `tabular` format has not yet been implemented in NDPP

``<scatt_order>`` Element
-------------------------

The ``<scatt_order>`` element has no attributes and contains a single integer. 
As discussed in the ``<scatt_type>`` section, the ``<scatt_order>`` element 
provides either the number of tabular bins or the scattering order to NDPP, 
depending on the value of ``<scatt_type>``.

  *Default*: 5
  
``<mu_bins>`` Element
---------------------

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

  *Default*: The :envvar:`CROSS_SECTIONS` environment variable will be used to 
  find the path to the XML cross section listing.

``<energy_bins>`` Element
-------------------------

The ``<energy_bins>`` element provides the energy group structure to NDPP.
``<energy_bins>`` simply contains a monotonically increasing list of 
bounding energies for a number of groups. For example, if this element is specified as
``<energy_bins> 0.0 1.0 20.0 </energy_bins>``, then two energy groups
will be created, one with energies between 0 and 1 MeV and the other with
energies between 1 and 20 MeV.

``<integrate_chi>`` Element
---------------------------

The ``<integrate_chi>`` element has no attributes and has an accepted value of
"true" or "false". If set to "true", all fissionable nuclides will have their
fission neutron spectrum (:math:`\chi\left(E\right)`) integrated over the 
provided energy group structure and writen to the output files.  
If "false", then the :math:`\chi\left(E\right)` integration will not be performed.

  *Default*: true

``<thinning_tol>`` Element
--------------------------

The ``<thinning_tol>`` element has no attributes and accepts a single
floating-point number.  This element is used to set the percent tolerance for 
thinning the energy grid of the calculated data (:math:`\chi\left(E\right)` and 
the scattering distributions). The larger this value is the smaller the memory 
footprint is of the resultant data, but with decreased inaccuracy.

  *Default*: 0.2%
  *NOTE*: This feature is not yet implemented in NDPP
  
``<print_tol>`` Element
-----------------------

The ``<print_tol>`` element has no attributes and accepts a single
floating-point number.  This element is used to set the minimum value of
group-to-group transfers that will be printed.  Increasing this value 
decreases the output file size but can reduce accuracy of the resultant
preprocessed data library.

  *Default*: 1.0E-8
  
``<output_format>`` Element
---------------------------

The ``<output_format>`` element determines what format the preprocessed data
libraries should use.  This element has no attributes and accepts a string.  
Valid options are "ascii", "binary", "hdf5", "human", and "none".  If "ascii" is
specified, an output library will be written for each entry in the 
cross_sections.xml file which contains the requested data in ASCII text. 
If "binary" is specified, the same will be written, but in a 
machine-readable binary format.  If "hdf5" is specified, a single binary HDF5 
library will be created which contains the data for all the cross_sections.xml
file entries. If "human" is specified, then a more verbose form of the "ascii" 
format will be written which is useful for manual inspection of results.  
Finally, if "none" is specified, then no library will be written.

  *Default*: "ascii"
  
---------------------------------------------------------
Cross-Section Library Specification -- cross_sections.xml
---------------------------------------------------------

The `cross_sections.xml` file uses the same format used in OpenMC_; its format
and generation strategies are discussed at cross_sections.xml_

.. _XML: http://www.w3.org/XML/
.. _OpenMC: https://github.com/mit-crpg/openmc
.. _cross_sections.xml: http://mit-crpg.github.io/openmc/usersguide/install.html#cross-section-configuration
