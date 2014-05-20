.. _usersguide_utilities:

==============
NDPP Utilities
==============

-----------------
NDPP_DATA Library
-----------------

NDPP is distributed with a python interface to read binary output libraries.
This library provides functionality to:

1)  Read in a binary library which was printed by NDPP

  a)  This is provided by the NDPP_lib class constructor

2)  Condense the outgoing angular moments to one-group values

  a)  This is provided by NDPP_lib.condense_outgoing_scatt method.

  b)  This method has an argument of `groups` which is a list of groups to which are to be condensed in to one.

3)  Expand the Legendre moments in to functions of angle

  a)  This is provided by NDPP_lib.expand_scatt

4)  A method to expand the library's distributions and report on incoming energy and outgoing group pairs who have negative probability distribution functions.

  a)  This is provided by NDPP_lib.test_scatt_positivity.

Examples of these methods being used can be seen in both ``test.py`` and
``validate_libary.py``.
