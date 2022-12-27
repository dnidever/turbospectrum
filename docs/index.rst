.. turbospectrum documentation master file, created by
   sphinx-quickstart on Tue Feb 16 13:03:42 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*************
Turbospectrum
*************

Introduction
============
|Turbospectrum| is a generical stellar spectral synthesis package that can be run from python.  It's essentially
a redistribution of the `Turbospectrum <https://github.com/bertrandplez/Turbospectrum_NLTE>`_
Fortran spectral synthesis code by Bertrand Plez and a python wrapper/driver software (mostly reused from Jon Holtzman's
code in the `APOGEE package <https://github.com/sdss/apogee>`_).  The setup.py file has also been modified to
automatically compile the Fortran code and copy them to the user's python scripts directory.

.. toctree::
   :maxdepth: 1

   install
   modules
	      

Description
===========
To run |Turbospectrum| you need 1) a model atmosphere, 2) a linelist (or multiple), and 3) the set of stellar parameters
and elemental abundances that you want to run.



Examples
========

.. toctree::
    :maxdepth: 1

    examples

*****
Index
*****

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
