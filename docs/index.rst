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

1) Model Atmospheres

   Turbospectrum can read both MARCS model atmospheres and Kurucz/ATLAS model atmosphers.

2) Linelists

   Turbospectrum requires a specific linelist format.  It is described in more detail `here <_static/Readme-Linelist_format>`_.
   
3) Stellar parameters and elemental abundances.

   The main stellar parameters are Teff, logg, [M/H], and [alpha/M].  These are the first four parameters in the
   main ``synthesis.synthesize()`` function.  The individual elements abundances are given in the ``elems`` parameters
   as a list of [element name, abundance] pairs, where abundance in the in [X/M] format relative to the overall metallicity.

More details on how the Fortran Turbospectrum code works can be found in the  `Turbospectrum_v20_Documentation.pdf <_static/Turbospectrum_v20_Documentation.pdf>`_ documentation.

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
