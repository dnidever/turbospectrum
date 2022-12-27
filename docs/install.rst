************
Installation
************


Important packages
==================
`turbospectrum` is package to run the `Turbospectrum <https://github.com/bertrandplez/Turbospectrum_NLTE>`_
spectral synthesis code by Bertrand Plez and python wrapper/driver software (mostly reused from Jon Holtzman's
code in the `APOGEE package <https://github.com/sdss/apogee>`_).

Installing Turbospectrum
========================

The easiest way to install the code is with pip.  This will both compile and install the Fortran code as
well as the Python wrapper code.

.. code-block:: bash

    pip install turbospectrum


Dependencies
============

- numpy
- scipy
- astropy
- matplotlib
- `dlnpyutils <https://github.com/dnidever/dlnpyutils>`_
