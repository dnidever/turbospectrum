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

Fortran code
------------
    
The pip install will attempt to automatically compile the Fortran code and copy the binaries to your
Python scripts directory (which should be in your path).  If this fails for some reason, then you'll
need to compile it yourself.  You'll likely want to do a full git clone of the repository for this.
To compile the code you need either the Intel Fortran compiler (``ifort``), which should be
freely available online, or the GNU Fortran compiler (``gfortran``).  The Fortran code lives in the `src/`
directory.  All you should need to do is to cd into that directory and type ``make``.  The binaries
``babsma_lu`` and ``bsyn_lu`` will be copied to the repository's `bin/` directory.  Copy these to
a directory in your path (e.g., ~/bin/ or /usr/local/bin/).  


Dependencies
============

- numpy
- scipy
- astropy
- matplotlib
- `dlnpyutils <https://github.com/dnidever/dlnpyutils>`_
