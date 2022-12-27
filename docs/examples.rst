********
Examples
********


Running Turbospectrum
=====================
You need to supply |Turbospectrum| with 1) a model atmosphere, 2) a linelist (or multiple linelists), and 3) the set of stellar
parameters and abundances.

Here's a simple example:

    >>> from turbospectrum import synthesis,utils
    >>> atmod = utils.testdir()+'s6000_g+1.0_m0.5.mod'
    >>> linelists = [utils.testdir()+'nlte_linelist_test.txt',utils.datadir()+'Hlinedata']
    >>> flux,cont,wave = synthesis.synthesize(6000.0,1.0,0.5,atmod=atmod,linelists=linelists,wrange=[5000,5200])

Now plot the spectrum:

    >>> import matplotlib.pyplot as plt
    >>> plt.plot(wave,flux)

It should look like this.

.. image:: spectrum_example.png
  :width: 600
  :alt: Example Turbospectrum synthetic spectrum

Abundances
----------
	
You can modify the global alpha abundance with `am` or individual abundances with `elems`.  The `elems` parameter
takes a list of [element name, abundance] pairs, where the abundance should be in the form [X/M], where M is the
overall metallicity that is used to scale the individual abundances.  For example, ``elems=[['Mg',0.55],['Ba',-0.15]]``
means a relative Magnesium abundance of +0.55 and a relative Barium abundance of -0.15.

Let's try it out:

    >>> flux2,cont2,wave2 = synthesis.synthesize(6000.0,1.0,0.5,atmod=atmod,linelists=linelists,wrange=[5000,5200],elems=[['Mg',0.55],['Ba',-0.15]])
    >>> plt.plot(wave,flux)
    >>> plt.plot(wave2,flux2)
    >>> plt.xlim(5150,5200)
