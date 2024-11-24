import os
import numpy as np

# The data directory
def datadir():
    """ Return the data/ directory."""
    fil = os.path.abspath(__file__)
    codedir = os.path.dirname(fil)
    datadir = codedir+'/data/'
    return datadir

# The test directory
def testdir():
    """ Return the test/ directory."""
    fil = os.path.abspath(__file__)
    codedir = os.path.dirname(fil)
    testdir = codedir+'/test/'
    return testdir

def identify_atmos(modelfile,verbose=False):
    """
    Idenfies the type of model atmosphere in an input file

    Valid options are kurucz, marcs, tlusty (.7) or phoenix

    Parameters
    ----------
    modelfile: str
      file with a model atmosphere

    Returns
    -------
    atmostype: str
       can take the value 'kurucz', 'marcs', 'tlusty' or 'phoenix' 

   """

    if ('PHOENIX' in modelfile and 'fits' in modelfile):
        atmostype = 'phoenix'
    else: 
        if modelfile[-3:] == '.gz':
            f = gzip.open(modelfile,'rt')
        else:
            f = open(modelfile,'r')
        line = f.readline()
        if verbose:
            print('modelfile / line=',modelfile,line)
        if ('TEFF' in line):
            atmostype = 'kurucz'
        else: 
            line = f.readline()
            if ('Teff' in line):
                atmostype = 'marcs'
            else:
                atmostype = 'tlusty'
        f.close()
   
    return atmostype

