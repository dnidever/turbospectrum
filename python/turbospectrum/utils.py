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
