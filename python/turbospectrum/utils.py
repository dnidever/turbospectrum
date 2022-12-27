import os
import numpy as np

# The data directory
def datadir():
    """ Return the  data/ directory."""
    fil = os.path.abspath(__file__)
    codedir = os.path.dirname(fil)
    datadir = codedir+'/data/'
    return datadir
