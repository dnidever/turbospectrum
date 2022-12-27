import os
import numpy as np
import shutil
import subprocess
import tempfile
import time
from . import utils, atomic, atmos


def synthesize(teff,logg,mh=0.0,am=0.0,cm=0.0,nm=0.0,vmicro=2.0,elems=None,
               wrange=[15000.0,17000.0],dw=0.1,atmod=None,atmos_type='marcs',
               dospherical=True,linelists=None,solarisotopes=False,workdir=None,
               save=False,verbose=False):
    """
    Code to synthesize a spectrum with Turbospectrum.
    
    Parameters
    ----------
    teff : float
       Effective temperature in K.
    logg : float
       Surface gravity.
    mh : float, optional
       Metallicity, [M/H].  Deftauls is 0.0 (solar).
    am : float, optional
       Alpha abundance, [alpha/M].  Default is 0.0 (solar).
    cm : float, optional
       Carbon abundance, [C/M].  Default is 0.0 (solar).
    nm : float, optional
       Nitrogen abundance, [N/M].  Default is 0.0 (solar).
    vmicro : float, optional
       Microturbulence in km/s.  Default is 2 km/s.
    solarisotopes : bool, optional
       Use solar isotope ratios, else "giant" isotope ratios ( default False ).
    elems : list, optional
       List of [element name, abundance] pairs.
    wrange : list, optional
       Two element wavelength range in A.  Default is [15000.0,17000.0].
    dw : float, optional
       Wavelength step.  Default is 0.1 A.
    atmod : str, optional
       Name of atmosphere model (default=None, model is determined from input parameters).
    atmos_type : str, optional
       Type of model atmosphere file.  Default is 'marcs'.
    dospherical : bool, optional
       Perform spherically-symmetric calculations (otherwise plane-parallel).  Default is True.
    linelists : list,
       List of linelists to use.
    save : bool, optional
       Save temporary directory and files for synthesis.  Default=False.
    workdir : str, optional
       Directory to perform the work in.  By default a temporary directory is
         created and deleted after the work is done (unless save=True).
    verbose : bool, optional
       Verbose output to the screen.

    Returns
    -------
    flux : numpy array
       The fluxed synthetic spectrum.
    continuum : numpy array
       The continuum of the spectrum.
    wave : numpy array
       Wavelength array in A.

    Example
    -------

    flux,cont,wave = synthesize(5000.0,2.5,-1.0)

    """

    t0 = time.time()
    
    # Default abundances
    abundances = atomic.solar()
    abundances[2:] += mh
    abundances[6-1] += cm
    abundances[7-1] += nm
    for i in [8,10,12,14,16,18,20,22]:
        abundances[i-1] += am
    # Abundance overrides from els, given as [X/M]
    if elems is not None :
        for el in elems:
            atomic_num = atomic.periodic(el[0])
            abundances[atomic_num-1] = atomic.solar(el[0]) + mh + el[1]

    # Change to temporary directory
    if workdir is None:
        workdir = tempfile.mkdtemp(prefix='turbo')
    cwd = os.getcwd()
    os.chdir(workdir)

    # Create the root name from the input parameters
    root = (atmos_type+'_t{:04d}g{:s}m{:s}a{:s}c{:s}n{:s}v{:s}').format(int(teff), atmos.cval(logg), 
                      atmos.cval(mh), atmos.cval(am), atmos.cval(cm), atmos.cval(nm),atmos.cval(vmicro))

    # Check that linelists and model atmosphere files exit
    if type(linelists) is str:
        linelists = [linelists]
    for l in linelists:
        if os.path.exists(l):
            raise FileNotFoundError(l)
    if os.path.exists(atmod)==False:
        raise FileNotFoundError(atmod)

    if dospherical and ('marcs' in atmos_type) and logg <= 3.001:
        spherical= True
    else:
        spherical = False
    flux,cont,wave = do_turbospec(root,atmod,linelists,mh,am,abundances,wrange,dw,
                                  save=save,solarisotopes=solarisotopes,
                                  babsma=None,atmos_type=atmos_type,
                                  spherical=spherical,tfactor=1.0)

    os.chdir(cwd)
    if not save:
        shutil.rmtree(workdir)

    if verbose:
        print('dt = {:.3f}s'.format(time.time()-t0))
        
    return flux,cont,wave

    
def do_turbospec(root,atmod,linelists,mh,am,abundances,wrange,dw,save=False,
                 solarisotopes=False,babsma=None,atmos_type='marcs',
                 spherical=True,vmicro=2.0,tfactor=1.0,verbose=False):
    """
    Runs Turbospectrum for specified input parameters.

    Parameters
    ----------
    root : str
       Root of filenames to use for this Turbospectrum run.
    atmod : str, optional
       Name of atmosphere model (default=None, model is determined from input parameters).
    linelists : list
       List of linelist filenames.
    mh : float, optional
       Metallicity, [M/H].  Default is 0.0 (solar).
    am : float, optional
       Alpha abundance, [alpha/M].  Default is 0.0 (solar).
    abundances : list
       List of abundances.
    wrange : list, optional
       Two element wavelength range in A.  Default is [15000.0,17000.0].
    dw : float, optional
       Wavelength step.  Default is 0.1 A.
    solarisotopes : bool, optional
       Use solar isotope ratios, else "giant" isotope ratios.  Default is False.
    babsma : bool, optional
       The name of hte babsma output file of opacities.  Default is None to run babmsa.
    atmos_type : str, optional
       Model atmosphere type.  Default is 'marcs'.
    spherical : bool, optional
       Spherical atmosphere.  Default is True.
    vmicro : float, optional
       Microturbulent velocity in km/s.  Default is 2.0 km/s.
    save : bool, optional
       Save temporary directory and files for synthesis.  Default=False.
    verbose : bool, optional
       Verbose output to the screen.

    Returns
    -------
    flux : numpy array
       The fluxed synthetic spectrum.
    continuum : numpy array
       The continuum of the spectrum.
    wave : numpy array
       Wavelength array in A.

    Example
    -------

    flux,cont,wave = do_turbospec(root,atmod,linelists,-0.1,0.2,abund,wrange=[15000.0,17000.0],dw=0.1)

    """

    # Turbospectrum setup
    datadir = utils.datadir()
    os.symlink(datadir,'./DATA')
    shutil.copy(atmod,'./'+os.path.basename(atmod))
    
    if verbose:
        stdout = None
    else:
        stdout = open(os.devnull, 'w')

    # Individual element grid?
    nels = len(abundances)

    welem = np.array(wrange)
    # Only compute opacities for a single nominal abundance
    if babsma is None:
        fout = open(root+'_babsma.csh','w')
        fout.write("#!/bin/csh -f\n")
        fout.write("babsma_lu << EOF\n")
        fout.write("'LAMBDA_MIN:'   '{:12.3f}'\n".format(welem.min()-dw))
        fout.write("'LAMBDA_MAX:'   '{:12.3f}'\n".format(welem.max()+dw))
        fout.write("'LAMBDA_STEP:'  '{:8.3f}'\n".format(dw))
        fout.write("'MODELINPUT:'  '{:s}'\n".format(os.path.basename(atmod)))
        if atmos_type != 'marcs':
            fout.write("'MARCS-FILE:'  '.false.'\n")
        fout.write("'MODELOPAC:'  '{:s}opac'\n".format(os.path.basename(root)))
        fout.write("'METALLICITY:'  '{:8.3f}'\n".format(mh))
        fout.write("'ALPHA/Fe:'  '{:8.3f}'\n".format(am))
        fout.write("'HELIUM:'  '{:8.3f}'\n".format(0.00))
        fout.write("'R-PROCESS:'  '{:8.3f}'\n".format(0.00))
        fout.write("'S-PROCESS:'  '{:8.3f}'\n".format(0.00))
        fout.write("'INDIVIDUAL ABUNDANCES:'  '{:2d}'\n".format(nels))
        for iel,abun in enumerate(abundances):
            fout.write("{:5d}  {:8.3f}\n".format(iel+1,abun))
        if not solarisotopes:
            fout.write("'ISOTOPES:'  '2'\n")
            # Adopt ratio of 12C/13C=15
            fout.write("   6.012 0.9375\n")
            fout.write("   6.013 0.0625\n")
        fout.write("'XIFIX:'  'T'\n")
        fout.write("{:8.3f}\n".format(vmicro))
        fout.write("EOF\n")
        fout.close()
        if run:
            os.chmod(root+'_babsma.csh', 0o777)
            ret = subprocess.check_output(['./'+os.path.basename(root)+'_babsma.csh'],stderr=subprocess.STDOUT)
            # Save the log file            
            if type(ret) is bytes:
                ret = ret.decode()
            with open(root+'_babsma.log','w') as f:
                f.write(ret)
        babsma = os.path.basename(root)+'opac'

    # Create bsyn control file
    bsynfile = root
    fout = open(bsynfile+'.inp','w')
    fout.write("'LAMBDA_STEP:'  '{:8.3f}'\n".format(dw))
    fout.write("'LAMBDA_MIN:'   '{:12.3f}'\n".format(welem.min()))
    fout.write("'LAMBDA_MAX:'   '{:12.3f}'\n".format(welem.max()))
    fout.write("'INTENSITY/FLUX:'  'Flux'\n")
    fout.write("'COS(THETA):'  '1.00'\n")
    fout.write("'ABFIND:'  '.false.'\n")
    fout.write("'MODELINPUT:'  '{:s}'\n".format(os.path.basename(atmod)))
    if atmos_type != 'marcs':
        fout.write("'MARCS-FILE:'  '.false.'\n")
    fout.write("'MODELOPAC:'  '{:s}'\n".format(babsma))
    fout.write("'RESULTFILE:'  '{:s}'\n".format(os.path.basename(root)))
    fout.write("'METALLICITY:'  '{:8.3f}'\n".format(mh))
    fout.write("'ALPHA/Fe:'  '{:8.3f}'\n".format(am))
    fout.write("'HELIUM:'  '{:8.3f}'\n".format(0.00))
    fout.write("'R-PROCESS:'  '{:8.3f}'\n".format(0.00))
    fout.write("'S-PROCESS:'  '{:8.3f}'\n".format(0.00))
    fout.write("'INDIVIDUAL ABUNDANCES:'  '{:2d}'\n".format(len(abundances)))
    for iel,abun in enumerate(abundances):
        fout.write("{:5d}  {:8.3f}\n".format(iel+1,abun))
    if not solarisotopes:
        fout.write("'ISOTOPES:'  '2'\n")
        # adopt ratio of 12C/13C=15
        fout.write("   6.012 0.9375\n")
        fout.write("   6.013 0.0625\n")
    fout.write("'NFILES:'  '{:4d}'\n".format(len(linelists)))
    for linelist in linelists: 
        fout.write(linelist+"\n")
    if spherical:
        fout.write("'SPHERICAL:'  'T'\n")
    else:
        fout.write("'SPHERICAL:'  'F'\n")
    fout.write("30\n")
    fout.write("300.00\n")
    fout.write("15\n")
    fout.write("1.3\n")
    fout.close()

    # Control file, with special handling in case bsyn goes into infinite loop ...
    fout = open(root+"_bsyn.csh",'w')
    fout.write("#!/bin/csh -f\n")
    fout.write("bsyn_lu < {:s} &\n".format(os.path.basename(bsynfile)+'.inp'))
    fout.close()

    # Run bsyn
    os.chmod(root+'_bsyn.csh', 0o777)
    ret = subprocess.check_output(['./'+os.path.basename(root)+'_bsyn.csh'],stderr=subprocess.STDOUT)
    # Save the log file
    if type(ret) is bytes:
        ret = ret.decode()
    with open(root+'_bsyn.log','w') as f:
        f.write(ret)
    try:
        out = np.loadtxt(root)
        wave = out[:,0]
        fluxnorm = out[:,1]
        flux = out[:,2]
        cont = flux/fluxnorm
    except :
        print('failed...',root,atmod,mh,am)
        return 0.,0.,0.
    return flux,cont,wave
