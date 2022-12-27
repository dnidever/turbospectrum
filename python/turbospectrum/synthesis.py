import os
import numpy as np
import shutil
import subprocess
import tempfile
import time
from astropy.table import Table
from . import atoms, atomic


def synthesize(teff,logg,mh=0.0,am=0.0,cm=0.0,nm=0.0,
               vmicro=2.0,elems=None,wrange=[15000.0,17000.0],dw=0.1,
               atoms=True,molec=True,save=False,run=True,verbose=False):
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
    elems : list, optional
       List of [element name, abundance] pairs.
    wrange : list, optional
       Two element wavelength range in A.  Default is [15000.0,17000.0].
    dw : float, optional
       Wavelength step.  Default is 0.1 A.
    atoms : bool, optional
       Include atomic lines in the calculation.  Default is True.
    molec : bool, optional
       Include molecular lines in the calculation.  Default is True.
    save : bool, optional
       Save temporary directory and files for synthesis.  Default=False.
    run : bool, optional
       Actually do the synthesis.  Default is True.
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
    cwd = os.getcwd()
    os.chdir(workdir)
            
    # Do the opacity calculations first, then synthesis for all abund
    out = do_turbospec(root,atmod,linelists,mh,am,abundances,wrange,dw,save=save,run=run,
                       solarisotopes=solarisotopes,bsyn=False,atmos_type=atmos_type,vmicro=vmicro)


    if dospherical and ('marcs' in atmos_type) and logg <= 3.001:
        spherical= True
    else:
        spherical = False
    if ielem == 0: 
        wave,flux,fluxnorm = do_turbospec(file,atmod,linelists,mh,am,abundances,wrange,dw,
                                          save=save,run=run,solarisotopes=solarisotopes,
                                          babsma=root+'opac',atmos_type=atmos_type,
                                          spherical=spherical,tfactor=tfactor)
    else:
        out = do_turbospec(file,atmod,linelists,mh,am,abundances,wrange,dw,
                           save=save,run=run,solarisotopes=solarisotopes,
                           babsma=root+'opac',atmos_type=atmos_type,spherical=spherical,
                           tfactor=tfactor)

    # Load into final arrays
    if ielem == 0:
        spec = flux
        specnorm = fluxnorm
    else:
        spec = np.vstack([spec,flux])
        specnorm = np.vstack([specnorm,fluxnorm])

    os.chdir(cwd)
    if run:
        if not save:
            shutil.rmtree(workdir)
        return spec, specnorm

    
def do_turbospec(filebase,atmod,linelists,mh,am,abundances,wrange,dw,save=False,run=True,
                 solarisotopes=False,babsma=None,bsyn=True,atmos_type='marcs',
                 spherical=True,vmicro=2.0,tfactor=1.0,verbose=False):
    """
    Runs Turbospectrum for specified input parameters.

    Parameters
    ----------
    filebase : str
       Base of filenames to use for this Turbospectrum run.
    atmod : str, optional
       Name of atmosphere model (default=None, model is determined from input parameters).
    linelists : list
       List of linelist filenames.
    mh : float, optional
       Metallicity, [M/H].  Deftauls is 0.0 (solar).
    am : float, optional
       Alpha abundance, [alpha/M].  Default is 0.0 (solar).
    abundances : list
       List of abundances.
    wrange : list, optional
       Two element wavelength range in A.  Default is [15000.0,17000.0].
    dw : float, optional
       Wavelength step.  Default is 0.1 A.
    save : bool, optional
       Save temporary directory and files for synthesis.  Default=False.
    run : bool, optional
       Actually do the synthesis.  Default is True.
    solarisotopes : bool, optional
       Use solar isotope ratios, else "giant" isotope ratios.  Default is False.
    babsma : bool, optional
       Run Turbospectrum "babsma" binary code to calculate opacities.  Default is None.
    bsyn : bool, optional
       Run Turbospectrum "bsyn" binary code to calculate spectral synthesis.  Default is True.
    atmos_type : str, optional
       Model atmosphere type.  Default is 'marcs'.
    spherical : bool, optional
       Spherical atmosphere.  Default is True.
    vmicro : float, optional
       Microturbulent velocity in km/s.  Default is 2.0 km/s.
    tfactor : float, optional
       Some factor.  Default is 1.0.
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

    flux,cont,wave = do_turbospec(filebase,atmod,linelists,-0.1,0.2,abund,wrange=[15000.0,17000.0],dw=0.1)

    """

    # Turbospectrum setup
    try:
        os.symlink(os.environ['APOGEE_DIR']+'/src/turbospec/DATA','./DATA')
    except:
        pass

    if verbose:
        stdout = None
    else:
        stdout = open(os.devnull, 'w')

    # Individual element grid?
    nels = len(abundances)

    welem = np.array(wrange)
    # Only compute opacities for a single nominal abundance
    if babsma is None:
        fout = open(filebase+'_babsma.csh','w')
        fout.write("#!/bin/csh -f\n")
        fout.write("{:s}/bin/babsma_lu << EOF\n".format(os.environ['APOGEE_DIR']))
        fout.write("'LAMBDA_MIN:'   '{:12.3f}'\n".format(welem.min()-dw))
        fout.write("'LAMBDA_MAX:'   '{:12.3f}'\n".format(welem.max()+dw))
        fout.write("'LAMBDA_STEP:'  '{:8.3f}'\n".format(dw))
        fout.write("'MODELINPUT:'  '{:s}'\n".format(os.path.basename(atmod)))
        if atmos_type != 'marcs' : fout.write("'MARCS-FILE:'  '.false.'\n")
        fout.write("'MODELOPAC:'  '{:s}opac'\n".format(os.path.basename(filebase)))
        fout.write("'METALLICITY:'  '{:8.3f}'\n".format(mh))
        fout.write("'ALPHA/Fe:'  '{:8.3f}'\n".format(am))
        fout.write("'HELIUM:'  '{:8.3f}'\n".format(0.00))
        fout.write("'R-PROCESS:'  '{:8.3f}'\n".format(0.00))
        fout.write("'S-PROCESS:'  '{:8.3f}'\n".format(0.00))
        fout.write("'INDIVIDUAL ABUNDANCES:'  '{:2d}'\n".format(nels))
        for iel,abun in enumerate(abundances) :
            fout.write("{:5d}  {:8.3f}\n".format(iel+1,abun))
        if not solarisotopes :
          fout.write("'ISOTOPES:'  '2'\n")
          # Adopt ratio of 12C/13C=15
          fout.write("   6.012 0.9375\n")
          fout.write("   6.013 0.0625\n")
        fout.write("'XIFIX:'  'T'\n")
        fout.write("{:8.3f}\n".format(vmicro))
        fout.write("EOF\n")
        fout.close()
        if run:
            os.chmod(filebase+'_babsma.csh', 0o777)
            subprocess.call(['time','./'+os.path.basename(filebase)+'_babsma.csh'],stdout=stdout,stderr=stdout)
        babsma = os.path.basename(filebase)+'opac'

    if not bsyn:
        return

    # Create bsyn control file
    bsynfile = filebase
    fout = open(bsynfile+'.inp','w')
    fout.write("'LAMBDA_STEP:'  '{:8.3f}'\n".format(dw))
    fout.write("'LAMBDA_MIN:'   '{:12.3f}'\n".format(welem.min()))
    fout.write("'LAMBDA_MAX:'   '{:12.3f}'\n".format(welem.max()))
    fout.write("'INTENSITY/FLUX:'  'Flux'\n")
    fout.write("'COS(THETA):'  '1.00'\n")
    fout.write("'ABFIND:'  '.false.'\n")
    fout.write("'MODELINPUT:'  '{:s}'\n".format(os.path.basename(atmod)))
    if atmos_type != 'marcs' : fout.write("'MARCS-FILE:'  '.false.'\n")
    fout.write("'MODELOPAC:'  '{:s}'\n".format(babsma))
    fout.write("'RESULTFILE:'  '{:s}'\n".format(os.path.basename(filebase)))
    fout.write("'METALLICITY:'  '{:8.3f}'\n".format(mh))
    fout.write("'ALPHA/Fe:'  '{:8.3f}'\n".format(am))
    fout.write("'HELIUM:'  '{:8.3f}'\n".format(0.00))
    fout.write("'R-PROCESS:'  '{:8.3f}'\n".format(0.00))
    fout.write("'S-PROCESS:'  '{:8.3f}'\n".format(0.00))
    fout.write("'INDIVIDUAL ABUNDANCES:'  '{:2d}'\n".format(len(abundances)))
    for iel,abun in enumerate(abundances):
        fout.write("{:5d}  {:8.3f}\n".format(iel+1,abun))
    if  not solarisotopes :
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
    fout = open(filebase+"_bsyn.csh",'w')
    fout.write("#!/bin/csh -f\n")
    fout.write("{:s}/bin/bsyn_lu < {:s} &\n".format(os.environ['APOGEE_DIR'],os.path.basename(bsynfile)+'.inp'))
    fout.write('set bsynjob = $!\n')
    fout.write("set ok = 0\n")
    fout.write("set runtime = `ps -p $bsynjob -o cputime | tail -1 | awk -F: '{print ($1*3600)+($2*60)+$3}'`\n")
    tmax = 120*int(0.05/min([0.05,dw]))*tfactor
    fout.write('while ( $runtime < {:d} && $ok == 0 )\n'.format(tmax))
    fout.write('  usleep 200000\n')
    fout.write("  set runtime = `ps -p $bsynjob -o cputime | tail -1 | awk -F: '{print ($1*3600)+($2*60)+$3}'`\n")
    fout.write('  if ( `ps -p $bsynjob -o comm=` == "" ) then\n')
    fout.write('    echo process done, exiting!\n')
    fout.write('    set ok = 1\n')
    fout.write('  endif\n')
    fout.write('end\n')
    fout.write('if ( $ok == 0 ) then\n')
    fout.write('  echo expired, killing job\n')
    fout.write('  kill $bsynjob\n')
    fout.write('endif\n')
    fout.close()
    if run:
        os.chmod(filebase+'_bsyn.csh', 0o777)
        subprocess.call(['time','./'+os.path.basename(filebase)+'_bsyn.csh'],stdout=stdout,stderr=stdout)
        try:
            out = np.loadtxt(filebase)
            wave = out[:,1]
            specnorm = out[:,1]
            spec = out[:,2]
        except :
            print('failed...',filebase,atmod,mh,am)
            return 0.,0.,0.
        return wave,spec,specnorm
