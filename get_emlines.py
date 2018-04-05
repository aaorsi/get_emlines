import numpy as np
import os.path
import read_photoion as rp
import pylab as pl
import matplotlib.gridspec as gridspec
from matplotlib.mlab import griddata
import scipy.stats
import warnings

def qZrelation(Zgas, q0,g0):# = 2.8e7, g0 = -1.3):
# Default parameters from Orsi+14.  
  return q0 * (Zgas/0.012)**g0


def get_lumlines(sfr, metal, LineProps, interp_func = 'interp2d', verbose = False):
  """
  LineProps contains Linesinfo, LinesArr, q(z) relation, line name and flux limit
  GalArr is the galaxy data dict.
  """
  Rootdir = os.path.dirname(rp.__file__)
  
  if verbose:
    print 'Rootdir %s' % Rootdir
  if Rootdir == '':
    warnings.warn('Rootdir is empty!')
  linesinfo, linesarr = rp.read_photoion(Rootdir)

  linename = LineProps['linename']

  ngals = len(sfr) if hasattr(sfr,"__len__") else 1
  qgals = qZrelation(metal, LineProps['q0'],LineProps['g0'])

  Nlyc = np.log10(1.35) + np.log10(sfr) + 53.0

  if verbose:
    print 'Computing emission lines for %d galaxies\n' % (ngals)
  lum_line = [] # np.zeros(ngals)
  linefunc = rp.get_2dfunc(linesinfo, linesarr, lname=linename, interp_func=interp_func) 
  hafunc   = rp.get_2dfunc(linesinfo, linesarr, lname='Halpha', interp_func=interp_func)
  
  n_alllines = len(linesinfo['Linename'])
  nlines = len(linename)

  if linename == 'All' or linename == 'all':
    nlines = n_alllines
    linename = linesinfo['linename']


  if verbose:
    print 'Running lines with g0=%f' % LineProps['g0']

#  import pdb ; pdb.set_trace()
  lum_line = np.zeros([ngals,nlines])

  ltype = []
  for i in range(nlines):
    _lt = tuple([linename[i],np.float32])
    ltype.append(_lt)
  
  lum_line = np.zeros(ngals,dtype=np.dtype(ltype))

  for i in range(ngals):
    lum = (rp.integ_line(linefunc,hafunc, qgals[i],metal[i],
             Nlyc[i],nlines, lname=linename) if hasattr(qgals,"__len__") 
             else rp.integ_line(linefunc,hafunc, qgals,metal,Nlyc,nlines, lname=linename,
             all_lines=all_lines))

    for j in range(nlines):
      jname = linesinfo['Linename'][j]
      lum_line[jname][i] = lum[j]

  return lum_line



def get_emlines(sfr,metals,q0=2.8e7,g0=-1.3,linename='All', interp_func ='interp2d',test_interp = False, verbose=False):
  
  # interp_func can be 'interp2d' or 'RectBivariateSpline'
  # test_interp is set to True for making plots that show the goodness of the interpolation within the grid. (TODO) 

  #if np.isnan(metals).any() or np.isnan(sfr).any():
  #  raise ValueError('\nget_emlines(): Input Zgas or sfr contain(s) NaN')

  ngal = len(sfr) if hasattr(sfr, "__len__") else 1

  LineProps = {'q0':float(q0), 'g0':float(g0),'linename':linename}
  ll = get_lumlines(sfr,metals,LineProps, interp_func= interp_func, verbose=verbose)
  lumline = np.array(ll)

  return lumline
