import numpy as np
import os.path
import read_photoion as rp
import pylab as pl
import matplotlib.gridspec as gridspec
from matplotlib.mlab import griddata
import scipy.stats

def qZrelation(Zgas, q0,g0):# = 2.8e7, g0 = -1.3):
# Default parameters from Orsi+14.  
  return q0 * (Zgas/0.012)**g0


def get_lumlines(sfr, metal, LineProps,all_lines=False):
  """
  LineProps contains Linesinfo, LinesArr, q(z) relation, line name and flux limit
  GalArr is the galaxy data dict.
  """
  Rootdir = os.path.dirname(rp.__file__)
  linesinfo, linesarr = rp.read_photoion(Rootdir)

  linename = LineProps['linename']

  ngals = len(sfr) if hasattr(sfr,"__len__") else 1
  qgals = qZrelation(metal, LineProps['q0'],LineProps['g0'])

  Nlyc = np.log10(1.35) + np.log10(sfr) + 53.0

  print 'Computing emission lines for %d galaxies\n' % (ngals)
  lum_line = [] # np.zeros(ngals)
  linefunc = rp.get_2dfunc(linesinfo, linesarr, lname=linename, all_lines=all_lines) 
  hafunc   = rp.get_2dfunc(linesinfo, linesarr, lname='Halpha',all_lines = False)
  nlines = len(linesinfo['Linename'])

  print 'Running lines with g0=%f' % LineProps['g0']

#  import pdb ; pdb.set_trace()
  lum_line = np.zeros([ngals,nlines])

  ltype = []
  for i in range(nlines):
    _lt = tuple([linesinfo['Linename'][i],np.float32])
    ltype.append(_lt)
  
  lum_line = np.zeros(ngals,dtype=np.dtype(ltype))

  for i in range(ngals):
    lum = (rp.integ_line(linefunc,hafunc, qgals[i],metal[i],
             Nlyc[i],nlines, lname=linename,all_lines=all_lines) if hasattr(qgals,"__len__") 
             else rp.integ_line(linefunc,hafunc, qgals,metal,Nlyc,nlines, lname=linename,
             all_lines=all_lines))

    for j in range(nlines):
      jname = linesinfo['Linename'][j]
      lum_line[jname][i] = lum[j]

  return lum_line



def get_emlines(linename,sfr,metals,q0=2.8e7,g0=-1.3,all_lines=False):

  if np.isnan(metals).any() or np.isnan(sfr).any():
    raise ValueError('\nget_emlines(): Input Zgas or sfr contain(s) NaN')

  ngal = len(sfr) if hasattr(sfr, "__len__") else 1

  LineProps = {'q0':float(q0), 'g0':float(g0),'linename':linename}
  print 'computing line luminosities...'
  ll = get_lumlines(sfr,metals,LineProps,all_lines=all_lines)
  lumline = np.array(ll)

  return lumline
