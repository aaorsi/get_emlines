import numpy as np
import os.path
from read_photoion import *
import pylab as pl
import matplotlib.gridspec as gridspec
from matplotlib.mlab import griddata
import scipy.stats

def qZrelation(Zgas, q0 = 2.8e7, g0 = -1.3):
# Default parameters from Orsi+14.  
  return q0 * (Zgas/0.012)**g0


def get_emlines(sfr, metal, LineProps):
  """
  LineProps contains Linesinfo, LinesArr, q(z) relation, line name and flux limit
  GalArr is the galaxy data dict.
  """
  linesinfo, linesarr = read_photoion()

  linename = LineProps['linename']

  ngals = len(sfr)
  qgals = qZrelation(metal, LineProps['q0'],LineProps['g0'])

  Nlyc = np.log10(1.35) + np.log10(sfr) + 53.0

  print 'Computing emission lines for %d galaxies\n' % (ngals)
  lum_line = np.zeros(ngals)
  for i in range(ngals):
    lum_line[i] = integ_line(linesinfo,linesarr,qgals[i],metal[i],
             Nlyc[i],lname=linename)

  return lum_line           



def read_eaglelines(linename,sfr,metals,q0=2.8e7,g0=-1.3):


  ngal = len(sfr)
  LineProps = {'q0':float(q0), 'g0':float(g0),'linename':linename}
  print 'computing line luminosities...'
  ll = get_emlines(sfr,metals,LineProps)
  lumline = np.asarray(ll)

  return lumline
