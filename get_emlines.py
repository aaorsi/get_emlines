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


def get_lumlines(sfr, metal, LineProps,all_lines=False):
  """
  LineProps contains Linesinfo, LinesArr, q(z) relation, line name and flux limit
  GalArr is the galaxy data dict.
  """
  linesinfo, linesarr = read_photoion()

  linename = LineProps['linename']

  ngals = len(sfr) if hasattr(sfr,"__len__") else 1
  qgals = qZrelation(metal, q0 = LineProps['q0'],g0 = LineProps['g0'])

  Nlyc = np.log10(1.35) + np.log10(sfr) + 53.0

  print 'Computing emission lines for %d galaxies\n' % (ngals)
  lum_line = [] # np.zeros(ngals)
  for i in range(ngals):
    lum_line.append(integ_line(linesinfo,linesarr,qgals[i],metal[i],
             Nlyc[i],lname=linename,all_lines=all_lines) if hasattr(qgals,"__len__") else integ_line(
             linesinfo,linesarr,qgals,metal,
             Nlyc,lname=linename,all_lines=all_lines))

  return lum_line           



def get_emlines(linename,sfr,metals,q0=2.8e7,g0=-1.3,all_lines=False):


  ngal = len(sfr) if hasattr(sfr, "__len__") else 1

  LineProps = {'q0':float(q0), 'g0':float(g0),'linename':linename}
  print 'computing line luminosities...'
  ll = get_lumlines(sfr,metals,LineProps,all_lines=all_lines)
  lumline = np.asarray(ll)

  return lumline
