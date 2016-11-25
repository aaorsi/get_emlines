#import h5py
import numpy as np
import os.path
from read_photoion import *
import pylab as pl
import matplotlib.gridspec as gridspec
from matplotlib.mlab import griddata
import scipy.stats

def make_grid(xdata,ydata,nbins):
  finx = np.where(np.isfinite(xdata))
  mindm = np.min(xdata[finx])
  maxdm = np.max(xdata[finx])
  
  finy = np.where(np.isfinite(ydata))
  minl = np.min(ydata[finy])
  maxl = np.max(ydata[finy])
  
  grid_x, grid_y = np.mgrid[mindm:maxdm:nbins, minl:maxl:nbins]
 
  grid_x = np.linspace(mindm,maxdm,num=nbins)
  grid_y = np.linspace(minl,maxl,num=nbins)
  xbin = grid_x[1] - grid_x[0]
  ybin = grid_y[1] - grid_y[0]

  griddata = np.zeros([nbins,nbins])

  for ix in range(nbins):
    for iy in range(nbins):
      cc = np.where( (xdata > grid_x[ix] - xbin/2.) & 
                     (xdata < grid_x[ix] + xbin/2.) &
                     (ydata > grid_y[iy] - ybin/2.) & 
                     (ydata < grid_y[iy] + ybin/2.))
      griddata[ix,iy] = len(cc[0])


  return grid_x, grid_y, griddata

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



def read_eaglelines(linename,snapshot,delta=False,ivbin=800,nvol=8000):

  outfile = '../out/'
  q0 = '2.8e7'
  g0 = '-1.3'
  #linename = 'Halpha'
  gtype = 'satellites'
  msrange = ['9.75','10.25']
  add_msrange = False 

  lumfile = outfile + linename + '_'+q0+'_'+g0+'.hdf5'

#  datadir = '/home/CEFCA/jchaves/Leiden/data/'
  datadir = '/home/CEFCA/jchaves/elg/data/'

  cat     = str(snapshot)+'.dat'
  datacat = np.loadtxt(datadir+cat,delimiter=',',usecols=(10,13,17),skiprows=31)
  ngal = len(datacat[:,0])
  sfr      = datacat[:,0]
  mstellar = datacat[:,1]
  metals   = datacat[:,2]


#  catfile = 'catCRalvaro.hdf5'
#  sfrfile = 'SFRalvaro.hdf5'

  deltadir = '/home/CEFCA/aaorsi/work/eagle_sham/out/delta/'
  

#  datacat = h5py.File(datadir + catfile,'r')
#  sfrcat  = h5py.File(datadir + sfrfile,'r')

  LineProps = {'q0':float(q0), 'g0':float(g0),'linename':linename}

#  ngal = len(sfrcat['SFR030kpc'])
#  metals = datacat['Metallicity'][0,0:]
#  sfr    = sfrcat['SFR030kpc'][0:]
#  mstellar = datacat['STMass30kpc']


#  if os.path.isfile(lumfile) is False: 
  print 'computing line luminosities...'
  ll = get_emlines(sfr,metals,LineProps)
#    lfile = h5py.File(lumfile,'w')
#    dset = lfile.create_dataset('lumline',(ngal,),dtype='f')
#    dset[0:ngal] = ll
#    lfile.close()
#  else:
#    lfile = h5py.File(lumfile,'r')
  #ll = lfile['lumline']
#
  
  import pdb ; pdb.set_trace()  
  lumline = np.asarray(ll)

  if delta:
    print 'reading overdensity files'
    nfiles = nvol / ivbin
    j = 0

    for i in range(nfiles):
      fn_i = "deltag_%d_%d.txt" % (i*ivbin, (i+1)*ivbin-1)
      print fn_i
      ifile = deltadir + fn_i
      
      f = open(ifile,'r')
      data_i = f.readlines() 
      f.close()
      ndata = len(data_i)
      data_fl = np.zeros(ndata)
      if i == 0:    
        strarr = data_i[0].strip().split()
        if len(strarr) == 1:
          idg = int(strarr[0][0])
          k = 0
          strarr = data_i[k+1].strip().split()
          while len(strarr) == 2:
            k += 1
            strarr = data_i[k+1].strip().split()
          nbins = k
        
        deltag = np.zeros([ngal,nbins])
      
      radg   = np.zeros(nbins)
      ngalsfile = ndata / (nbins + 1)
      kk = 0
      for ig in range(ngalsfile):
        id0 = (nbins+1)*ig
        strarr = data_i[id0].strip().split()
        if len(strarr) >1:
          print 'problem reading file. check'
          import pdb; pdb.set_trace()

        idg = int(strarr[0])
        k = 0
        id0 += 1
        while (k < nbins):
          strarr = data_i[k+id0].strip().split()
          if ig == 0:
            radg[k] = float(strarr[0])
          deltag[idg-1,k] = float(strarr[1])
          k += 1
          
    deltadata = {'deltag':deltag,'r':radg}


  if delta:
    return datacat, sfrcat, lumline, deltadata
  else:  
    return datacat, sfrcat, lumline

def get_nm(idhalos, datamhalos,mass_arr):
  j = 0
  mbin = mass_arr[1] - mass_arr[0]
  nghalo = np.zeros(len(mass_arr))
  ndata = len(idhalos)

  while (j < ndata):
    
    idj = idhalos[j]
    mj  = datamhalos[idj]
    j0 = j
    while (idj == idhalos[j]):
      j   += 1
      if j == ndata:
        break
    
    if ((mj > mass_arr[0]-mbin/2.) &
       (mj < mass_arr[-1] + mbin/2.)):
      im = np.round((mj - mass_arr[0])/mbin)
      nghalo[im] += j - j0
  
  return nghalo

