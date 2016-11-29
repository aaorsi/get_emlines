import numpy as np
import sys

def read_photoion(debug=0,MappingsModel='Levesque10'):
  LineDataDir     = 'HIImodels/Lines/'
  LineData_file   = LineDataDir + 'LineData_'+MappingsModel+'dens_1e1'
  LineInfo_file   = LineDataDir + 'LineInfo_'+MappingsModel
  l_info = np.genfromtxt(LineInfo_file,dtype=[('id','i8'),('Linename','S10'),('lambda0','f8')],skip_header=2)
  print 'Lines available:',l_info['Linename']

  fd = open(LineData_file,'rb')
  nlines = np.fromfile(fd,dtype='i4',count=1)

  if (nlines != len(l_info['id'])):
    print 'Warning read_photoion.py: check number of lines read'
    print 'nlines=',nlines,' len(l_info):',len(l_info['id'])
  
  lambda0 = np.fromfile(fd,dtype='f4',count=nlines)
  nz      = np.fromfile(fd,dtype='i4',count=1)
  nq      = np.fromfile(fd,dtype='i4',count=1)
  nage    = np.fromfile(fd,dtype='i4',count=1)
  
#  if debug:
  print 'checking lambda0:',lambda0,l_info['lambda0']
  print 'nz, nq, nage,',nz,nq,nage

  ZArray  = np.fromfile(fd,dtype='f4',count=nz)
  QArray  = np.fromfile(fd,dtype='f4',count=nq)
  AgeArray= np.fromfile(fd,dtype='f4',count=nage)

  LinesArr = np.fromfile(fd,dtype='f8',count=nz*nq*nage*nlines)

  Linesinfo = {'Linename':l_info['Linename'], 'lambda0':l_info['lambda0'],
  'nlines':nlines,'ZArray':ZArray,'QArray':QArray,'AgeArray':AgeArray,
  'nz':nz,'nq':nq,'nage':nage}
# return l_info,nlines,LinesArr,ZArray,QArray,AgeArray,nz,nq,nage
  return Linesinfo,LinesArr


def calc_emlines(Linesinfo,LinesArr,qgas,zgas,lname='Halpha',all_lines=False):

  """
  This function computes the emission line luminosity of line lname as a function
  of ionization parameter, qgas, and metallicity, zgas.

  arguments:
    qgas: ionisation parameter, [cm/s]
    zgas: metallicity of the gas
    lname: String containing the name of the line.
    all_lines : returns all line fluxes
  """

  ZArray = Linesinfo['ZArray']
  QArray = Linesinfo['QArray']

  nz = Linesinfo['nz']
  nq = Linesinfo['nq']

    
  nlines = len(Linesinfo['Linename'])

  if all_lines:
    idl = range(nlines)
  else:
    idl = np.where(Linesinfo['Linename'] == lname)
    idl = idl[0]

  if len(idl) == 0:
    print 'calc_emlines: line name ',lname,' not found/recognised'
    print 'possible linenames are: ',Linesinfo['Linename']
    sys.exit()

  i = 0

  zgas = ZArray[0] if zgas < ZArray[0] else ZArray[-1] if zgas > ZArray[-1] else zgas

  if (zgas > 1 or zgas < 1e-4):
    line = 1e-30
    return line
  
  if (zgas <= ZArray[0]):
    iz0 = 0
    iz1 = 1
    z0  = ZArray[iz0]
    z1  = ZArray[iz1]
  else:
    i = 0
    while (zgas >= ZArray[i] and i < nz-1):
      z0 = ZArray[i]
      iz0= i
      i +=1

    iz1 = i
    z1 = ZArray[iz1]

  if (qgas <= QArray[0]):
    iq0 = 0
    iq1 = 1
    q0 = QArray[0]
    q1 = QArray[1]
  else:
    i = 0
    while(qgas >= QArray[i] and i < nq-1):
      q0 = QArray[i]
      iq0 = i
      i += 1
    
    iq1 = i
    q1 = QArray[iq1]

  if all_lines:
    idlarr = range(nlines)
    line = np.zeros(nlines)
    for idl in idlarr:
      f00 = LinesArr[idl + nlines*iz0 + nz*nlines*iq0]
      f10 = LinesArr[idl + nlines*iz1 + nz*nlines*iq0]
      f01 = LinesArr[idl + nlines*iz0 + nz*nlines*iq1]
      f11 = LinesArr[idl + nlines*iz1 + nz*nlines*iq1]

      den = (z1 - z0)*(q1 - q0)
      t1 = f00 * (z1 - zgas)*(q1-qgas)/den
      t2 = f10 * (zgas - z0)*(q1-qgas)/den
      t3 = f01 * (z1 - zgas)*(qgas-q0)/den
      t4 = f11 * (zgas - z0)*(qgas-q0)/den

      line_ = t1+t2+t3+t4
      if line_ < 0:
        line_ = 1e-30
      line[idl] = line_
  else:
    f00 = LinesArr[idl + nlines*iz0 + nz*nlines*iq0]
    f10 = LinesArr[idl + nlines*iz1 + nz*nlines*iq0]
    f01 = LinesArr[idl + nlines*iz0 + nz*nlines*iq1]
    f11 = LinesArr[idl + nlines*iz1 + nz*nlines*iq1]

    den = (z1 - z0)*(q1 - q0)
    t1 = f00 * (z1 - zgas)*(q1-qgas)/den
    t2 = f10 * (zgas - z0)*(q1-qgas)/den
    t3 = f01 * (z1 - zgas)*(qgas-q0)/den
    t4 = f11 * (zgas - z0)*(qgas-q0)/den

    line = t1+t2+t3+t4
    
    
  return line


def integ_line(lineinfo,LinesArr,qgas,zgas,nlyc,lname='Halpha',all_lines=False):

  qmin = 1.0e7
  qmax = 1.0e8

  qgas = qmin if qgas < qmin else qmax if qgas > qmax else qgas


  frachyda = calc_emlines(lineinfo,LinesArr,qgas,zgas,'Halpha')
  alpha = np.log10(1.37) - 12
  if all_lines: 
    cel      = calc_emlines(lineinfo,LinesArr,qgas,zgas,all_lines=True)
  else:
    cel      = calc_emlines(lineinfo,LinesArr,qgas,zgas,lname=lname)

  nlines = len(lineinfo['Linename'])
  
  if frachyda == 1e-30:
    if all_lines == True:
      line = np.zeros(nlines) + 1e-30
    else:
      line = 1e-30
  
  else:
    line = nlyc + alpha + np.log10(cel) - np.log10(frachyda)

  if np.isnan(line):
    print 'line is nan!'
    line = -999 

#  if line[0] == line[1]:
  return line





#if __name__ == "__read_photoion__":
#  read_photoion()

