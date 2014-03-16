#############################################################################
# Convert Eni data from segy to dat

from imports import *

# Samplings
ss = Sampling(3661,0.0125,1.0) # shot
sr = Sampling(99,0.0125,0.0) # receiver (correct: dx=-0.0125, fx=0.0)
st = Sampling(2501,0.002,0.0) # time
sz = Sampling(401,0.010,0.0) # depth
sx = Sampling(7416,0.00625,0.0) # offset
#ns,ds,fs = ss.count,ss.delta,ss.first
#nr,dr,fr = sr.count,sr.delta,sr.first
#nt,dt,ft = st.count,st.delta,st.first
#nz,dz,fz = sz.count,sz.delta,sz.first
#nx,dx,fx = sx.count,sx.delta,sx.first
ns,nr,nt,nz,nx = ss.count,sr.count,st.count,sz.count,sx.count
sgyDir = '/data/sluo/eni/segy/'

datDir = None
datDir = '/data/sluo/eni/dat/'

#############################################################################

def main(args):
  #showData()
  #showVelocity()
  #writeSubsetA()
  #writeSubsetB()
  #writeSubsetC()
  #writeSubsetD()
  #writeSubsetE()
  #writeSubsetF()
  writeSubsetG()

def showVelocity():
  segy = SegyImage(sgyDir+'velo.segy')
  #segy.printSummaryInfo()
  segy.printAllInfo()
  segy.setFormat(1) # IBM format
  n1 = segy.getN1(); print 'n1=%d'%n1
  n2 = segy.getN2(); print 'n2=%d'%n2
  v = zerofloat(n1,n2)
  for i2 in range(n2):
    segy.getTrace(i2,v[i2])
  plot(v,cmap=jet)
  #SimplePlot.asPixels(sz,sx,v)
    
def showData():
  segy = SegyImage(sgyDir+'shot.segy')
  #segy.printSummaryInfo()
  segy.printAllInfo()
  n1 = segy.getN1(); print 'n1=%d'%n1
  #n2 = segy.getN2(); print 'n2=%d'%n2
  n2 = 10*99 # show some gathers
  d = zerofloat(n1,n2)
  for i2 in range(n2):
    segy.getTrace(i2,d[i2])
  plot(d,cmap=gray,sperc=98.0)

def writeSubsetA():
  fix = 1600
  nix = 1101
  #niz = 161
  niz = 181
  #t = 2.2
  tt = 1.5
  dt = 0.0004
  ddir = 'suba/'
  writeSubset(fix,nix,niz,tt,dt,ddir)

def writeSubsetB():
  fix = 3520
  nix = 1101
  #niz = 161
  niz = 181
  #tt = 2.2
  tt = 1.5
  dt = 0.0004
  ddir = 'subb/'
  writeSubset(fix,nix,niz,tt,dt,ddir)

def writeSubsetC():
  #fix = 1120
  fix = 1200
  nix = 3201
  niz = 181
  #tt = 1.8
  tt = 1.5
  dt = 0.0004
  ddir = 'subc/'
  writeSubset(fix,nix,niz,tt,dt,ddir)

def writeSubsetD():
  fix = 2080
  nix = 1121
  niz = 181
  tt = 1.5
  dt = 0.0004
  ddir = 'subd/'
  writeSubset(fix,nix,niz,tt,dt,ddir)

def writeSubsetE():
  fix = 1600
  nix = 1121
  niz = 181
  tt = 1.5
  dt = 0.0004
  ddir = 'sube/'
  writeSubset(fix,nix,niz,tt,dt,ddir)

def writeSubsetF():
  fix = 1120
  nix = 2241
  niz = 181
  tt = 1.5
  dt = 0.0004
  ddir = 'subf/'
  writeSubset(fix,nix,niz,tt,dt,ddir)

def writeSubsetG():
  fix = 1600
  nix = 1921
  niz = 181
  tt = 1.5
  dt = 0.0004
  ddir = 'subg/'
  writeSubset(fix,nix,niz,tt,dt,ddir)

def writeSubset(fix,nix,niz,tt,dt,ddir):
  """
  Writes a velocity and data subset.
  Parameters:
    fix - first x index
    nix - number of samples in x
    niz - number of samples in z (for dz = 0.00625)
    tt - total time
    dt - time resampling interval
    ddir - dat directory
  """
  Check.argument(fix%2==0,'fix%2==0')
  writeVelocitySubset(fix,nix,niz,ddir)
  writeDataSubset(fix,nix,tt,dt,ddir)

def writeVelocitySubset(fix,nix,niz,ddir):
  """
  Writes a velocity subset interpolated onto a uniformly sampled
  grid in both x and z (dz = dx = 0.00625 km).
  """
  segy = SegyImage(sgyDir+'velo.segy')
  li = LinearInterpolator()
  li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT)
  li.setUniformSampling(401,0.010,0.0)
  v = zerofloat(niz,nix)
  for ix in range(fix,fix+nix):
    vi = zerofloat(401)
    segy.getTrace(ix,vi)
    li.setUniformSamples(vi)
    for iz in range(niz):
      v[ix-fix][iz] = li.interpolate(iz*0.00625)
  plot(v,cmap=jet)
  if datDir is not None:
    write(datDir+ddir+'v.dat',v)
  print 'nx =',nix
  print 'nz =',niz

def writeDataSubset(fix,nix,tt,dt,ddir):
  """
  Writes data subsets interpolated to dr = 0.00625 (half original spacing)
  """
  #fsou = int((0.225+fix*0.00625)/0.0125)
  #lsou = int(((fix+nix-1)*0.00625-1.0)/0.0125)
  fsou = fix/2+18
  lsou = (fix+nix-1)/2-80
  nsou = 1+lsou-fsou
  nt = 1+int(tt/dt)
  print 'nt =',nt
  print 'ns =',nsou
  print 'fs =',fsou
  print 'ls =',lsou
  si = SincInterp()
  #for isou in range(fsou,lsou+1):
  #for isou in range(fsou,fsou+nsou/2):
  for isou in range(fsou+nsou/2,lsou+1):
    print '  isou = %d'%isou
    d = getGather(isou)
    e = zerofloat(nt,2*nr-1)
    #fname = 'd_'+str(isou+722)+'.dat'
    fname = 'd_'+str(isou-fsou)+'.dat'
    for ir in range(2*nr-1):
      r = 0.5*ir*0.0125
      for it in range(nt):
        t = it*0.0004
        e[ir][it] = si.interpolate(st,sr,d,t,r)
    if datDir is not None:
      write(datDir+ddir+fname,e)
  print 'done!'

def getGather(isou,d=None):
  segy = SegyImage(sgyDir+'shot.segy')
  fr = isou*nr # first trace
  if not d:
    d = zerofloat(nt,nr)
  for ir in range(nr):
    #segy.getTrace(fr+ir,d[ir])
    segy.getTrace(fr+nr-1-ir,d[ir]) # reverse order
  return d

def write(fname,image):
  aos = ArrayOutputStream(fname)
  aos.writeFloats(image)
  aos.close()

def read(name,image):
  fileName = name
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  

#############################################################################
# plotting

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,cmap=gray,perc=100.0,sperc=None,cmin=0.0,cmax=0.0,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  cb = sp.addColorBar()
  cb.setWidthMinimum(100)
  sp.setSize(1010,740)
  if title:
    sp.addTitle(title)
  pv = sp.addPixels(x)
  pv.setColorModel(cmap)
  if perc<100.0:
    pv.setPercentiles(100.0-perc,perc)
  if sperc is not None: # symmetric percentile clip (for plotting gradients)
    clips = Clips(100-sperc,sperc,x)
    clip = max(abs(clips.getClipMin()),abs(clips.getClipMax()))
    pv.setClips(-clip,clip)
  if cmin<cmax:
    pv.setClips(cmin,cmax)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
