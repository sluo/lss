#############################################################################
# Migration demo

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util.ArrayMath import *
from lss.mod import *
from fah import FakeData
from jarray import zeros

#############################################################################

tmax = 2.5 # maximum time
fpeak = 10.0 # peak frequency of Ricker wavelet
nabsorb = 20 # number of samples in absorbing boundary
dx = 0.020 # spatial sampling interval in x direction
dz = 0.020 # spatial sampling interval in z direction
dt = 0.002 # time sampling interval
#nx = 102 # number of samples in x direction
#nz = 101 # number of samples in z direction
nt = 1+int(tmax/dt) # number of time samples
nxp = 102+2*nabsorb # total number of samples in x direction
nzp = 101+2*nabsorb # total number of samples in z direction
niter 1 = # number of iterations when computing the image
ds = 10 # shot spacing in samples
ns = 1+101/ds # number of shots
np = min(16,ns) # number of shots to migrate in parallel
print 'nshot =',ns
print 'niter =',niter

def getSlowness():
  f = FakeData.seismicAndSlopes3d2014A('OA',3,False,False,True,False,0.0)
  r = f[0][50]
  sub(r,min(r),r)
  mul(1.0/max(r),r,r)
  mul(2.0,r,r)
  add(2.0,r,r)
  div(1.0,r,r) # convert to slowness
  return transpose(r)

timer = Timer()
def main(args):
  s = getSlowness()
  #checkStability(s)
  #goModeling(s)
  #r = goMigration(s)
  r = goMig(s)
  #pixels(transpose(s),cmap=jet)
  #pixels(transpose(r),cmap=gray)

def goMig(s):

  # Allocate wavefields and sources/receivers
  timer.start('setup')
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)
  src,rco = getSourcesAndReceivers()
  timer.stop('setup')

  # Compute observed data
  wave = WaveOperatorS(s,dx,dt,nabsorb,u,a)
  timer.start('observed data')
  wave.applyForward(src,rco)
  timer.stop('observed data')

  # Compute the image
  born = BornOperatorS(s,dx,dt,nabsorb,u,a)
  ref = RecursiveExponentialFilter(4.0)
  bs = BornSolver(born,src,rcp,rco,ref)
  r = bs.solve(niter)
  
  return r

def getSourcesAndReceivers():
  src = zeros(ns,Source) # source
  rco = zeros(ns,Receiver) # observed data
  xr = rampint(0,1,102)
  zr = zeroint(102)
  ns = 1+100/ds # nx=102
  for isou in range(ns):
    xs = isou*ds
    src[isou] = Source.RickerSource(xs,0,dt,fpeak)
    rco[isou] = Receiver(xr,zr,nt)
  return src,rco

def goMigration(s):
  nabsorb = 20
  nx,nz,nt = len(s[0]),len(s),1+int(tmax/dt)
  u = zerofloat(nx+2*nabsorb,nz+2*nabsorb,nt)
  a = zerofloat(nx+2*nabsorb,nz+2*nabsorb,nt)
  wave = WaveOperator(s,dx,dt,nabsorb)
  xs,zs = (nx/2)*dx,0.0
  src = Source.RickerSource(nx/2,0,dt,fpeak)
  rec = Receiver(rampint(0,1,nx),fillint(0,nx),nt)
  wave.applyForward(src,rec,u)
  wave.applyAdjoint(Source.ReceiverSource(rec),a)
  r = WaveOperator.collapse(u,a,nabsorb)
  pixels(r,perc=99.5)

def goModeling(s):
  nx,nz,nt = len(s[0]),len(s),1+int(tmax/dt)
  wave = WaveOperator(s,dx,dt,20)
  src = Source.RickerSource(nx/2,0,dt,fpeak)
  rec = Receiver(rampint(0,1,nx),fillint(0,nx),nt)
  wave.applyForward(src,rec)
  data = rec.getData()
  print sum(data)
  pixels(data,perc=99.5)

def goModelingSinc(s):
  nabsorb = 20
  nx,nz,nt = len(s[0]),len(s),1+int(tmax/dt)
  wave = WaveOperator(s,dx,dt,nabsorb)
  #xs,zs = (nx/2)*dx,0.0
  xs,zs = (nx/2)*dx+dx/2,0.010
  print xs
  print zs
  src = SincSource.RickerWavelet(xs,zs,fpeak)
  src.setupForDomain(nx,dx,nz,dz,nt,dt,nabsorb)
  #rec = Receiver(rampint(0,1,nx),fillint(0,nx),nt)
  rec = SincReceiver(xs+1.123,zs+1.123,nt)
  rec.setupForDomain(nx,dx,nz,dz,nt,dt,nabsorb)
  wave.applyForward(src,rec)
  data = rec.getData()
  print sum(data)
  #pixels(data,perc=99.5)
  SimplePlot.asPoints(data[0])

def checkStability(s):
  vmax = 1.0/min(s)
  vmin = 1.0/max(s)
  ppw = vmin/(dx*fpeak) # points per wavelength
  cfl = dx/(dt*vmax) # cfl
  print 'vmin =',vmin
  print 'vmax =',vmax
  print 'ppw =',ppw
  print 'cfl =',cfl

import socket
def getMarmousi(stride=3):
  """Marmousi model.
  Parameters:
    sigma - half-width of smoothing window for scale separation
    econst - constant error in background slowness s0
    egauss - gaussian error in background slowness s0
    erand - random error in background slowness s0
  Returns:
    tt - true model
    t0 - true background slowness
    t1 - true reflectivity
    s0 - true background slowness
    s1 - true reflectivity
  """
  p = zerofloat(751,2301)
  if socket.gethostname()=='backus.Mines.EDU':
    readImage("/data/sluo/marmousi/marmousi.dat",p)
  else:
    readImage("/data/seis/marmousi/marmousi.dat",p)
  p = copy(743,2301,8,0,p)
  div(1000.0,p,p) # slowness (s/km)
  wa = iceil(0.2/(0.004*stride))
  nz = iceil(743.0/stride)+wa
  nx = iceil(2301.0/stride)
  #print "wa =",wa
  #print "nz =",nz
  #print "nx =",nx
  t = fillfloat(2.0/3.0,nz,nx)
  copy(nz-wa,nx,0,0,stride,stride,p,wa,0,1,1,t)
  return transpose(t)

def iceil(x):
  return int(ceil(x))

def readImage(name,image=None):
  if not image:
    image = zerofloat(nz,nx)
  fileName = name
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(fname,image):
  aos = ArrayOutputStream(fname)
  aos.writeFloats(image)
  aos.close()

#############################################################################
# plotting

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def pixels(x,cmap=gray,perc=100.0,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  sp.setSize(600,600)
  if title:
    sp.addTitle(title)
  pv = sp.addPixels(x)
  pv.setColorModel(cmap)
  if perc<100.0:
    pv.setPercentiles(100.0-perc,perc)

#############################################################################
# Do everything on Swing thread.

import sys
from java.lang import Runnable
from javax.swing import SwingUtilities
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
