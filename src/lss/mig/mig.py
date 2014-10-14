#############################################################################
# Migration demo
#############################################################################

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util.ArrayMath import *
from lss.mod import *
from lss.mig import *
from lss.util import *
from fah import FakeData
from jarray import zeros

#############################################################################
# Parameters

# Might want to change these
ds = 10 # shot spacing in samples
niter = 2 # number of iterations when computing the image
npmax = 18 # max number of shots migrated in parallel (reduce to save memory)

# Probably should not change these
tmax = 2.5 # maximum recording time
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
ns = 1+101/ds # number of shots
np = min(npmax,ns) # number of shots migrated in parallel
print 'nshot =',ns
print 'nparallel =',np

#############################################################################

timer = Timer()
def main(args):
  s = getSlowness()
  #r = goLsMigration(s)
  pixels(transpose(s),cmap=jet)
  #pixels(transpose(r),cmap=gray,perc=98.0)

def goLsMigration(s):

  # Allocate wavefields and sources/receivers
  timer.start('setup')
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)
  src,rco = getSourcesAndReceivers()
  timer.stop('setup')

  # Compute observed data
  wave = WaveOperatorS(s,dx,dt,nabsorb,u,a)
  timer.start('observed data')
  wave.applyForward(src,rco) # probably should use different wavelet
  timer.stop('observed data')

  # Compute the image
  born = BornOperatorS(s,dx,dt,nabsorb,u,a)
  ref = RecursiveExponentialFilter(8.0)
  mp = getMask()
  bs = BornSolver(born,src,rco,ref,mp)
  r = bs.solve(niter)
  return r

def getMask():
  m = zerofloat(102,101)
  ref = RecursiveGaussianFilter(24.0)
  for ix in range(102):
    m[0][ix] = 1.0
  ref.applyX0(m,m)
  mul(1.0/max(m),m,m)
  sub(1.0,m,m)
  """
  for ix in range(102):
    for iz in range(4):
      m[iz][ix] = 0.0 # mask water layer
  """
  #pixels(m)
  return m

def getSourcesAndReceivers():
  src = zeros(ns,Source) # source
  rco = zeros(ns,Receiver) # observed data
  xr = rampint(0,1,102)
  zr = zeroint(102)
  for isou in range(ns):
    xs = isou*ds
    src[isou] = Source.RickerSource(xs,0,dt,fpeak)
    rco[isou] = Receiver(xr,zr,nt)
  return src,rco

def getSlowness():
  f = FakeData.seismicAndSlopes3d2014A('OA',3,False,False,True,False,0.0)
  r = f[0][50]
  sub(r,min(r),r)
  mul(1.0/max(r),r,r)
  mul(2.0,r,r)
  add(2.0,r,r)
  for ix in range(102):
    for iz in range(4):
      r[ix][iz] = 2.0 # add water layer
  div(1.0,r,r) # convert to slowness
  return transpose(r)

#############################################################################
# Plotting

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def pixels(x,cmap=gray,perc=100.0,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  cb = sp.addColorBar()
  cb.setWidthMinimum(130)
  sp.setSize(750,600)
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
