#############################################################################
# Farhad's code

from imports import *
from lss.mod import WaveOperator

#############################################################################

sz = Sampling(201,16.0,0.0)
sx = Sampling(202,16.0,0.0)
st = Sampling(203,0.0012,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
xs,zs = dx*(nx/2),0.0 # source location
xr,zr = rampfloat(0.0,dx,nx),fillfloat(0.0,nx)
fpeak = 10.0 # Ricker wavelet peak frequency
#nabsorb = 12 # absorbing boundary size

def main(args):
  #adjointTest()
  SimplePlot.asPixels(zerofloat(101,102))
  print 'done'

def adjointTest():
  v = filldouble(4000.0,nz,nx)
  #add(mul(add(randdouble(Random(0),nz,nx),-0.5),100.0),v,v)
  random = Random(01234)
  #random = Random()
  sa = sub(randdouble(random,nz,nx,nt),0.0) # wavefield source
  sb = sub(randdouble(random,nz,nx,nt),0.0) # wavefield source
  print sum(sa)
  print sum(sb)
  ua = zerodouble(nz,nx,nt)
  ub = zerodouble(nz,nx,nt)
  va = copy(sa)
  vb = copy(sb)
  wave = FarhadWavefield(sz,sx,st)
  wave.setAdjoint(True)
  sw = Stopwatch(); sw.start()
  wave.modelAcousticWavefield(
    FarhadWavefield.WavefieldSource(sa),v,ua)
  wave.modelAcousticWavefield(
    FarhadWavefield.WavefieldAdjointSource(sb),v,ub)
  sw.stop(); print 'time:',sw.time()
  sum1 = dot(ua,vb)
  sum2 = dot(ub,va)
  print "adjoint test:",FarhadWavefield.compareDigits(sum1,sum2)
  print sum1
  print sum2

def dot(u,a):
  return FarhadWavefield.dot(u,a)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
