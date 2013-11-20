##############################################################################
# Sinc interpolation and accumulation

from imports import *
from dnp import InverseInterpolator

##############################################################################

nt = 101
st = Sampling(nt)

savDir = None
#savDir = '/home/sluo/Desktop/pngdat/'
#savDir = '/home/sluo/Desktop/pngdat2/'
#savDir = '/home/sluo/Desktop/pngdat3/'

##############################################################################

def main(args):
  goInterpolateAndAccumulate()

def points(f,cmin=0.0,cmax=0.0):
  sp = SimplePlot()
  if cmin<cmax:
    sp.setVLimits(cmin,cmax)
  pv = sp.addPoints(f)

def goInterpolateAndAccumulate():

  si = SincInterp()
  #si = SincInterp.fromErrorAndFrequency(0.001,0.499)

  u = makeShifts() # shifts
  SimplePlot.asPoints(u)
  s = add(rampfloat(0.0,1.0,nt),u) # shifted coordinates
  t = zerofloat(nt); InverseInterpolator(nt,nt).invert(s,t) # inverse
  SimplePlot.asPoints(t)

  f = zerofloat(nt,nt) # input
  for it in range(nt):
    f[it][it] = 1.0
  RecursiveGaussianFilter(nt/50.0).apply00(f,f); mul(1.0/max(f),f,f)
  g = zerofloat(nt,nt) # interpolation
  h = zerofloat(nt,nt) # inverse interpolation
  p = zerofloat(nt,nt) # accumulation
  class Loop(Parallel.LoopInt):
    #def xcompute(self,it):
    #  f[it][it] = 1.0
    #  for jt in range(nt):
    #    g[it][jt] = si.interpolate(st,f[it],jt+u[jt])
    #    #h[it][jt] = si.accumulate(st,f[it],jt+u[jt])
    def compute(self,it):
      si.interpolate(nt,1.0,0.0,f[it],nt,s,g[it])
      si.interpolate(nt,1.0,0.0,f[it],nt,t,h[it])
      si.accumulate(nt,s,f[it],nt,1.0,0.0,p[it])
      #si.interpolate(nt,1.0,0.0,g[it],nt,t,h[it])
      #si.accumulate(nt,s,g[it],nt,1.0,0.0,p[it])
  Parallel.loop(nt,Loop())
  pixels(f,cmap=jet,cmin=-0.2*max(f),cmax=1.0*max(f),title='input')
  pixels(g,cmap=jet,cmin=-0.2*max(g),cmax=1.0*max(g),title='interpolate')
  pixels(h,cmap=jet,cmin=-0.2*max(h),cmax=1.0*max(h),title='inverse')
  pixels(p,cmap=jet,cmin=-0.2*max(p),cmax=1.0*max(p),title='accumulate')

  # Compare to identity matrix.
  #r = zerofloat(nt,nt)
  #for it in range(nt):
  #  r[it][it] = 1.0
  #pixels(sub(r,mmul(g,h)),cmap=jet,cmin=-1.0,cmax=1.0,
  #  title='(identity) - (interpolation) x (inverse)')
  #pixels(sub(r,mmul(g,p)),cmap=jet,cmin=-1.0,cmax=1.0,
  #  title='(identity) - (interpolation) x (accumulation)')

  # Compare inverse interpolation and accumulation.
  d = sub(mmul(g,h),mmul(g,p))
  pixels(d,cmap=jet,cmin=-1.0*max(abs(d)),cmax=1.0*max(abs(d)),
    title='(interpolate) x (inverse) - (interpolate) x (accumulate)')

  # Dot-product test
  random = Random(0123)
  va = randfloat(random,nt)
  vb = randfloat(random,nt)
  print dot(vmul(g,va),vb)
  print dot(vmul(p,vb),va)

def dot(f,g):
  """Dot-product between two vectors."""
  s = 0.0
  for it in range(nt):
    s += f[it]*g[it]
  return s

def mmul(f,g):
  """Matrix-matrix multiplication."""
  h = zerofloat(nt,nt)
  class Loop(Parallel.LoopInt):
    def compute(self,i):
      for j in range(nt):
        hij = 0.0
        for k in range(nt):
          #hij += f[i][k]*g[k][j]
          hij += f[k][i]*g[j][k]
        h[i][j] = hij
  Parallel.loop(nt,Loop())
  return h

def vmul(f,g):
  """Matrix-vector multiplication."""
  h = zerofloat(nt)
  class Loop(Parallel.LoopInt):
    def compute(self,i):
      for k in range(nt):
        h[i] += f[k][i]*g[k]
  Parallel.loop(nt,Loop())
  return h

def applyShifts(u,f,g=None):
  n = len(f)
  if g is None:
    g = zerofloat(n)
  si = SincInterp()
  for i in range(n):
    g[i] = si.interpolate(st,f,i+u[i])
  return g

def makeShifts():
  u = zerofloat(nt)
  nhp = 2.0 # number of half periods
  umax = 10.0 # max shift
  for it in range(nt):
    k = nhp*it/(nt-1)
    u[it] = umax*sin(k*FLT_PI)
  return u

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def pixels(f,cmap=gray,cmin=0.0,cmax=0.0,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(800,700)
  cb = sp.addColorBar()
  cb.setWidthMinimum(100)
  if title:
    sp.setTitle(title)
  pv = sp.addPixels(f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cmap)
  if cmin<cmax:
    pv.setClips(cmin,cmax)

##############################################################################
# Do everything on Swing thread.
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
if __name__=='__main__':
  SwingUtilities.invokeLater(RunMain())
