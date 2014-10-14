"""
Combined waveform and traveltime inversion
"""
from imports import *

#############################################################################
# Parameters

sz = Sampling(701,4.0,0.0)
sx = Sampling(401,4.0,0.0)
st = Sampling(1300,0.0005,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
kxs,kzs = [0],[nz/2]
#kxs,kzs = [0,0],[nz/3,2*nz/3]
#kxs,kzs = fillint(0,69),rampint(0,10,71) 
#kxs,kzs = fillint(0,29),rampint(0,25,29) 
#kxs,kzs = fillint(0,15),rampint(0,50,15) 
kxr,kzr = fillint(nx-1,nz),rampint(0,1,nz) 
ns,nr = len(kxs),len(kxr)
fpeak = 40.0 # Ricker wavelet peak frequency

pngDir = None
#pngDir = '/Users/sluo/Desktop/png/'
#############################################################################

parallel = False
def main(args):
  setModel()
  #showWeights()
  goWaveform()
  #goTraveltime()
  #goCombinedA()
  #goCombinedB()

def setModel():
  global v,c
  #v = getMarmousi()
  v = getGaussian(0.02,0.02)
  #v = getGaussian(0.20,0.20)
  #v = getGaussian(0.20,0.02)
  c = fillfloat(4000,nz,nx)
  plot(v,cmap=jet,title='true')
  plot(c,cmap=jet,cmin=min(v),cmax=max(v),title='initial')

class Modeler():
  def __init__(self,jzs,jxs,v,c,aw=None):
    self.jzs = jzs
    self.jxs = jxs
    self.aw = AcousticWavefield(sz,sx,st) if aw==None else aw
  def modelDataAndSourceWavefield(self):
    source = AcousticWavefield.RickerSource(fpeak,self.jzs,self.jxs)
    self.aw.forwardPropagate(source,v)
    do = self.aw.getWavefield(kzr,kxr) # observed data
    self.aw.forwardPropagate(source,c)
    ds = self.aw.getWavefield(kzr,kxr) # simulated data
    u = self.aw.getWavefield() # simulated wavefield
    return do,ds,u
  def modelReceiverWavefield(self,r):
    source = AcousticWavefield.AdjointSource(dt,kzr,kxr,r)
    self.aw.backPropagate(source,c)
    return self.aw.getWavefield()

#############################################################################
# Waveform

def goWaveform():
  if parallel:
    g = Parallel.reduce(ns,WaveformParallel())
  else:
    g = zerofloat(nz,nx)
    aw = AcousticWavefield(sz,sx,st)
    for isou in range(ns):
      print "isou =",isou
      add(computeWaveformGradient(isou,aw),g,g)
  div(g,max(abs(g)),g)
  plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='gradient (all shots)')

def computeWaveformGradient(isou=0,aw=None):
  modeler = Modeler(kzs[isou],kxs[isou],v,c,aw)
  do,ds,u = modeler.modelDataAndSourceWavefield()
  r = sub(ds,do) # residual
  a = modeler.modelReceiverWavefield(r) # receiver wavefield
  g = makeGradient(isou,u,a) # gradient
  plot(ds,perc=99.5,title='observed')
  plot(do,perc=99.5,title='simulated')
  plot(r,perc=99.5,title='residual')
  return g

class WaveformParallel(Parallel.ReduceInt):
  def compute(self,isou):
    return computeWaveformGradient(isou)
  def combine(self,v1,v2):
    return add(v1,v2)

#############################################################################
# Traveltime

def goTraveltime():
  g = zerofloat(nz,nx)
  aw = AcousticWavefield(sz,sx,st)
  for isou in range(ns):
    print "isou =",isou
    add(computeTraveltimeGradient(isou,aw),g,g)
  div(g,max(abs(g)),g)
  plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='gradient (all shots)')

def computeTraveltimeGradient(isou=0,aw=None):
  modeler = Modeler(kzs[isou],kxs[isou],v,c,aw)
  do,ds,u = modeler.modelDataAndSourceWavefield()

  # Residual
  dw,s = warp(ds,do) # order matters!
  #r = timeDerivative(dw)
  r = timeDerivative(ds) # XXX use ds instead of warped do
  mul(s,r,r)

  a = modeler.modelReceiverWavefield(r) # receiver wavefield
  g = makeGradient(isou,u,a) # gradient
  #plot(ds,perc=99.5,title='observed')
  #plot(do,perc=99.5,title='simulated')
  #plot(dw,perc=99.5,title='warped')
  #plot(r,perc=99.5,title='residual')
  #plot(s,cmap=jet,title='shifts')
  return g

#############################################################################
# Combined (Version A)

def goCombinedA():
  g = zerofloat(nz,nx)
  aw = AcousticWavefield(sz,sx,st)
  for isou in range(ns):
    print "isou =",isou
    add(computeCombinedAGradient(isou,aw),g,g)
  div(g,max(abs(g)),g)
  plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='gradient (all shots)')

def computeCombinedAGradient(isou=0,aw=None):
  modeler = Modeler(kzs[isou],kxs[isou],v,c,aw)
  do,ds,u = modeler.modelDataAndSourceWavefield()

  # Residual
  p = sub(ds,do) # amplitude residual
  dw,s = warp(ds,do) # order matters!
  #q = mul(s,timeDerivative(dw))
  q = mul(s,timeDerivative(ds)) # traveltime residual
  #plot(p,perc=99.5,title='amplitude residual')
  #plot(q,perc=99.5,title='traveltime residual')
  r,w = like(p),like(p)
  for ir in range(nr):
    for it in range(nt):
      t = s[ir][it]*dt
      wi = rcos(t)
      w[ir][it] = wi
      r[ir][it] = wi*p[ir][it]+(1.0-wi)*q[ir][it]

  a = modeler.modelReceiverWavefield(r) # receiver wavefield
  g = makeGradient(isou,u,a) # gradient
  #plot(do,perc=99.5,title='observed')
  #plot(ds,perc=99.5,title='simulated')
  #plot(dw,perc=99.5,title='warped')
  #plot(r,perc=99.5,title='residual')
  #plot(s,cmap=jet,title='shifts')
  #plot(w,cmap=jet,title='weights')
  return g


#############################################################################
# Combined (Version B)

def goCombinedB():
  do,ds,u = modelDataAndSourceWavefield() # data and source wavefield

  # Residual
  dw,s = warp(ds,do) # order matters!
  #z = timeDerivative(dw)
  z = timeDerivative(ds) # XXX use ds instead of warped do
  y = sub(ds,dw)
  plot(z,perc=99.5,title='traveltime residual')
  plot(y,perc=99.5,title='amplitude residual (after shifting)')

#############################################################################

def warp(p,q):
  print "warping..."
  d = 2 # decimate by this factor
  a = addRandomNoise(1.0e1,p,sigma=2.0)
  b = addRandomNoise(1.0e1,q,sigma=2.0)
  f = copy(nt/d,nr/d,0,0,d,d,a)
  g = copy(nt/d,nr/d,0,0,d,d,b)
  shiftMax = int(128/d)
  strainMax1 = 1.0
  strainMax2 = 1.0
  dw = DynamicWarping(-shiftMax,shiftMax)
  dw.setStrainMax(strainMax1,strainMax2)
  dw.setShiftSmoothing(8.0/d,8.0/d) # shift smoothing
  dw.setErrorSmoothing(2) # number of smoothings of alignment errors
  s = dw.findShifts(f,g)
  mul(d,s,s) # scale shifts to compensate for decimation
  li = LinearInterpolator()
  li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT)
  li.setUniform(nt/d,d*dt,0.0,nr/d,d*dz,0.0,s)
  r = like(p)
  for ir in range(nr):
    z = ir*dz
    for it in range(nt):
      t = it*dt
      r[ir][it] = li.interpolate(t,z)
  h = dw.applyShifts(r,q)
  #plot(s,cmap=jet,title='shifts')
  #plot(r,cmap=jet,title='shifts (interpolated)')
  #plot(f,perc=99.5,title='f')
  #plot(g,perc=99.5,title='g')
  #plot(h,perc=99.5,title='h')
  print "done"
  return h,r

def timeDerivative(f):
  """Derivative in time, assumed to be the 1st dimension."""
  n = len(f)
  #odt = 0.5/dt
  odt = 0.5
  g = like(f)
  for i in range(n):
    for it in range(1,nt-1):
      g[i][it] = odt*(f[i][it+1]-f[i][it-1])
  return g

def addRandomNoise(r,x,sigma=1.0):
  n1,n2 = len(x[0]),len(x)
  xrms = sqrt(sum(mul(x,x))/n1/n2) # rms of signal
  s = sub(randfloat(Random(3),n1,n2),0.5)
  #s = sub(randfloat(n1,n2),0.5)
  #sigma = 1.0
  RecursiveGaussianFilter(sigma).apply00(s,s); # bandlimited noise
  srms = sqrt(sum(mul(s,s))/n1/n2) # rms of noise
  return add(mul(xrms/(srms*r),s),x); # r = rms-signal / rms-noise

def rcos(t):
  """Amplitude response of a (modified) raised-cosine
     filter, used as a weight function."""
  oott = 300.0*dt/fpeak # 1/2T is the center of the transition zone
  beta = 0.25 # 0<beta<1 (beta=0 gives boxcar, beta=0 gives cosine)
  t = abs(t)
  p = 0.5/oott # period T
  a = (1.0-beta)*oott
  b = (1.0+beta)*oott
  if t<a:
    return 1.0
  elif t<b:
    return 0.5*(1.0+cos(PI*p*(t-a)/beta))
  else:
    return 0.0

def showWeights():
  r = zerofloat(200)
  for it in range(200):
    t = it*0.01
    r[it] = rcos(t)
  SimplePlot.asPoints(r)

#############################################################################

def modelDataAndSourceWavefield(jzs,jxs):
  aw = AcousticWavefield(sz,sx,st)

  # Observed data
  source = AcousticWavefield.RickerSource(fpeak,jzs,jxs)
  aw.forwardPropagate(source,v)
  do = aw.getWavefield(kzr,kxr)

  # Simulated data
  aw.forwardPropagate(source,c)
  ds = aw.getWavefield(kzr,kxr)
  u = aw.getWavefield()

  return do,ds,u

def modelReceiverWavefield(r):
  source = AcousticWavefield.AdjointSource(dt,kzr,kxr,r)
  aw.backPropagate(source,c)
  return aw.getWavefield()

def makeGradient(isou,u,a):
  g = correlate(u,a)
  m = zerofloat(nz,nx)
  m[kxs[isou]][kzs[isou]] = 1.0
  RecursiveGaussianFilter(32.0/dx).apply00(m,m)
  div(m,max(m),m)
  sub(1.0,m,m)
  mul(m,g,g)
  #div(g,max(abs(g)),g)
  #plot(m,jet)
  return g

def correlate(u,a):
  """Zero-lag correlation."""
  g = zerofloat(nz,nx)
  for it in range(1,nt-1):
    #t = add(add(u[it-1],u[it+1]),mul(-2.0,u[it])) # 2nd time derivative
    #add(mul(a[it],t),g,g)
    add(mul(a[it],u[it]),g,g)
  return g

def getMarmousi():
  v = zerofloat(751,2301)
  read("/data/seis/marmousi/marmousi.dat",v);
  c = copy(701,401,50,1550,v)
  RecursiveExponentialFilter(4.0).apply(c,c)
  RecursiveExponentialFilter(4.0).apply(c,c)
  #plot(v,cmap=jet)
  #plot(c,cmap=jet)
  return c

def getGaussian(pupper=0.05,plower=0.05):
  """ pupper - percent slowness perturbation of upper anomaly """
  """ plower - percent slowness perturbation of lower anomaly """
  v = fillfloat(4000,701,401);
  div(1.0,v,v); s0 = v[0][0]
  t = zerofloat(701,401)
  t[200][175] = -1.0
  t[200][525] = plower/pupper
  RecursiveGaussianFilter(60.0).apply00(t,t)
  div(t,max(abs(t)),t)
  mul(t,pupper*s0,t)
  add(t,v,v)
  div(1.0,v,v)
  return v

def like(x):
  return zerofloat(len(x[0]),len(x))

def read(name,image=None):
  if not image:
    image = zerofloat(n1,n2,n3)
  #fileName = dataDir+name+".dat"
  fileName = name
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def write(name,image):
  fileName = dataDir+name+'.dat'
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,cmap=gray,cmin=0,cmax=0,perc=100,cbar=None,title=None):
  panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  cb = panel.addColorBar()
  if cbar:
    cb.setLabel(cbar)
    cb.setWidthMinimum(120)
  else:
    cb.setWidthMinimum(80)
  pixel = panel.addPixels(x)
  pixel.setColorModel(cmap)
  if cmin<cmax:
    pixel.setClips(cmin,cmax)
  if perc<100:
    pixel.setPercentiles(100-perc,perc)
  pixel.setInterpolation(PixelsView.Interpolation.LINEAR)
  frame = PlotFrame(panel)
  if (len(x[0])>nz):
    frame.setSize(800,1000)
  else:
    frame.setSize(600,800)
  if title:
    frame.setTitle(title)
  frame.setVisible(True)
  if title and pngDir:
    frame.paintToPng(360,3.0,pngDir+title+'.png')

def xplot(x,cmap=gray,perc=100,cmin=0,cmax=0,title=None):
  pp = PlotPanel
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if (len(x[0])==nt):
    sp.setSize(800,1000)
  else:
    sp.setSize(600,800)
  sp.addColorBar()
  if title:
    sp.addTitle(title)
  pv = sp.addPixels(x)
  pv.setColorModel(cmap)
  #pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if perc<100:
    pv.setPercentiles(100-perc,perc)
  if cmin<cmax:
    pv.setClips(cmin,cmax)

#############################################################################
# Do everything on Swing thread.
import sys,time
class RunMain(Runnable):
  def run(self):
    start = time.time()
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02ds'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
