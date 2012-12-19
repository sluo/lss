"""
Combined waveform and traveltime inversion
"""
from imports import *

#############################################################################
# Parameters

sz = Sampling(401,0.008,0.0)
sx = Sampling(701,0.008,0.0)
st = Sampling(3000,0.0010,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
kxs,kzs = [nx/2],[0]
#kxs,kzs = [nx/3,2*nx/3],[0]
#kxs,kzs = rampint(0,10,71),fillint(0,69)
#kxs,kzs = rampint(0,25,29),fillint(0,29)
#kxs,kzs = rampint(0,50,15),fillint(0,15)
kxr,kzr = rampint(0,1,nx),fillint(0,nx) 
ns,nr = len(kxs),len(kxr)
fpeak = 12.0 # Ricker wavelet peak frequency

pngDir = None
#pngDir = '/Users/sluo/Desktop/png/'
#############################################################################

parallel = False
def main(args):
  #go()
  #old()
  #new()

  #vmin,vmax = min(v),max(v)
  #v0 = lowpass(v) 
  #plot(v,cmap=jet)
  #plot(v0,cmap=jet,cmin=vmin,cmax=vmax)

  #acoustic()
  #acousticX()
  born()

def acousticX():
  setModel()
  div(1.0,v,v)
  aw = AcousticWavefieldY(sz,sx,st)
  source = AcousticWavefieldY.RickerSource(fpeak,kzs[0],kxs[0])
  aw.forwardPropagate(source,v)
  d = aw.getWavefield(kzr,kxr)
  for ir in range(nr):
    for it in range(2000):
      d[ir][it] = 0.0
  plot(d,perc=99.5)


def acoustic():
  setModel()
  w = Wavefield(sz,sx,st)
  source = Wavefield.RickerSource(fpeak,kzs[0],kxs[0])
  receiver = Wavefield.Receiver(kzr,kxr)
  d = zerofloat(nt,nr)
  w.modelAcousticData(source,receiver,v,d)
  #plot(d,perc=99.5)

  for ir in range(nr):
    for it in range(2000):
      d[ir][it] = 0.0
  #plot(d,perc=99.5)
  plot(d,cmin=-0.04,cmax=0.07)

def born():
  setModel()
  s = v

  ref = RecursiveExponentialFilter(16.0)
  smin,smax = min(s),max(s)
  s0 = like(s)
  ref.apply(s,s0)
  s1 = sub(s,s0)
  plot(s0,cmap=jet,title="background slowness")
  plot(s1,cmap=jet,title="perturbation slowness")

  w = Wavefield(sz,sx,st)
  source = Wavefield.RickerSource(fpeak,kzs[0],kxs[0])
  receiver = Wavefield.Receiver(kzr,kxr)
  d = zerofloat(nt,nr)
  w.modelBornData(source,receiver,s0,s1,d)
  #plot(d,perc=99.5)
  plot(d,cmin=-0.04,cmax=0.07)


def lowpass(x):
  bp = BandPassFilter(0.005,0.5,0.02,0.01)
  y = like(x)
  bp.apply(x,y)
  return y

count = 1
def old():
  setModel()
  aw = AcousticWavefieldY(sz,sx,st)
  source = AcousticWavefieldY.RickerSource(fpeak,kzs[0],kxs[0])
  print "start"
  sw = Stopwatch(); sw.restart()
  for i in range(count):
    aw.forwardPropagate(source,v)
  sw.stop()
  print "stop", sw.time()
  d = aw.getWavefield(kzr,kxr)
  plot(d,perc=99.5)

def new():
  setModel()
  w = Wavefield(sz,sx,st)
  source = Wavefield.RickerSource(fpeak,kzs[0],kxs[0])
  receiver = Wavefield.Receiver(kzr,kxr)
  d = zerofloat(nt,nr)
  print "start"
  sw = Stopwatch(); sw.restart()
  for i in range(count):
    w.modelAcousticData(source,receiver,v,d)
  sw.stop()
  print "stop", sw.time()
  plot(d,perc=99.5)

def go():
  setModel()
  #showWeights()
  #goWaveform()
  #goTraveltime()
  #goCombinedA()
  #goCombinedB()

  modeler = ReflectionModeler(kzs[0],kxs[0],v,c,b)
  do,ds,u = modeler.modelDataAndSourceWavefield()
  cmax = max(max(do),max(ds))
  cmin = min(min(do),min(ds))
  plot(do,cmin=cmin,cmax=cmax,title="observed data")
  plot(ds,cmin=cmin,cmax=cmax,title="simulated data")


  # Residual
  p = sub(ds,do) # amplitude residual
  dw,s = warp(ds,do) # order matters!
  #q = mul(s,timeDerivative(dw))
  q = mul(s,timeDerivative(ds)) # traveltime residual
  plot(p,cmin=cmin,cmax=cmax,title='amplitude residual')
  plot(q,cmin=cmin,cmax=cmax,title='traveltime residual')
  plot(dw,cmin=cmin,cmax=cmax,title='shifted observed data')
  plot(sub(ds,dw),cmin=cmin,cmax=cmax,title='shifted amplitude residual')
  """
  r,w = like(p),like(p)
  for ir in range(nr):
    for it in range(nt):
      t = s[ir][it]*dt
      wi = rcos(t)
      w[ir][it] = wi
      r[ir][it] = wi*p[ir][it]+(1.0-wi)*q[ir][it]
  """
  dg = gain(dw,ds)
  plot(dg,cmin=cmin,cmax=cmax,title='gained simulated data')

  #a = modeler.modelReceiverWavefield(r) # receiver wavefield
  #g = makeGradient(isou,u,a) # gradient
  #plot(dw,cmin=cmin,cmax=cmax,title='warped')
  #plot(r,cmin=cmin,cmax=cmax,title='residual')
  plot(s,cmap=jet,title='shifts')
  #plot(w,cmap=jet,title='weights')

def gain(f,g):
  sigma1 = 32.0
  sigma2 = 8.0
  eps = 1.0e-4 # stabilize division
  p = mul(f,f)
  q = mul(g,g)
  #ref1 = RecursiveGaussianFilter(sigma1)
  #ref2 = RecursiveGaussianFilter(sigma2)
  #ref1.apply0X(p,p); ref2.applyX0(p,p)
  #ref1.apply0X(q,q); ref2.applyX0(q,q)
  ref1 = RecursiveExponentialFilter(sigma1)
  ref2 = RecursiveExponentialFilter(sigma2)
  ref1.apply1(p,p); ref2.apply2(p,p)
  ref1.apply1(q,q); ref2.apply2(q,q)
  r = div(p,add(q,fillfloat(eps,len(f[0]),len(f))))
  r = sqrt(r)
  h = mul(r,g)
  plot(r,cmap=jet,title="gain")
  return h

def setModel():
  global v,c,b
  #v = getMarmousi()
  #v = getGaussian(0.02,0.02)
  #v = getGaussian(0.20,0.20)
  """
  v = getGaussian(0.20,0.02)
  c = fillfloat(2500,nz,nx)
  b = copy(c)
  for ix in range(nx):
    for iz in range(4*nz/5,nz):
      c[ix][iz] = 5000.0
  """
  v,c,b = getLayered(0.05)
  #v,c,b = getGaussian1(0.05)
  #v,c,b = getGaussian2(0.02,0.02)
  #v,c,b = getGaussian2(0.20,0.20)
  #v,c,b = getGaussian2(0.20,0.02)
  div(1.0,v,v),div(1.0,c,c),div(1.0,b,b)
  #plot(v,cmap=jet,title='true')
  #plot(c,cmap=jet,cmin=min(v),cmax=max(v),title='initial')
  #plot(b,cmap=jet,cmin=min(v),cmax=max(v),title='background')
  return v,c,b

class ReflectionModeler():
  def __init__(self,jzs,jxs,v,c,b,aw=None):
    self.jzs = jzs
    self.jxs = jxs
    self.aw = AcousticWavefieldY(sz,sx,st) if aw==None else aw
  def modelDataAndSourceWavefield(self):
    source = AcousticWavefieldY.RickerSource(fpeak,self.jzs,self.jxs)
    self.aw.forwardPropagate(source,v)
    do = self.aw.getWavefield(kzr,kxr) # observed data
    self.aw.forwardPropagate(source,c)
    ds = self.aw.getWavefield(kzr,kxr) # simulated data
    u = self.aw.getWavefield() # simulated wavefield TODO: remove direct wave?
    """
    self.aw.forwardPropagate(source,b)
    da = self.aw.getWavefield(kzr,kxr) # direct arrival
    sub(ds,da,ds) # remove direct arrival
    sub(do,da,do)
    """
    for ix in range(nx):
      for iz in range(1500):
        do[ix][iz] = 0.0
        ds[ix][iz] = 0.0
    return do,ds,u
  def modelReceiverWavefield(self,r):
    source = AcousticWavefieldY.AdjointSource(dt,kzr,kxr,r)
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
  #plot(ds,perc=99.5,title='observed')
  #plot(do,perc=99.5,title='simulated')
  #plot(r,perc=99.5,title='residual')
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

def getLayered(p=0.05):
  vfill = 2.5
  v = fillfloat(vfill,nz,nx)
  c,b = copy(v),copy(v)
  div(1.0,v,v)
  mul(1.0+p,v,v)
  div(1.0,v,v)
  for ix in range(nx):
    for iz in range(4*nz/5,nz):
      c[ix][iz] = 4.0
      v[ix][iz] = 6.0
  return v,c,b

def getGaussian1(p=0.05):
  """ p - percent slowness perturbation """
  vfill = 2500.0
  v = fillfloat(vfill,nz,nx)
  c,b = copy(v),copy(v)
  div(1.0,v,v); s0 = v[0][0]
  t = zerofloat(nz,nx)
  t[nx/2][2*nz/5] = -1.0
  RecursiveGaussianFilter(60.0).apply00(t,t)
  div(t,max(abs(t)),t)
  mul(t,p*s0,t)
  add(t,v,v)
  div(1.0,v,v)
  for ix in range(nx):
    for iz in range(4*nz/5,nz):
      c[ix][iz] = 5000.0
      v[ix][iz] = 6000.0
  return v,c,b

def getGaussian2(pleft=0.05,pright=0.05):
  """ pleft - percent slowness perturbation of left anomaly """
  """ pright - percent slowness perturbation of right anomaly """
  vfill = 2500.0
  v = fillfloat(vfill,nz,nx)
  c,b = copy(v),copy(v)
  div(1.0,v,v); s0 = v[0][0]
  t = zerofloat(nz,nx)
  t[  nx/4][2*nz/5] = -1.0
  t[3*nx/4][2*nz/5] = pright/pleft
  RecursiveGaussianFilter(60.0).apply00(t,t)
  div(t,max(abs(t)),t)
  mul(t,pleft*s0,t)
  add(t,v,v)
  div(1.0,v,v)
  for ix in range(nx):
    for iz in range(4*nz/5,nz):
      c[ix][iz] = 5000.0
      v[ix][iz] = 6000.0
  return v,c,b

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
    frame.setSize(1000,500)
  else:
    frame.setSize(800,500)
  if title:
    frame.setTitle(title)
  frame.setVisible(True)
  if title and pngDir:
    frame.paintToPng(360,3.0,pngDir+title+'.png')

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
