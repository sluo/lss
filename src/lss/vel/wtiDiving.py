"""
Combined waveform and traveltime inversion
"""
from imports import *

#############################################################################
sz = Sampling(265,0.012,0.0)
sx = Sampling(767,0.012,0.0)
st = Sampling(4000,0.0015,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta

#kxs,kzs = [0],[0]
kxs,kzs = [nx/2],[0]
#kxs,kzs = [nx-1],[0]
#kxs,kzs = rampint(8,50,16),fillint(0,16)
#kxs,kzs = rampint(8,30,26),fillint(0,26)
#kxs,kzs = rampint(3,20,39),fillint(0,39)
#kxs,kzs = rampint(1,15,52),fillint(0,52)
#kxs,kzs = rampint(3,10,77),fillint(0,77)
kxr,kzr = rampint(0,1,nx),fillint(0,nx)
ns,nr = len(kxs),len(kxr)
fpeak = 5.0
niter = 1

pngDir = None
datDir = None
#pngDir = os.getenv('HOME')+'/Desktop/png/'
#datDir = os.getenv('HOME')+'/Desktop/dat/'

sfile = None
#sfile = '/home/sluo/Desktop/s_iter9.dat'

gfile = None
#gfile = '/home/sluo/Desktop/g_iter0.dat'

#############################################################################

# TODO: use traveltime misfit for split inversion
# TODO: better line search?
# TODO: better agc
def main(args):
  setModel()
  showData()
  #goMigration()
  #WaveformInversion()
  #TraveltimeInversion()
  #SplitInversion()

def setModel():
  global _t,_s
  #_t,_s = getGaussian()
  #_t,_s = getMarmousi(0.2,1.0)
  #_t,_s = getMarmousi(0.5,1.0)
  #_t,_s = getMarmousi(2.0,1.0)
  _t,_s = getMarmousi(2.0,0.95)
  if sfile is not None:
    print 'reading sfile'
    _s = read(sfile)
  g = sub(_s,_t); div(g,max(abs(g)),g)
  plot(_t,cmap=jet,cbar='Slowness (s/km)',title='s_true')
  plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
       title='s_init')
  plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='g_true')

#############################################################################
# Migration

def goMigration():
  #sigma = 0.10
  #s = esmooth(sigma,_t)

  #s = read('/home/sluo/Desktop/fwi/dat/s_iter9.dat')
  s = read('/home/sluo/Desktop/combined/dat/s_iter9.dat')

  plot(s,cmap=jet,cbar='Smoothed slowness (s/km)',title='s')

  r = migration(s)
  m = zerofloat(nz,nx)
  for ix in range(nx):
    m[ix][0] = 1.0
  RecursiveGaussianFilter(16.0).apply0X(m,m)
  div(m,max(m),m)
  sub(1.0,m,m)
  plot(m)
  
  mul(m,r,r)

  plot(r,title='image')
  r = laplacian(r)
  for ix in range(nx):
    for iz in range(20):
      #r[ix][iz] = 0.0
      pass
  plot(r,title='laplacian')

def migration(s):
  #maskWaterLayer(s,2.0/3.0)
  do = modelData(_t) # observed data
  da = modelData(fillfloat(2.0/3.0,nz,nx)) # direct arrival
  sub(do,da,do)
  r = zerofloat(nz,nx) # image
  p = zerofloat(nz,nx) # preconditioner
  u,a = zerofloat(nz,nx,nt),zerofloat(nz,nx,nt)
  wave = Wavefield(sz,sx,st)
  for isou in range(ns):
    print 'isou =',isou

    source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
    receiver = Wavefield.Receiver(kzr,kxr)
    wave.modelAcousticWavefield(source,s,u)

    source = Wavefield.AdjointSource(dt,kzr,kxr,do[isou])
    wave.modelAcousticWavefield(source,s,a)
    num,den = makeGradient(u,a,d2=False)

    add(num,r,r)
    add(den,p,p)
  div(r,p,r)
  #maskWaterLayer(r)
  div(r,max(abs(r)),r)
  return r

def esmooth(sigma,x):
  ref = RecursiveExponentialFilter(sigma/dx)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  y = like(x)
  ref.apply(x,y)
  return y

def laplacian(x):
  n1,n2 = len(x[0]),len(x)
  y = like(x)
  class Loop(Parallel.LoopInt):
    def compute(self,i2):
      for i1 in range(1,n1-1):
        y[i2][i1] = (1.0/6.0)*(-20.0*x[i2][i1]+
          4.0*(x[i2-1][i1  ]+x[i2  ][i1+1]+x[i2  ][i1-1]+x[i2+1][i1  ])+
               x[i2+1][i1+1]+x[i2+1][i1-1]+x[i2-1][i1+1]+x[i2-1][i1-1]);
  Parallel.loop(1,n2-1,Loop())
  return y

#############################################################################
# Data

def modelData(s,isou=None):
  sw = Stopwatch(); sw.start()
  if isou is None:
    d = zerofloat(nt,nr,ns)
    Parallel.loop(ns,DataP(s,d))
  else:
    d = zerofloat(nt,nr)
    wave = Wavefield(sz,sx,st)
    source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
    receiver = Wavefield.Receiver(kzr,kxr)
    wave.modelAcousticData(source,receiver,s,d)
  sw.stop(); print 'data: %.2fs'%sw.time()
  return d

def modelDirectArrival(s,isou=None):
  t = copy(s)
  # Find water bottom then flood below
  t00 = t[0][0]
  zbot = 0
  while t[0][zbot]==t00:
    zbot += 1
  for ix in range(nx):
    for iz in range(zbot,nz):
      t[ix][iz] = t[ix][zbot] # flood
  # Scale and smooth model
  t = transpose(t)
  sigma = 20.0
  for iz in range(zbot+3,nz):
    r = 1.0+0.1*(iz-zbot-3.0)/(nz-zbot-10.0)
    mul(r,t[iz],t[iz])
    ref = RecursiveExponentialFilter(sigma)
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
    ref.apply(t[iz],t[iz])
    sigma += 1.0
  t = transpose(t)
  #plot(t,cmap=jet)
  return modelData(t,isou)

class DataP(Parallel.LoopInt):
  def __init__(self,s,d):
    self.s = s # slowness
    self.d = d # output array
  def compute(self,isou):
    wave = Wavefield(sz,sx,st)
    source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
    receiver = Wavefield.Receiver(kzr,kxr)
    wave.modelAcousticData(source,receiver,self.s,self.d[isou])

def showData():
  sw = Stopwatch(); sw.restart()
  do = modelData(_t) # observed data
  ds = modelData(_s) # simulated data
  eo = modelDirectArrival(_t)  # direct arrival for observed data
  es = modelDirectArrival(_s)  # direct arrival for simulated data
  sw.stop(); print 'data total: %.2fs'%sw.time()
  do = do[ns/2]
  ds = ds[ns/2]
  eo = eo[ns/2]
  es = es[ns/2]
  sub(do,eo,do) # observed
  sub(ds,es,ds) # simulated
  ra = sub(ds,do) # amplitude residual
  s,dw = like(do),like(do)
  rt,rw,rc = makeWarpedResiduals(do,ds,s,dw) # warped residuals
  smax = 0.9*max(abs(s))
  rmin,rmax = 0.5*min(min(ra),min(rc)),0.5*max(max(ra),max(rc))
  dmin,dmax = 0.5*min(min(do),min(ds)),0.5*max(max(do),max(ds))
  #emin,emax = 0.5*min(min(eo),min(es)),0.5*max(max(eo),max(es))
  emin,emax = 2.0*dmin,2.0*dmax
  plot(eo,cmin=emin,cmax=emax,title='direct arrival (observed)')
  plot(es,cmin=emin,cmax=emax,title='direct arrival (simulated)')
  plot(do,cmin=dmin,cmax=dmax,title='observed')
  plot(ds,cmin=dmin,cmax=dmax,title='simulated')
  plot(dw,cmin=dmin,cmax=dmax,title='warped')
  plot(s,cmap=rwb,cmin=-smax,cmax=smax,title='shifts')
  plot(ra,cmin=rmin,cmax=rmax,title='residual')
  plot(rw,cmin=rmin,cmax=rmax,title='warped residual')
  plot(rt,cmin=rmin,cmax=rmax,title='traveltime residual')
  plot(rc,cmin=rmin,cmax=rmax,title='traveltime+warped residual')

#############################################################################
# Line Search

def updateModel(misfitFunction):
  print 'searching for step length...'
  a = -0.5*max(abs(sub(_t,_s)))
  #a = -0.5*(rms(sub(_t,_s))/rms(misfitFunction.g))
  tol = 0.2*(-a)
  sw = Stopwatch(); sw.restart()
  step = BrentMinFinder(misfitFunction).findMin(a,0,tol)
  sw.stop()
  print 'a =',a
  print 'step =',step
  print 'line search: %.2fs'%sw.time()
  add(mul(step,misfitFunction.g),_s,_s)

class MisfitFunction(BrentMinFinder.Function):
  """Abstract misfit function."""
  def __init__(self,g,isou,do):
    self.g = g
    self.isou = isou
    self.do = do[isou]
    self.wave = Wavefield(sz,sx,st)
  def evaluate(self,a):
    print 'evaluating'
    s = add(_s,mul(a,self.g))
    ds = modelData(s,self.isou)
    es = modelDirectArrival(s,self.isou)
    sub(ds,es,ds)
    r = self.residual(self.do,ds)
    return sum(mul(r,r))
  def residual(self,do,ds):
    pass

class WaveformMisfitFunction(MisfitFunction):
  def residual(self,do,ds):
    return sub(ds,do)

class TraveltimeMisfitFunction(MisfitFunction):
  def residual(self,do,ds):
    rt,_,_ = makeWarpedResiduals(do,ds)
    return rt

class SplitMisfitFunction(MisfitFunction):
  def residual(self,do,ds):
    _,_,rc = makeWarpedResiduals(do,ds)
    return rc

#############################################################################
# Inversion

class Inversion():
  """Abstract inversion class."""
  def __init__(self):
    self.ds = zerofloat(nt,nr)
    self.u = zerofloat(nz,nx,nt)
    self.a = zerofloat(nz,nx,nt)
    self.wave = Wavefield(sz,sx,st)
    self.do = modelData(_t) # observed data
    eo = modelDirectArrival(_t) # direct arrival
    sub(self.do,eo,self.do)
    if gfile is not None:
      print 'reading gfile'
      g = read(gfile)
      updateModel(self.getMisfitFunction(g,ns/2,self.do))
      plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='g_iter')
      plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
           title='s_iter')
    self.invert()
  def invert(self):
    for iter in range(niter):
      print '\niteration',iter
      sw = Stopwatch(); sw.start()
      g = self.gradient()
      print 'gradient: %.2fm'%(sw.time()/60.0)
      plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='g_iter'+str(iter))
      if niter>1:
        updateModel(self.getMisfitFunction(g,ns/2,self.do))
        plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
             title='s_iter'+str(iter))
      sw.stop(); print 'iteration: %.2fm'%(sw.time()/60.0)
  def gradient(self):
    g = zerofloat(nz,nx)
    p = zerofloat(nz,nx) # preconditioner
    es = modelDirectArrival(_s) # direct arrival for simulated data
    for isou in range(ns):
      num,den = self.gradientForOneSource(isou,es)
      add(num,g,g)
      add(den,p,p)
    div(g,p,g)
    maskWaterLayer(g)
    div(g,max(abs(g)),g)
    return g
  def gradientForOneSource(self,isou,es):
    sw0 = Stopwatch(); sw0.start()
    source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
    receiver = Wavefield.Receiver(kzr,kxr)
    sw = Stopwatch(); sw.start()
    self.wave.modelAcousticDataAndWavefield(source,receiver,_s,self.ds,self.u)
    print 'forward: %.2fs'%sw.time()
    sub(self.ds,es[isou],self.ds) # subtract direct arrival
    checkForNaN(self.ds) # throws an exception if NaN
    r = self.residual(self.do[isou],self.ds) # residual
    source = Wavefield.AdjointSource(dt,kzr,kxr,r)
    sw.restart()
    self.wave.modelAcousticWavefield(source,_s,self.a)
    print 'reverse: %.2fs'%sw.time()
    g,p = makeGradient(self.u,self.a)
    sw0.stop(); print 'source %d: %.2fs'%(isou,sw0.time())
    return g,p
  def getMisfitFunction(g,isou,do):
    pass
  def residual(do,ds):
    pass

class WaveformInversion(Inversion):
  def residual(self,do,ds):
    return sub(ds,do)
  def getMisfitFunction(self,g,isou,do):
    return WaveformMisfitFunction(g,isou,do)

class TraveltimeInversion(Inversion):
  def residual(self,do,ds):
    rt,_,_ = makeWarpedResiduals(do,ds)
    return rt
  def getMisfitFunction(self,g,isou,do):
    return TraveltimMisfitFunction(g,isou,do)

class SplitInversion(Inversion):
  def residual(self,do,ds):
    _,_,rc = makeWarpedResiduals(do,ds)
    return rc
  def getMisfitFunction(self,g,isou,do):
    return SplitMisfitFunction(g,isou,do)

#############################################################################
# Dynamic Warping

def makeWarpedResiduals(do,ds,u=None,wd=None):
  reverseOrder = True
  if reverseOrder:
    v,dw,_ = warp(do,ds) # wrong order
    mul(-1.0,v,v)
    rt = mul(v,timeDerivative(ds)) # traveltime residual
    rw = sub(dw,do) # warped residual
  else:
    v,dw,rt = warp(ds,do) # right order
    rt = mul(v,timeDerivative(ds)) # traveltime residual
    rw = sub(ds,dw) # warped residual
  #print '  max shift =',max(abs(v))
  #print '  max combined residual =',max(abs(add(rt,rw)))
  if u is not None:
    copy(v,u)
  if wd is not None:
    copy(dw,wd)
  return rt,rw,add(rt,rw)

def warp(p,q):
  td = 5 # time decimation
  rd = 1 # receiver decimation
  qc = copy(q)
  p,q = agc(p,q)
  a,b = addRandomNoise(10.0,p,q,sigma=1.0)
  f = copy(nt/td,nr/rd,0,0,td,rd,a)
  g = copy(nt/td,nr/rd,0,0,td,rd,b)
  shiftMax = int(400/td)
  #strainMax1 = 0.25
  #strainMax2 = 0.10
  strainMax1 = 0.10
  strainMax2 = 0.05
  dw = DynamicWarping(-shiftMax,shiftMax)
  dw.setErrorExponent(1.0)
  dw.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.AVERAGE)
  dw.setStrainMax(strainMax1,strainMax2)
  dw.setShiftSmoothing(32.0/td,8.0/rd) # shift smoothing
  dw.setErrorSmoothing(2) # number of smoothings of alignment errors
  sw = Stopwatch(); sw.start()
  s = dw.findShifts(f,g)
  sw.stop(); print 'warping: %.2fs'%sw.time()
  mul(td,s,s) # scale shifts to compensate for decimation
  li = LinearInterpolator()
  li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT)
  li.setUniform(nt/td,td*dt,0.0,nr/rd,rd*dx,0.0,s)
  r = like(p)
  for ir in range(nr):
    z = ir*dz
    for it in range(nt):
      t = it*dt
      r[ir][it] = li.interpolate(t,z) # interpolate shifts
  h = dw.applyShifts(r,qc)
  hr = dw.applyShifts(r,mul(r,timeDerivative(qc)))
  #plot(s,cmap=jet,title='shifts')
  #plot(r,cmap=jet,title='shifts (interpolated)')
  #plot(f,title='f')
  #plot(g,title='g')
  #plot(a,title='a')
  #plot(b,title='b')
  #plot(h,title='h')
  #plot(hr,title='hr')
  return r,h,hr

def agc(do,ds):
  # AGC observed data
  ro = mul(do,do)
  ref = RecursiveExponentialFilter(50.0)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref.apply(ro,ro); div(ro,max(ro),ro); add(0.05,ro,ro)
  eo = div(do,ro)

  # AGC simulated data
  rs = mul(ds,ds)
  ref = RecursiveExponentialFilter(1000.0) # larger window
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref.apply(rs,rs); div(rs,max(rs),rs); add(0.05,rs,rs)
  es = div(ds,rs)

  # Divide by (scaled) rms
  """
  scale = 2.5
  print 'rms(eo) =',rms(eo)
  print 'rms(es) =',rms(es)
  div(eo,rms(eo),eo)
  div(es,rms(es)*scale,es)
  """

  #plot(ro,cmap=jet,title='ro')
  #plot(rs,cmap=jet,title='rs')
  #plot(do,title='do')
  #plot(ds,title='ds')
  #plot(eo,title='eo')
  #plot(es,title='es')
  return eo,es

def rms(x):
  return sqrt(sum(mul(x,x))/len(x[0])/len(x))

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

def addRandomNoise(snr,f,g,sigma=1.0):
  n1,n2 = len(f[0]),len(f)
  frms = sqrt(sum(mul(f,f))/n1/n2)
  grms = sqrt(sum(mul(g,g))/n1/n2)
  xrms = 0.5*(frms+grms)  # (average) rms of signal
  s = sub(randfloat(n1,n2),0.5)
  RecursiveGaussianFilter(sigma).apply00(s,s) # bandlimited noise
  srms = sqrt(sum(mul(s,s))/n1/n2) # rms of noise
  mul(xrms/(srms*snr),s,s)
  return add(f,s),add(g,s)

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

#############################################################################

def makeGradient(u,a,d2=True):
  g = correlate(u,a,d2) # gradient
  p = correlate(u,u) # illumination preconditioner
  return g,p

def correlate(u,a,d2=False):
  """Zero-lag correlation."""
  class ReduceInt(Parallel.ReduceInt):
    def compute(self,it):
      if d2:
        t = add(add(u[it-1],u[it+1]),mul(-2.0,u[it])) # 2nd time derivative
        return mul(-1.0,mul(a[it],t))
      else:
        return mul(a[it],u[it])
    def combine(self,g1,g2):
      return add(g1,g2)
  return Parallel.reduce(1,nt-1,ReduceInt())

def maskWaterLayer(g,value=0.0):
  t0 = _t[0][0]
  for ix in range(nx):
    iz = 0
    while _t[ix][iz]==t0:
        g[ix][iz] = value
        iz += 1

import socket
def getMarmousi(sigmaC=0.5,smul=1.0):
  p = zerofloat(751,2301)
  if socket.gethostname()=='backus.Mines.EDU':
    read("/data/sluo/marmousi/marmousi.dat",p)
  else:
    read("/data/seis/marmousi/marmousi.dat",p)
  p = copy(743,2301,8,0,p)
  div(1000.0,p,p) # convert to slowness
  q = copy(p)
  refC = RecursiveExponentialFilter(sigmaC/0.004)
  refC.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  refC.apply(q,q)
  mul(smul,q,q) # scale slowness, e.g., so that simulated data arrives earlier
  v = fillfloat(2.0/3.0,nz,nx)
  c = copy(v)
  copy(248,767,0,0,3,3,p,17,0,1,1,v)
  copy(248,767,0,0,3,3,q,17,0,1,1,c)
  return v,c

def getGaussian(splus=0.01):
  v = rampfloat(0.5,-0.001,0.0,nz,nx)
  c = copy(v)
  t = zerofloat(nz,nx); t[nx/2][nz/2] = 1.0
  RecursiveGaussianFilter(40.0).apply00(t,t)
  div(t,max(t),t)
  mul(t,splus,t)
  add(t,v,v)
  for ix in range(nx):
    for iz in range(17):
      v[ix][iz] = 2.0/3.0
      c[ix][iz] = 2.0/3.0
  return v,c

def like(x):
  return zerofloat(len(x[0]),len(x))

def checkForNaN(x):
  n1,n2 = len(x[0]),len(x)
  class Loop(Parallel.LoopInt):
    def compute(self,i2):
      for i1 in range(n1):
        if Float.isNaN(x[i2][i1]):
          raise RuntimeError('found NaN')
  Parallel.loop(n2,Loop())

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,cmap=gray,cmin=0,cmax=0,perc=100,cbar=None,title=None):
  panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  cb = panel.addColorBar()
  if cbar:
    cb.setLabel(cbar)
  cb.setWidthMinimum(160)
  if len(x[0])==nz and len(x)==nx:
    pixel = panel.addPixels(sz,sx,x)
    panel.setHLabel('Offset (km)')
    panel.setVLabel('Depth (km)')
  elif len(x[0])==nt and len(x)==nr:
    pixel = panel.addPixels(st,Sampling(len(x)),x)
    panel.setHLabel('Receiver')
    panel.setVLabel('Time (s)')
  else:
    pixel = panel.addPixels(x)
  pixel.setColorModel(cmap)
  if cmin<cmax:
    pixel.setClips(cmin,cmax)
  if perc<100:
    pixel.setPercentiles(100-perc,perc)
  pixel.setInterpolation(PixelsView.Interpolation.LINEAR)
  frame = PlotFrame(panel)
  frame.setFontSizeForSlide(1.5,1.5)
  if (len(x[0])==nz):
    frame.setSize(1200,500)
  else:
    frame.setSize(1200,1000)
  if title:
    frame.setTitle(title)
  frame.setVisible(True)
  if title and pngDir:
    frame.paintToPng(360,3.0,pngDir+title+'.png')
  if title and datDir:
    write(datDir+title+'.dat',x)

def read(name,image=None):
  if not image:
    image = zerofloat(nz,nx)
  fileName = name
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def write(fname,image):
  aos = ArrayOutputStream(fname)
  aos.writeFloats(image)
  aos.close()

def cleanDir(dir):
  os.chdir(dir)
  for f in os.listdir(dir):
    os.remove(f)

#############################################################################
# Do everything on Swing thread.
import sys,time
class RunMain(Runnable):
  def run(self):
    start = time.time()
    if pngDir is not None:
      print 'cleaning pngDir'
      cleanDir(pngDir)
    if datDir is not None:
      print 'cleaning datDir'
      cleanDir(datDir)
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02d'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
