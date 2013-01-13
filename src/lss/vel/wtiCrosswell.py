"""
Combined waveform and traveltime inversion
"""
from imports import *

#############################################################################

sz = Sampling(701,0.004,0.0)
sx = Sampling(401,0.004,0.0)
st = Sampling(1750,0.0005,0.0) # for gaussian
##st = Sampling(3000,0.0006,0.0) # for marmousi
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
#kxs,kzs = [0],[0]
#kxs,kzs = [0],[nz/2]
#kxs,kzs = [0,0],[nz/3,2*nz/3]
#kxs,kzs = fillint(0,15),rampint(0,50,15) 
#kxs,kzs = fillint(0,29),rampint(0,25,29) 
kxs,kzs = fillint(0,71),rampint(0,10,71) 
kxr,kzr = fillint(nx-1,nz),rampint(0,1,nz) 
ns,nr = len(kxs),len(kxr)
fpeak = 50.0 # Ricker wavelet peak frequency
niter = 2

pngDir = None
datDir = None
pngDir = os.getenv('HOME')+'/Desktop/png/'
datDir = os.getenv('HOME')+'/Desktop/dat/'

sfile = None
#sfile = '/Users/sluo/Desktop/s_iter4.dat'

#############################################################################

def main(args):
  setModel(sfile)
  #showWeights()
  #goWaveform()
  #goTraveltime()
  #goCombined()
  #goCombinedX()

  #showData()
  WaveformInversion()
  SplitInversion()
      
def setModel(ffile=None):
  global _t,_s
  #_t,_s = getGaussian(0.02,0.02),fillfloat(0.25,nz,nx)
  _t,_s = getGaussian(0.20,0.20),fillfloat(0.25,nz,nx)
  #_t,_s = getGaussian(0.20,0.02),fillfloat(0.25,nz,nx)
  #_t,_s = getMarmousi(),fillfloat(0.35,nz,nx)
  if ffile is not None:
    _s = read(ffile)
  g = sub(_s,_t); div(g,max(abs(g)),g)
  plot(_t,cmap=jet,cbar='Slowness (s/km)',title='s_true')
  plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
       title='s_init')
  plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='g_true')

#############################################################################
# Data

def modelData(s):
  d = zerofloat(nt,nr,ns)
  sw = Stopwatch(); sw.restart()
  Parallel.loop(ns,DataP(s,d))
  sw.stop(); print 'data:',sw.time(),'s'
  return d

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
  do = zerofloat(nt,nr,ns) # observed data
  ds = zerofloat(nt,nr,ns) # simulated data
  sw = Stopwatch(); sw.restart()
  Parallel.loop(ns,DataP(_t,do))
  Parallel.loop(ns,DataP(_s,ds))
  sw.stop(); print 'data:',sw.time(),'s'
  do = do[ns/2]
  ds = ds[ns/2]
  ra = sub(ds,do) # amplitude residual

  #dw,v,_ = warp(ds,do) # right order
  #w,m = makeWeights(v) # weighting function
  #rt = mul(v,timeDerivative(ds)) # traveltime residual
  #rw = sub(ds,dw) # warped residual
  #"""
  #dw,v, = warp(do,ds) # wrong order
  #w,m = makeWeights(v) # weighting function
  #rt = mul(mul(-1.0,v),timeDerivative(ds)) # traveltime residual
  #rw = sub(dw,do) # warped residual
  #"""
  #rd = add(mul(w,ra),mul(m,rt)) # combined residual (traveltime+amplitude)

  v,dw = like(do),like(do)
  rt,rw,rc = makeWarpedResiduals(do,ds,v,dw) # residual
  w,m = makeWeights(v) # weighting function
  rd = add(mul(w,ra),mul(m,rt)) # combined residual (traveltime+amplitude)

  vmax = 0.9*max(abs(v))
  rmin,rmax = 0.8*min(rc),0.8*max(rc)
  dmin,dmax = 0.8*min(do),0.8*max(do)
  plot(do,cmin=dmin,cmax=dmax,title='observed')
  plot(ds,cmin=dmin,cmax=dmax,title='simulated')
  plot(dw,cmin=dmin,cmax=dmax,title='warped')
  plot(v,cmap=rwb,cmin=-vmax,cmax=vmax,title='shifts')
  #plot(w,cmap=jet,cmin=0.0,cmax=1.0,title='weights')
  #plot(ra,cmin=rmin,cmax=rmax,title='amplitude_residual')
  #plot(rw,cmin=rmin,cmax=rmax,title='warped_residual')
  #plot(rt,cmin=rmin,cmax=rmax,title='traveltime_residual')
  #plot(rc,cmin=rmin,cmax=rmax,title='traveltime_plus_warped_residual')
  #plot(rd,cmin=rmin,cmax=rmax,title='weighted_traveltime_plus_amplitude')

#############################################################################
# Line Search

def findStepLength(g,isou,do,da=None):
  print 'searching for step length...'
  a = -max(abs(sub(_t,_s)))
  tol = 0.20*(-a)
  sw = Stopwatch(); sw.restart()
  step = BrentMinFinder(WaveformMisfitFunction(g,isou,do,da)).findMin(a,0,tol)
  sw.stop()
  print 'a =',a
  print 'step =',step
  print 'line search:',sw.time(),'s'
  return step

class WaveformMisfitFunction(BrentMinFinder.Function):
  def __init__(self,g,isou,do,da=None):
    self.g = g
    self.isou = isou
    self.do = do[isou]
    self.da = None if da is None else da[isou]
    self.wave = Wavefield(sz,sx,st)
  def evaluate(self,a):
    print '  evaluating'
    s = add(_s,mul(a,self.g))
    ds = zerofloat(nt,nr)
    source = Wavefield.RickerSource(fpeak,kzs[self.isou],kxs[self.isou])
    receiver = Wavefield.Receiver(kzr,kxr)
    self.wave.modelAcousticData(source,receiver,s,ds) # simulated data
    if self.da is not None:
      sub(ds,self.da,ds) # subtract direct arrival
    dif = sub(ds,self.do)
    return sum(mul(dif,dif))

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

class WaveformMisfitFunction(BrentMinFinder.Function):
  def __init__(self,g,isou,do):
  #def __init__(self,g,isou,do,da=None):
    self.g = g
    self.isou = isou
    self.do = do[isou]
    #self.da = None if da is None else da[isou]
    self.wave = Wavefield(sz,sx,st)
  def evaluate(self,a):
    print 'evaluating'
    s = add(_s,mul(a,self.g))
    ds = modelData(s,self.isou)
    dif = sub(ds,self.do)
    return sum(mul(dif,dif))

class SplitMisfitFunction(BrentMinFinder.Function):
  def __init__(self,g,isou,do):
    self.g = g
    self.isou = isou
    self.do = do[isou]
    self.wave = Wavefield(sz,sx,st)
  def evaluate(self,a):
    print 'evaluating'
    s = add(_s,mul(a,self.g))
    ds = modelData(s,self.isou)
    rt,rw,rc = makeWarpedResiduals(self.do,ds) # residual
    return sum(mul(rc,rc))

class TraveltimeMisfitFunction(BrentMinFinder.Function):
  def __init__(self,g,isou,do):
    self.g = g
    self.isou = isou
    self.do = do[isou]
    self.wave = Wavefield(sz,sx,st)
  def evaluate(self,a):
    print 'evaluating'
    s = add(_s,mul(a,self.g))
    ds = modelData(s,self.isou)
    rt,rw,rc = makeWarpedResiduals(self.do,ds) # residual
    return sum(mul(rt,rt))

#############################################################################
# Waveform

class WaveformInversion():
  def __init__(self):
    self.ds = zerofloat(nt,nr)
    self.u = zerofloat(nz,nx,nt)
    self.a = zerofloat(nz,nx,nt)
    self.wave = Wavefield(sz,sx,st)
    self.do = modelData(_t) # observed data
    for iter in range(niter):
      print '\niteration',iter
      sw = Stopwatch(); sw.start()
      g = self.gradient()
      print 'gradient: %.2fm'%(sw.time()/60.0)
      plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='g_iter'+str(iter))
      if niter>1:
        updateModel(WaveformMisfitFunction(g,ns/2,self.do))
        plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
             title='s_iter'+str(iter))
      sw.stop(); print 'iteration: %.2fm'%(sw.time()/60.0)
  def gradient(self):
    g = zerofloat(nz,nx)
    p = zerofloat(nz,nx) # preconditioner
    for isou in range(ns):
      num,den = self.gradientForOneSource(isou)
      add(num,g,g)
      add(den,p,p)
    div(g,p,g)
    maskWaterLayer(g)
    div(g,max(abs(g)),g)
    return g
  def gradientForOneSource(self,isou):
    sw0 = Stopwatch(); sw0.start()
    source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
    receiver = Wavefield.Receiver(kzr,kxr)
    sw = Stopwatch(); sw.start()
    self.wave.modelAcousticDataAndWavefield(source,receiver,_s,self.ds,self.u)
    print 'forward: %.2fs'%sw.time()
    checkForNaN(self.ds) # throws an exception if NaN
    r = sub(self.ds,self.do[isou]) # residual
    source = Wavefield.AdjointSource(dt,kzr,kxr,r)
    sw.restart()
    self.wave.modelAcousticWavefield(source,_s,self.a)
    print 'reverse: %.2fs'%sw.time()
    g,p = makeGradient(self.u,self.a)
    sw0.stop(); print 'source %d: %.2fs'%(isou,sw0.time())
    return g,p

def goWaveform():
  d = modelData(_t) # observed data
  for iter in range(niter):
    print '\niteration',iter
    sw = Stopwatch(); sw.restart()
    g = waveformGradientS(d)
    sw.stop(); print 'gradient:',sw.time(),'s'
    plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='g_iter'+str(iter))
    if niter>1:
      step = findStepLength(g,ns/2,d)
      add(mul(step,g),_s,_s)
      plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
           title='s_iter'+str(iter))

def waveformGradientS(d):
  g = zerofloat(nz,nx)
  p = zerofloat(nz,nx) # preconditioner
  ds,u,a = zerofloat(nt,nr),zerofloat(nz,nx,nt),zerofloat(nz,nx,nt)
  for isou in range(ns):
    print 'isou =',isou
    num,den = waveformGradient(isou,d[isou],ds,u,a)
    add(num,g,g)
    add(den,p,p)
  div(g,p,g)
  smoothSourceLocations(g)
  div(g,max(abs(g)),g)
  return g

def smoothSourceLocations(g):
  for ix in range(20):
    RecursiveExponentialFilter(20-ix).apply(g[ix],g[ix])

def waveformGradient(isou,do,ds,u,a):
  wave = Wavefield(sz,sx,st)
  source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
  receiver = Wavefield.Receiver(kzr,kxr)
  wave.modelAcousticDataAndWavefield(source,receiver,_s,ds,u)
  checkForNaN(ds)
  r = sub(ds,do) # residual
  source = Wavefield.AdjointSource(dt,kzr,kxr,r)
  wave.modelAcousticWavefield(source,_s,a)
  return makeGradient(u,a)

#############################################################################
# Traveltime

def goTraveltime():
  g = zerofloat(nz,nx)
  for isou in range(ns):
    print "isou =",isou
    add(computeTraveltimeGradient(isou),g,g)
  smoothSourceLocations(g)
  div(g,max(abs(g)),g)
  plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='gradient (all shots)')

def computeTraveltimeGradient(isou=0):
  modeler = Modeler(kzs[isou],kxs[isou],v,c)
  do,ds,u = modeler.modelDataAndSourceWavefield()

  # Residual
  dw,s,_ = warp(ds,do) # order matters!
  #r = timeDerivative(dw)
  r = timeDerivative(ds) # XXX use ds instead of warped do
  mul(s,r,r)

  a = modeler.modelReceiverWavefield(r) # receiver wavefield
  g = makeGradient(u,a,isou) # gradient
  #plot(do,perc=99.5,title='observed')
  #plot(ds,perc=99.5,title='simulated')
  #plot(dw,perc=99.5,title='warped')
  #plot(r,perc=99.5,title='residual')
  #plot(s,cmap=jet,title='shifts')
  return g
  #return zerofloat(nz,nx)

#############################################################################
# Split (Combined)

class SplitInversion():
  def __init__(self):
    self.ds = zerofloat(nt,nr)
    self.u = zerofloat(nz,nx,nt)
    self.a = zerofloat(nz,nx,nt)
    self.wave = Wavefield(sz,sx,st)
    self.do = modelData(_t) # observed data
    for iter in range(niter):
      print '\niteration',iter
      sw = Stopwatch(); sw.start()
      g = self.gradient()
      print 'gradient: %.2fm'%(sw.time()/60.0)
      plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='g_iter'+str(iter))
      if niter>1:
        #updateModel(SplitMisfitFunction(g,ns/2,self.do))
        updateModel(TraveltimeMisfitFunction(g,ns/2,self.do))
        plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
             title='s_iter'+str(iter))
      sw.stop(); print 'iteration: %.2fm'%(sw.time()/60.0)
  def gradient(self):
    g = zerofloat(nz,nx)
    p = zerofloat(nz,nx) # preconditioner
    for isou in range(ns):
      num,den = self.gradientForOneSource(isou)
      add(num,g,g)
      add(den,p,p)
    div(g,p,g)
    smoothSourceLocations(g)
    div(g,max(abs(g)),g)
    return g
  def gradientForOneSource(self,isou):
    sw0 = Stopwatch(); sw0.start()
    source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
    receiver = Wavefield.Receiver(kzr,kxr)
    sw = Stopwatch(); sw.start()
    self.wave.modelAcousticDataAndWavefield(source,receiver,_s,self.ds,self.u)
    print 'forward: %.2fs'%sw.time()
    checkForNaN(self.ds) # throws an exception if NaN
    _,_,rc = makeWarpedResiduals(self.do[isou],self.ds) # residual
    source = Wavefield.AdjointSource(dt,kzr,kxr,rc)
    sw.restart()
    self.wave.modelAcousticWavefield(source,_s,self.a)
    print 'reverse: %.2fs'%sw.time()
    g,p = makeGradient(self.u,self.a)
    sw0.stop(); print 'source %d: %.2fs'%(isou,sw0.time())
    return g,p

def goCombined():
  do = modelData(_t) # observed data
  for iter in range(niter):
    print '\niteration',iter
    sw = Stopwatch(); sw.restart()
    g = combinedGradientS(do)
    sw.stop(); print 'gradient:',sw.time(),'s'
    plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='g_iter'+str(iter))
    if niter>1:
      step = findStepLength(g,ns/2,do)
      add(mul(step,g),_s,_s)
      plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
           title='s_iter'+str(iter))

def combinedGradientS(do):
  g = zerofloat(nz,nx)
  p = zerofloat(nz,nx) # preconditioner
  ds,u,a = zerofloat(nt,nr),zerofloat(nz,nx,nt),zerofloat(nz,nx,nt)
  for isou in range(ns):
    print 'isou =',isou
    num,den = combinedGradient(isou,do[isou],ds,u,a)
    add(num,g,g)
    add(den,p,p)
  div(g,p,g)
  smoothSourceLocations(g)
  div(g,max(abs(g)),g)
  return g

def combinedGradient(isou,do,ds,u,a):
  wave = Wavefield(sz,sx,st)
  source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
  receiver = Wavefield.Receiver(kzr,kxr)
  wave.modelAcousticDataAndWavefield(source,receiver,_s,ds,u)
  checkForNaN(ds)

  # Residual
  dw,s,_ = warp(do,ds) # wrong order
  rt = mul(mul(-1.0,s),timeDerivative(ds)) # traveltime residual
  rw = sub(dw,do) # warped residual
  rc = add(rt,rw) # combined

  source = Wavefield.AdjointSource(dt,kzr,kxr,rc)
  wave.modelAcousticWavefield(source,_s,a)
  return makeGradient(u,a)

#############################################################################
# Combined (using weighting function)

def goCombinedX():
  do = modelData(_t) # observed data
  for iter in range(niter):
    print '\niteration',iter
    sw = Stopwatch(); sw.restart()
    g = combinedGradientXS(do)
    sw.stop(); print 'gradient:',sw.time(),'s'
    if niter>1:
      step = findStepLength(g,ns/2,do)
      add(mul(step,g),_s,_s)
      plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
           title='s_iter'+str(iter))

def combinedGradientXS(do):
  g = zerofloat(nz,nx)
  p = zerofloat(nz,nx) # preconditioner
  ds,u,a = zerofloat(nt,nr),zerofloat(nz,nx,nt),zerofloat(nz,nx,nt)
  for isou in range(ns):
    print 'isou =',isou
    num,den = combinedGradientX(isou,do[isou],ds,u,a)
    add(num,g,g)
    add(den,p,p)
  div(g,p,g)
  div(g,max(abs(g)),g)
  return g

def combinedGradientX(isou,do,ds,u,a):
  wave = Wavefield(sz,sx,st)
  source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
  receiver = Wavefield.Receiver(kzr,kxr)
  wave.modelAcousticDataAndWavefield(source,receiver,_s,ds,u)

  # Residual
  dw,v,_ = warp(do,ds) # wrong order
  w,m = makeWeights(v)
  ra = sub(ds,do) # amplitude residual
  rt = mul(mul(-1.0,v),timeDerivative(ds)) # traveltime residual
  rc = add(mul(w,ra),mul(m,rt)) # combined residual

  source = Wavefield.AdjointSource(dt,kzr,kxr,rc)
  wave.modelAcousticWavefield(source,_s,a)
  return makeGradient(u,a)

def makeWeights(v):
  w = zerofloat(nt,nr)
  m = zerofloat(nt,nr)
  for ir in range(nr):
    for it in range(nt):
      t = v[ir][it]*dt
      wi = rcos(t)
      w[ir][it] = wi
      m[ir][it] = 1.0-wi
  return w,m

def rcos(t):
  """Amplitude response of a (modified) raised-cosine
     filter, used as a weighting function."""
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
  #p,q = agc(p,q)
  a,b = addRandomNoise(10.0,p,q,sigma=1.0)
  f = copy(nt/td,nr/rd,0,0,td,rd,a)
  g = copy(nt/td,nr/rd,0,0,td,rd,b)
  shiftMax = int(400/td)
  #strainMax1 = 0.25
  #strainMax2 = 0.10
  strainMax1 = 1.00
  strainMax2 = 0.50
  dw = DynamicWarping(-shiftMax,shiftMax)
  #dw.setErrorExponent(1.0)
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

def xwarp(p,q):
  td = 1 # time decimation
  rd = 2 # receiver decimation
  a,b = addRandomNoise(2.0,p,q,sigma=1.0)
  f = copy(nt/td,nr/rd,0,0,td,rd,a)
  g = copy(nt/td,nr/rd,0,0,td,rd,b)
  #shiftMax = int(400/td)
  shiftMax = int(200/td)
  strainMax1 = 1.00
  strainMax2 = 1.00
  dw = DynamicWarping(-shiftMax,shiftMax)
  #dw.setErrorExponent(2.0)
  dw.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.AVERAGE)
  dw.setStrainMax(strainMax1,strainMax2)
  dw.setShiftSmoothing(32.0/td,8.0/rd) # shift smoothing
  dw.setErrorSmoothing(2) # number of smoothings of alignment errors
  sw = Stopwatch(); sw.restart()
  s = dw.findShifts(f,g)
  sw.stop(); print "warping:",sw.time(),"s"
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
  h = dw.applyShifts(r,q)
  #plot(s,cmap=jet,title='shifts')
  #plot(r,cmap=jet,title='shifts (interpolated)')
  #plot(f,perc=99.5,title='f')
  #plot(g,perc=99.5,title='g')
  #plot(a,perc=99.5,title='a')
  #plot(b,perc=99.5,title='b')
  #plot(h,perc=99.5,title='h')
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

#############################################################################

def makeGradient(u,a,d2=True):
  g = correlate(u,a,d2) # gradient
  #p = correlate(u,u) # illumination preconditioner
  p = fillfloat(1.0,nz,nx)
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

import socket
def getMarmousi():
  p = zerofloat(751,2301)
  if socket.gethostname()=='backus.Mines.EDU':
    read("/data/sluo/marmousi/marmousi.dat",p)
  else:
    read("/data/seis/marmousi/marmousi.dat",p)
  div(1000.0,p,p) # slowness
  #s = copy(701,401,50,1550,p)
  #s = copy(701,401,50,200,p)
  s = copy(701,401,50,1200,p)
  return s

def getGaussian(pupper=0.05,plower=0.05):
  """ pupper - percent slowness perturbation of upper anomaly """
  """ plower - percent slowness perturbation of lower anomaly """
  s = fillfloat(4.0,701,401);
  div(1.0,s,s); s0 = s[0][0]
  t = zerofloat(701,401)
  t[200][175] = -1.0
  t[200][525] = plower/pupper
  #RecursiveGaussianFilter(60.0).apply00(t,t)
  RecursiveGaussianFilter(0.20/dx).apply00(t,t)
  div(t,max(abs(t)),t)
  mul(t,pupper*s0,t)
  add(t,s,s)
  return s

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
    frame.setSize(850,1000)
  else:
    frame.setSize(1000,1000)
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
  #fileName = dataDir+name+".dat"
  fileName = name
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def write(fname,image):
  #fname = datDir+name+'.dat'
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
    print '%02d:%02d:%02ds'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
