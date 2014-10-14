"""
Reflectivity inversion
"""
from imports import *

#############################################################################

sz = Sampling(265,0.012,0.0)
sx = Sampling(767,0.012,0.0)
st = Sampling(2500,0.0015,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
#kxs,kzs = [0],[0]
kxs,kzs = [nx/2],[0]
#kxs,kzs = [nx-1],[0]
#kxs,kzs = rampint(8,30,26),fillint(0,26)
#kxs,kzs = rampint(3,20,39),fillint(0,39)
#kxs,kzs = rampint(1,15,52),fillint(0,52)
#kxs,kzs = rampint(3,10,77),fillint(0,77)
kxr,kzr = rampint(0,1,nx),fillint(0,nx)
ns,nr = len(kxs),len(kxr)
fpeak = 10.0 # Ricker wavelet peak frequency
niter = 1

pngDir = None
datDir = None
pngDir = os.getenv('HOME')+'/Desktop/png/'
datDir = os.getenv('HOME')+'/Desktop/dat/'

sfile = None
#sfile = '/home/sluo/Desktop/10hz_new/cmb2/dat/s_iter1.dat'

parallel = False
#############################################################################

def main(args):
  setModel(sfile)
  #showData()
  #goWaveform()
  #goTraveltime()
  #goCombined()
  goInversion()

def goInversion():
  warping = True
  tb = fillfloat(_t[0][0],nz,nx) # background velocity for true model
  sb = fillfloat(_s[0][0],nz,nx) # background velocity for init model
  wave = Wavefield(sz,sx,st)

  # Observed and simulated data
  print 'data...'
  do = zerofloat(nt,nr)
  ds = zerofloat(nt,nr)
  da = zerofloat(nt,nr)
  source = Wavefield.RickerSource(fpeak,kzs[0],kxs[0])
  receiver = Wavefield.Receiver(kzr,kxr)
  wave.modelAcousticData(source,receiver,_t,do)
  wave.modelAcousticData(source,receiver,tb,da)
  sub(do,da,do)
  wave.modelAcousticData(source,receiver,_s,ds)
  wave.modelAcousticData(source,receiver,sb,da)
  sub(ds,da,ds)

  if warping:
    v,dw,_ = warp(ds,do) # right order
    r = sub(ds,dw)
  else:
    r = sub(ds,do) # residual
  #GaussianTaper.apply(r,r)

  # Compute gradient in smooth background model
  print 'gradient...'
  u,a = zerofloat(nz,nx,nt),zerofloat(nz,nx,nt)
  source = Wavefield.RickerSource(fpeak,kzs[0],kxs[0])
  receiver = Wavefield.Receiver(kzr,kxr)
  wave.modelAcousticWavefield(source,sb,u)
  source = Wavefield.AdjointSource(dt,kzr,kxr,r)
  wave.modelAcousticWavefield(source,sb,a)
  g = correlate(u,a,d2=True)
  for ix in range(1,nx):
    add(g[0],g[ix],g[0])
  for ix in range(1,nx):
    copy(g[0],g[ix])

  #div(g,max(abs(g)),g)
  #mul(0.1,g,g)
  #plot(sub(div(g,sb),1.0),cmap=jet,title='reflectivity')

  dmin,dmax = min(min(do),min(ds)),max(max(do),max(ds))
  plot(do,cmin=dmin,cmax=dmax,title='do')
  plot(ds,cmin=dmin,cmax=dmax,title='ds')
  if warping:
    plot(dw,cmin=dmin,cmax=dmax,title='dw')
    plot(v,cmap=rwb,cmin=-max(abs(v)),cmax=max(abs(v)),title='shifts')
  plot(r,title='residual')
  plot(g,cmap=jet,title='gradient')

#############################################################################

def setModel(ffile=None):
  global _t,_s
  #_t,_s = getMarmousi(2.000)
  #_t,_s = getLayered(0.0,0.05)
  _t,_s = getLayered(0.05,0.0)
  if ffile is not None:
    _s = read(ffile)
  plot(_t,cmap=jet,cbar='Slowness (s/km)',title='s_true')
  plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
       title='s_init')

def showData():
  do = zerofloat(nt,nr,ns) # observed data
  ds = zerofloat(nt,nr,ns) # simulated data
  da = zerofloat(nt,nr,ns) # direct arrival
  rr = zerofloat(nt,nr) # residual
  sw = Stopwatch(); sw.restart()
  Parallel.loop(ns,DataP(_t,do))
  Parallel.loop(ns,DataP(fillfloat(_t[0][0],nz,nx),da))
  Parallel.loop(ns,DataP(_s,ds))
  sw.stop(); print 'data:',sw.time(),'s'
  do = do[ns/2]
  ds = ds[ns/2]
  da = da[ns/2]
  sub(do,da,do) # observed
  sub(ds,da,ds) # simulated
  sub(ds,do,rr) # residual
  agc(do,ds)

#  s,dw,_ = warp(ds,do) # right order
#  rt = mul(s,timeDerivative(ds)) # traveltime residual
#  rw = sub(ds,dw) # warped residual
#  rc = add(rt,rw) # combined residual (traveltime+warped)
#  """
#  s,dw,_ = warp(do,ds) # wrong order
#  rt = mul(mul(-1.0,s),timeDerivative(ds)) # traveltime residual
#  rw = sub(dw,do) # warped residual
#  rc = add(rt,rw) # combined residual (traveltime+warped)
#  """

  s,dw = like(do),like(do)
  rt,rw,rc = makeWarpedResiduals(do,ds,s,dw)

  smax = 0.9*max(abs(s))
  rmin,rmax = 0.9*min(rr),0.9*max(rr)
  dmin,dmax = 0.9*min(min(do),min(ds)),0.9*max(max(do),max(ds))
  plot(do,cmin=dmin,cmax=dmax,title='observed')
  plot(ds,cmin=dmin,cmax=dmax,title='simulated')
  plot(dw,cmin=dmin,cmax=dmax,title='warped')
  plot(s,cmap=rwb,cmin=-smax,cmax=smax,title='shifts')
  plot(rr,cmin=rmin,cmax=rmax,title='residual')
  plot(rw,cmin=rmin,cmax=rmax,title='warped residual')
  plot(rt,cmin=rmin,cmax=rmax,title='traveltime residual')
  plot(rc,title='traveltime+warped residual')

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

#############################################################################
# Line Search

def findStepLength(g,isou,do,da=None):
  print 'searching for step length...'
  a = -0.5*max(abs(sub(_t,_s)))
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

# TODO: traveltime misfit function

#############################################################################
# Waveform

def goWaveform():
  d = modelData(_t) # observed data
  for iter in range(niter):
    print '\niteration',iter
    sw = Stopwatch(); sw.restart()
    g = waveformGradientP(d) if parallel else waveformGradientS(d)

    #gg = zerofloat(nz,nx,ns)
    #Parallel.loop(ns,WaveformParallelLoop(gg))
    #li = LinearInterpolator()
    #li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT)
    #li.setUniform(nz,1.0,0.0,nx,1.0,0.0,ns,kzs[1]-kzs[0],kzs[0],gg)
    #print 'interpolating...'
    #class InterpolateParallel(Parallel.ReduceInt):
    #  def compute(self,js):
    #    g = zerofloat(nz,nx)
    #    for ix in range(nx):
    #      for iz in range(nz):
    #        g[ix][iz] = li.interpolate(iz,ix,js)
    #    return g
    #  def combine(self,g1,g2):
    #    return add(g1,g2)
    #g = Parallel.reduce(nz,InterpolateParallel())

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
  maskWaterLayer(g)
  div(g,max(abs(g)),g)
  return g

def waveformGradientP(d):
  class Reducer(Parallel.ReduceInt):
    def __init__(self,d):
      self.d = d # observed data (precomputed)
      self.ds = Parallel.Unsafe(); # simulated data
      self.u = Parallel.Unsafe(); # source wavefield
      self.a = Parallel.Unsafe(); # receiver wavefield
    def compute(self,isou):
      ds,u,a = self.initializeUnsafe()
      return waveformGradient(isou,d[isou],ds,u,a)
    def xcompute(self,isou):
      sw = Stopwatch(); sw.restart()
      wave = Wavefield(sz,sx,st)
      ds,u,a = self.initializeUnsafe()
      #print sw.time(),'(1)'
      source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
      receiver = Wavefield.Receiver(kzr,kxr)
      wave.modelAcousticDataAndWavefield(source,receiver,_s,ds,u)
      #print sw.time(),'(2)'
      r = sub(ds,self.d[isou]) # residual
      source = Wavefield.AdjointSource(dt,kzr,kxr,r)
      wave.modelAcousticWavefield(source,_s,a)
      #print sw.time(),'(3)'
      g = makeGradient(u,a)
      #print sw.time(),'(4)'
      #return makeGradient(u,a)
      return g
    def initializeUnsafe(self):
      ds,u,a = self.ds.get(),self.u.get(),self.a.get()
      if ds is None:
        ds = zerofloat(nt,nr)
        self.ds.set(ds)
      if u is None:
        u = zerofloat(nz,nx,nt)
        self.u.set(u)
      if a is None:
        a = zerofloat(nz,nx,nt)
        self.a.set(a)
      return ds,u,a
    def combine(self,g1,g2):
      return add(g1,g2)
  g = Parallel.reduce(ns,Reducer(d))
  maskWaterLayer(g)
  div(g,max(abs(g)),g)
  return g


# TODO
class WaveformGradientChunker():
  def __init__(self,chunksize):
    pass

class WaveformGradientChunkP(Parallel.ReduceInt):
  def __init__(self):
    pass
  def compute(self,ichunk):
    pass
  def combine(self,chunk1,chunk2):
    pass


def waveformGradient(isou,do,ds,u,a):
  wave = Wavefield(sz,sx,st)
  #ds,u,a = zerofloat(nt,nr),zerofloat(nz,nx,nt),zerofloat(nz,nx,nt)
  source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
  receiver = Wavefield.Receiver(kzr,kxr)
  wave.modelAcousticDataAndWavefield(source,receiver,_s,ds,u)
  r = sub(ds,do) # residual
  source = Wavefield.AdjointSource(dt,kzr,kxr,r)
  wave.modelAcousticWavefield(source,_s,a)
  return makeGradient(u,a)

def maskWaterLayer(g,value=0.0):
  t0 = _t[0][0]
  for ix in range(nx):
    iz = 0
    while _t[ix][iz]==t0:
        g[ix][iz] = value
        iz += 1

#############################################################################
# Traveltime

def goTraveltime():
  do = modelData(_t) # observed data
  da = modelData(fillfloat(_t[0][0],nz,nx)) # direct arrival
  sub(do,da,do)
  g = zerofloat(nz,nx)
  for iter in range(niter):
    print '\niteration',iter
    sw = Stopwatch(); sw.restart()
    g = traveltimeGradientS(do,da)
    sw.stop(); print 'gradient:',sw.time(),'s'
    plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='g_iter'+str(iter))
    if niter>1:
      step = findStepLength(g,ns/2,do,da)
      add(mul(step,g),_s,_s)
      plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
           title='s_iter'+str(iter))

def traveltimeGradientS(do,da):
  g = zerofloat(nz,nx)
  ds,u,a = zerofloat(nt,nr),zerofloat(nz,nx,nt),zerofloat(nz,nx,nt)
  for isou in range(ns):
    print 'isou =',isou
    num,den = traveltimeGradient(isou,do[isou],da[isou],ds,u,a)
    add(num,g,g)
    add(den,p,p)
  div(g,p,g)
  maskWaterLayer(g)
  div(g,max(abs(g)),g)
  return g

def traveltimeGradient(isou,do,da,ds,u,a):
  wave = Wavefield(sz,sx,st)
  source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
  receiver = Wavefield.Receiver(kzr,kxr)
  wave.modelAcousticDataAndWavefield(source,receiver,_s,ds,u)
  sub(ds,da,ds) # subtract direct arrival

  # Residual
  """
  s,dw,_ = warp(ds,do) # order matters!
  #r = mul(s,timeDerivative(dw))
  r = mul(s,timeDerivative(ds))
  """
  s,dw,_ = warp(do,ds) # wrong order
  r = mul(mul(-1.0,s),timeDerivative(ds))

  source = Wavefield.AdjointSource(dt,kzr,kxr,r)
  wave.modelAcousticWavefield(source,_s,a)
  return makeGradient(u,a)

#############################################################################
# Combined

def goCombined():
  do = modelData(_t) # observed data
  da = modelData(fillfloat(_t[0][0],nz,nx)) # direct arrival
  sub(do,da,do) # subtract direct arrival

  # XXX
  #g = read('/home/sluo/Desktop/10hz_new/cmb/dat/g_iter0.dat')
  #step = findStepLength(g,ns/2,do,da)
  #add(mul(step,g),_s,_s)
  #plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='g_init0')
  #plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
  #     title='s_init0')

  g = zerofloat(nz,nx)
  for iter in range(niter):
    print '\niteration',iter
    sw = Stopwatch(); sw.restart()
    g = combinedGradientS(do,da)
    sw.stop(); print 'gradient:',sw.time(),'s'
    plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='g_iter'+str(iter))
    if niter>1:
      step = findStepLength(g,ns/2,do,da)
      add(mul(step,g),_s,_s)
      plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
           title='s_iter'+str(iter))

def combinedGradientS(do,da):
  g = zerofloat(nz,nx)
  p = zerofloat(nz,nx) # preconditioner
  ds,u,a = zerofloat(nt,nr),zerofloat(nz,nx,nt),zerofloat(nz,nx,nt)
  for isou in range(ns):
    print 'isou =',isou
    num,den = combinedGradient(isou,do[isou],da[isou],ds,u,a)
    add(num,g,g)
    add(den,p,p)
  div(g,p,g)
  maskWaterLayer(g)
  div(g,max(abs(g)),g)
  return g

def combinedGradient(isou,do,da,ds,u,a):
  wave = Wavefield(sz,sx,st)
  source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
  receiver = Wavefield.Receiver(kzr,kxr)
  wave.modelAcousticDataAndWavefield(source,receiver,_s,ds,u)
  sub(ds,da,ds) # subtract direct arrival

#  # Residual
#  s,dw,_ = warp(do,ds) # wrong order
#  rt = mul(mul(-1.0,s),timeDerivative(ds)) # traveltime residual
#  rw = sub(dw,do) # warped residual
#  rc = add(rt,rw) # combined

  rt,rw,rc = makeWarpedResiduals(do,ds)

  source = Wavefield.AdjointSource(dt,kzr,kxr,rc)
  wave.modelAcousticWavefield(source,_s,a)
  return makeGradient(u,a)

#############################################################################
# Dynamic Warping

def agc(do,ds):
  average = False
  n1,n2 = len(do[0]),len(do)
  ro = mul(do,do)
  rs = mul(ds,ds)
  #rs = mul(do,do)

  rgf = RecursiveGaussianFilter(80.0)
  rgf.apply00(ro,ro)
  rgf.apply00(rs,rs)
  """
  ref = RecursiveExponentialFilter(50.0)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  for i in range(4):
    ref.apply(ro,ro)
    ref.apply(rs,rs)
  """

  div(ro,max(ro),ro)
  div(rs,max(rs),rs)
  add(0.01,ro,ro)
  add(0.01,rs,rs)
  if average:
    t = mul(0.5,add(ro,rs))
    copy(t,ro)
    copy(t,rs)
  eo = div(do,ro)
  es = div(ds,rs)
  #div(eo,sqrt(sum(mul(eo,eo))/n1/n2),eo) # divide
  #div(es,sqrt(sum(mul(es,es))/n1/n2),es) # by rms

  #plot(ro,cmap=jet,title='ro')
  #plot(rs,cmap=jet,title='rs')
  #plot(do,title='do')
  #plot(ds,title='ds')
  #plot(eo,title='eo')
  #plot(es,title='es')
  return eo,es

def makeWarpedResiduals(do,ds,u=None,wd=None):
  v,dw,_ = warp(do,ds) # wrong order
  mul(-1.0,v,v)
  rt = mul(v,timeDerivative(ds)) # traveltime residual
  rw = sub(dw,do) # warped residual
  """
  v,dw,rt = warp(ds,do) # right order
  rw = sub(ds,dw)
  """

  print 'max shift =',max(abs(v))
  print 'max combined residual =',max(abs(add(rt,rw)))
  if u is not None:
    copy(v,u)
  if wd is not None:
    copy(dw,wd)
  return rt,rw,add(rt,rw)

def warp(p,q):
  qc = copy(q)
  p,q = agc(p,q)
  td = 1 # time decimation
  rd = 1 # receiver decimation
  a,b = addRandomNoise(2.0,p,q,sigma=1.0)
  f = copy(nt/td,nr/rd,0,0,td,rd,a)
  g = copy(nt/td,nr/rd,0,0,td,rd,b)
  #shiftMax = int(400/td)
  shiftMax = int(200/td)
  strainMax1 = 0.25
  strainMax2 = 0.25
  dw = DynamicWarping(-shiftMax,shiftMax)
  dw.setErrorExponent(1.0)
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
  d = 2 # decimate by this factor
  a,b = addRandomNoise(1.0e1,p,q,sigma=2.0)
  f = copy(nt/d,nr/d,0,0,d,d,a)
  g = copy(nt/d,nr/d,0,0,d,d,b)
  shiftMax = int(600/d)
  strainMax1 = 0.50
  strainMax2 = 0.25
  dw = DynamicWarping(-shiftMax,shiftMax)
  dw.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.REFLECT)
  dw.setStrainMax(strainMax1,strainMax2)
  dw.setShiftSmoothing(32.0/d) # shift smoothing
  #dw.setShiftSmoothing(1.0) # shift smoothing
  dw.setErrorSmoothing(2) # number of smoothings of alignment errors
  sw = Stopwatch(); sw.restart()
  s = dw.findShifts(f,g)
  sw.stop(); print "warping:",sw.time(),"s"
  mul(d,s,s) # scale shifts to compensate for decimation
  li = LinearInterpolator()
  li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT)
  li.setUniform(nt/d,d*dt,0.0,nr/d,d*dz,0.0,s)
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

def warpLsf(p,q):
  d = 2 # decimate
  sigma1 = 50.0 # local correlation window
  sigma2 = 10.0 # local correlation window
  maxShift = 400 # maximum shift
  lsf = LocalShiftFinder(sigma1,sigma2)
  f = copy(nt/d,nr,0,0,d,1,p)
  g = copy(nt/d,nr,0,0,d,1,q)
  u = like(f)
  lsf.find1(-maxShift,maxShift,f,g,u)
  mul(d,u,u)
  li = LinearInterpolator()
  li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT)
  li.setUniform(nt/d,d*dt,0.0,nr,dx,0.0,u)
  r = like(p)
  for ir in range(nr):
    z = ir*dz
    for it in range(nt):
      t = it*dt
      r[ir][it] = li.interpolate(t,z) # interpolate shifts
  h = DynamicWarping(-1,1).applyShifts(r,q)
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

import socket
def getMarmousi(sigmaC=0.5,sigmaV=0.0):
  p = zerofloat(751,2301)
  if socket.gethostname()=='backus.Mines.EDU':
    read("/data/sluo/marmousi/marmousi.dat",p)
  else:
    read("/data/seis/marmousi/marmousi.dat",p)
  p = copy(743,2301,8,0,p)
  div(1000.0,p,p) # slowness
  q = copy(p)
  refC = RecursiveExponentialFilter(sigmaC/0.004)
  refC.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  refC.apply(q,q)
  if sigmaV>0.0:
    refV = RecursiveExponentialFilter(sigmaV/0.004)
    refV.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
    refV.apply(p,p)
  v = fillfloat(2.0/3.0,nz,nx)
  c = copy(v)
  #copy(372,1151,0,0,2,2,p,13,0,1,1,v)
  #copy(372,1151,0,0,2,2,q,13,0,1,1,c)
  copy(248,767,0,0,3,3,p,17,0,1,1,v)
  copy(248,767,0,0,3,3,q,17,0,1,1,c)
  return v,c

def getLayered(errorTop=0.0,errorBottom=0.01):
  stop = 0.4
  sbot = 0.2
  t = fillfloat(stop,nz,nx)
  s = fillfloat(stop*(1.0-errorTop),nz,nx)
  for ix in range(nx):
    for iz in range(2*nz/3,nz):
      t[ix][iz] = sbot
      s[ix][iz] = sbot*(1.0-errorBottom)
  return t,s

def like(x):
  return zerofloat(len(x[0]),len(x))

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
      cleanDir(pngDir)
    if datDir is not None:
      cleanDir(datDir)
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02ds'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
