"""
Combined waveform and traveltime inversion
"""
from imports import *

#############################################################################

sz = Sampling(265,0.012,0.0)
sx = Sampling(767,0.012,0.0)
#st = Sampling(4000,0.0015,0.0)
st = Sampling(4000,0.0016,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
#kxs,kzs = [0],[0]
#kxs,kzs = [nx/2],[0]
#kxs,kzs = [nx-1],[0]
#kxs,kzs = rampint(8,30,26),fillint(0,26)
#kxs,kzs = rampint(3,20,39),fillint(0,39)
#kxs,kzs = rampint(1,15,52),fillint(0,52)
kxs,kzs = rampint(3,10,77),fillint(0,77)
kxr,kzr = rampint(0,1,nx),fillint(0,nx)
ns,nr = len(kxs),len(kxr)
fpeak = 10.0 # Ricker wavelet peak frequency
niter = 10

pngDir = None
datDir = None
pngDir = os.getenv('HOME')+'/Desktop/png/'
datDir = os.getenv('HOME')+'/Desktop/dat/'
parallel = False

def main(args):
  setModel()
  #showData()
  #goWaveform()
  #goTraveltime()
  goCombined()
  #goMigration()

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

def setModel():
  global _t,_s
  #_t,_s = getMarmousi(0.200)
  #_t,_s = getMarmousi(1.000)
  _t,_s = getMarmousi(2.000)
  #_t,_s = getMarmousi(2.000,0.200)
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

  smax = 0.8*max(abs(s))
  rmin,rmax = 0.8*min(rr),0.8*max(rr)
  dmin,dmax = 0.8*min(min(do),min(ds)),0.8*max(max(do),max(ds))
  plot(do,cmin=dmin,cmax=dmax,title='observed')
  plot(ds,cmin=dmin,cmax=dmax,title='simulated')
  plot(dw,cmin=dmin,cmax=dmax,title='warped')
  plot(s,cmap=rwb,cmin=-smax,cmax=smax,title='shifts')
  plot(rr,cmin=rmin,cmax=rmax,title='residual')
  plot(rw,cmin=rmin,cmax=rmax,title='warped residual')
  plot(rt,cmin=rmin,cmax=rmax,title='traveltime residual')
  plot(rc,title='traveltime+warped residual')
  print 'diff =',sum(abs(sub(rc,rr)))

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
  #a = -1.0
  a = -2.0*max(abs(sub(_t,_s)))
  b =  0.0
  tol = 0.20*(b-a)
  sw = Stopwatch(); sw.restart()
  step = BrentMinFinder(WaveformMisfitFunction(g,isou,do,da)).findMin(a,b,tol)
  sw.stop()
  print 'a =',a
  #print 'b =',b
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

def makeWarpedResiduals(do,ds,u=None,wd=None):
  v,dw,_ = warp(do,ds) # wrong order
  mul(-1.0,v,v)
  rt = mul(v,timeDerivative(ds)) # traveltime residual
  rw = sub(dw,do) # warped residual
  """
  v,dw,rt = warp(ds,do) # right order
  rw = sub(ds,dw)
  """

  if u is not None:
    copy(v,u)
  if wd is not None:
    copy(dw,wd)
  return rt,rw,add(rt,rw)

def warp(p,q):
  td = 1 # time decimation
  rd = 1 # receiver decimation
  a,b = addRandomNoise(2.0,p,q,sigma=1.0)
  f = copy(nt/td,nr/rd,0,0,td,rd,a)
  g = copy(nt/td,nr/rd,0,0,td,rd,b)
  shiftMax = int(400/td)
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
  hr = dw.applyShifts(r,mul(r,timeDerivative(q)))
  #plot(s,cmap=jet,title='shifts')
  #plot(r,cmap=jet,title='shifts (interpolated)')
  #plot(f,perc=99.5,title='f')
  #plot(g,perc=99.5,title='g')
  #plot(a,perc=99.5,title='a')
  #plot(b,perc=99.5,title='b')
  #plot(h,perc=99.5,title='h')
  #plot(hr,perc=99.5,title='hr')
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
