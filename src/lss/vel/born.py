"""
Born modeling
"""
from imports import *

#############################################################################

sz = Sampling(265,0.012,0.0)
sx = Sampling(767,0.012,0.0)
st = Sampling(2500,0.0015,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
#kxs,kzs = [0],[0]
#kxs,kzs = [nx/4],[0]
kxs,kzs = [nx/2],[0]
#kxs,kzs = [nx-1],[0]
#kxs,kzs = rampint(8,30,26),fillint(0,26)
#kxs,kzs = rampint(3,20,39),fillint(0,39)
#kxs,kzs = rampint(1,15,52),fillint(0,52)
#kxs,kzs = rampint(3,10,77),fillint(0,77)
kxr,kzr = rampint(0,1,nx),fillint(0,nx)
#kxr,kzr = [3*nx/4],[0]
ns,nr = len(kxs),len(kxr)
fpeak = 10.0 # Ricker wavelet peak frequency
niter = 1

pngDir = None
datDir = None
#pngDir = os.getenv('HOME')+'/Desktop/png/'
#datDir = os.getenv('HOME')+'/Desktop/dat/'

sfile = None
#sfile = '/home/sluo/Desktop/10hz_new/cmb2/dat/s_iter1.dat'

parallel = False
#############################################################################

def main(args):
  setModel()
  #showData()

  #goData()
  #goBorn()
  #goMigration()
  #goMigrationFakeBorn()
  #goBornMigration()
  goWarping()

#############################################################################
# Dynamic Warping

def goWarping():
  isou = ns/2
  wave = Wavefield(sz,sx,st)
  do = sub(modelData(_t),modelData(fillfloat(_t[0][0],nz,nx)))
  ds = sub(modelData(_s),modelData(fillfloat(_s[0][0],nz,nx)))
  do = do[isou]
  ds = ds[isou]
  ra = sub(ds,do)
  s,dw = like(do),like(do)
  rt,rw,rc = makeWarpedResiduals(do,ds,s,dw)
  plot(do,title='do')
  plot(ds,title='ds')
  plot(dw,title='dw')
  plot(s,cmap=jet,title='shifts')
  rmin,rmax = 1.0*min(min(rt),min(rw)),1.0*max(max(rw),max(rt))
  plot(ra,cmin=rmin,cmax=rmax,title='residual')
  plot(rw,cmin=rmin,cmax=rmax,title='warped residual')
  plot(rt,cmin=rmin,cmax=rmax,title='traveltime residual')

def makeWarpedResiduals(do,ds,u=None,wd=None):
  reverseOrder = False
  if reverseOrder:
    v,dw,_ = warp(do,ds) # wrong order
    mul(-1.0,v,v)
    rt = mul(v,timeDerivative(ds)) # traveltime residual
    rw = sub(dw,do) # warped residual
  else:
    v,dw,rt = warp(ds,do) # right order
    #rt = mul(v,timeDerivative(ds)) # traveltime residual
    rw = sub(ds,dw) # warped residual
  #print '  max shift =',max(abs(v))
  #print '  max combined residual =',max(abs(add(rt,rw)))
  if u is not None:
    copy(v,u)
  if wd is not None:
    copy(dw,wd)
  return rt,rw,add(rt,rw)

def warp(p,q):
  td = 4 # time decimation
  rd = 1 # receiver decimation
  qc = copy(q)
  p,q = agc(p,q)
  a,b = addRandomNoise(10.0,p,q,sigma=1.0)
  f = copy(nt/td,nr/rd,0,0,td,rd,a)
  g = copy(nt/td,nr/rd,0,0,td,rd,b)
  shiftMax = int(400/td)
  #strainMax1 = 1.00
  #strainMax2 = 0.50
  strainMax1 = 0.50
  strainMax2 = 0.10
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

#############################################################################

def goData():
  do = modelData(_t); do = do[ns/2]
  da = modelData(_s); da = da[ns/2]
  sub(do,da,do)
  plot(do,title='data')

def goMigration():
  isou = ns/2
  u = zerofloat(nz,nx,nt)
  a = zerofloat(nz,nx,nt)

  wave = Wavefield(sz,sx,st)
  #do = sub(modelData(_t),modelData(_s)); do = do[isou]; plot(do,title='data')
  do = sub(modelData(_t),modelData(fillfloat(_t[0][0],nz,nx)))
  ds = sub(modelData(_s),modelData(fillfloat(_s[0][0],nz,nx)))
  SimplePlot.asPoints(do[isou][0])
  SimplePlot.asPoints(ds[isou][0])
  do = do[isou]
  ds = ds[isou]
  sub(do,ds,do)

  # Forward
  print 'forward'
  source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
  wave.modelAcousticWavefield(source,_s,u)

  # Reverse
  print 'reverse'
  source = Wavefield.AdjointSource(dt,kzr,kxr,do)
  wave.modelAcousticWavefield(source,_s,a)

  # Image
  r = correlate(u,a)

  #plot(r,title='r')
  div(r,max(abs(r)),r)
  plot(r,cmap=rwb,cmin=-0.95,cmax=0.95,title='r')

def goMigrationFakeBorn():
  isou = ns/2
  u0 = zerofloat(nz,nx,nt)
  u1 = zerofloat(nz,nx,nt)
  a0 = zerofloat(nz,nx,nt)
  a1 = zerofloat(nz,nx,nt)
  s0,s1 = makeBornModel(_t)

  wave = Wavefield(sz,sx,st)
  #do = sub(modelData(_t),modelData(_s)); do = do[isou]; plot(do,title='data')
  do = sub(modelData(_t),modelData(fillfloat(_t[0][0],nz,nx)))
  ds = sub(modelData(_s),modelData(fillfloat(_s[0][0],nz,nx)))
  do = do[isou]
  ds = ds[isou]
  sub(do,ds,do)

  # Forward
  print 'forward'
  source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
  wave.modelAcousticWavefield(source,_s,u0) # u0
  wave.modelAcousticWavefield(source,_t,u1); sub(u1,u0,u1) # u1

  # Reverse
  print 'reverse'
  source = Wavefield.AdjointSource(dt,kzr,kxr,do)
  wave.modelAcousticWavefield(source,_s,a0) # a0
  wave.modelAcousticWavefield(source,_t,a1); sub(a1,a0,a1) # a1

  # Images
  r00 = correlate(u0,a0) # down/up
  r01 = correlate(u0,a1) # down/down
  r10 = correlate(u1,a0) # up/up
  r11 = correlate(u1,a1) # up/down

  plot(s0,cmap=jet,title='s0')
  plot(s1,cmap=jet,title='s1')
  #plot(r00,title='r00')
  #plot(r01,title='r10')
  #plot(r10,title='r01')
  #plot(r11,title='r11')
  rmax = max(max(max(abs(r00)),max(abs(r01))),max(max(abs(r10)),max(abs(r11))))
  div(r00,rmax,r00)
  div(r01,rmax,r01)
  div(r10,rmax,r10)
  div(r11,rmax,r11)
  plot(r00,cmap=rwb,cmin=-0.95,cmax=0.95,title='r00')
  plot(r01,cmap=rwb,cmin=-0.95,cmax=0.95,title='r10')
  plot(r10,cmap=rwb,cmin=-0.95,cmax=0.95,title='r01')
  plot(r11,cmap=rwb,cmin=-0.95,cmax=0.95,title='r11')


def goBorn():
  s0,s1 = makeBornModel(_t)
  do = modelBornData(s0,s1); do = do[ns/2]
  plot(s0,title='s0')
  plot(s1,title='s1')
  plot(do,title='born_data')

def goBornMigration():
  isou = ns/2
  u0 = zerofloat(nz,nx,nt)
  u1 = zerofloat(nz,nx,nt)
  a0 = zerofloat(nz,nx,nt)
  a1 = zerofloat(nz,nx,nt)
  t0,t1 = makeBornModel(_t)
  s0,s1 = makeBornModel(_s)

  wave = Wavefield(sz,sx,st)
  do = modelBornData(t0,t1); do = do[isou]; plot(do,title='observed_born_data')
  ds = modelBornData(s0,s1); ds = ds[isou]; plot(ds,title='simulated_born_data')
  sub(do,ds,do)

  # Forward
  print 'forward'
  source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
  wave.modelBornWavefield(source,s0,s1,u0,u1)

  # Reverse
  print 'reverse'
  source = Wavefield.AdjointSource(dt,kzr,kxr,do)
  wave.modelBornWavefield(source,s0,s1,a0,a1)

  # Images
  r00 = correlate(u0,a0) # down/up
  r01 = correlate(u0,a1) # down/down
  r10 = correlate(u1,a0) # up/up
  r11 = correlate(u1,a1) # up/down

  smin,smax = min(min(_t),min(_s)),max(max(_t),max(_s))
  plot(_t,cmap=jet,cmin=smin,cmax=smax,title='t')
  plot(t0,cmap=jet,cmin=smin,cmax=smax,title='t0')
  plot(t1,cmap=jet,title='t1')
  plot(_s,cmap=jet,cmin=smin,cmax=smax,title='s')
  plot(s0,cmap=jet,cmin=smin,cmax=smax,title='s0')
  plot(s1,cmap=jet,title='s1')
  #plot(r00,title='r00')
  #plot(r01,title='r10')
  #plot(r10,title='r01')
  #plot(r11,title='r11')
  rmax = max(max(max(abs(r00)),max(abs(r01))),max(max(abs(r10)),max(abs(r11))))
  div(r00,rmax,r00)
  div(r01,rmax,r01)
  div(r10,rmax,r10)
  div(r11,rmax,r11)
  plot(r00,cmap=rwb,cmin=-0.95,cmax=0.95,title='r00')
  plot(r01,cmap=rwb,cmin=-0.95,cmax=0.95,title='r10')
  plot(r10,cmap=rwb,cmin=-0.95,cmax=0.95,title='r01')
  plot(r11,cmap=rwb,cmin=-0.95,cmax=0.95,title='r11')

def makeBornModel(s):
  sigma = 0.1
  s0,s1 = like(s),like(s)
  ref = RecursiveExponentialFilter(sigma/dx)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref.apply(s,s0)
  sub(s,s0,s1)
  for i in range(5): # extra smoothing
    ref.apply(s0,s0)
    pass
  s0 = fillfloat(s[0][0],nz,nx)
  return s0,s1

#############################################################################

def setModel():
  global _t,_s
  #_t,_s = getMarmousi(2.000)
  _t,_s = getLayered()
  if sfile is not None:
    _s = read(sfile)
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

def modelBornData(s0,s1):
  d = zerofloat(nt,nr,ns)
  sw = Stopwatch(); sw.restart()
  Parallel.loop(ns,BornDataP(s0,s1,d))
  sw.stop(); print 'born data:',sw.time(),'s'
  return d

class BornDataP(Parallel.LoopInt):
  def __init__(self,s0,s1,d):
    self.s0 = s0 # background slowness
    self.s1 = s1 # slowness perturbation
    self.d = d # output array
  def compute(self,isou):
    wave = Wavefield(sz,sx,st)
    source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
    receiver = Wavefield.Receiver(kzr,kxr)
    wave.modelBornData(source,receiver,self.s0,self.s1,self.d[isou])

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
  rd = 2 # receiver decimation
  a,b = addRandomNoise(5.0,p,q,sigma=1.0)
  f = copy(nt/td,nr/rd,0,0,td,rd,a)
  g = copy(nt/td,nr/rd,0,0,td,rd,b)
  #shiftMax = int(400/td)
  shiftMax = int(200/td)
  strainMax1 = 0.50
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

def getLayered():
  t = fillfloat(0.4,nz,nx)
  s = copy(t)
  for ix in range(nx):
    for iz in range(2*nz/3,nz):
      t[ix][iz] = 0.2
      s[ix][iz] = 0.3
  """
  t = fillfloat(0.4,nz,nx)
  s = mul(0.95,t)
  for ix in range(nx):
    for iz in range(2*nz/3,nz):
      t[ix][iz] = 0.2
      s[ix][iz] = 0.2
  """

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
    #frame.setSize(1200,1000)
    frame.setSize(1200,600)
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
