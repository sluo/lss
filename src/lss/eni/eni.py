##############################################################################
# Eni data

from imports import *
from dnp import *

#subDir = '/data/sluo/eni/dat/suba/'
#subDir = '/data/sluo/eni/dat/subb/'
subDir = '/data/sluo/eni/dat/subc/'

savDir = None
#savDir = '/home/sluo/Desktop/pngdat/'
#savDir = '/home/sluo/Desktop/pngdat2/'
#savDir = '/home/sluo/Desktop/pngdat3/'

##############################################################################

def main(args):
  #showFiles()
  #readFiles()
  #goBornData()
  #goAcousticData()
  #resimulateData()
  #compareWavelets()
  #estimateWavelet(toFile=False)
  goAmplitudeInversionQs()

def setGlobals():
  global sx,sz,st#,ss,sr
  global nx,nz,nt,ns,nr
  global dx,dz,dt,ds,dr
  global fx,fz,ft,fs,fr
  global nxp,nzp,np
  global fmax,nabsorb,stride
  subset = subDir.split('/')[-2]
  if subset=='suba' or subset=='subb':
    ss = Sampling(453,0.0125,1.225) # shot (relative offset)
    sr = Sampling(197,0.00625,-1.225) # receiver
    sx = Sampling(1101,0.00625,0.0) # distance
    sz = Sampling(181,0.00625,0.0) # depth
    st = Sampling(3751,0.0004,0.0) # time
    npmax = 16 # max number of parallel shots
  elif subset=='subc':
    ss = Sampling(1503,0.0125,1.225) # shot (relative offset)
    sr = Sampling(197,0.00625,-1.225) # receiver
    sx = Sampling(3201,0.00625,0.0) # distance
    sz = Sampling(181,0.00625,0.0) # depth
    st = Sampling(3751,0.0004,0.0) # time
    npmax = 8 # max number of parallel shots
  #stride = 1
  stride = 2
  #stride = 500
  ns,ds,fs = int((ss.count+stride-1)/stride),ss.delta,ss.first
  nr,dr,fr = sr.count,sr.delta,sr.first
  nt,dt,ft = st.count,st.delta,st.first
  nz,dz,fz = sz.count,sz.delta,sz.first
  nx,dx,fx = sx.count,sx.delta,sx.first
  np = min(ns,npmax) # number of parallel shots
  fmax = 30.0 # max passband frequency
  nabsorb = 22 # absorbing boundary size
  nxp,nzp = nx+2*nabsorb,nz+2*nabsorb
  print 'subset=%r'%subset
  print 'ns=%r'%ns

def getSourceAndReceiver():
  w = getWavelet() # wavelet
  src = BornOperatorS.getSourceArray(ns) # source
  rcp = BornOperatorS.getReceiverArray(ns) # predicted data
  rco = BornOperatorS.getReceiverArray(ns) # observed data
  for isou in range(ns):
    xs = fs+isou*stride*ds
    kxs = sx.indexOf(xs)
    #print 'kxs=%d'%kxs
    src[isou] = Source.WaveletSource(kxs,0,w)
    kxr = zeroint(nr)
    kzr = zeroint(nr)
    for ir in range(nr):
      xr = xs+fr+ir*dr
      #print ' ',sx.indexOf(xr)
      kxr[ir] = sx.indexOf(xr)
    d = getGather(isou*stride)
    bandpass1(d,d,10.0,fmax) # bandpass filter
    #GaussianTaper.apply2(0.5,d,d) # taper
    e = timeDelay(2.0/fmax,d) # time delay to match first arrivals
    rco[isou] = Receiver(kxr,kzr,e)
    rcp[isou] = Receiver(kxr,kzr,len(e[0]))
  return src,rcp,rco

def getWavelet():
  return readWavelet()
  #return estimateWavelet()
  #return makeRickerWavelet() # Ricker wavelet

def goAmplitudeInversionQs():
  vz = True # 1D velocity?
  warp3d = True # 3D warping?
  nouter,ninner,nfinal = 5,2,5 # outer, inner, final after last outer
  #nouter,ninner,nfinal = 0,0,5 # outer, inner, final after last outer
  print 'vz=%r'%vz
  print 'warp3d=%r'%warp3d
  print 'nouter=%r'%nouter
  print 'ninner=%r'%ninner
  print 'nfinal=%r'%nfinal

  # BornSolver
  print "allocating..."
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)
  s,m = getBackgroundAndMask(vz)
  src,rcp,rco = getSourceAndReceiver()
  #sigma = 0.25*nx*nz/(sum(s)*fmax*dx*sqrt(2.0)) # quarter wavelength
  sigma = 0.25*averageWavelength()/(sqrt(2.0)*dx) # quarter wavelength
  ref = RecursiveExponentialFilter(sigma)
  born = BornOperatorS(s,dx,dt,nabsorb,u,a)
  bs = BornSolver(born,src,rcp,rco,ref,m)

  # DataWarping
  td = 4 # time decimation
  maxShift = 0.1 # max shift (seconds)
  strainT,strainR,strainS = 0.20,0.20,0.50*stride if warp3d else -1.0
  smoothT,smoothR,smoothS = 32.0,4.0,4.0
  warping = DataWarping(
    strainT,strainR,strainS,smoothT,smoothR,smoothS,maxShift,dt,td)
  u = zerofloat(nt,nr,ns) # warping shifts

  for iouter in range(nouter+1):
    sw = Stopwatch(); sw.start()
    r = bs.solve(nfinal if iouter==nouter else ninner);
    pixels(r,cmap=gray,title='r%d'%iouter)
    if iouter<nouter and nouter>0:
      born.applyForward(src,r,rcp)
      print 'warping (3d=%r)'%warp3d
      rcw = warping.warp(rcp,rco,u)
      bs.setObservedData(rcw)
      pixels(rcp[ns/2].getData(),title='rcp%d'%iouter)
      pixels(rcw[ns/2].getData(),title='rcw%d'%iouter)
      pixels(u[ns/2],cmap=rwb,sperc=100.0,title='shifts%d'%iouter)
      if iouter==(nouter-1) and savDir is not None:
        dp = zerofloat(nt,nr,ns)
        for isou in range(ns):
          copy(rcp[isou].getData(),dp[isou])
        write(savDir+'dp.dat',dp)
        write(savDir+'u.dat',u)
    sw.stop(); print 'outer: %.1f minutes'%(sw.time()/60.0)

  pixels(rco[ns/2].getData(),title='rco')
  pixels(s,cmap=jet,title='s')

def showFiles():
  w = readWavelet(); points(w)
  d = getGather(ns/2); pixels(d,perc=99.8); points(d[196])
  s = getSlowness(False); pixels(s,cmap=jet)
  s0,_ = getBackgroundAndMask(False); pixels(t0,cmap=jet)
  t0,m = getBackgroundAndMask(True); pixels(s0,cmap=jet); pixels(m,cmap=gray)
  pixels(sub(s0,t0),sperc=100.0)

def readFiles():
  ra = zerofloat(nz,nx)
  rb = zerofloat(nz,nx)
  read('/home/sluo/Desktop/save/eni/subc/iter525/r5.dat',ra);
  read('/home/sluo/Desktop/save/eni/subc/iter005/r0.dat',rb);
  pixels(ra,cmap=gray,sperc=98.0)
  pixels(rb,cmap=gray,sperc=98.0)
  """
  dp = zerofloat(nt,nr,ns)
  read('/home/sluo/Desktop/eni_save/less_time_delay/525_vz/dp.dat',dp)
  dp = dp[100]
  do = timeDelay(2.0/fmax,getGather(200))
  bandpass1(do,do,10.0,fmax)
  clip = max(max(abs(dp)),max(abs(do)))
  pixels(dp,cmin=-clip,cmax=clip)
  pixels(do,cmin=-clip,cmax=clip)
  """

##############################################################################
# wavelet

def compareWavelets():
  w = readWavelet()
  #w = estimateWavelet(rotate=0.25*FLT_PI)
  v = estimateWaterBottomWavelet()
  #mul(1.0/max(abs(w)),w,w)
  mul(1.0/max(abs(v)),v,v)
  wm = zeroint(1); max(w,wm); wm = wm[0]
  vm = zeroint(1); max(v,vm); vm = vm[0]
  t = copy(v); zero(v); copy(nt-(vm-wm),vm-wm,t,0,v)
  points(copy(400,w),cmin=-1.0,cmax=1.0)
  points(copy(400,v),cmin=-1.0,cmax=1.0)

def readWavelet():
  print 'reading wavelet'
  w = zerofloat(nt)
  read(subDir+'w.dat',w)
  #points(w)
  return w

def estimateWavelet(toFile=False,rotate=0.25*FLT_PI):
  print 'estimating wavelet'
  w = estimateZeroPhaseWavelet()
  #w = estimateMinimumPhaseWavelet()

  # Phase rotation and bandpass filter.
  w = Frequency(nt).rotatePhase(w,rotate)
  bandpass(w,w,10.0,fmax)

#  # Scale.
#  mul(1.0*max(abs(getGather(0)[196]))/max(abs(w)),w,w)
#  #mul(6.0*max(abs(getGather(0)[196]))/max(abs(w)),w,w)
#  #mul(1.0e-6*max(abs(getGather(0)[196]))/max(abs(w)),w,w)
#  print 1.0*max(abs(getGather(0)[196]))/max(abs(w))
#
#  # Additional scale.
#  mul(0.0337413899989,w,w)

#  # Scale
#  #mul(0.115,w,w)
#  mul(0.12,w,w)

  # Normalize
  mul(1.0/max(abs(w)),w,w)

  points(amplitudeSpectrum(w))
  points(w)
  if toFile:
    print 'writing wavelet'
    write(subDir+'w.dat',w)
  return w

def estimateWaterBottomWavelet():
  showPlots = False

  # Read data.
  irec = 196
  d = zerofloat(nt,ns)
  e = zerofloat(nt,ns)
  for isou in range(ns):
    copy(getGather(isou)[irec],d[isou])

  # Pick first breaks.
  #p = FirstBreaks(nt/2).pick(d)
  p = zerofloat(ns)
  for isou in range(ns):
    t = zeroint(1)
    max(d[isou],t)
    p[isou] = t[0]
  if showPlots:
    pixels(d,perc=100.0).\
      addPoints(p,rampfloat(0.0,1.0,ns)).\
      setLineColor(Color.YELLOW)

  # Flatten on first breaks.
  pmax = max(p)
  add(500-pmax,p,p)
  for isou in range(ns):
    pi = int(p[isou])
    copy(nt-pi,pi,d[isou],0,e[isou])
  if showPlots:
    pixels(e,perc=100.0)

  # Stack traces.
  w = zerofloat(nt)
  for isou in range(ns):
    add(e[isou],w,w)
  if showPlots:
    points(w)

  # Process wavelet.
  bandpass(w,w,10.0,fmax) # bandpass
  s = zerofloat(nt)
  i = zeroint(1)
  max(w,i)
  #for it in range(i[0]):
  #  s[it] = 1.0
  s[i[0]] = 1.0
  sigma = 1.0/fmax/dt
  RecursiveGaussianFilter(sigma).apply0(s,s)
  mul(1.0/max(s),s,s)
  if showPlots:
    points(s)
  mul(s,w,w)
  if showPlots:
    points(w)
  w = Frequency(nt).rotatePhase(w,0.25*FLT_PI)
  if showPlots:
    points(w)

  # Spectra.
  ts = Sampling(nt,dt,-i[0]*dt)
  ap = Spectrum.computeAmplitudeAndPhase(ts,w)
  a,p = ap[0],ap[1] # amplitude and phase
  if showPlots:
    points(a)
    points(p)

  return mul(1.0/max(abs(w)),w)

def pickFirstBreaks(d):
  p = zerofloat(ns)
  for isou in range(ns):
    t = zeroint(1)
    max(d[isou],t)
    p[isou] = t[0]
  return p

def estimateZeroPhaseWavelet(): # average spectra
  freq = Frequency(nt,4)
  c = freq.amplitudeSpectrum(getGather(0)[196])
  count = 1.0
  for isou in range(1,ns):
    #for irec in range(nr):
    for irec in range(nr-1,nr):
      add(freq.amplitudeSpectrum(getGather(isou)[irec]),c,c)
      count += 1.0
  #print count
  mul(1.0/count,c,c)
  #RecursiveGaussianFilter(32.0).apply0(c,c)
  RecursiveGaussianFilter(64.0).apply0(c,c)
  #points(d)
  #points(c)
  w = freq.findZeroPhase(c)
  return w

def xestimateZeroPhaseWavelet(): # average wavelets
  w = zerofloat(nt)
  count = 0.0
  for isou in range(ns):
    d = getGather(isou)
    t = zerofloat(nt)
    #bandpass(d[196],t)
    copy(d[196],t)
    c = amplitudeSpectrum(t)
    r = findZeroPhase(c)
    add(r,w,w)
    count += 1.0
  mul(1.0/count,w,w)
  z = zerofloat(nt)
  kt = int(2.0/fmax/dt) # offset (backwards from it=0)
  copy(kt,nt-1-kt,w,0,z)
  copy(nt-kt,0,w,kt,z)
  return w

def xestimateMinimumPhaseWavelet(): # average wavelets
  t = zerofloat(nt)
  count = 0.0
  for isou in range(ns):
    #print 'isou=%d'%isou
    d = getGather(isou)
    for ir in range(196,197):
      #bandpass(d[ir],t)
      copy(d[ir],t)
      c = squaredSpectrum(t)
      r = findMinimumPhase(c)
      add(r,w,w); count += 1.0
  mul(1.0/count,w,w)
  #bandpass(w,w)
  #points(squaredSpectrum(t))
  #points(squaredSpectrum(w))
  #points(w)
  return w

def estimateMinimumPhaseWavelet(): # average spectra
  freq = Frequency(nt,1)
  c = freq.squaredSpectrum(getGather(0)[196])
  count = 1.0
  for isou in range(1,ns):
    #print 'isou=%d'%isou
    add(freq.squaredSpectrum(getGather(isou)[196]),c,c)
    count += 1.0
  mul(1.0/count,c,c)
  points(c)
  RecursiveGaussianFilter(20.0).apply0(c,c)
  w = freq.findMinimumPhase(c)
  #points(squaredSpectrum(t))
  #points(squaredSpectrum(w))
  #points(w)
  return w

def makeRickerWavelet():
  print 'making Ricker wavelet'
  def ricker(t):
    x = FLT_PI*fmax*(t-2.0/fmax)
    xx = x*x
    return (1.0-2.0*xx)*exp(-xx);
  w = zerofloat(nt)
  for it in range(nt):
    t = ft+it*dt
    w[it] = ricker(t)
  w = Frequency(nt).rotatePhase(w,0.25*FLT_PI)
  mul(1.0/max(w),w,w)
  return w

class Frequency:
  def __init__(self,nt,pad=4):
    self.rx = zerofloat(pad*nt)
    self.fft = Fft(pad*nt)
    self.nt = nt
  def amplitudeSpectrum(self,rx):
    self.set(rx)
    cx = self.fft.applyForward(self.rx)
    return cabs(cx)
  def phaseSpectrum(self,rx):
    self.set(rx)
    cx = self.fft.applyForward(self.rx)
    return carg(cx)
  def squaredSpectrum(self,rx):
    cx = self.amplitudeSpectrum(rx)
    mul(cx,cx,cx)
    #add(1.0e-6*max(cx),cx,cx) # stabilize for logarithm
    return cx
  def findZeroPhase(self,cx):
    nc = len(cx)
    cx = cmplx(cx,zerofloat(nc)) # amplitude spectrum
    ry = self.fft.applyInverse(cx)
    nr = len(ry)
    #points(ry)
    z = zerofloat(self.nt)
    kt = int(2.0/fmax/dt) # offset (backwards from it=0)
    #copy(kt,nr-kt,ry,0,z)
    copy(self.nt-kt,0,ry,kt,z)
    for it in range(0,kt+1):
      z[kt-it] = z[kt+it]
    return z
  def findMinimumPhase(self,cx):
    cx = cmplx(cx,zerofloat(len(cx))) # spectrum
    cx = clog(cx)
    ry = self.fft.applyInverse(cx)
    nr = len(ry)
    ry[0   ] *= 0.5
    ry[nr/2] *= 0.5
    for it in range(nr/2+1,nr):
      ry[it] = 0.0
    cz = self.fft.applyForward(ry)
    cz = cexp(cz)
    rz = self.fft.applyInverse(cz)
    return copy(self.nt,rz)
  def rotatePhase(self,rx,p=0.25*FLT_PI):
    self.set(rx)
    sf = self.fft.getFrequencySampling1()
    nf = sf.count
    cy = self.fft.applyForward(self.rx) # forward FFT
    t = zerofloat(2*nf)
    for i in range(nf):
      w = sf.getValue(i)
      #t[2*i  ] = w*w*cos(p) # phase rotation and negative 2nd time derivative
      #t[2*i+1] = w*w*sin(p) # phase rotation and negative 2nd time derivative
      t[2*i  ] = cos(p) # phase rotation only
      t[2*i+1] = sin(p) # phase rotation only
    cmul(t,cy,cy)
    ry = self.fft.applyInverse(cy) # inverse FFT
    return copy(self.nt,ry)
  def set(self,rx):
    zero(self.rx)
    copy(rx,self.rx)

##############################################################################
# data

def goBornData(s0,s1):
  bwo = BornOperator(s0,dx,dt,nabsorb)
  b = zerofloat(nxp,nzp,nt)
  u = zerofloat(nxp,nzp,nt)
  src,rcp,rco = getSourceAndReceiver()
  sou = src[ns/2]
  rec = rco[ns/2]
  d = rec.getData()
  dobs = copy(d)
  pixels(dobs,cmap=gray,sperc=99.8,title="observed data")
  points(dobs[196])
  sw = Stopwatch(); sw.start()
  bwo.applyForward(sou,b,s1,rec,u)
  print sw.time()
  dpre = rec.getData()
  pixels(s0,cmap=jet,title="background slowness (s/km)")
  pixels(s1,sperc=100.0,title="reflectivity")
  pixels(dpre,cmap=gray,sperc=99.8,title="predicted data")
  points(dpre[196])
  return dpre,dobs

def goAcousticData(s=None):
  if s is None:
    s = getSlowness()
  awo = WaveOperator(s,dx,dt,nabsorb)
  src,rcp,rco = getSourceAndReceiver()
  sou = src[ns/2]
  rec = rcp[ns/2]
  u = zerofloat(nxp,nzp,nt)
  awo.applyForward(sou,rec,u)
  dpre = rec.getData()
  pixels(dpre,cmap=gray,sperc=99.8,title="predicted data")
  points(dpre[196])
  points(dpre[98])
  points(dpre[0])
  #pixels(u[1000])
  #pixels(u[1500])
  #pixels(u[2000])
  #pixels(u[2500])

def resimulateData():
  s0 = getBackground()
  r = zerofloat(nz,nx)
  read('/home/sluo/Desktop/eni_save/less_time_delay/5iter/r0.dat',r)
  s1 = transpose(r)
  m = getMask()
  mul(m,s1,s1)
  dp,do = goBornData(s0,s1)
  rmsFilter(dp,do,2.0/fmax)
  pixels(dp,cmap=gray,sperc=99.8,title="predicted data (rms filtered)")
  pixels(do,cmap=gray,sperc=99.8,title="observed data (rms filtered)")
  print "rms scale factor =",rms(do)/rms(dp)

def rmsFilter(ds,do,sigma):
  x,y = copy(ds),copy(do)
  rmsx = rms(x)
  rmsy = rms(y)
  # equalize rms
  if rmsx>rmsy:
    mul(rms(y)/rms(x),x,x)
  else:
    mul(rms(x)/rms(y),y,y)
  xx = mul(x,x)
  yy = mul(y,y)
  rgf = RecursiveGaussianFilter(sigma/dt)
  rgf.apply00(xx,xx)
  rgf.apply00(yy,yy)
  num = mul(mul(2.0,xx),yy)
  den = add(mul(xx,xx),mul(yy,yy))
  add(1.0e-6,den,den)
  div(num,den,den)
  #mul(den,x,x)
  #mul(den,y,y)
  mul(den,ds,ds)
  mul(den,do,do)
  #plot(den,cmap=jet,title='rms_weights')

##############################################################################
# dat

def getGather(isou):
  d = zerofloat(nt,nr)
  fname = subDir+'d_'+str(isou)+'.dat'
  read(fname,d)
  return d

def getVelocity(vz=False):
  v = zerofloat(nz,nx)
  read(subDir+'v.dat',v)
  #if vz: # 1D velocity model
  #  v00 = v[0][0]
  #  total = zerofloat(nz)
  #  count = zerofloat(nz)
  #  bottom = zeroint(nx)
  #  for ix in range(nx):
  #    bz = 0
  #    while v[ix][bz]==v00:
  #      bz += 1
  #    bz += 2 # extra
  #    for iz in range(bz,nz):
  #      total[iz] += v[ix][iz]
  #      count[iz] += 1.0
  #    bottom[ix] = bz
  #  div(total,count,total)
  #  for iz in range(min(bottom)):
  #    total[iz] = total[min(bottom)]
  #  RecursiveExponentialFilter(4.0).apply(total,total)
  #  for ix in range(nx):
  #    for iz in range(bottom[ix],nz):
  #      v[ix][iz] = total[iz]
  return transpose(v)

def getSlowness(vz=False):
  v = getVelocity(vz)
  div(1000.0,v,v)
  if vz: # linear slowness model
    v = transpose(v)
    v00 = v[0][0]
    a,b,c = 0.0,0.0,0.0
    for ix in range(int(0.31*nx),nx): # 0.31 just because
      x,y,xy,xx,n = 0.0,0.0,0.0,0.0,0.0
      for iz in range(nz):
        if v[ix][iz]!=v00:
          x += iz
          y += v[ix][iz]
          xy += iz*v[ix][iz]
          xx += iz*iz
          n += 1.0
      s = (xy-x*y/n)/(xx-x*x/n)
      b += s
      a += y/n-s*x/n
      c += 1.0
    a = a/c
    b = b/c
    for ix in range(nx):
      for iz in range(nz):
        if v[ix][iz]!=v00:
          v[ix][iz] = a+b*iz
    v = transpose(v)
  return v

def getBackgroundAndMask(vz=False):
  return getBackground(vz),getMask()

def getBackground(vz=False):
  s = getSlowness(vz)
  sigma = 0.5*averageWavelength() # half wavelength
  print 'sigma0=%f'%sigma
  #esmooth(sigma/dx,s,s)
  gsmooth(sigma/dx,s,s)
  return s

def getMask():
  v = zerofloat(nz,nx)
  read(subDir+'v.dat',v)
  m = fillfloat(1.0,nz,nx)
  v00 = v[0][0]
  sigma = 2
  for ix in range(nx):
    for iz in range(nz-6*sigma):
      if v[ix][iz+6*sigma]==v00:
        m[ix][iz] = 0.0
  #esmooth(sigma,m,m)
  gsmooth(sigma,m,m)
  return transpose(m)

def averageWavelength():
  ss = getSlowness(vz=False)
  s,n = 0.0,0.0
  s00 = ss[0][0]
  for iz in range(nz):
    for ix in range(nx):
      if ss[iz][ix]!=s00:
        s += ss[iz][ix]
        n += 1.0
  #print n/s/fmax
  return n/s/fmax

def esmooth(sigma,x,y):
  ref = RecursiveExponentialFilter(sigma)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref.apply(x,y)

def gsmooth(sigma,x,y):
  extend = int(6*sigma)
  n1,n2 = len(x[0]),len(x)
  e = zerofloat(n1+2*extend,n2+2*extend)
  copy(n1,n2,0,0,x,extend,extend,e)
  for i2 in range(extend,n2+extend):
    for i1 in range(extend):
      e[i2][i1] = e[i2][extend]
      e[i2][n1+extend+i1] = e[i2][n1+extend-1]
  for i2 in range(extend):
    copy(e[extend],e[i2])
    copy(e[n2+extend-1],e[n2+extend+i2])
  #pixels(e,cmap=jet)
  RecursiveGaussianFilter(2.0*sigma).apply0X(e,e) # more horizontal smoothing
  RecursiveGaussianFilter(    sigma).applyX0(e,e)
  copy(n1,n2,extend,extend,e,0,0,y)

##############################################################################
# util

def amplitudeSpectrum(x):
  return cabs(Fft(x).applyForward(x))

def bandpass(x,y,fmin=None,fmax=None):
  fnyq = 0.5/dt # Nyquist
  if fmin is None:
    klower = 0.01
  else:
    klower = fmin/fnyq
  if fmax is None:
    kupper = 40.0/fnyq
  else:
    kupper = fmax/fnyq
  bp = BandPassFilter(klower,kupper,0.01,0.01)
  #points(spectrum(bp.getCoefficients1())) # amplitude spectrum
  bp.apply(x,y)

def bandpass1(x,y,fmin=None,fmax=None):
  n2 = len(x)
  for i2 in range(n2):
    bandpass(x[i2],y[i2],fmin,fmax)

def timeDelay(t,d):
  e = zerofloat(nt,nr)
  kt = int(t/dt)
  copy(nt-kt,nr,0,0,d,kt,0,e)
  return e

def rms(x):
  return sqrt(sum(mul(x,x))/len(x[0])/len(x))

def like(x):
  return zerofloat(len(x[0]),len(x))

def read(name,image):
  fileName = name
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()

def write(fname,image):
  aos = ArrayOutputStream(fname)
  aos.writeFloats(image)
  aos.close()

def cleanDir(dir):
  os.chdir(dir)
  for f in os.listdir(dir):
    os.remove(f)

##############################################################################
# plots

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
gyr = ColorMap.GRAY_YELLOW_RED
def pixels(x,cmap=gray,perc=100.0,sperc=None,cmin=0.0,cmax=0.0,title=None):
  if (len(x)==nz):
    x = transpose(x)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  cb = sp.addColorBar()
  cb.setWidthMinimum(100)
  sp.setSize(1010,740)
  if title:
    sp.addTitle(title)
  pv = sp.addPixels(x)
  pv.setColorModel(cmap)
  #pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if perc<100.0:
    pv.setPercentiles(100.0-perc,perc)
  if sperc is not None: # symmetric percentile clip (for plotting gradients)
    clips = Clips(100-sperc,sperc,x)
    clip = max(abs(clips.getClipMin()),abs(clips.getClipMax()))
    pv.setClips(-clip,clip)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if title and savDir:
    sp.paintToPng(360,3.33,savDir+title+'.png')
    write(savDir+title+'.dat',x)
  return sp

def points(x,cmin=0.0,cmax=0.0):
  sp = SimplePlot.asPoints(x)
  if cmin<cmax:
    sp.setVLimits(cmin,cmax)

##############################################################################
# Do everything on Swing thread.
class RunMain(Runnable):
  def run(self):
    start = time.time()
    if savDir is not None:
      print 'cleaning '+savDir.split('/')[-2]
      cleanDir(savDir)
    setGlobals()
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02d'%(h,m,s)
if __name__=='__main__':
  SwingUtilities.invokeLater(RunMain())
