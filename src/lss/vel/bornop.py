#############################################################################
# Wavefield modeling using BornWaveOperator

from imports import *
from dnp import *

#############################################################################

sz = Sampling(201,0.016,0.0)
sx = Sampling(202,0.016,0.0)
st = Sampling(2003,0.0010,0.0)
#sz = Sampling(401,0.0025,0.0) # for 40 Hz
#sx = Sampling(402,0.0025,0.0)
#st = Sampling(5003,0.0001,0.0)
#sz = Sampling(201,0.00625,0.0) # for 30 Hz
#sx = Sampling(202,0.00625,0.0)
#sz = Sampling(313,0.004,0.0) # for 30 Hz
#sx = Sampling(314,0.004,0.0)
#st = Sampling(3003,0.0004,0.0)
#sz = Sampling(265,0.012,0.0); stride = 3
#sx = Sampling(767,0.012,0.0)
#st = Sampling(5501,0.0010,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
fz,fx,ft = sz.first,sx.first,st.first
#xs,zs = [0],[0]
#xs,zs = [nx/2],[0]
#xs,zs = [nx/2],[nz/2]
xs,zs = [nx/4,nx/2,3*nx/4],[0,0,0]
#xs,zs = [nx/5,2*nx/5,3*nx/5,4*nx/5],[0,0,0,0]
#xs,zs = rampint(1,15,52),fillint(0,52) # for marmousi
#xs,zs = rampint(3,10,77),fillint(0,77) # for marmousi
#xs,zs = rampint(3,5,153),fillint(0,153) # for marmousi
xr,zr = rampint(0,1,nx),fillint(0,nx)
ns,nr = len(xs),len(xr)
fpeak = 10.0 # Ricker wavelet peak frequency
nabsorb = 22 # absorbing boundary size
nxp,nzp = nx+2*nabsorb,nz+2*nabsorb
np = min(16,ns) # number of parallel sources

pngdatDir = None
#pngdatDir = os.getenv('HOME')+'/Desktop/pngdat/'
#pngdatDir = os.getenv('HOME')+'/Desktop/pngdat2/'
#pngdatDir = os.getenv('HOME')+'/Desktop/pngdat3/'

# TODO: Test adjoint with TransformQuadratic

def main(args):
  #goAcousticData()
  #goBornData()
  #makeSource()
  #adjointTest()
  #adjointTestMultiSource()
  #adjointTestMultiSourceParallel()
  #goAmplitudeInversion()
  goAmplitudeInversionQs()
  #goInversionCg()
  #goInversionQs()
  #goInversionOfMultiplesQs()

def getModelAndMask():
  return getLayeredModelAndMask()
  #return getMarmousiModelAndMask()

def getMarmousiModelAndMask():
  s = getMarmousi(stride)
  s0,s1 = makeBornModel(s)
  m = fillfloat(1.0,nx,nz)
  ref = RecursiveExponentialFilter(1.0)
  for iz in range(13):
    for ix in range(nx):
      m[iz][ix] = 0.0
  ref.apply2(m,m)
  return s0,s1,m

def getLayeredModelAndMask():
  s = getLayeredModel()
  s0,s1 = makeBornModel(s)
  s0 = fillfloat(0.25,nx,nz)
  m = None
  return s0,s1,m

def goAmplitudeInversion():
  nouter = 1 # outer iterations
  ninner = 3 # inner iterations
  
  # Slowness and reflectivity models.
  s = getMarmousi(stride); s0,s1 = makeBornModel(s)
  #s = getLayeredModel(); s0,s1 = makeBornModel(s); s0 = fillfloat(0.25,nx,nz)

  # Allocate wavefields.
  print "allocating..."
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)

  # Born operator.
  born = BornOperatorS(s0,dx,dt,nabsorb,u,a)

  # Sources and receivers.
  src = BornOperatorS.getSourceArray(ns)
  rcp = BornOperatorS.getReceiverArray(ns) # predicted
  rco = BornOperatorS.getReceiverArray(ns) # observed
  for isou in range(ns):
    src[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
    rcp[isou] = Receiver(xr,zr,nt)
    rco[isou] = Receiver(xr,zr,nt)

  # Observed data.
  born.applyForward(src,s1,rco)

  # Slowness error.
  mul(0.95,s0,s0) 

  # CG solver.
  ref = RecursiveExponentialFilter(0.5/fpeak/dx)
  mm = PreconditionOperator(born,src,ref)
  ma = HessianOperator(born,src,rcp)
  cg = CgSolver(0.0,ninner);

  # Data warping.
  td = 4 # time decimation
  maxShift = 0.1 # max shift (seconds)
  strainT,strainR,strainS = 0.50,0.20,0.20 # 3d warping
  #strainT,strainR,strainS = 0.50,0.20,-1.0 # 2d warping
  smoothT,smoothR,smoothS = 16.0,4.0,4.0
  warp = DataWarping(
    strainT,strainR,strainS,smoothT,smoothR,smoothS,maxShift,dt,td)

  for iout in range(nouter):

    # RHS vector.
    print "computing RHS..."
    vb = VecArrayFloat2(zerofloat(nx,nz))
    born.applyAdjoint(src,rco,vb.getArray())

    # Solution vector.
    vx = VecArrayFloat2(zerofloat(nx,nz))

    # Solve.
    print "solving..."
    cg.solve(ma,mm,vb,vx);

    pixels(rco[ns/2].getData())
    pixels(rcp[ns/2].getData())
    
    # Warp.
    if nouter>1:
      print "warping..."
      rco = warp.warp(rcp,rco)

    pixels(rco[ns/2].getData())

  pixels(s0,cmap=jet)
  pixels(s1,sperc=99.5)
  pixels(vx.getArray(),sperc=99.5)


def goInversionCg():
  #s = getMarmousi(stride); s0,s1 = makeBornModel(s)
  s = getLayeredModel(); s0,s1 = makeBornModel(s); s0 = fillfloat(0.25,nx,nz)

  # Allocate wavefields.
  print "allocating"
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)

  # Sources and receivers.
  source = BornOperatorS.getSourceArray(ns)
  receiver = BornOperatorS.getReceiverArray(ns)
  for isou in range(ns):
    source[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
    receiver[isou] = Receiver(xr,zr,nt)

  # Born operator.
  born = BornOperatorS(s0,dx,dt,nabsorb,u,a)

  # RHS vector.
  print "computing RHS"
  rb = like(s1)
  vb = VecArrayFloat2(rb)
  born.applyHessian(source,receiver,s1,rb)

  # Solution vector.
  rx = zerofloat(nx,nz)
  vx = VecArrayFloat2(rx)

  # CG solver.
  niter = 1
  #if niter>1:
  if True:
    print "starting CgSolver"
    ref = RecursiveExponentialFilter(0.5/fpeak/dx)
    mm = PreconditionOperator(born,source,ref)
    ma = HessianOperator(born,source,receiver)
    cg = CgSolver(0.0,niter);
    cg.solve(ma,mm,vb,vx);
  else:
    copy(rb,rx)
    ref = RecursiveExponentialFilter(0.5/fpeak/dx)
    mm = PreconditionOperator(born,source,ref)
    mm.apply(vx,vx)

  pixels(s0,cmap=jet)
  pixels(s1,sperc=99.5)
  pixels(rx,sperc=99.5)

class HessianOperator(CgSolver.A):
  def __init__(self,born,source,receiver):
    self.born = born
    self.source = source
    self.receiver = receiver
  def apply(self,vx,vy):
    print "applying LHS"
    x = vx.getArray()
    y = vy.getArray()
    self.born.applyHessian(self.source,self.receiver,x,y)

class PreconditionOperator(CgSolver.A):
  def __init__(self,born,source,ref):
    self.born = born
    self.source = source
    self.ref = ref
    self.illum = None
  def apply(self,vx,vy):
    print "applying preconditioner"
    x = vx.getArray()
    y = vy.getArray()
    self.born.applyAdjointRoughen(self.ref,x,y)
    self.compensateIllum(y,y)
    self.born.applyForwardRoughen(self.ref,y,y)
  def compensateIllum(self,x,y):
    if self.illum is None:
      print "computing preconditioner"
      nx,nz = len(y[0]),len(y)
      self.illum = m = zerofloat(nx,nz)
      self.born.applyForIllumination(self.source,m)
      div(1.0,m,m) # 1/illumination
      mul(m,m,m) # squared
      mul(1.0/max(abs(m)),m,m) # normalized
      pixels(m,cmap=jet)
    mul(self.illum,x,y)

def goAmplitudeInversionQs():
  nouter,ninner,nfinal = 4,2,2 # outer, inner, inner for last outer
  #nouter,ninner,nfinal = 5,2,10 # outer, inner, inner for last outer
  #nouter,ninner,nfinal = 0,0,10 # outer, inner, inner for last outer
  warp3d = False # use 3D warping
  s,r,m = getModelAndMask()
  #e = mul(0.95,s) # erroneous background slowness
  e = mul(0.85,s) # erroneous background slowness

  # Wavefields
  print "allocating"
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)

  # BornOperator
  born = BornOperatorS(e,dx,dt,nabsorb,u,a) # erroneous slowness
  bornt = BornOperatorS(s,dx,dt,nabsorb,u,a) # true slowness

  # Sources and receivers
  src = BornOperatorS.getSourceArray(ns) # sources
  rco = BornOperatorS.getReceiverArray(ns) # receivers
  rcp = BornOperatorS.getReceiverArray(ns) # receivers
  for isou in range(ns):
    src[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
    rco[isou] = Receiver(xr,zr,nt)
  bornt.applyForward(src,r,rco) # observed data
  for isou in range(ns):
    rcp[isou] = Receiver(rco[isou])

  # BornSolver
  ref = RecursiveExponentialFilter(0.5/(fpeak*dx))
  bs = BornSolver(born,src,rcp,rco,ref,m)

  # DataWarping.
  td = 4 # time decimation
  maxShift = 0.1 # max shift (seconds)
  strainT,strainR,strainS = 0.50,0.20,0.20 if warp3d else -1.0
  smoothT,smoothR,smoothS = 16.0,4.0,4.0
  warp = DataWarping(
    strainT,strainR,strainS,smoothT,smoothR,smoothS,maxShift,dt,td)

  w = zerofloat(nt,nr,ns) # warping shifts
  for iouter in range(nouter+1):
    rx = bs.solve(nfinal if iouter==nouter else ninner);
    born.applyForward(src,rx,rcp)
    pixels(rx,sperc=99.5,title='r%d'%iouter)
    pixels(rcp[ns/2].getData(),title='rcp%d'%iouter)
    if nouter>1 and iouter<nouter:
      print 'warping (3d=%r)'%warp3d
      rcw = warp.warp(rcp,rco,w)
      bs.setObservedData(rcw)
      pixels(rcw[ns/2].getData(),title='rcw%d'%iouter)
      pixels(w[ns/2],cmap=rwb,sperc=100.0,title='shifts%d'%iouter)

  pixels(rco[ns/2].getData(),title='rco')
  pixels(s,cmap=jet,title='s')
  pixels(e,cmap=jet,title='e')
  pixels(r,sperc=99.5,title='r')

def goInversionQs():
  niter = 2
  #s,r = makeBornModel(getMarmousi(stride))
  s,r,m = getModelAndMask()
  e = mul(0.95,s) # erroneous background slowness

  src = BornOperatorS.getSourceArray(ns) # sources
  rco = BornOperatorS.getReceiverArray(ns) # receivers
  for isou in range(ns):
    src[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
    rco[isou] = Receiver(xr,zr,nt)
  print "allocating"
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)
  born = BornOperatorS(e,dx,dt,nabsorb,u,a)
  bornt = BornOperatorS(s,dx,dt,nabsorb,u,a) # true slowness

  ref = RecursiveExponentialFilter(0.5/(fpeak*dx));
  bs = BornSolver(born,src,rco,ref)
  bs.setTrueBornOperator(bornt)
  bs.setTrueReflectivity(r)
  """
  ref = RecursiveExponentialFilter(0.5/(fpeak*dx));
  bornt.applyForward(src,r,rco)
  bs = BornSolver(born,src,rco,ref)
  """

  rx = bs.solve(niter);
  pixels(s,cmap=jet,title="s")
  pixels(e,cmap=jet,title="e")
  pixels(r,sperc=99.5,title="r")
  pixels(rx,sperc=99.5,title="r%d"%niter)

def xgoInversionQs():
  niter = 2
  #s = getMarmousi(stride); s0,s1 = makeBornModel(s)
  s = getLayeredModel(); s0,s1 = makeBornModel(s); s0 = fillfloat(0.25,nx,nz)
  ################
  print "allocating"
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)
  born = BornOperatorS(mul(0.95,s0),dx,dt,nabsorb,u,a) # XXX
  #born = BornOperatorS(mul(0.85,s0),dx,dt,nabsorb,u,a) # XXX
  bornt = BornOperatorS(s0,dx,dt,nabsorb,u,a) # true slowness
  src = BornOperatorS.getSourceArray(ns)
  rcp = BornOperatorS.getReceiverArray(ns) # predicted
  rco = BornOperatorS.getReceiverArray(ns) # predicted
  for isou in range(ns):
    src[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
    rcp[isou] = Receiver(xr,zr,nt)
    rco[isou] = Receiver(xr,zr,nt)
  #qt = QuadraticTransform(0.5/fpeak/dx,born,src,rcp,refl=s1)
  qt = QuadraticTransformX(0.5/fpeak/dx,born,bornt,src,rcp,rco,s1) # XXX
  rx = QuadraticSolver(qt).solve(niter,None).getData()
  pixels(s0,cmap=jet)
  pixels(s1,sperc=99.5)
  pixels(rx,sperc=99.5)

def goInversionOfMultiplesQs():
  niter = 5
  s = getShallowLayeredModel();
  s0,_ = makeBornModel(s)
  print "allocating"
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)
  born = BornOperatorS(s0,dx,dt,nabsorb,u,a)
  src = BornOperatorS.getSourceArray(ns)
  rcp = BornOperatorS.getReceiverArray(ns) # predicted
  rco = BornOperatorS.getReceiverArray(ns) # observed
  for isou in range(ns):
    src[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
    rcp[isou] = Receiver(xr,zr,nt)
    rco[isou] = Receiver(xr,zr,nt)

  wave = WaveOperatorS(s,dx,dt,nabsorb,u,a)
  wave.applyForward(src,rco)
  wave = WaveOperatorS(fillfloat(s[0][0],nx,nz),dx,dt,nabsorb,u,a)
  wave.applyForward(src,rcp)
  for isou in range(ns):
    sub(rco[isou].getData(),rcp[isou].getData(),rco[isou].getData())
  pixels(rco[ns/2].getData(),perc=99.8)

  qt = QuadraticTransform(0.5/fpeak/dx,born,src,rcp,rco)
  rx = QuadraticSolver(qt).solve(niter,None).getData()
  pixels(s0,cmap=jet)
  pixels(rx,sperc=99.5)

class QuadraticTransformX(Quadratic):
  def __init__(self,sigma,born,bornt,src,rcp,rco,refl,mask=None):
    self.ref = RecursiveExponentialFilter(sigma)
    nxz = born.getDimensions()
    self.nx = nxz[0]
    self.nz = nxz[1]
    self.born = born # BornOperatorS with incorrect background
    self.bornt = bornt # BornOperatorS with true background
    self.src = src # source
    self.rcp = rcp # predicted data
    self.rco = rco # observed data
    self.rr = refl # reflectivity
    self.mask = mask # optional mask
    self.illum = None # inverse illumination
  def getB(self):
    #print "! computing RHS"
    rb = zerofloat(self.nx,self.nz)
    vb = ArrayVect2f(rb,1.0)

    #self.bornt.applyHessian(self.src,self.rcp,self.rr,rb)
    self.bornt.applyForward(self.src,self.rr,self.rco)
    self.born.applyAdjoint(self.src,self.rco,rb)

    mul(-1.0,rb,rb) # gradient = adjoint applied to negative observed data?
    return vb
  def multiplyHessian(self,vx):
    #print "! applying LHS"
    rx = vx.getData()
    self.born.applyHessian(self.src,self.rcp,rx,rx)
  def inverseHessian(self,vx):
    #print "applying preconditioner"
    rx,ry = vx.getData(),zerofloat(nx,nz)
    self.born.applyAdjointRoughen(self.ref,rx,ry)
    self.applyIllum(ry,ry)
    self.applyMask(ry,ry)
    self.born.applyForwardRoughen(self.ref,ry,rx)
  def applyIllum(self,rx,ry):
    if self.illum is None:
      #print "! computing illumination"
      self.illum = m = zerofloat(self.nx,self.nz)
      self.born.applyForIllumination(self.src,m)
      mul(m,m,m) # square
      pixels(m,cmap=jet)
      add(1.0e-8*max(m),m,m) # stabilize
      div(1.0,m,m) # reciprocal
      mul(1.0/max(m),m,m) # normalize
      pixels(m,cmap=jet)
    mul(self.illum,rx,ry)
  def applyMask(self,rx,ry):
    if self.mask is not None:
      mul(self.mask,rx,ry)

class QuadraticTransform(Quadratic):
  def __init__(self,sigma,born,src,rcp,rco=None,refl=None,mask=None):
    self.ref = RecursiveExponentialFilter(sigma)
    nxz = born.getDimensions()
    self.nx = nxz[0]
    self.nz = nxz[1]
    self.born = born # BornOperatorS
    self.src = src # source
    self.rcp = rcp # predicted data
    self.rco = rco # observed data
    self.rr = refl # reflectivity
    self.mask = mask # optional mask
    self.illum = None # inverse illumination
  def getB(self):
    #print "! computing RHS"
    rb = zerofloat(self.nx,self.nz)
    vb = ArrayVect2f(rb,1.0)
    if self.rr is not None:
      self.born.applyHessian(self.src,self.rcp,self.rr,rb)
    else:
      self.born.applyAdjoint(self.src,self.rco,rb)
    mul(-1.0,rb,rb) # gradient = adjoint applied to negative observed data?
    return vb
  def multiplyHessian(self,vx):
    #print "! applying LHS"
    rx = vx.getData()
    self.born.applyHessian(self.src,self.rcp,rx,rx)
  def inverseHessian(self,vx):
    #print "applying preconditioner"
    rx,ry = vx.getData(),zerofloat(nx,nz)
    self.born.applyAdjointRoughen(self.ref,rx,ry)
    self.applyIllum(ry,ry)
    self.applyMask(ry,ry)
    self.born.applyForwardRoughen(self.ref,ry,rx)
  def applyIllum(self,rx,ry):
    if self.illum is None:
      #print "! computing illumination"
      self.illum = m = zerofloat(self.nx,self.nz)
      self.born.applyForIllumination(self.src,m)
      mul(m,m,m) # square
      pixels(m,cmap=jet)
      add(1.0e-8*max(m),m,m) # stabilize
      div(1.0,m,m) # reciprocal
      mul(1.0/max(m),m,m) # normalize
      pixels(m,cmap=jet)
    mul(self.illum,rx,ry)
  def applyMask(self,rx,ry):
    if self.mask is not None:
      mul(self.mask,rx,ry)

def goAcousticData():
  #s = getLayeredModel()
  s = getMarmousi(stride)
  s,_ = makeBornModel(s)
  awo = WaveOperator(s,dx,dt,nabsorb)
  source = Source.RickerSource(xs[0],zs[0],dt,fpeak)
  receiver = Receiver(xr,zr,nt)
  u = zerofloat(nxp,nzp,nt)
  sw = Stopwatch(); sw.start()
  awo.applyForward(source,receiver,u)
  print sw.time()
  d = receiver.getData()
  for it in range(1300): # mute direct arrival
    for ir in range(nr):
      d[ir][it] = 0.0
  points(d[nr/2])
  pixels(d,cmap=gray,sperc=99.9,title="data")
  #pixels(d,cmap=gray,cmin=-0.35,cmax=0.35,title="data")
  #pixels(s,cmap=jet,title="slowness (s/km)")

def goBornData():
  #s = getLayeredModel()
  s = getMarmousi(stride)
  s0,s1 = makeBornModel(s)
  #zero(s1)
  #for ix in range(nx):
  #  s1[nz/3][ix] = 1.00
  bwo = BornWaveOperator(
    s,dx,dt,nabsorb)
  receiver = Receiver(xr,zr,nt)
  u = zerofloat(nxp,nzp,nt)
  source = Source.RickerSource(xs[0],zs[0],dt,fpeak)
  #source = Source.Gaussian4Source(xs[0],zs[0],dt,fpeak)
  #source = makeSource()
  sw = Stopwatch(); sw.start()
  bwo.applyForward(source,u,s1,receiver)
  print sw.time()
  d = receiver.getData()
  #pixels(d,cmap=gray,cmin=-0.35,cmax=0.35,title="data")
  points(d[nr/2])
  pixels(d,cmap=gray,sperc=99.9,title="data")
  pixels(s0,cmap=jet,title="background slowness (s/km)")
  pixels(s1,sperc=100.0,title="reflectivity")

# Phase rotated Ricker
def makeSource():
  def ricker(t):
    x = FLT_PI*fpeak*(t-1.5/fpeak)
    xx = x*x
    return (1.0-2.0*xx)*exp(-xx);
  w = zerofloat(nt)
  for it in range(nt):
    t = ft+it*dt
    w[it] = ricker(t)
  #points(w)
  w = rotateAndDifferentiate(w)
  points(w)
  return Source.WaveletSource(xs[0],zs[0],w)
def rotateAndDifferentiate(rx):
  fft = Fft(rx)
  sf = fft.getFrequencySampling1()
  nf = sf.count
  cy = fft.applyForward(rx) # forward FFT
  p = 0.25*FLT_PI # phase rotation angle
  t = zerofloat(2*nf)
  for i in range(nf):
    w = sf.getValue(i)
    t[2*i  ] = cos(p) # phase rotation
    t[2*i+1] = sin(p) # phase rotation
  cmul(t,cy,cy)
  #t = zerofloat(2*nf)
  #for i in range(nf):
  #  w = sf.getValue(i)
  #  t[2*i  ] = w*w # negative 2nd time derivative
  #  #t[2*i+1] = w   # 1st time derivative
  #cmul(t,cy,cy)
  ry = fft.applyInverse(cy) # inverse FFT
  return ry

def makeBornModel(s):
  """
  sigma0: smoothing for background model
  sigma1: smoothing for perturbation
  """
  #sigma0 = 0.5*nx*nz/sum(s)/fpeak # half wavelength
  #sigma1 = 0.5*sigma0 # quarter wavelength
  sigma0 = 1.0/fpeak
  sigma1 = 0.5*sigma0
  print 'sigma0=%f'%sigma0
  print 'sigma1=%f'%sigma1
  s0,s1 = like(s),like(s)
  ref0 = RecursiveExponentialFilter(sigma0/dx)
  ref0.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref0.apply(s,s0)
  t = copy(s)
  ref1 = RecursiveExponentialFilter(sigma1/dx)
  ref1.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref1.apply(s,t)
  sub(s,t,s1)
  r = sub(mul(s,s),mul(t,t))
  GaussianTaper.apply1(0.25,r,r)
  return s0,r

def adjointTest():
  s = fillfloat(0.25,nx,nz)
  add(mul(sub(randfloat(Random(0),nx,nz),0.5),0.05),s,s)
  random = Random(01234)
  #random = Random()
  b = randfloat(random,nxp,nzp,nt); sub(b,0.5,b) # background wavefield
  bwo = BornWaveOperator(s,dx,dt,nabsorb)
  bwo.setIlluminationCompensation(True)
  #bwo.setReflectivityRoughening(1.0)
  ra = randfloat(random,nx,nz) # random reflectivity
  rb = zerofloat(nx,nz) 
  db = randfloat(random,nt,nr) # random data
  reca = Receiver(xr,zr,nt)
  recb = Receiver(xr,zr,db)
  bwo.applyForward(b,ra,reca)
  bwo.applyAdjoint(b,zerofloat(nxp,nzp,nt),recb,rb)
  sum1 = dot(ra,rb)
  sum2 = dot(reca.getData(),recb.getData())
  print "adjoint test:",WaveOperator.compareDigits(sum1,sum2)
  print sum1
  print sum2

def adjointTestMultiSource():
  nsou = 1
  #random = Random()
  random = Random(01234)
  s = fillfloat(0.25,nx,nz)
  add(mul(sub(randfloat(Random(0),nx,nz),0.5),0.05),s,s)
  b = WaveOperator.randfloat(random,nxp,nzp,nt,nsou)
  born = BornWaveOperator(s,dx,dt,nabsorb)

  # Forward
  ra = randfloat(random,nx,nz) # random reflectivity
  da = zerofloat(nt,nr,nsou)
  for isou in range(nsou):
    bi = b[isou]
    rec = Receiver(xr,zr,da[isou])
    born.applyForward(bi,ra,rec)

  # Adjoint
  rb = zerofloat(nx,nz)
  db = randfloat(random,nt,nr,nsou) # random data
  for isou in range(nsou):
    bi = b[isou]
    rec = Receiver(xr,zr,nt)
    copy(db[isou],rec.getData())
    rt = zerofloat(nx,nz)
    born.applyAdjoint(bi,zerofloat(nxp,nzp,nt),rec,rt)
    add(rt,rb,rb)

  # Dot product
  print 'sum(ra)=%f'%sum(ra)
  print 'sum(rb)=%f'%sum(rb)
  print 'sum(da)=%f'%sum(da)
  print 'sum(db)=%f'%sum(db)
  sum1 = dot(ra,rb)
  sum2 = dot(da,db)
  print "adjoint test:",WaveOperator.compareDigits(sum1,sum2)
  print sum1
  print sum2

def adjointTestMultiSourceParallel():
  nsou = 4 # number of sources
  psou = 2 # number of sources to compute in parallel
  #random = Random()
  random = Random(01234)
  s = fillfloat(0.25,nx,nz)
  add(mul(sub(randfloat(Random(0),nx,nz),0.5),0.05),s,s)
  u = SharedFloat4(WaveOperator.randfloat(random,nxp,nzp,nt,psou))
  a = SharedFloat4(nxp,nzp,nt,psou)
  born = BornOperatorS(s,dx,dt,nabsorb,u,a)
  #born.setReflectivityRoughening(1.0)

  # Forward
  ra = randfloat(random,nx,nz) # random reflectivity
  da = zerofloat(nt,nr,nsou)
  reca = BornOperatorS.getReceiverArray(nsou)
  for isou in range(nsou):
    reca[isou] = Receiver(xr,zr,da[isou])
  born.applyForward(ra,reca)

  # Adjoint
  rb = zerofloat(nx,nz)
  db = randfloat(random,nt,nr,nsou) # random data
  recb = BornOperatorS.getReceiverArray(nsou)
  for isou in range(nsou):
    recb[isou] = Receiver(xr,zr,db[isou])
  born.applyAdjoint(recb,rb)

  # Dot product
  sum1 = dot(ra,rb)
  sum2 = dot(da,db)
  print "adjoint test:",WaveOperator.compareDigits(sum1,sum2)
  print sum1
  print sum2

  # Preconditioner test
  ra = randfloat(random,nx,nz)
  rb = randfloat(random,nx,nz)
  ca = zerofloat(nx,nz)
  cb = zerofloat(nx,nz)
  ref = RecursiveExponentialFilter(1.0)
  born.applyForwardRoughen(ref,ra,ca)
  born.applyAdjointRoughen(ref,rb,cb)
  sum1 = dot(ra,cb)
  sum2 = dot(rb,ca)
  print "preconditioner adjoint test:",\
    WaveOperator.compareDigits(sum1,sum2)
  print sum1
  print sum2

def dot(u,a):
  return WaveOperator.dot(u,a)

def getLayeredModel():
  """Make slowness (s/km) model."""
  s = fillfloat(1.0/1.5,nx,nz) # water velocity
  #s = fillfloat(0.5,nx,nz)
  for iz in range(nz/3,2*nz/3):
    for ix in range(nx):
      s[iz][ix] = 0.50
  for iz in range(2*nz/3,nz):
    for ix in range(nx):
      s[iz][ix] = 0.2
      #s[iz][ix] = 0.5
  return s

def getShallowLayeredModel():
  s = fillfloat(0.5,nx,nz)
  for iz in range(nz/5,2*nz/5):
    for ix in range(nx):
      s[iz][ix] = 0.2
  for iz in range(2*nz/5,nz):
    for ix in range(nx):
      s[iz][ix] = 0.5
  return s

def like(x):
  return zerofloat(len(x[0]),len(x))

import socket
def getMarmousi(stride=3):
  """Marmousi model.
  Parameters:
    sigma - half-width of smoothing window for scale separation
    econst - constant error in background slowness s0
    egauss - gaussian error in background slowness s0
    erand - random error in background slowness s0
  Returns:
    tt - true model
    t0 - true background slowness
    t1 - true reflectivity
    s0 - true background slowness
    s1 - true reflectivity
  """
  p = zerofloat(751,2301)
  if socket.gethostname()=='backus.Mines.EDU':
    read("/data/sluo/marmousi/marmousi.dat",p)
  else:
    read("/data/seis/marmousi/marmousi.dat",p)
  p = copy(743,2301,8,0,p)
  div(1000.0,p,p) # slowness (s/km)
  wa = iceil(0.2/(0.004*stride))
  nz = iceil(743.0/stride)+wa
  nx = iceil(2301.0/stride)
  #print "wa =",wa
  #print "nz =",nz
  #print "nx =",nx
  t = fillfloat(2.0/3.0,nz,nx)
  copy(nz-wa,nx,0,0,stride,stride,p,wa,0,1,1,t)
  return transpose(t)

def iceil(x):
  return int(ceil(x))

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
# plotting

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def pixels(x,cmap=gray,perc=100.0,sperc=None,cmin=0.0,cmax=0.0,title=None):
  if len(x)==nz:
    x = transpose(x)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  cb = sp.addColorBar()
  cb.setWidthMinimum(100)
  sp.setSize(1010,740)
  if title:
    sp.addTitle(title)
  pv = sp.addPixels(x)
  pv.setColorModel(cmap)
  if perc<100.0:
    pv.setPercentiles(100.0-perc,perc)
  if sperc is not None: # symmetric percentile clip (for plotting gradients)
    clips = Clips(100-sperc,sperc,x)
    clip = max(abs(clips.getClipMin()),abs(clips.getClipMax()))
    pv.setClips(-clip,clip)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if title and pngdatDir:
    sp.paintToPng(360,3.33,pngdatDir+title+'.png')
    write(pngdatDir+title+'.dat',x)

def points(x):
  SimplePlot.asPoints(x)

#############################################################################
# Do everything on Swing thread.
import sys,time
class RunMain(Runnable):
  def run(self):
    start = time.time()
    if pngdatDir is not None:
      print 'cleaning '+pngdatDir.split('/')[-2]
      cleanDir(pngdatDir)
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02d'%(h,m,s)
if __name__=='__main__':
  SwingUtilities.invokeLater(RunMain())
