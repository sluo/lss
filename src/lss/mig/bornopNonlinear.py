#############################################################################
# Wavefield modeling using BornWaveOperator

from imports import *

#############################################################################

savDir = None
savDir = os.getenv('HOME')+'/Desktop/pngdat/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat2/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat3/'

#############################################################################

def main(args):
  #readFiles()
  #getModelAndMask()
  #goNonlinearInversionQs() # inversion using line search
  goNonlinearAmplitudeInversionQs() # amplitude inversion using line search
  #goInversionQs() # inversion without line search
  #goAmplitudeInversionQs() # amplitude inversion without line search
  #goNewAmplitudeInversionQs() # shift predicted data instead

  #compareData()
  #test()

def compareData():
  born,bs,src,rcp,rco,warp = getInputs()
  rf = zerofloat(nx,nz)
  #read('/home/sluo/Desktop/new/rand/0-0-10/r0.dat',rf)
  read('/home/sluo/Desktop/pngdat3/r0.dat',rf)
  _,_,_,m,_ = getMarmousiModelAndMask()
  mul(m,rf,rf)
  print "computing predicted data..."
  born.applyForward(src,rf,rcp)
  pixels(rcp[ns/2].getData(),sperc=99.9,title='rcp')
  pixels(rco[ns/2].getData(),sperc=99.9,title='rco')

def test():
  """ Amplitude inversion using shifts computed from lsm.py result."""
  niter = 10
  born,bs,src,rcp,rco,warp = getInputs()

  # Time shifts
  rf = zerofloat(nz,nx)
  read('/home/sluo/Desktop/save/lsm/marmousi/95p/ares5/s1_19.dat',rf)
  rf = transpose(rf)
  print "computing predicted data..."
  born.applyForward(src,rf,rcp)
  w = zerofloat(nt,nr,ns) # warping shifts
  print 'warping...'
  rcw = warp.warp(rco,rcp,w)
  bs.setTimeShifts(w)

  # Solve
  rx = bs.solve(niter)
  pixels(rx,cmap=gray,sperc=100.0,title='r0')

  dmin,dmax = min(rco[ns-1].getData()),max(rco[ns-1].getData())
  pixels(w[ns-1],cmap=rwb,sperc=100.0,title='w0')
  pixels(rcp[ns-1].getData(),cmin=dmin,cmax=dmax,title='rcp0')
  pixels(rcw[ns-1].getData(),cmin=dmin,cmax=dmax,title='rcw0')
  pixels(rco[ns-1].getData(),cmin=dmin,cmax=dmax,title='rco')

def readFiles():
  setupForMarmousi()
  #ra = zerofloat(nz,nx)
  ra = zerofloat(nx,nz)
  rb = zerofloat(nx,nz)
  #read('/home/sluo/Desktop/save/lsm/marmousi/95p/ares5/s1_19.dat',ra)
  #read('/home/sluo/Desktop/pngdat3/r0.dat',rb); rb = transpose(rb)
  read('/home/sluo/Desktop/new/rand/0-0-10/r0.dat',ra)
  read('/home/sluo/Desktop/new/rand/0-0-10-A-256/r0.dat',rb)
  pixels(ra,sperc=99.9)
  pixels(rb,sperc=99.9)

def setupForMarmousi():
  global sz,sx,st,nz,nx,nt,nxp,nzp,dz,dx,dt
  global zs,xs,zr,xr,ns,nr,fpeak,nabsorb,np
  global sigma0,sigma1
  sz = Sampling(265,0.012,0.0)
  sx = Sampling(767,0.012,0.0)
  st = Sampling(5501,0.0010,0.0)
  nz,nx,nt = sz.count,sx.count,st.count
  dz,dx,dt = sz.delta,sx.delta,st.delta
  fz,fx,ft = sz.first,sx.first,st.first
  #xs,zs = [0],[0]
  #xs,zs = [nx/2],[0]
  #xs,zs = rampint(1,15,52),fillint(0,52)
  xs,zs = rampint(3,10,77),fillint(0,77)
  #xs,zs = rampint(3,8,96),fillint(0,96)
  #xs,zs = rampint(3,5,153),fillint(0,153)
  #xs,zs = rampint(1,3,256),fillint(0,256)
  xr,zr = rampint(0,1,nx),fillint(0,nx)
  ns,nr = len(xs),len(xr)
  fpeak = 10.0 # Ricker wavelet peak frequency
  nabsorb = 22 # absorbing boundary size
  nxp,nzp = nx+2*nabsorb,nz+2*nabsorb
  np = min(8,ns) # number of parallel sources
  sigma0 = 1.0/fpeak
  sigma1 = 1.0/fpeak
  print 'ns=%d'%ns
  return getMarmousiModelAndMask()

def setupForLayered():
  global sz,sx,st,nz,nx,nt,nxp,nzp,dz,dx,dt
  global zs,xs,zr,xr,ns,nr,fpeak,nabsorb,np
  global sigma0,sigma1
  #sz = Sampling(301,0.012,0.0) # for layered model
  #sx = Sampling(501,0.012,0.0)
  #st = Sampling(2750,0.0015,0.0)
  sz = Sampling(201,0.016,0.0)
  sx = Sampling(202,0.016,0.0)
  st = Sampling(2003,0.0010,0.0)
  nz,nx,nt = sz.count,sx.count,st.count
  dz,dx,dt = sz.delta,sx.delta,st.delta
  fz,fx,ft = sz.first,sx.first,st.first
  #xs,zs = [0],[0]
  #xs,zs = [nx/2],[0]
  #xs,zs = [nx/4,nx/2,3*nx/4],[0,0,0]
  xs,zs = [nx/5,2*nx/5,3*nx/5,4*nx/5],[0,0,0,0]
  #xs,zs = rampint(1,10,51),fillint(0,51) # for layered
  xr,zr = rampint(0,1,nx),fillint(0,nx)
  ns,nr = len(xs),len(xr)
  fpeak = 20.0 # Ricker wavelet peak frequency
  nabsorb = 22 # absorbing boundary size
  nxp,nzp = nx+2*nabsorb,nz+2*nabsorb
  np = min(16,ns) # number of parallel sources
  sigma0 = 4.0/fpeak
  sigma1 = 1.0/fpeak
  print 'ns=%d'%ns
  return getLayeredModelAndMask()

def getMarmousiModelAndMask():
  s = getMarmousi()
  s0,s1 = makeBornModel(s)

  # Water-layer mask.
  m = fillfloat(1.0,nx,nz)
  for iz in range(14):
    for ix in range(nx):
      m[iz][ix] = 0.0
  RecursiveExponentialFilter(1.0).apply2(m,m)

  # Preconditioning by v^2.
  p = div(1.0,s)
  mul(p,p,p)
  sigma = 1.0 # half-width in km
  niter = 8 # number of times to apply filter
  sigma = sigma/(sqrt(niter)*dx)
  ref = RecursiveExponentialFilter(sigma)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  for i in range(niter):
    ref.apply(p,p)
  mul(1.0/max(p),p,p)
  mul(m,p,p)

  #pixels(s,cmap=jet)
  #pixels(s0,cmap=jet)
  #pixels(s1,sperc=100.0)
  #pixels(p,cmap=gray)
  return s,s0,mul(s1,m),m,p

def getLayeredModelAndMask():
  constantBackground = True
  tb = 0.25 # background slowness
  t = getLayered2()
  t0,t1 = makeBornModel(t)
  GaussianTaper.apply1(t1,t1)
  if constantBackground:
    fill(tb,t0)
  m = fillfloat(1.0,nx,nz)
  for iz in range(nz/3-6):
    for ix in range(nx):
      m[iz][ix] = 0.0
  RecursiveExponentialFilter(1.0).apply2(m,m)
  #return t,t0,t1,m
  return t,t0,t1,None,None

def getLayered2(s0mul=1.0):
  t = fillfloat(1.0/1.5,nz,nx)
  for ix in range(nx):
    for iz in range(nz/3,2*nz/3):
      t[ix][iz] = 0.5
    for iz in range(2*nz/3,nz):
      t[ix][iz] = 0.2
  #t = fillfloat(0.35,nz,nx)
  #for ix in range(nx):
  #  for iz in range(nz/3,2*nz/3):
  #    t[ix][iz] = 0.3
  #  for iz in range(2*nz/3,nz):
  #    t[ix][iz] = 0.2
  return transpose(t)

#############################################################################
# Inversion without line search

def addRandomError(s,emax=0.010,m=None):
  random = Random(0)
  r = sub(randfloat(random,nz,nx),0.5); r = transpose(r)
  RecursiveGaussianFilter(62.5).apply00(r,r)
  mul(emax/max(abs(r)),r,r)
  if m is not None:
    mul(m,r,r)
  return add(s,r)

def addGaussianError(s,emax=0.010,m=None):
  g = zerofloat(nx,nz)
  g[nz/2][nx/2] = 1.0
  RecursiveGaussianFilter(nz/5.0).apply00(g,g)
  mul(emax/max(g),g,g)
  if m is not None:
    mul(m,g,g)
  return add(s,g)

def getModelAndMask():
  t,s,r,m,p = setupForMarmousi()
  #t,s,r,m,p = setupForLayered()
  e = mul(0.95,s) # erroneous background slowness
  #e = mul(0.85,s) # erroneous background slowness
  #e = addGaussianError(s,emax=0.025,m=m)
  #e = addRandomError(s,emax=0.025,m=m)
  #e = addRandomError(s,emax=0.050,m=m)
  pixels(t,cmap=jet,title='t')
  pixels(s,cmap=jet,title='s')
  pixels(e,cmap=jet,title='e')
  if m is not None:
    pixels(m,cmap=gray,title='m')
  pixels(sub(s,e),cmap=rwb,sperc=100.0,title='s-e')
  pixels(r,cmap=gray,sperc=100.0,title='r')
  return t,s,e,r,p

def getInputs():
  useAcoustic = False # use WaveOperator for observed data
  warp3d = False # use 3D warping
  print 'warp3d=%r'%warp3d
  t,s,e,r,m = getModelAndMask()

  # Wavefields.
  print "allocating..."
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)

  # Sources and receivers.
  print "computing observed data..."
  src = zeros(ns,Source)
  rco = zeros(ns,Receiver)
  rcp = zeros(ns,Receiver)
  for isou in range(ns):
    rco[isou] = Receiver(xr,zr,nt)
    rcp[isou] = Receiver(xr,zr,nt)
  if useAcoustic:
    v = fillfloat(t[0][0],nx,nz)
    sro = BornOperatorS.getSourceArray(ns) # source for observed data
    rct = BornOperatorS.getReceiverArray(ns) # temp receivers
    for isou in range(ns):
      sro[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
      src[isou] = Source.WaveletSource(xs[isou],zs[isou],rotatedRicker())
      rct[isou] = Receiver(xr,zr,nt)
    wavet = WaveOperatorS(t,dx,dt,nabsorb,u,a) # true slowness
    wavev = WaveOperatorS(v,dx,dt,nabsorb,u,a) # for direct arrival
    wavet.applyForward(sro,rco) # observed data
    wavev.applyForward(sro,rct) # direct arrival
    for isou in range(ns):
      sub(rco[isou].getData(),rct[isou].getData(),rco[isou].getData())
  else:
    for isou in range(ns):
      src[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
    bornt = BornOperatorS(s,dx,dt,nabsorb,u,a) # true slowness
    bornt.applyForward(src,r,rco) # observed data

  # Born modeling and solver.
  born = BornOperatorS(e,dx,dt,nabsorb,u,a)
  ref = RecursiveExponentialFilter(1.0/(fpeak*dx*sqrt(2.0)))
  bs = BornSolver(born,src,rcp,rco,ref,m) # solver

  # Warping.
  td = 4 # time decimation
  maxShift = 0.1 # max shift (seconds)
  strainT,strainR,strainS = 0.50,0.20,1.00 if warp3d else -1.0
  #smoothT,smoothR,smoothS = 16.0,4.0,4.0
  smoothT,smoothR,smoothS = 32.0,8.0,8.0
  warp = ImageWarping(
    strainT,strainR,strainS,smoothT,smoothR,smoothS,maxShift,dt,td)
  w = zerofloat(nt,nr,ns) # warping shifts

  return born,bs,src,rcp,rco,warp

def goAmplitudeInversionQs():
  nouter,ninner,nfinal = 5,2,5 # outer, inner, inner for last outer
  #nouter,ninner,nfinal = 0,0,5 # outer, inner, inner for last outer
  zeroReflectivity = True # zero reflectivity after each outer iteration
  """""" 
  born,bs,src,rcp,rco,warp = getInputs()
  rx = zerofloat(nx,nz) # reflectivity image
  w = zerofloat(nt,nr,ns) # warping shifts
  dmin,dmax = min(rco[ns-1].getData()),max(rco[ns-1].getData())
  for iouter in range(nouter+1):
    bs.solve(nfinal if iouter==nouter else ninner,rx);
    pixels(rx,cmap=gray,sperc=100.0,title='r%d'%iouter)
    if nouter>0 and iouter<nouter:
      born.applyForward(src,rx,rcp)
      print 'warping'
      rcw = warp.warp(rcp,rco,w)
      bs.setObservedData(rcw)
      pixels(rcp[ns-1].getData(),cmin=dmin,cmax=dmax,title='rcp%d'%iouter)
      pixels(rcw[ns-1].getData(),cmin=dmin,cmax=dmax,title='rcw%d'%iouter)
      pixels(w[ns-1],cmap=rwb,sperc=100.0,title='w%d'%iouter)
    if zeroReflectivity:
      zero(rx)
  pixels(rco[ns-1].getData(),cmin=dmin,cmax=dmax,title='rco')

def goNewAmplitudeInversionQs():
  #nouter,ninner,nfinal = 5,2,5 # outer, inner, inner for last outer
  nouter,ninner,nfinal = 0,0,10 # outer, inner, inner for last outer
  saveReflectivity = False # save reflectivity between outer iterations
  shiftsFromFile = True # compute shifts from reflectivity from file
  """""" 
  print 'nouter=%d'%nouter
  print 'ninner=%d'%ninner
  print 'nfinal=%d'%nfinal
  print 'saveReflectivity=%r'%saveReflectivity
  print 'shiftsFromFile=%r'%shiftsFromFile
  born,bs,src,rcp,rco,warp = getInputs()
  w = zerofloat(nt,nr,ns) # warping shifts
  if shiftsFromFile:
    rf = zerofloat(nx,nz)
    read('/home/sluo/Desktop/new/const/0-0-10-96s-2/r0.dat',rf)
    _,_,_,m,_ = getMarmousiModelAndMask()
    mul(m,rf,rf)
    print "computing predicted data..."
    born.applyForward(src,rf,rcp)
    print "warping..."
    warp.warp(rco,rcp,w)
    bs.setTimeShifts(w)
  rx = zerofloat(nx,nz) # reflectivity image
  dmin,dmax = min(rco[ns-1].getData()),max(rco[ns-1].getData())
  for iouter in range(nouter+1):
    bs.solve(nfinal if iouter==nouter else ninner,rx);
    pixels(rx,cmap=gray,sperc=100.0,title='r%d'%iouter)
    if nouter>0 and iouter<nouter:
      print "computing predicted data..."
      born.applyForward(src,rx,rcp) # simulate predicted data
      print 'warping...'
      rcw = warp.warp(rco,rcp,w) # estimate new time shifts
      bs.setTimeShifts(w) # set new time shifts
      pixels(rcp[ns/2].getData(),cmin=dmin,cmax=dmax,title='rcp%d'%iouter)
      pixels(rcw[ns/2].getData(),cmin=dmin,cmax=dmax,title='rcw%d'%iouter)
      pixels(w[ns/2],cmap=rwb,sperc=100.0,title='w%d'%iouter)
    if not saveReflectivity:
      zero(rx)
  pixels(rco[ns/2].getData(),cmin=dmin,cmax=dmax,title='rco')

def rotatedRicker():
  """Phase-rotated and twice-differentiated Ricker wavelet."""
  def ricker(t):
    x = FLT_PI*fpeak*(t-1.5/fpeak)
    xx = x*x
    return (1.0-2.0*xx)*exp(-xx);
  w = zerofloat(nt)
  for it in range(nt):
    t = it*dt
    w[it] = ricker(t)
  w = rotateAndDifferentiate(w)
  mul(1.0/max(abs(w)),w,w)
  #points(w)
  return w
def rotateAndDifferentiate(rx):
  fft = Fft(rx)
  sf = fft.getFrequencySampling1()
  nf = sf.count
  cy = fft.applyForward(rx) # forward FFT
  p = 0.25*FLT_PI # phase rotation angle
  t = zerofloat(2*nf)
  for i in range(nf):
    w = sf.getValue(i)
    #t[2*i  ] = cos(p) # phase rotation
    #t[2*i+1] = sin(p) # phase rotation
    #t[2*i  ] = w*w # negative 2nd time derivative
    #t[2*i+1] = 0.0 # negative 2nd time derivative
    t[2*i  ] = w*w*cos(p) # phase rotation and negative 2nd time derivative
    t[2*i+1] = w*w*sin(p) # phase rotation and negative 2nd time derivative
  cmul(t,cy,cy)
  ry = fft.applyInverse(cy) # inverse FFT
  return ry

def goInversionQs():
  niter = 2
  #s,r = makeBornModel(getMarmousi(stride))
  _,s,_,r,m = getModelAndMask()
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

#############################################################################
# Nonlinear inversion using line search

def goNonlinearInversionQs():
  niter = 1
  #s0scale = 0.95
  s0scale = 1.00
  amplitudeResidual = False

  # Slowness models
  _,t0,_,t1,m = getModelAndMask()
  s0 = mul(s0scale,t0) # background slowness for migration

  # Allocate Wavefields.
  print "allocating"
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)

  # Wave operators.
  born = BornOperatorS(t0,dx,dt,nabsorb,u,a) # Born with true slowness
  wave = WaveOperator(s0,dx,dt,nabsorb) # Wave with erroneous slowness
  wave.setAdjoint(False) # use correct adjoint?

  # Sources and receivers.
  src = BornOperatorS.getSourceArray(ns) # sources
  rco = BornOperatorS.getReceiverArray(ns) # receivers
  rcp = BornOperatorS.getReceiverArray(ns) # receivers
  for isou in range(ns):
    src[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
    rco[isou] = Receiver(xr,zr,nt)
    rcp[isou] = Receiver(xr,zr,nt)
  born.applyForward(src,t1,rco) # observed data

  # Preconditioning by v^2.
  w = div(1.0,mul(s0,s0)); mul(1.0/max(w)/ns,w,w)

  # Gradient computation.
  class GradientParallelReduce(Parallel.ReduceInt):
    def compute(self,isou):
      ui,ai = u.get(isou),a.get(isou)
      wave.applyForward(src[isou],ui);
      #wave.applyLaplacian(ui) # 2nd time derivative
      if amplitudeResidual:
        rcr = makeAmplitudeResidual(rcp[isou],rco[isou])
      else:
        rcr = makeWaveformResidual(rcp[isou],rco[isou])
      #twiceIntegrate(rcr); # XXX
      wave.applyAdjoint(Source.ReceiverSource(rcr),ai)
      gi = applyImagingCondition(ui,ai)
      #roughen(gi) # roughen
      #mul(w,gi,gi) # precondition
      return gi
    def combine(self,ga,gb):
      return add(ga,gb)

  # Iterate...
  s1 = like(t1)
  gm,pm = None,None
  for iiter in range(niter):
    sw = Stopwatch(); sw.start()
    g = PartialParallel(np).reduce(ns,GradientParallelReduce())
    p = conjugateDirection(g,gm,pm)
    report('gradient',sw)
    if niter>1:
      amf = AmplitudeMisfitFunction(s1,p,src[ns/2],rco[ns/2],born)
      updateModel(amf,s1) # line search
    gm,pm = g,p
    pixels(gm,cmap=rwb,sperc=100.0,title='g_'+str(iiter))
    pixels(pm,cmap=rwb,sperc=100.0,title='p_'+str(iiter))
    pixels(s1,cmap=rwb,sperc=100.0,title='s1_'+str(iiter))

  #pixels(rco[ns/2].getData(),title='rco')
  pixels(t1,cmap=rwb,sperc=100.0,title='t1')
  pixels(t0,cmap=jet,title='t0')
  pixels(s0,cmap=jet,title='s0')

def applyImagingCondition(u,a):
  return WaveOperator.collapse(u,a,nabsorb)

def goNonlinearAmplitudeInversionQs():
  niter = 5
  _,t0,_,t1,m = getModelAndMask()
  s0 = mul(0.95,t0) # erroneous background slowness
  #s0 = mul(1.00,t0)

  # Wavefields
  print "allocating"
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)

  # BornOperator
  born = BornOperatorS(s0,dx,dt,nabsorb,u,a) # erroneous slowness
  bornt = BornOperatorS(t0,dx,dt,nabsorb,u,a) # true slowness

  # Sources and receivers
  src = BornOperatorS.getSourceArray(ns) # sources
  rco = BornOperatorS.getReceiverArray(ns) # receivers
  rcp = BornOperatorS.getReceiverArray(ns) # receivers
  for isou in range(ns):
    src[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
    rco[isou] = Receiver(xr,zr,nt)
    rcp[isou] = Receiver(xr,zr,nt)
  bornt.applyForward(src,t1,rco) # observed data

  w = div(1.0,mul(s0,s0)); mul(1.0/max(w)/ns,w,w) # v^2 preconditioning
  s1 = like(t1) # computed reflectivity
  gm,pm = None,None # previous gradient, conjugate gradient
  for iiter in range(niter):
    if iiter>0:
      born.applyForward(src,s1,rcp)
      rcr = makeAmplitudeResiduals(rcp,rco)
    else:
      rcr = BornOperatorS.getReceiverArray(ns)
      for isou in range(ns):
        rcr[isou] = Receiver(xr,zr,mul(-1.0,rco[isou].getData()))
    g = like(s1)
    #twiceIntegrate(rcr); # XXX
    born.applyAdjoint(src,rcr,g) # gradient
    roughen(g) # roughen
    mul(w,g,g) # precondition
    if m is not None:
      mul(m,g,g) # mask
    p = conjugateDirection(g,gm,pm) # conjugate gradient
    if niter>1:
      amf = AmplitudeMisfitFunction(s1,p,src[ns/2],rco[ns/2],born)
      updateModel(amf,s1) # line search
    pixels(g,cmap=rwb,sperc=100.0,title='g_'+str(iiter))
    pixels(p,cmap=rwb,sperc=100.0,title='p_'+str(iiter))
    pixels(s1,cmap=rwb,sperc=100.0,title='s1_'+str(iiter))
    gm,pm = g,p

  #pixels(rco[ns/2].getData(),title='rco')
  pixels(t1,cmap=rwb,sperc=100.0,title='t1')
  pixels(t0,cmap=jet,title='t0')
  pixels(s0,cmap=jet,title='s0')

def twiceIntegrate(rec):
  try:
    ns = len(rec)
    for isou in range(ns):
      d = rec[isou].getData()
      Util.integrate1(d,d)
      Util.integrate1(d,d)
      mul(-1.0,d,d)
  except:
    d = rec.getData()
    Util.integrate1(d,d)
    Util.integrate1(d,d)
    mul(-1.0,d,d)

#############################################################################
# Line search

def updateModel(misfitFunction,s1):
  print 'searching for step length...'
  _,_,_,t1,_ = getModelAndMask()
  #a,b = -1.0*max(abs(t1)),0.1*max(abs(t1)); tol = 0.10*abs(b-a)
  #a,b = -1.0*max(abs(t1)),0.1*max(abs(t1)); tol = 0.20*abs(b-a)
  #a,b = -1.0*max(abs(t1)),0.0*max(abs(t1)); tol = 0.20*abs(b-a)
  a,b = -0.5*max(abs(t1)),0.1*max(abs(t1)); tol = 0.20*abs(b-a)
  #a,b = -0.5*max(abs(t1)),0.0*max(abs(t1)); tol = 0.20*abs(b-a)
  sw = Stopwatch(); sw.restart()
  step = BrentMinFinder(misfitFunction).findMin(a,b,tol)
  print 'a =',a
  #print 'b =',b
  #print 'tol =',tol
  print 'step =',step
  report('line search',sw)
  add(mul(step,misfitFunction.p),s1,s1)

class AmplitudeMisfitFunction(BrentMinFinder.Function):
  def __init__(self,s1,p,src,rco,born):
    self.s1 = s1
    self.p = div(p,max(abs(p))) # normalized conjugate ascent direction
    self.src = BornOperatorS.getSourceArray(1)
    self.rco = BornOperatorS.getReceiverArray(1)
    self.rcp = BornOperatorS.getReceiverArray(1)
    self.src[0] = src
    self.rco[0] = rco
    self.rcp[0] = Receiver(xr,zr,nt)
    self.born = born
  def evaluate(self,a):
    print 'evaluating'
    s1p = add(self.s1,mul(a,self.p))
    self.born.applyForward(self.src,s1p,self.rcp)
    r = self.residual(self.rcp,self.rco)
    return sum(mul(r,r))
  def residual(self,rcp,rco):
    rcr = makeAmplitudeResiduals(rcp,rco)
    return rcr[0].getData()

class MisfitFunction(BrentMinFinder.Function):
  def __init__(self,p,isou,do):
    self.p = div(p,max(abs(p))) # normalized conjugate ascent direction
    self.isou = isou
    self.do = do[isou]
    #plot(self.do,title='do')
  def evaluate(self,a):
    print 'evaluating'
    s1p = add(s1,mul(a,self.p))
    ds = modelBornData(s0,s1p,self.isou)
    #plot(ds,title='a='+str(a))
    r = self.residual(ds,self.do)
    return sum(mul(r,r))

class WaveformMisfitFunction(MisfitFunction):
  def residual(self,ds,do,rd=None):
    if rd is None:
      rd = like(ds)
    sub(ds,do,rd)
    return rd

class xAmplitudeMisfitFunction(MisfitFunction):
  def residual(self,ds,do,ra=None):
    if ra is None:
      ra = like(ds)
    makeWarpedResidual(ds,do,ra=ra)
    return ra

#############################################################################
# Warping and more

def roughen(g):
  h = copy(g)
  #sigma = 0.100
  sigma = 1.0/fpeak
  ref = RecursiveExponentialFilter(sigma/dx)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref.apply(g,h)
  sub(g,h,h)
  copy(h,g)

def mask(g,l2=0):
  n1,n2 = len(g[0]),len(g)
  for i2 in range(l2):
    for i1 in range(n1):
      g[i2][i1] = 0.0

def makeAmplitudeResiduals(rcp,rco):
  ns = len(rcp)
  rcr = zeros(ns,Receiver) # residuals
  for isou in range(ns):
    dp = rcp[isou].getData()
    do = rco[isou].getData()
    ra = like(dp) # amplitude residual
    makeWarpedResidual(dp,do,ra=ra)
    rcr[isou] = Receiver(xr,zr,ra)
  return rcr

def makeAmplitudeResidual(rcp,rco):
  dp = rcp.getData()
  do = rco.getData()
  ra = like(dp) # amplitude residual
  makeWarpedResidual(dp,do,ra=ra)
  return Receiver(xr,zr,ra)

def makeWaveformResidual(rcp,rco):
  dp = rcp.getData()
  do = rco.getData()
  return Receiver(xr,zr,sub(dp,do))

def makeWarpedResidual(ds,do,ra=None,rt=None,rc=None,dw=None,v=None):
  if dw is None:
    dw = like(ds)
  if v is None:
    v = like(ds)
  if (rt is None) and (rc is not None):
    rt = like(ds)
  warp(ds,do,dw,v,rt) # right order
  if ra is not None:
    sub(ds,dw,ra)
  if rc is not None:
    sub(ds,dw,rc)
    add(rt,rc,rc)

def warp(ds,do,dw=None,v=None,rt=None):
  doRmsFilter = True # filter by local rms amplitudes
  td = 5 # time decimation
  rd = 1 # receiver decimation
  maxShift = 0.5 # max shift in seconds
  sw = Stopwatch(); sw.start()
  doc = copy(do)
  if doRmsFilter:
    ds,do = rmsFilter(ds,do,sigma=0.5*maxShift)
  ds,do = addRandomNoise(10.0,ds,do,sigma=1.0)
  ds = copy(nt/td,nr/rd,0,0,td,rd,ds)
  do = copy(nt/td,nr/rd,0,0,td,rd,do)
  #strainMax1,strainMax2 = 1.00,0.50
  strainMax1,strainMax2 = 0.50,0.20
  #strainMax1,strainMax2 = 0.25,0.10
  #strainMax1,strainMax2 = 0.10,0.05
  shiftMax = int(maxShift/(td*dt))
  warp = DynamicWarping(-shiftMax,shiftMax)
  warp.setErrorExponent(1.0)
  warp.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.AVERAGE)
  warp.setStrainMax(strainMax1,strainMax2)
  warp.setShiftSmoothing(32.0/td,8.0/rd) # shift smoothing
  warp.setErrorSmoothing(2) # number of smoothings of alignment errors
  u = warp.findShifts(ds,do)
  mul(td,u,u) # scale shifts to compensate for decimation
  li = LinearInterpolator()
  li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT)
  li.setUniform(nt/td,td*dt,0.0,nr/rd,rd*dx,0.0,u)
  if v is None:
    v = zerofloat(nt,nr)
  for ir in range(nr):
    z = ir*dz
    for it in range(nt):
      t = it*dt
      v[ir][it] = li.interpolate(t,z) # interpolate shifts
  if dw is None:
    dw = zerofloat(nt,nr)
  warp.applyShifts(v,doc,dw)
  if rt is not None:
    warp.applyShifts(v,timeDerivative(doc),rt)
    mul(v,rt,rt)
  #report('warping',sw)
  #plot(ds,title='f') # used for dynamic warping
  #plot(do,title='g') # used for dynamic warping
  #plot(dw,title='h') # warped
  #plot(v,sperc=100.0,title='shifts')
  return dw,v
  
def rmsFilter(ds,do,sigma=0.1):
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
  mul(den,x,x)
  mul(den,y,y)
  #plot(den,cmap=jet,title='rms_weights')
  return x,y

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
  random = Random(314159)
  s = sub(randfloat(random,n1,n2),0.5)
  RecursiveGaussianFilter(sigma).apply00(s,s) # bandlimited noise
  srms = sqrt(sum(mul(s,s))/n1/n2) # rms of noise
  mul(xrms/(srms*snr),s,s)
  return add(f,s),add(g,s)

def conjugateDirection(g,gm=None,pm=None):
  """Polak-Ribiere nonlinear conjugate gradient method.
  Parameters:
    g - gradient ascent direction for current iteration.
    gm - gradient ascent direction for previous iteration.
    pm - conjugate ascent direction for previous iteration.
  Returns:
    the conjugate ascent direction for the current iteration.
  """
  if gm is None and pm is None:
    # For first iteration, use steepest descent direction.
    return g
  else:
    b = sum(mul(g,sub(g,gm)))/sum(mul(gm,gm))
    if b<0.0:
      b = 0.0
      print "  CG DIRECTION RESET"
    return add(g,mul(b,pm))

def report(str,stopwatch):
  s = stopwatch.time()
  h = int(s/3600); s -= h*3600
  m = int(s/60); s -= m*60
  if h>0:
    print str+': %02d:%02d:%02d'%(h,m,s)
  else:
    print str+': %02d:%02d'%(m,s)

#############################################################################

def makeBornModel(s):
  """
  sigma0: smoothing for background model
  sigma1: smoothing for perturbation
  """
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
  if nx!=767:
    GaussianTaper.apply1(0.25,r,r)
  return s0,r

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
# Plots

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def pixels(x,cmap=gray,perc=100.0,sperc=None,cmin=0.0,cmax=0.0,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #sp.setSize(1010,740)
  sp.setSize(1000,650)
  #sp.setFontSizeForSlide(1.0,1.0)
  cb = sp.addColorBar()
  #cb.setWidthMinimum(100)
  cb.setWidthMinimum(150)
  if title:
    sp.addTitle(title)
    pass
  if len(x)==nz:
    pv = sp.addPixels(sz,sx,transpose(x))
    sp.setHLabel("Distance (km)")
    sp.setVLabel("Depth (km)")
  elif len(x[0])==nt:
    pv = sp.addPixels(st,sx,x)
    sp.setHLabel("Distance (km)")
    sp.setVLabel("Time (s)")
  else:
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
  if title and savDir:
    sp.paintToPng(360,3.33,savDir+title+'.png')
    write(savDir+title+'.dat',x)

def points(x):
  SimplePlot.asPoints(x)

#############################################################################
# Do everything on Swing thread.
import sys,time
class RunMain(Runnable):
  def run(self):
    start = time.time()
    if savDir is not None:
      print 'cleaning '+savDir.split('/')[-2]
      cleanDir(savDir)
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02d'%(h,m,s)
if __name__=='__main__':
  SwingUtilities.invokeLater(RunMain())
