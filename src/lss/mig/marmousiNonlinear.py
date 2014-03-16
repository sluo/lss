#############################################################################
# Marmousi nonlinear inversion. Testing shift predicted vs observed.

from imports import *

#############################################################################

savDir = None
#savDir = os.getenv('HOME')+'/Desktop/pngdat/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat2/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat3/'

timer = Timer()
#############################################################################

def main(args):
  #goNonlinearAmplitudeInversionO() # shift observed
  goNonlinearAmplitudeInversionP() # shift predicted

#############################################################################
# Setup

def getInputs():
  t,s,r,p = getMarmousiModelAndMask()
  e = mul(0.95,s) # erroneous background slowness

  # Wavefields.
  timer.start("allocating")
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)
  timer.stop("allocating")

  # Sources and receivers.
  src = zeros(ns,Source)
  rco = zeros(ns,Receiver)
  rcp = zeros(ns,Receiver)
  for isou in range(ns):
    rco[isou] = Receiver(xr,zr,nt)
    rcp[isou] = Receiver(xr,zr,nt)
  for isou in range(ns):
    src[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
  bornt = BornOperatorS(s,dx,dt,nabsorb,u,a) # true slowness

  # Observed data.
  timer.start("observed data")
  bornt.applyForward(src,r,rco) # observed data
  timer.stop("observed data")

  # Born modeling and solver.
  born = BornOperatorS(e,dx,dt,nabsorb,u,a)
  ref = RecursiveExponentialFilter(1.0/(fpeak*dx*sqrt(2.0)))
  bs = BornSolver(born,src,rcp,rco,ref,p) # solver

  # Warping.
  td = 4 # time decimation
  maxShift = 0.1 # max shift (seconds)
  #strainT,strainR,strainS = 0.50,0.20,-1.00 # 2D warping
  strainT,strainR,strainS = 0.50,0.20,1.00
  #smoothT,smoothR,smoothS = 16.0,4.0,4.0
  smoothT,smoothR,smoothS = 32.0,8.0,8.0
  warp = ImageWarping(
    strainT,strainR,strainS,smoothT,smoothR,smoothS,maxShift,dt,td)

  pixels(t,cmap=jet,title='t')
  pixels(s,cmap=jet,title='s')
  pixels(e,cmap=jet,title='e')
  pixels(sub(s,e),cmap=rwb,sperc=100.0,title='s-e')
  pixels(r,cmap=gray,sperc=100.0,title='r')
  pixels(p,cmap=gray,title='p')
  return born,bs,src,rcp,rco,warp,p

def getMarmousiModelAndMask():
  setupForMarmousi()
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

  return s,s0,mul(s1,m),p

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
  #xs,zs = [nx/4,nx/2,3*nx/4],[0,0,0]
  #xs,zs = rampint(1,15,52),fillint(0,52)
  #xs,zs = rampint(3,10,77),fillint(0,77) #
  xs,zs = rampint(3,8,96),fillint(0,96)
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

def goNonlinearAmplitudeInversionP():
  goNonlinearAmplitudeInversion(shiftPredicted=True)
def goNonlinearAmplitudeInversionO():
  goNonlinearAmplitudeInversion(shiftPredicted=False)
def goNonlinearAmplitudeInversion(shiftPredicted):
  niter = 2 # nonlinear iterations
  born,bs,src,rcp,rco,warping,mp = getInputs()
  r = zerofloat(nx,nz) # reflectivity image
  gm,pm = None,None # previous gradient, previous conjugate gradient
  for iiter in range(niter):
    g = computeGradientP(r,mp,born,src,rcp,rco,warping) if shiftPredicted\
      else computeGradientO(r,mp,born,src,rcp,rco,warping)
    p = conjugateDirection(g,gm,pm) # conjugate gradient direction
    add(mul(-0.05/max(abs(p)),p),r,r) # fixed step length XXX
    pixels(g,cmap=rwb,sperc=100.0,title='g_'+str(iiter))
    pixels(p,cmap=rwb,sperc=100.0,title='p_'+str(iiter))
    pixels(r,cmap=rwb,sperc=100.0,title='r_'+str(iiter))
    gm,pm = g,p

def computeGradientO(r,mp,born,src,rcp,rco,warping):
  g = zerofloat(nx,nz) # gradient
  w = zerofloat(nt,nr,ns) # warping shifts
  if sum(r)!=0.0:
    timer.start("predicted data")
    born.applyForward(src,r,rcp)
    timer.stop("predicted data")
    timer.start("warping")
    rcw = warping.warp(rcp,rco,w) # warp observed to predicted
    timer.stop("warping")
  else:
    rcw = rco
  rcr = rsub(rcp,rcw) # residual
  timer.start("gradient")
  born.applyAdjoint(src,rcr,g) # gradient
  timer.stop("gradient")
  roughen(g) # roughen
  mul(mp,g,g) # precondition
  return g

def computeGradientP(r,mp,born,src,rcp,rco,warping):
  g = zerofloat(nx,nz) # gradient
  w = zerofloat(nt,nr,ns) # warping shifts
  if sum(r)!=0.0:
    timer.start("predicted data")
    born.applyForward(src,r,rcp)
    timer.stop("predicted data")
    timer.start("warping")
    rcw = warping.warp(rco,rcp,w) # warp predicted to observed
    timer.stop("warping")
  else:
    rcw = rcp
  rcr = rsub(rcw,rco) # residual
  timer.start("gradient")
  #born.applyAdjoint(src,rcr,w,g) # gradient, including warping shifts
  born.applyAdjoint(src,rcr,g) # gradient, including warping shifts
  timer.stop("gradient")
  roughen(g) # roughen
  mul(mp,g,g) # precondition
  return g

def rsub(rcx,rcy,rcz=None):
  if rcz is None:
    rcz = zeros(ns,Receiver)
    for isou in range(ns):
      rcz[isou] = Receiver(xr,zr,nt)
  for isou in range(ns):
    sub(rcx[isou].getData(),rcy[isou].getData(),rcz[isou].getData())
  return rcz

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
