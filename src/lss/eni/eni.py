##############################################################################
# Eni data

from imports import *
from dnp import *

#subset = 'suba'
#subset = 'subb'
#subset = 'subc'
#subset = 'subd'
#subset = 'sube'
#subset = 'subf'
subset = 'subg'
subDir = '/data/sluo/eni/dat/'+subset+'/'

savDir = None
savDir = '/Users/sluo/Dropbox/png/'
#savDir = '/Users/sluo/Desktop/png/'
#savDir = '/home/sluo/Desktop/pngdat/'
#savDir = '/home/sluo/Desktop/pngdat2/'
#savDir = '/home/sluo/Desktop/pngdat3/'

timer = Timer()
##############################################################################

def main(args):
  #showFiles()
  #readFiles()
  #goBornData()
  #goAcousticData()
  #resimulateData()
  #compareWavelets()
  #estimateWavelet(toFile=False,rotate=0.25*FLT_PI,d2=False)
  #estimateWavelet(toFile=False,rotate=0.50*FLT_PI,d2=True)
  #goAmplitudeInversionO() # shift observed data
  #goAmplitudeInversionP() # shift predicted data
  #goNonlinearWaveformInversion() # nonlinear inversion using data residual
  #goNonlinearAmplitudeInversionO() # nonlinear inversion, shift observed
  #goNonlinearAmplitudeInversionP() # nonlinear inversion, shift predicted
  #goPredictedData() # resimulate predicted data
  goObservedData() # read and plot observed data
  #testSmoothingError()

def testSmoothingError():
  b = getBackground(); div(1.0,b,b)
  #s = getSlowness(); div(1.0,s,s)
  s = getBackground(vz=True,smin=0.72,sder=-0.0023); div(1.0,s,s)
  cmin = min(min(b),min(s))
  cmax = max(max(b),max(s))
  pixels(s,cmap=jet,cmin=cmin,cmax=cmax,title='provided')
  pixels(b,cmap=jet,cmin=cmin,cmax=cmax,title='smoothed')
  pixels(sub(b,s),cmap=jet,sperc=100.0,title='smoothed-provided') 

def getWavelet():
  return readWavelet()
  #return estimateWavelet(rotate=0.50*FLT_PI,d2=True)
  #return makeRickerWavelet() # Ricker wavelet

def setGlobals():
  global sx,sz,st#,ss,sr
  global nx,nz,nt,ns,nr
  global dx,dz,dt,ds,dr
  global fx,fz,ft,fs,fr
  global nxp,nzp,np
  global fmax,nabsorb,stride
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
    #npmax = 6 # max number of parallel shots
  elif subset=='subd' or subset=='sube':
    ss = Sampling(463,0.0125,1.225) # shot (relative offset)
    sr = Sampling(197,0.00625,-1.225) # receiver
    sx = Sampling(1121,0.00625,0.0) # distance
    sz = Sampling(181,0.00625,0.0) # depth
    st = Sampling(3751,0.0004,0.0) # time
    npmax = 16 # max number of parallel shots
  elif subset=='subf':
    ss = Sampling(1023,0.0125,1.225) # shot (relative offset)
    sr = Sampling(197,0.00625,-1.225) # receiver
    sx = Sampling(2241,0.00625,0.0) # distance
    sz = Sampling(181,0.00625,0.0) # depth
    st = Sampling(3751,0.0004,0.0) # time
    npmax = 16 # max number of parallel shots
  elif subset=='subg':
    ss = Sampling(863,0.0125,1.225) # shot (relative offset)
    sr = Sampling(197,0.00625,-1.225) # receiver
    sx = Sampling(1921,0.00625,0.0) # distance
    sz = Sampling(181,0.00625,0.0) # depth
    st = Sampling(3751,0.0004,0.0) # time
    npmax = 12 # max number of parallel shots
  #stride = 1
  stride = 2
  #stride = 4
  #stride = 1000
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
  src = zeros(ns,Source) # source
  rcp = zeros(ns,Receiver) # predicted data
  rco = zeros(ns,Receiver) # observed data
  timer.start('reading observed data')
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
  timer.stop('reading observed data')
  return src,rcp,rco

def getInputs():
  vz = True # 1D velocity
  if vz:
    if subset=='subc':
      smin,sder = None,None
    elif subset=='sube':
      smin,sder = 0.76,-0.0027
    elif subset=='subg':
      smin,sder = 0.72,-0.0023
    else:
      smin,sder = 0.75,-0.0025
  else:
    smin,sder = None,None
  print 'vz=%r'%vz

  # ImageWarping
  td = 4 # time decimation
  maxShift = 0.1 # max shift (seconds)
  #maxShift = 0.05 # max shift (seconds)
  #strainT,strainR,strainS = 0.2,0.2,0.2
  #strainT,strainR,strainS = 0.2,0.2,0.5
  #strainT,strainR,strainS = 0.3,0.3,0.3
  strainT,strainR,strainS = 0.4,0.4,0.3 # 3D warping (best)
  #strainT,strainR,strainS = 0.4,0.4,-1.0 # 2D warping
  #strainT,strainR,strainS = 0.4,-1.0,-1.0 # 1D warping
  #strainT,strainR,strainS = 0.4,0.4,0.4
  #strainT,strainR,strainS = 0.4,0.4,0.5
  #strainT,strainR,strainS = 0.5,0.5,0.5
  #strainT,strainR,strainS = 0.4,0.4,1.0
  #strainT,strainR,strainS = 0.5,0.5,1.0
  #strainT,strainR,strainS = 1.0,1.0,1.0
  smoothT,smoothR,smoothS = 32.0,4.0,4.0
  warping = ImageWarping(
    strainT,strainR,strainS,smoothT,smoothR,smoothS,maxShift,dt,td)
  #warping.setNoise(False) # XXX
  print 'bstrainT=%d'%int(ceil(1.0/strainT))
  print 'bstrainR=%d'%int(ceil(1.0/strainR))
  print 'bstrainS=%d'%int(ceil(1.0/strainS))

  # BornSolver
  src,rcp,rco = getSourceAndReceiver()
  s,m = getBackgroundAndMask(vz,smin,sder)
  timer.start('allocating')
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)
  timer.stop('allocating')
  #sigma = 0.25*nx*nz/(sum(s)*fmax*dx*sqrt(2.0)) # quarter wavelength
  sigma = 0.25*averageWavelength()/(sqrt(2.0)*dx) # quarter wavelength
  ref = RecursiveExponentialFilter(sigma)
  born = BornOperatorS(s,dx,dt,nabsorb,u,a)
  bs = BornSolver(born,src,rcp,rco,ref,m)

  pixels(s,cmap=jet,title='s')
  pixels(m,cmap=gray,title='m')
  return born,bs,src,rcp,rco,warping,s,m,ref

def goAmplitudeInversionO():
  """Shift observed data."""
  nouter,ninner,nfinal = 5,2,5 # outer, inner, final after last outer
  #nouter,ninner,nfinal = 0,0,5 # outer, inner, final after last outer
  zeroReflectivity = True # zero reflectivity between outer iterations
  born,bs,src,rcp,rco,warping,_,_,_ = getInputs()
  print 'nouter=%r'%nouter
  print 'ninner=%r'%ninner
  print 'nfinal=%r'%nfinal
  r = zerofloat(nx,nz) # reflectivity
  w = zerofloat(nt,nr,ns) # warping shifts
  for iouter in range(nouter+1):
    sw = Stopwatch(); sw.start()
    bs.solve(nfinal if iouter==nouter else ninner,r,r)
    pixels(r,cmap=gray,sperc=100.0,title='r%d'%iouter)
    if iouter<nouter and nouter>0:
      print "computing predicted data..."
      born.applyForward(src,r,rcp)
      print 'warping...'
      rcw = warping.warp(rcp,rco,w)
      bs.setObservedData(rcw)
      pixels(rcp[ns/2].getData(),title='rcp%d'%iouter)
      pixels(rcw[ns/2].getData(),title='rcw%d'%iouter)
      pixels(w[ns/2],cmap=rwb,sperc=100.0,title='w%d'%iouter)
      if iouter==(nouter-1) and savDir is not None:
        dp = zerofloat(nt,nr,ns)
        for isou in range(ns):
          copy(rcp[isou].getData(),dp[isou])
        write(savDir+'dp.dat',dp)
        write(savDir+'w.dat',w)
    if zeroReflectivity:
      zero(r)
    sw.stop(); print 'outer: %.1f minutes'%(sw.time()/60.0)
  pixels(rco[ns/2].getData(),title='rco')

def goAmplitudeInversionP():
  """Shift predicted data."""
  #nouter,ninner,nfinal = 5,1,5 # outer, inner, final after last outer
  nouter,ninner,nfinal = 5,2,5 # outer, inner, final after last outer
  #nouter,ninner,nfinal = 0,0,5 # outer, inner, final after last outer
  zeroReflectivity = True # zero reflectivity between outer iterations
  born,bs,src,rcp,rco,warping,_,_,_ = getInputs()
  print 'nouter=%r'%nouter
  print 'ninner=%r'%ninner
  print 'nfinal=%r'%nfinal
  r = zerofloat(nx,nz) # reflectivity
  w = zerofloat(nt,nr,ns) # warping shifts
  for iouter in range(nouter+1):
    sw = Stopwatch(); sw.start()
    bs.solve(nfinal if iouter==nouter else ninner,r,r);
    pixels(r,cmap=gray,sperc=100.0,title='r%d'%iouter)
    if iouter<nouter and nouter>0:
      print "computing predicted data..."
      born.applyForward(src,r,rcp) # simulate predicted data
      print 'warping...'
      rcw = warping.warp(rco,rcp,w) # estimate new time shifts
      bs.setTimeShifts(w) # set new time shifts
      pixels(rcp[ns/2].getData(),title='rcp%d'%iouter)
      pixels(rcw[ns/2].getData(),title='rcw%d'%iouter)
      pixels(w[ns/2],cmap=rwb,sperc=100.0,title='w%d'%iouter)
      if iouter==(nouter-1) and savDir is not None:
        dp = zerofloat(nt,nr,ns)
        for isou in range(ns):
          copy(rcp[isou].getData(),dp[isou])
        write(savDir+'dp.dat',dp)
        write(savDir+'w.dat',w)
    if zeroReflectivity:
      zero(r)
    sw.stop(); print 'outer: %.1f minutes'%(sw.time()/60.0)
  pixels(rco[ns/2].getData(),title='rco')

##############################################################################
# Nonlinear inversion with line search

def xgoObservedData():
  src,rcp,rco = getSourceAndReceiver()
  do = zerofloat(nt,nr,ns)
  for isou in range(ns):
    copy(rco[isou].getData(),do[isou])
  #do = randfloat(nt,nr,ns)
  ss = Sampling(ns,0.0125*stride,1.225) # shot (relative offset)
  sr = Sampling(197,0.00625,-1.225) # receiver
  st = Sampling(3751,0.0004,0.0) # time
  pixels3(do,s1=st,s2=sr,s3=ss,cmin=-8.0,cmax=8.0,title='ppp3')

def goObservedData():
  rdir = '/Users/sluo/Desktop/subg/nonlinear/vz/ares05/warp1d/'
  #rdir = '/Users/sluo/Desktop/subg/nonlinear/vz/ares05/warp2d/'
  #rdir = '/Users/sluo/Desktop/subg/nonlinear/vz/ares05/warp3d/'
  do = zerofloat(nt,nr,ns); read(rdir+'do.dat',do)
  dp = zerofloat(nt,nr,ns); read(rdir+'dp.dat',dp)
  dw = zerofloat(nt,nr,ns); read(rdir+'dw.dat',dw)
  w = zerofloat(nt,nr,ns); read(rdir+'w.dat',w); mul(1000.0*dt,w,w)
  ss = Sampling(ns,0.0125*stride,1.225) # shot (relative offset)
  sr = Sampling(197,0.00625,-1.225) # receiver
  st = Sampling(3751,0.0004,0.0) # time
  k1,k2,k3 = nt/2,nr-1,25 # x=1.85 (shown in CWP report)
  #k1,k2,k3 = nt/2,nr-1,173 # x=5.55 (smallest shifts)
  pixels3(do,s1=st,s2=sr,s3=ss,cmin=-8.0,cmax=8.0,
    k1=k1,k2=k2,k3=k3,title='do')
  pixels3(dp,s1=st,s2=sr,s3=ss,cmin=-8.0,cmax=8.0,
    k1=k1,k2=k2,k3=k3,title='dp')
  pixels3(dw,s1=st,s2=sr,s3=ss,cmin=-8.0,cmax=8.0,
    k1=k1,k2=k2,k3=k3,title='dw')
  pixels3(w,s1=st,s2=sr,s3=ss,cmap=jet,lineColor=Color.BLACK,
    k1=k1,k2=k2,k3=k3,cbar='Time shift (ms)',title='w',
    cmin=-17.5,cmax=17.5)
    #sperc=100.0)

def goPredictedData():
  #rdir = '/home/sluo/Desktop/subg/nonlinear/vz/dres00/'
  rdir = '/home/sluo/Desktop/subg/nonlinear/vz/ares05/'
  rfile = rdir+'r9.dat'
  sfile = rdir+'s.dat'
  r = zerofloat(nz,nx); read(rfile,r); r = transpose(r)
  s = zerofloat(nz,nx); read(sfile,s); s = transpose(s)
  pixels(r,sperc=99,title='r')
  pixels(s,cmap=jet,title='s')
  born,bs,src,rcp,rco,warping,ss,m,ref = getInputs()
  Check.argument(sum(sub(s,ss))==0.0,'sum(sub(s,ss))==0.0')
  timer.start('predicted data')
  born.applyForward(src,r,rcp) # simulate predicted data
  timer.stop('predicted data')
  w = zerofloat(nt,nr,ns) # warping shifts
  timer.start('warping')
  rcw = warping.warp(rcp,rco,w) # warping
  timer.stop('warping')
  ksou = [0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425]
  clip = getSymmetricClip(0.98,rco[ns/2].getData())
  for isou in ksou:
    x = fs+isou*stride*ds
    pixels(rco[isou].getData(),cmin=-clip,cmax=clip,title='do_x=%f'%x)
    pixels(rcp[isou].getData(),cmin=-clip,cmax=clip,title='dp_x=%f'%x)
    pixels(rcw[isou].getData(),cmin=-clip,cmax=clip,title='dw_x=%f'%x)
    pixels(w[isou],cmap=rwb,sperc=99.9,title='w_x=%f'%x)
  do = zerofloat(nt,nr,ns)
  dp = zerofloat(nt,nr,ns)
  dw = zerofloat(nt,nr,ns)
  for isou in range(ns):
    copy(rco[isou].getData(),do[isou])
    copy(rcp[isou].getData(),dp[isou])
    copy(rcw[isou].getData(),dw[isou])
  write(savDir+'do.dat',do)
  write(savDir+'dp.dat',dp)
  write(savDir+'dw.dat',dw)
  write(savDir+'w.dat',w)

def goNonlinearWaveformInversion():
  goNonlinearInversion(0)
def goNonlinearAmplitudeInversionO():
  goNonlinearInversion(1)
def goNonlinearAmplitudeInversionP():
  goNonlinearInversion(2)
def goNonlinearInversion(residualType=0):
  """Nonlinear inversion via a line search.
  residualType = 0: data residual
  residualType = 1: amplitude residual with shifted observed data
  residualType = 2: amplitude residual with shifted predicted data
  """
  #niter = 2 # number of iterations
  niter = 10 # number of iterations
  print 'niter=%d'%niter
  born,bs,src,rcp,rco,warping,s,m,ref = getInputs()
  if residualType==0:
    res = WaveformResidual()
  elif residualType==1:
    res = AmplitudeResidualO(warping)
  elif residualType==2:
    res = AmplitudeResidualP(warping)
  r = zerofloat(nx,nz) # reflectivity
  gm,pm = None,None # gradient & cg directions
  for iiter in range(niter):
    print ''
    timer.start('ITERATION %d'%iiter)
    if residualType==0:
      g = computeGradient(iiter,r,born,src,rcp,rco)
    elif residualType==1:
      g = computeGradientO(iiter,r,born,src,rcp,rco,warping)
    elif residualType==2:
      g = computeGradientP(iiter,r,born,src,rcp,rco,warping)
    preconditionGradient(s,m,g)
    p = g # steepest descent
    #p = conjugateDirection(g,gm,pm)
    if niter>1:
      lineSearchUpdate(p,r,src,rco,born,res)
    gm,pm = g,p
    pixels(g,cmap=rwb,sperc=100.0,title='g'+str(iiter))
    pixels(p,cmap=rwb,sperc=100.0,title='p'+str(iiter))
    pixels(r,cmap=gray,sperc=100.0,title='r'+str(iiter))
    timer.stop('ITERATION %d'%iiter)
  pixels(rco[ns/2].getData(),title='rco')

def preconditionGradient(s,m,g):
  sigma = 0.25*averageWavelength()/(sqrt(2.0)*dx) # quarter wavelength
  ref = RecursiveExponentialFilter(sigma)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  for i in range(2):
    roughen(g,ref)
  p = div(1.0,mul(s,s)) # velocity-squared preconditioner
  mul(1.0/max(p),p,p) # normalized
  mul(p,g,g) # apply preconditioner
  mul(m,g,g) # apply mask

def xpreconditionGradient(s,m,g):
  h = like(g)
  sigma = 0.25*averageWavelength()/dx # quarter wavelength
  RecursiveGaussianFilter(sigma).apply00(g,h) # smooth...
  sub(g,h,g) # ...and subtract
  p = div(1.0,mul(s,s)) # velocity-squared preconditioner
  mul(1.0/max(p),p,p) # normalized
  mul(p,g,g) # apply preconditioner
  mul(m,g,g) # apply mask

def computeGradient(iiter,r,born,src,rcp,rco):
  """Gradient for waveform residual."""
  g = zerofloat(nx,nz) # gradient
  if iiter>0:
    timer.start('predicted data')
    born.applyForward(src,r,rcp) # simulate predicted data
    timer.stop('predicted data')
    pixels(rcp[ns/2].getData(),title='rcp%d'%iiter)
  rsub(rcp,rco,rcp) # residual
  timer.start('gradient')
  born.applyAdjoint(src,rcp,g) # gradient
  timer.stop('gradient')
  return g

def computeGradientO(iiter,r,born,src,rcp,rco,warping):
  """Gradient for amplitude residual with shifted observed data."""
  g = zerofloat(nx,nz) # gradient
  w = zerofloat(nt,nr,ns) # warping shifts
  if iiter>0:
    timer.start('predicted data')
    born.applyForward(src,r,rcp) # simulate predicted data
    timer.stop('predicted data')
    timer.start("warping")
    rcw = warping.warp(rcp,rco,w) # warp observed to predicted
    timer.stop("warping")
    pixels(rcp[ns/2].getData(),title='rcp%d'%iiter)
    pixels(rcw[ns/2].getData(),title='rcw%d'%iiter)
    pixels(w[ns/2],cmap=rwb,sperc=100.0,title='w%d'%iiter)
  else:
    rcw = rco
  rsub(rcp,rcw,rcp) # residual
  timer.start('gradient')
  born.applyAdjoint(src,rcp,g) # gradient
  timer.stop('gradient')
  return g

def computeGradientP(iiter,r,born,src,rcp,rco,warping):
  """Gradient for amplitude residual with shifted predicted data."""
  g = zerofloat(nx,nz) # gradient
  w = zerofloat(nt,nr,ns) # warping shifts
  if iiter>0:
    timer.start('predicted data')
    born.applyForward(src,r,rcp) # simulate predicted data
    timer.stop('predicted data')
    timer.start("warping")
    rcw = warping.warp(rco,rcp,w) # warp predicted to observed
    timer.stop("warping")
    pixels(rcp[ns/2].getData(),title='rcp%d'%iiter)
    pixels(rcw[ns/2].getData(),title='rcw%d'%iiter)
    pixels(w[ns/2],cmap=rwb,sperc=100.0,title='w%d'%iiter)
  else:
    rcw = rcp
  #rsub(rcw,rco,rcw) # residual
  #timer.start('gradient')
  #born.applyAdjoint(src,rcw,w,g) # gradient, including warping shifts
  #timer.stop('gradient')
  rsub(rcw,rco,rcp) # residual
  timer.start('gradient')
  born.applyAdjoint(src,rcp,w,g) # gradient, including warping shifts
  timer.stop('gradient')
  return g

def lineSearchUpdate(p,r,src,rco,born,res):
  amin,amax,atol = -0.10,0.01,0.005 # step length bounds and tolerance
  #amin,amax,atol = -0.10,0.01,0.01 # step length bounds and tolerance
  #amin,amax,atol = -0.05,0.01,0.01 # step length bounds and tolerance
  #amin,amax,atol = -0.02,0.01,0.005 # step length bounds and tolerance
  nsou = 16 # number of shots used to evaluate misfit function
  #nsou = 4 # number of shots used to evaluate misfit function
  misfit = MisfitFunction(p,r,src,rco,born,res,nsou)
  timer.start('line search')
  aopt = BrentMinFinder(misfit).findMin(amin,amax,atol)
  timer.stop('line search')
  print 'neval=%d'%misfit.neval
  print 'amin=%f'%amin
  #print 'amax=%f'%amax
  #print 'atol=%f'%atol
  print 'aopt=%f'%aopt
  add(mul(aopt,misfit.p),r,r)

class MisfitFunction(BrentMinFinder.Function):
  def __init__(self,p,r,src,rco,born,res,nsou):
    self.r = r
    self.p = div(p,max(abs(p))) # normalized conjugate ascent direction
    self.src = zeros(nsou,Source) # source
    self.rcp = zeros(nsou,Receiver) # predicted data
    self.rco = zeros(nsou,Receiver) # observed data
    for isou in range(nsou):
      ksou = (1+isou)*ns/(1+nsou)
      #print 'ksou=%d'%ksou
      self.src[isou] = src[ksou]
      self.rcp[isou] = Receiver(rco[ksou])
      self.rco[isou] = rco[ksou]
    self.born = born
    self.res = res
    self.nsou = nsou
    self.neval = 0
  def evaluate(self,atry):
    self.neval += 1
    rr = add(self.r,mul(atry,self.p))
    self.born.applyForward(self.src,rr,self.rcp)
    misfit = 0.0
    for isou in range(self.nsou):
      dp,do = self.rcp[isou].getData(),self.rco[isou].getData()
      dr = self.res.compute(dp,do)
      mul(dr,dr,dr)
      misfit += sum(dr)
    return misfit
    #class Reduce(Parallel.ReduceInt):
    #  def compute(self,isou):
    #    dp,do = self.rcp[isou].getData(),self.rco[isou].getData()
    #    dr = self.res.compute(dp,do)
    #    mul(dr,dr,dr)
    #    return sum(dr)
    #  def combine(self,misfit1,misfit2):
    #    return misfit1+misfit2
    #return Parallel.reduce(self.nsou,Reduce())

class WaveformResidual():
  def compute(self,dp,do):
    return sub(dp,do)

class AmplitudeResidualO():
  def __init__(self,warping):
    self.warping = warping
  def compute(self,dp,do):
    u = self.warping.findShifts(dp,do)
    dw = self.warping.applyShifts(u,do)
    sub(dp,dw,dw)
    return dw

class AmplitudeResidualP():
  def __init__(self,warping):
    self.warping = warping
  def compute(self,dp,do):
    u = self.warping.findShifts(do,dp)
    dw = self.warping.applyShifts(u,dp)
    sub(dw,do,dw)
    return dw

def roughen(g,ref):
  h = copy(g)
  ref.apply(g,h)
  sub(g,h,h)
  copy(h,g)

def rsub(rcx,rcy,rcz):
  nsou = len(rcz)
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      sub(rcx[isou].getData(),rcy[isou].getData(),rcz[isou].getData())
  Parallel.loop(nsou,Loop())

def conjugateDirection(gi,gm=None,pm=None):
  """Polak-Ribiere nonlinear conjugate gradient method.
  Parameters:
    gi - gradient ascent direction for current iteration.
    gm - gradient ascent direction for previous iteration.
    pm - conjugate ascent direction for previous iteration.
  Returns:
    the conjugate ascent direction for the current iteration.
  """
  if gm is None and pm is None:
    return gi
  else:
    beta = sum(mul(gi,sub(gi,gm)))/sum(mul(gm,gm))
    if beta<0.0:
      print "restarting CG (beta<0.0)"
      beta = 0.0
    return add(gi,mul(beta,pm))

##############################################################################

def showFiles():
  #vz,smin,sder = True,0.75,-0.0025 # 1D velocity
  #vz,smin,sder = True,0.73,-0.0024 # 1D velocity
  vz,smin,sder = True,0.72,-0.0023 # 1D velocity
  #vz,smin,sder = True,0.76,-0.0027 # 1D velocity
  w = readWavelet(); points(w)
  d = getGather(ns/2); pixels(d,perc=99.8); points(d[196])
  s = getSlowness(False); pixels(s,cmap=jet)
  s0,_ = getBackgroundAndMask(False,None,None)
  e0,m = getBackgroundAndMask(vz,smin,sder)
  pixels(s0,cmap=jet)
  pixels(e0,cmap=jet)
  pixels(m,cmap=gray)
  pixels(sub(s0,e0),cmap=rwb,sperc=100.0)
  ts = transpose(s0)
  te = transpose(e0)
  SimplePlot.asPoints(ts[nx/2])
  SimplePlot.asPoints(te[nx/2])

def readFiles():
  ra,rb,rc = zerofloat(nz,nx),zerofloat(nz,nx),zerofloat(nz,nx)
  #read('/home/sluo/Desktop/save/eni/subc/iter005vz/r0.dat',ra);
  #read('/home/sluo/Desktop/save/eni/subc/iter525vz/r5.dat',rb);
  #read('/home/sluo/Desktop/save/eni/subc/n/iter005vz/r0.dat',ra);
  #read('/home/sluo/Desktop/save/eni/subc/n/O525vz/r5.dat',rb);
  #read('/home/sluo/Desktop/save/eni/subc/n/P525vz/r5.dat',rc);
  #read('/home/sluo/Desktop/sube/nonlinear/aresVz6/r9.dat',ra);
  #read('/home/sluo/Desktop/sube/nonlinear/aresVz7/r9.dat',rb);
  #read('/home/sluo/Desktop/sube/nonlinear/vfile/dres/r9.dat',ra);
  #read('/home/sluo/Desktop/sube/nonlinear/vz/ares08/r9.dat',rb);
  #read('/home/sluo/Desktop/sube/nonlinear/vz/ares12/r19.dat',rc);
  read('/home/sluo/Desktop/subg/nonlinear/vfile/dres00/r9.dat',ra);
  read('/home/sluo/Desktop/subg/nonlinear/vz/dres00/r9.dat',rb);
  read('/home/sluo/Desktop/subg/nonlinear/vz/ares05/r9.dat',rc);
  pixels(ra,cmap=gray,sperc=98.0)
  pixels(rb,cmap=gray,sperc=98.0)
  pixels(rc,cmap=gray,sperc=98.0)
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
  #w = readWavelet()
  w = estimateWavelet(toFile=False,rotate=0.0*FLT_PI,d2=False)
  v = estimateWaterBottomWavelet(rotate=0.0*FLT_PI,d2=False)
  mul(1.0/max(abs(w)),w,w)
  mul(1.0/max(abs(v)),v,v)
  wm = zeroint(1); max(w,wm); wm = wm[0]
  vm = zeroint(1); max(v,vm); vm = vm[0]
  t = copy(v); zero(v); copy(nt-(vm-wm),vm-wm,t,0,v)
  sst = Sampling(301,dt,0.0)
  wc = copy(301,30,w)
  vc = copy(301,30,v)
  points(wc,sst,cmin=-1.1,cmax=1.1,label='Time (s)',title='w')
  points(vc,sst,cmin=-1.1,cmax=1.1,label='Time (s)',title='v')
  points(wc,sst,x2=vc,cmin=-1.1,cmax=1.1,label='Time (s)',title='wv')

def readWavelet():
  print 'reading wavelet'
  w = zerofloat(nt)
  #read(subDir+'w.dat',w)
  read('/data/sluo/eni/dat/wavelet.dat',w)
  #points(w)
  return w

def estimateWavelet(toFile=False,rotate=0.25*FLT_PI,d2=False):
  print 'estimating wavelet'
  print '  rotate=%f'%rotate
  print '  d2=%r'%d2
  w = estimateZeroPhaseWavelet()
  #w = estimateMinimumPhaseWavelet()

  # Phase rotation and bandpass filter.
  w = Frequency(nt).rotateAndDifferentiate(w,rotate,d2)
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

  amplitudeSpectrum(w,plot=True)
  points(w)
  if toFile:
    print 'writing wavelet'
    #write(subDir+'w.dat',w)
    write('/data/sluo/eni/dat/wavelet.dat',w)
  return w

def estimateWaterBottomWavelet(rotate=0.25*FLT_PI,d2=False):
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
  w = Frequency(nt).rotateAndDifferentiate(w,rotate,d2)
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
  #points(c)
  RecursiveGaussianFilter(64.0).apply0(c,c)
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
  w = Frequency(nt).rotateAndDifferentiate(w,0.25*FLT_PI,d2=False)
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
  def rotateAndDifferentiate(self,rx,p=0.25*FLT_PI,d2=False):
    self.set(rx)
    sf = self.fft.getFrequencySampling1()
    nf = sf.count
    cy = self.fft.applyForward(self.rx) # forward FFT
    t = zerofloat(2*nf)
    for i in range(nf):
      w = sf.getValue(i)
      t[2*i  ] = w*w*cos(p) if d2 else cos(p)
      t[2*i+1] = w*w*sin(p) if d2 else sin(p)
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

def getSlowness(vz=False,smin=None,sder=None):
  v = getVelocity(vz)
  div(1000.0,v,v)
  if vz: # linear slowness model
    if smin is not None and sder is not None:
      v00 = v[0][0]
      c = copy(v)
      for iz in range(nz):
        siz = smin+iz*sder
        fill(siz,v[iz])
      for iz in range(nz):
        for ix in range(nx):
          if c[iz][ix]==v00:
            v[iz][ix] = v00;
    else:
      v = transpose(v)
      v00 = v[0][0]
      a,b,c = 0.0,0.0,0.0
      for ix in range(int(0.31*nx),nx): # 0.31 just because
      #for ix in range(nx):
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

def getBackgroundAndMask(vz=False,smin=None,sder=None):
  return getBackground(vz,smin,sder),getMask()

def getBackground(vz=False,smin=None,sder=None):
  s = getSlowness(vz,smin,sder)
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

def amplitudeSpectrum(x,plot=False):
  fft = Fft(x)
  y = cabs(fft.applyForward(x))
  if plot:
    sw = fft.getFrequencySampling1()
    nw,dw,fw = sw.count,sw.delta,sw.first
    sf = Sampling(nw,dw/dt,fw)
    points(y,sf,label='Frequency (Hz)',title='amp')
  return y

def bandpass(x,y,fmin=None,fmax=None):
  fnyq = 0.5/dt # Nyquist
  if fmin is None:
    klower = 0.01
  else:
    klower = fmin/fnyq
    #klower = 0.5*(fmin/fnyq)
  if fmax is None:
    kupper = 40.0/fnyq
    #kupper = 0.5*(40.0/fnyq)
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
    clip = getSymmetricClip(sperc,x)
    pv.setClips(-clip,clip)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if title and savDir:
    write(savDir+title+'.dat',x)
    sp.paintToPng(360,3.33,savDir+title+'.png')
  return sp

def pixels3(x,s1=None,s2=None,s3=None,cmap=gray,lineColor=Color.YELLOW,
    perc=100.0,sperc=None,cmin=0.0,cmax=0.0,cbar=None,
    k1=None,k2=None,k3=None,title=None):
  if s1 is None: s1 = Sampling(len(x[0][0]))
  if s2 is None: s2 = Sampling(len(x[0]))
  if s3 is None: s3 = Sampling(len(x))
  o = PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT
  a = PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM
  panel = PlotPanelPixels3(o,a,s1,s2,s3,x)
  panel.setColorModel(cmap)
  panel.setSlices(
    k1 if k1 else s1.count/2,
    k2 if k2 else s2.count/2,
    k3 if k3 else s3.count/2)
  panel.setLabel1("Time (s)")
  panel.setLabel2("Offset (km)")
  panel.setLabel3("Position (km)")
  panel.setLineColor(lineColor)
  if cbar:
    panel.addColorBar(cbar)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(160)
  panel.setBackground(Color.WHITE)
  #panel.setInterval1(0.1)
  #panel.setInterval2(1.0)
  #panel.setInterval3(1.0)
  if cmin<cmax:
    panel.setClips(cmin,cmax)
  if perc<100.0:
    panel.setPercentiles(100.0-perc,perc)
  if sperc is not None:
    clip = getSymmetricClip(sperc,x)
    panel.setClips(-clip,clip)
  panel.mosaic.setHeightMinimum(0,300)
  panel.mosaic.setWidthMinimum(0,330)
  panel.mosaic.setWidthElastic(0,0)
  panel.mosaic.setHeightElastic(0,0)
  panel.setVLimits(1,0.3,1.4)
  frame = PlotFrame(panel)
  frame.setSize(980,800)
  #frame.setFontSizeForSlide(1.0,1.0,16.0/9.0)
  frame.setFontSize(45.0)
  frame.setVisible(True)
  if title and savDir:
    frame.paintToPng(1920,1.0,savDir+title+'.png')

def points(x,s1=None,x2=None,cmin=0.0,cmax=0.0,label=None,title=None):
  mul(1.0/max(abs(x)),x,x) # XXX
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  if x2 is not None:
    pv2 = sp.addPoints(s1,x2) if s1 else sp.addPoints(x2)
    pv2.setLineWidth(4.0)
    pv2.setLineColor(Color.RED)
    #pv2.setLineStyle(PointsView.Line.DASH)
  pv = sp.addPoints(s1,x) if s1 else sp.addPoints(x)
  pv.setLineWidth(4.0)
  if cmin<cmax:
    sp.setVLimits(cmin,cmax)
  #sp.setVLimits(0.0,1.1) # XXX
  #sp.setHLimits(0.0,100.0) # XXX
  if label:
    sp.setHLabel(label)
  sp.setFontSize(45.0)
  sp.setSize(800,805)
  if title and savDir:
    sp.paintToPng(1000,1.0,savDir+title+'.png')

def getSymmetricClip(sperc,x):
  clips = Clips(100-sperc,sperc,x)
  return max(abs(clips.getClipMin()),abs(clips.getClipMax()))

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
