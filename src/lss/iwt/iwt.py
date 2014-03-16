#############################################################################
# Image-warping tomography

from imports import *

#############################################################################

savDir = None
#savDir = os.getenv('HOME')+'/Desktop/pngdat/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat2/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat3/'
savDir = os.getenv('HOME')+'/Desktop/pngdat4/'
#savDir = '/Users/sluo/Dropbox/pngdat/'
#savDir = '/Users/sluo/Dropbox/pngdat2/'

datDir = None
#datDir = './dat/'

plots = True
#plots = False

#############################################################################

# TODO: change Ricker to Gaussian derivative
# TODO: second time derivative in adjoint source

timer = Timer()
def main(args):
  #getModelAndMask()
  #goMigration()
  #goGradient()
  goInversionX()

def setupForLayered():
  #global zs,xs,zr,xr,ns,nr
  global ns,nr
  global sz,sx,st,nz,nx,nt,nxp,nzp,dz,dx,dt
  global fpeak,nabsorb,np,sigma0,sigma1

  #sz = Sampling(301,0.012,0.0) # for layered model
  #sx = Sampling(501,0.012,0.0)
  #st = Sampling(2750,0.0015,0.0)
  sz = Sampling(126,0.016,0.0)
  sx = Sampling(189,0.016,0.0)
  #sx = Sampling(251,0.016,0.0)
  #st = Sampling(2003,0.0010,0.0)
  st = Sampling(1003,0.0015,0.0)
  nz,nx,nt = sz.count,sx.count,st.count
  dz,dx,dt = sz.delta,sx.delta,st.delta
  fz,fx,ft = sz.first,sx.first,st.first
  #xs,zs = rampint(0,1,nx),fillint(0,nx)
  #xr,zr = rampint(0,1,nx),fillint(0,nx)
  #ns,nr = len(xs),len(xr)
  ns = nx
  fpeak = 20.0 # Ricker wavelet peak frequency
  nabsorb = 22 # absorbing boundary size
  nxp,nzp = nx+2*nabsorb,nz+2*nabsorb
  np = min(16,ns) # number of parallel sources
  sigma0 = 4.0/fpeak
  sigma1 = 1.0/fpeak
  print 'ns=%d'%ns

  #return getSingleLayeredModelAndMask()
  return getMultiLayeredModelAndMask()

def getSourcesAndReceivers():
  #maxOffset = 80
  maxOffset = nx # full offsets, fixed spread
  src = zeros(nx,Source)
  rco = zeros(nx,Receiver)
  rcp = zeros(nx,Receiver)
  for ix in range(nx):
    fxr,lxr = max(ix-maxOffset,0),min(ix+maxOffset,nx-1)
    nxr = 1+lxr-fxr
    xr,zr = rampint(fxr,1,nxr),fillint(0,nxr)
    rco[ix] = Receiver(xr,zr,nt)
    rcp[ix] = Receiver(xr,zr,nt)
    src[ix] = Source.RickerSource(ix,0,dt,fpeak)
  ns = nx
  return src,rcp,rco

def getSingleLayeredModelAndMask():
  constantBackground = True
  muteFirstLayer = True
  tb = 0.25 # background slowness
  t = getLayered2()
  t0,t1 = makeBornModel(t)
  if muteFirstLayer:
    for iz in range(nz/2):
      for ix in range(nx):
        t1[iz][ix] = 0.0
  #GaussianTaper.apply1(0.25,t1,t1)
  if constantBackground:
    fill(tb,t0)
  m = fillfloat(1.0,nx,nz)
  for iz in range(nz/3-6):
    for ix in range(nx):
      m[iz][ix] = 0.0
  RecursiveExponentialFilter(1.0).apply2(m,m)
  return t,t0,t1,m

def getMultiLayeredModelAndMask():
  tb = 0.25 # background slowness
  nlayer = 10
  #vstart,vstep = 2.0,1.0
  #vstart,vstep = float(nlayer),-1.0
  vstart,vstep = 100.0,-1.0
  t = fillfloat(vstart,nx,nz)
  for ilayer in range(2,nlayer):
    vstart += vstep
    for iz in range(ilayer*nz/nlayer,(1+ilayer)*nz/nlayer):
      fill(vstart,t[iz])
  t0,t1 = makeBornModel(t)
  mul(1.0/max(abs(t1)),t1,t1) # normalize reflectivity
  fill(tb,t0)
  m = fillfloat(1.0,nx,nz)
  sigma = 3
  for iz in range(nz-3*sigma):
    for ix in range(nx):
      if t[iz+2*sigma][ix]==t[0][0]:
        m[iz][ix] = 0.0
  ref = RecursiveExponentialFilter(sigma)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref.apply2(m,m)
  return t,t0,t1,m

def getLayered2(s0mul=1.0):
  sfill = 1.0/1.5
  t = fillfloat(sfill,nz,nx)
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

def addGaussianError(s,emax=0.010,sigma=None):
  g = zerofloat(nx,nz)
  #g[nz/3][nx/2] = 1.0
  #g[2*nz/5][nx/2] = 1.0
  g[nz/2][nx/2] = 1.0
  if sigma is None:
    RecursiveGaussianFilter(nz/8.0).apply00(g,g)
  else:
    RecursiveGaussianFilter(sigma).apply00(g,g)
  mul(emax/max(g),g,g)
  return add(s,g)

def getModelAndMask():
  t,s,r,m = setupForLayered()
  e = copy(s)
  #s = mul(0.90,e) # erroneous background slowness
  #s = mul(1.10,e) # erroneous background slowness
  #s = addGaussianError(e,emax=0.05)
  #s = addGaussianError(e,emax=-0.05)
  s = addGaussianError(e,emax=0.05,sigma=nz/10.0) ###
  #s = addGaussianError(e,emax=-0.05,sigma=nz/10.0) ###
  #s = addGaussianError(e,emax=-0.05,sigma=nz/10.0)
  #s = addGaussianError(e,emax=0.04)
  #s = addGaussianError(e,emax=-0.04)
  #s = addGaussianError(e,emax=-0.04)
  pixels(t,cmap=jet,title='t')
  pixels(s,cmap=jet,title='s')
  pixels(e,cmap=jet,title='e')
  if m is not None:
    pixels(m,cmap=gray,title='m')
  pixels(sub(s,e),cmap=rwb,sperc=100.0,title='s-e')
  pixels(r,cmap=gray,sperc=100.0,title='r')
  return t,s,e,r,m

def getInputs(modelChanged=True):
  t,s,e,r,m = getModelAndMask()

  # Wavefields.
  timer.start('allocating')
  u4 = SharedFloat4(nxp,nzp,nt,np)
  a4 = SharedFloat4(nxp,nzp,nt,np)
  b4 = SharedFloat4(nxp,nzp,nt,np)
  timer.stop('allocating')

  # Sources and receivers.
  src,rcp,rco = getSourcesAndReceivers()
  #src = zeros(ns,Source)
  #rco = zeros(ns,Receiver)
  #rcp = zeros(ns,Receiver)
  #for isou in range(ns):
  #  rco[isou] = Receiver(xr,zr,nt)
  #  rcp[isou] = Receiver(xr,zr,nt)
  #  src[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
  if modelChanged: # if model has changed, resimulated data...
    bornt = BornOperatorS(s,dx,dt,nabsorb,u4,a4) # true slowness
    timer.start('observed data')
    bornt.applyForward(src,r,rco) # observed data
    timer.stop('observed data')
    timer.start('writing observed data')
    for isou in range(ns):
      write('./dat/d%d.dat'%isou,rco[isou].getData())
    timer.stop('writing observed data')
  else: # ...else read from disk.
    timer.start('reading observed data')
    for isou in range(ns):
      read('./dat/d%d.dat'%isou,rco[isou].getData())
    timer.stop('reading observed data')

  # Born modeling and solver.
  born = BornOperatorS(e,dx,dt,nabsorb,u4,a4)
  ref = RecursiveExponentialFilter(1.0/(fpeak*dx*sqrt(2.0)))
  bornsolver = BornSolver(born,src,rcp,rco,ref,m) # solver

  # Wave operator
  wave = WaveOperator(e,dx,dt,nabsorb)

  # Warping
  maxShift = 0.1 # max shift in km
  #strain1,strain2,strain3 = 0.1,0.1,-1.0
  #strain1,strain2,strain3 = 0.1,0.1,0.5
  strain1,strain2,strain3 = 0.1,0.1,1.0
  #strain1,strain2,strain3 = 0.2,0.2,-1.0
  #strain1,strain2,strain3 = 0.2,0.2,0.5
  #strain1,strain2,strain3 = 0.2,0.2,1.0
  #strain1,strain2,strain3 = 0.4,0.4,1.0
  #strain1,strain2,strain3 = 0.4,0.4,1.0
  #strain1,strain2,strain3 = 0.5,0.5,-1.0
  #strain1,strain2,strain3 = 0.5,0.5,0.5
  #strain1,strain2,strain3 = 0.5,0.5,1.0
  #strain1,strain2,strain3 = 1.0,1.0,1.0
  #strain1,strain2,strain3 = 1.0,1.0,-1.0
  smooth1,smooth2,smooth3 = 4.0,4.0,4.0
  warping = ImageWarping(
    strain1,strain2,strain3,smooth1,smooth2,smooth3,maxShift,dz)
  print 'bstrain1=%d'%int(ceil(1.0/strain1))
  print 'bstrain2=%d'%int(ceil(1.0/strain2))
  print 'bstrain3=%d'%int(ceil(1.0/strain3))

  return e,m,wave,born,bornsolver,src,rco,u4,a4,b4,warping

#############################################################################

def goInversionX():
  modelChanged = False # if true, resimulate and remigrated observed data
  doLineSearch = True
  ninv = 2 # inversion iterations
  nmig = 5 # migration iterations within each inversion iteration
  nlsr = 2 # migration iterations within each inversion iteration
  #offset,stride = 10,5 # shot offset for warping, stride for gradient
  #offset,stride = 20,5 # shot offset for warping, stride for gradient
  #offset,stride = 20,10 # shot offset for warping, stride for gradient
  #offset,stride = 30,5 # shot offset for warping, stride for gradient
  offset,stride = 40,5 # shot offset for warping, stride for gradient
  
  s,m,wave,born,bs,src,rco,u4,a4,b4,warping = getInputs(modelChanged)
  m = transpose(m) # mask (transposed)
  r = zerofloat(nx,nz,ns) # reflectivity images
  rr = zerofloat(nz,nx,ns) # reflectivity images (transposed)
  rn = zerofloat(nz,nx,ns) # normalized images (transposed)
  rm = zerofloat(nz,nx,ns) # shifted images (transposed)
  rp = zerofloat(nz,nx,ns) # shifted images (transposed)
  rmw = zerofloat(nz,nx,ns) # warped images (transposed)
  rpw = zerofloat(nz,nx,ns) # warped images (transposed)
  um = zerofloat(nz,nx,ns) # shifts (transposed)
  up = zerofloat(nz,nx,ns) # shifts (transposed)
  sm = zerofloat(nz,nx,ns) # weights for shifts (transposed)
  sp = zerofloat(nz,nx,ns) # weights for shifts (transposed)
  wm = zerofloat(nz,nx,ns) # weighted shifts (transposed)
  wp = zerofloat(nz,nx,ns) # weighted shifts (transposed)
  ru = zerofloat(nx,nz,ns) # adjoint source images
  gm,pm = None,None # gradient & conjugate-gradient directions

  for i in range(ninv):
    print ''
    timer.start('ITERATION %d'%i)

    # Migration
    if modelChanged or i>0:
      timer.start('migration')
      bs.solve(nmig,r)
      timer.stop('migration')
      if i==0:
        write('./dat/r.dat',r)
    else:
      read('./dat/r.dat',r)
    
    # Transpose images
    transposeImages(r,rr)

    # Model preconditioner TODO: here or below? (shouldn't matter)
    #mp = getModelPrecondition(rr) # preconditioner
    #applyModelPrecondition(mp,rr)
    #pixels(mp,title='mp')

    # Find warping shifts
    findShifts(offset,warping,m,rr,rm,rp,um,up)
  
    useNewScale = False
    useOldScale = True
    if useNewScale:

      # New denominator for normalizing image amplitudes
      makeScaleN(m,rr,sm,sp)

      # Apply warping shifts (QC only)
      normalizeImages(m,rr,rn)
      offsetImages(offset,rn,rm,rp)
      taperImages(offset,rn,rn)
      taperImages(offset,rm,rm)
      taperImages(offset,rp,rp)
      applyShifts(warping,um,up,rm,rp,rmw,rpw)

    elif useOldScale:

      # Model preconditioner TODO: here or above? (shouldn't matter)
      #mp = getModelPrecondition(rr) # preconditioner
      #applyModelPrecondition(mp,rr)
      #pixels(mp,title='mp')

      # Apply warping shifts
      offsetImages(offset,rr,rm,rp) # use original images
      applyShifts(warping,um,up,rm,rp,rmw,rpw)

      # Denominator
      makeScale(offset,m,rr,rmw,rpw,sm,sp)

      # Normalize images (QC only)
      normalizeImages(m,rr,rn)
      taperImages(offset,rn,rn)
      taperImages(offset,rm,rm)
      taperImages(offset,rp,rp)

    else:

      # Constant denominator
      sm = sp = like3(um); fill(1.0,sm)

      # Apply warping shifts (QC only)
      normalizeImages(m,rr,rn)
      offsetImages(offset,rn,rm,rp)
      taperImages(offset,rn,rn)
      taperImages(offset,rm,rm)
      taperImages(offset,rp,rp)
      applyShifts(warping,um,up,rm,rp,rmw,rpw)

    # Scale warping shifts
    scaleShifts(um,up,sm,sp,wm,wp)

    # Taper images and scaled shifts XXX
    taperImages(offset,rr,rr)
    taperImages(offset,wm,wm)
    taperImages(offset,wp,wp)

    # Adjoint sources tapered
    makeAdjointSources(wm,wp,rr,ru)
    #makeAdjointSourcesX(wm,wp,rmw,rpw,ru)
    #taperAdjointSources(ru)

    # Gradient and conjugate gradient
    g = computeGradient(offset,stride,ru,wave,src,rco,u4,a4,b4)
    #RecursiveExponentialFilter(1.0/(fpeak*dx*sqrt(2.0))).apply(g,g) # XXX
    RecursiveGaussianFilter(1.0/(fpeak*dx*sqrt(2.0))).apply00(g,g) # XXX
    mul(transpose(m),g,g) # mask XXX
    p = getConjugateDirection(g,gm,pm)
    print 'sum(g)=%f'%sum(g)
    gm,pm = g,p # rotate

    # Line search
    if doLineSearch:
      timer.start('line search')
      aopt = findStepLength(s,m,p,nlsr,offset,rco,u4,a4,warping)
      timer.stop('line search')
    else:
      aopt = 0.0

    # Update slowness model
    add(s,mul(aopt,p),s)
    born.setSlowness(s)
    wave.setSlowness(s)

    if ninv>1:
      pixels(rr[ns/2],sperc=99.9,title='rr_%d'%i)
      pixels(getStackedImage(rr),title='r_%d'%i)
      pixels(g,cmap=rwb,sperc=99.9,title='g_%d'%i)
      pixels(p,cmap=rwb,sperc=99.9,title='p_%d'%i)
      pixels(s,cmap=jet,title='s_%d'%i)
    else:
      plots(i,rr,rn,rm,rp,rmw,rpw,um,up,sm,sp,wm,wp,ru,g,p,s)
    timer.stop('ITERATION %d'%i)

def plots(iiter,rr,rn,rm,rp,rmw,rpw,um,up,sm,sp,wm,wp,ru,g,p,s):
  #ks = 20 # shot to plot
  ks = nx/2 # shot to plot
  pixels(rr[ks],sperc=99.9,title='rr_%d'%iiter)
  pixels(rn[ks],sperc=99.9,title='rn_%d'%iiter)
  pixels(rm[ks],sperc=99.9,title='rm_%d'%iiter)
  pixels(rp[ks],sperc=99.9,title='rp_%d'%iiter)
  pixels(rmw[ks],sperc=99.9,title='rmw_%d'%iiter)
  pixels(rpw[ks],sperc=99.9,title='rpw_%d'%iiter)
  #pixels(um[ks],cmap=rwb,sperc=99.9,title='um_%d'%iiter)
  pixels(up[ks],cmap=rwb,sperc=99.9,title='up_%d'%iiter)
  pixels(add(um[ks],up[ks]),cmap=rwb,sperc=99.9,title='um+up_%i'%iiter)
  #pixels(sm[ks],cmap=rwb,sperc=99.9,title='sm_%d'%iiter)
  #pixels(sp[ks],cmap=rwb,sperc=99.9,title='sp_%d'%iiter)
  pixels(wm[ks],cmap=rwb,sperc=99.9,title='wm_%d'%iiter)
  pixels(wp[ks],cmap=rwb,sperc=99.9,title='wp_%d'%iiter)
  #pixels(add(wm[ks],wp[ks]),cmap=rwb,sperc=99.9,title='wm+wp_%d'%iiter)
  pixels(ru[ks],sperc=99.9,title='ru_%d'%iiter)
  pixels(g,cmap=rwb,sperc=99.9,title='g_%d'%iiter)
  #pixels(p,cmap=rwb,sperc=99.9,title='p_%d'%iiter)
  pixels(s,cmap=jet,title='s_%d'%iiter)

def transposeImages(r,rr):
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      transpose12(r[isou],rr[isou])
  Parallel.loop(ns,Loop())

def offsetImages(offset,rr,rm,rp):
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      if isou>=offset:
        copy(rr[isou-offset],rm[isou])
      if isou<=ns-1-offset:
        copy(rr[isou+offset],rp[isou])
  Parallel.loop(ns,Loop())

def applyModelPrecondition(mp,rr):
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      mul(mp,rr[isou],rr[isou])
  Parallel.loop(ns,Loop())

addRandomNoise = False
def findShifts(offset,warping,m,rr,rm,rp,um,up):
  rn = like3(rr)
  normalizeImages(m,rr,rn) # normalize images before warping
  offsetImages(offset,rn,rm,rp) # offset normalized images
  taperImages(offset,rn,rn) # taper images
  taperImages(offset,rm,rm) # taper images
  taperImages(offset,rp,rp) # taper images
  #pixels(rn[ns/2],sperc=99.9,title='rn')
  #pixels(rm[ns/2],sperc=99.9,title='rm')
  #pixels(rp[ns/2],sperc=99.9,title='rp')
  timer.start('warping')
  if addRandomNoise:
    warping.filterAndFindShifts(copy(rn),copy(rm),um)
    warping.filterAndFindShifts(copy(rn),copy(rp),up)
  else:
    warping.findShifts(rn,rm,um)
    warping.findShifts(rn,rp,up)
    #warping.findShifts(rm,rn,um); mul(-1.0,um,um)
    #warping.findShifts(rp,rn,up); mul(-1.0,up,up)
  taperImages(offset,um,um) # taper shifts
  taperImages(offset,up,up) # taper shifts
  timer.stop('warping')

def applyShifts(warping,um,up,rm,rp,rmw,rpw):
  warping.applyShifts(um,rm,rmw) # shift image
  warping.applyShifts(up,rp,rpw) # shift image

def makeScale(offset,m,rr,rmw,rpw,dm,dp):
  niter = 4 # number of times to apply ref
  #sigma = 0.064 # sigma in km
  #sigma = 0.125 # sigma in km
  sigma = 0.250 # sigma in km
  #sigma = 0.375 # sigma in km
  refSmooth = RecursiveExponentialFilter(sigma/(dx*sqrt(niter)))
  refSmooth.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  #rgfSmooth = RecursiveGaussianFilter(0.125/dx)
  #rgfSmooth = RecursiveGaussianFilter(0.250/dx)
  #rgfSmooth = RecursiveGaussianFilter(0.500/dx)
  rgfDeriv = RecursiveGaussianFilter(1.0)
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      dma,dmb = copy(rr[isou]),copy(rr[isou]) # shifted image
      dpa,dpb = copy(rr[isou]),copy(rr[isou]) # shifted image
      #dma,dmb = copy(rmw[isou]),copy(rmw[isou]) # shifted image
      #dpa,dpb = copy(rpw[isou]),copy(rpw[isou]) # shifted image
      rgfDeriv.apply1X(dma,dma) # 1st derivative of shifted image
      rgfDeriv.apply1X(dpa,dpa) # 1st derivative of shifted image
      rgfDeriv.apply2X(dmb,dmb) # 2nd derivative of shifted image
      rgfDeriv.apply2X(dpb,dpb) # 2nd derivative of shifted image
      mul(dma,dma,dma) # 1st derivative squared, first term in denom
      mul(dpa,dpa,dpa) # 1st derivative squared, first term in denom
      mul(sub(rr[isou],rmw[isou]),dmb,dmb) # second term in denom
      mul(sub(rr[isou],rpw[isou]),dpb,dpb) # second term in denom
      #mul(sub(rmw[isou],rr[isou]),dmb,dmb) # second term in denom
      #mul(sub(rpw[isou],rr[isou]),dpb,dpb) # second term in denom

      #if isou==ns/2:
      #  clip = 0.95*max(max(abs(dpa)),max(abs(dpb)))
      #  pixels(dpa,cmap=rwb,cmin=-clip,cmax=clip,title='dpa')
      #  pixels(dpb,cmap=rwb,cmin=-clip,cmax=clip,title='dpb')
      #  pixels(sub(dpa,dpb),cmap=rwb,cmin=-clip,cmax=clip,title='dpa-dpb')
      #  pixels(add(dpa,dpb),cmap=rwb,cmin=-clip,cmax=clip,title='dpa+dpb')

      #sub(dma,dmb,dm[isou]) # denominator
      #sub(dpa,dpb,dp[isou]) # denominator
      copy(dma,dm[isou]) # denominator ignoring second term
      copy(dpa,dp[isou]) # denominator ignoring second term
      for i in range(niter):
        refSmooth.apply(dm[isou],dm[isou]) # smooth
        refSmooth.apply(dp[isou],dp[isou]) # smooth
      #rgfSmooth.apply00(dm[isou],dm[isou]) # smooth
      #rgfSmooth.apply00(dp[isou],dp[isou]) # smooth

      #if isou==ns/2:
      #  pixels(dm[isou],cmap=rwb,sperc=100.0,title='dm_smooth')
      #  pixels(dp[isou],cmap=rwb,sperc=100.0,title='dp_smooth')

  Parallel.loop(offset,ns-offset,Loop())
  mul(1.0/max(dm),dm,dm) # normalize
  mul(1.0/max(dp),dp,dp) # normalize
  add(0.01,dm,dm) # stabilize for division
  add(0.01,dp,dp) # stabilize for division
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      div(m,dm[isou],dm[isou]) # reciprocal
      div(m,dp[isou],dp[isou]) # reciprocal
  Parallel.loop(ns,Loop())

def makeDenominatorX(offset,rr,dd):
  """Simple denominator for amplitude normalization."""
  rgfDeriv = RecursiveGaussianFilter(1.0)
  rgfSmooth = RecursiveGaussianFilter(0.125/dx)
  #rgfSmooth = RecursiveGaussianFilter(0.250/dx)
  #rgfSmooth = RecursiveGaussianFilter(0.500/dx)
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      di = copy(rr[isou])
      rgfDeriv.apply1X(di,di) # 1st derivative
      abs(di,di) # 1st derivative absolute value
      #mul(di,di,di) # 1st derivative squared
      rgfSmooth.apply00(di,dd[isou]) # denominator
      #mul(1.0/max(dd[isou]),dd[isou],dd[isou]) # normalize TODO: outside?
      #add(0.01,dd[isou],dd[isou]) # stabilize for division
  Parallel.loop(offset,ns-offset,Loop())
  mul(1.0/max(dd),dd,dd) # normalize TODO: inside or outside loop?
  add(0.01,dd,dd) # stabilize for division

def normalizeImages(m,rr,ry):
  ss = like3(rr)
  makeScaleN(m,rr,ss,sp=None,useAbs=True)
  mul(ss,rr,ry)

def makeScaleN(m,rr,sm,sp=None,useAbs=False):
  niter = 4 # number of times to apply ref
  #sigma = 0.125 # sigma in km
  sigma = 0.250 # sigma in km
  nsou = len(rr)
  ref = RecursiveExponentialFilter(sigma/(dx*sqrt(niter)))
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      smi = sm[isou]
      if useAbs:
        abs(rr[isou],smi)
      else:
        mul(rr[isou],rr[isou],smi)
      for i in range(niter):
        ref.apply(smi,smi)
  Parallel.loop(nsou,Loop())
  mul(1.0/max(sm),sm,sm)
  add(0.01,sm,sm)
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      div(m,sm[isou],sm[isou])
  Parallel.loop(nsou,Loop())
  if sp is not None:
    copy(sm,sp)

def taperImages(offset,rr,ry=None):
  if ry is None:
    ry = like3(rr)
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      taperImage(isou,rr[isou],ry[isou])
      #w = zerofloat(nx)
      #fx,lx = max(isou-hx,0),min(isou+hx,nx-1)
      #for ix in range(fx,lx+1):
      #  w[ix] = 1.0
      #rgf.apply0(w,w)
      ##if isou==ns/2:
      ##  points(w)
      #for ix in range(nx):
      #  mul(w[ix],rr[isou][ix],ry[isou][ix])
  Parallel.loop(ns,Loop())

def taperImage(ix,rx,ry,imageIsTransposed=True):
  #hx = 10 # half taper width
  hx = 20 # half taper width
  #hx = 40 # half taper width
  #hx = 50 # half taper width
  #hx = 100 # half taper width
  #hx = 3*offset # half taper width
  #hx = 4*offset # half taper width
  #hx = 6*offset # half taper width
  #hx = nx # no taper
  #sigma = 0.125/dx # half width for taper transition
  sigma = 8 # half width for taper transition
  w = zerofloat(nx)
  fx,lx = max(ix-hx,0),min(ix+hx,nx-1)
  for ix in range(fx,lx+1):
    w[ix] = 1.0
  RecursiveGaussianFilter(sigma).apply0(w,w)
  if imageIsTransposed:
    for ix in range(nx):
      mul(w[ix],rx[ix],ry[ix])
  else:
    for iz in range(nz):
      mul(w,rx[iz],ry[iz])

def scaleShifts(um,up,sm,sp,wm,wp):
  mul(sm,um,wm)
  mul(sp,up,wp)

def makeAdjointSources(wm,wp,rr,ru):
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      for ix in range(nx):
        for iz in range(1,nz-1):
          ru[isou][iz][ix] =\
            (wm[isou][ix][iz]+wp[isou][ix][iz])*\
            (rr[isou][ix][iz+1]-rr[isou][ix][iz-1])
  Parallel.loop(ns,Loop())

def makeAdjointSourcesX(wm,wp,rmw,rpw,ru):
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      for ix in range(nx):
        for iz in range(1,nz-1):
          ru[isou][iz][ix] =\
            wm[isou][ix][iz]*(rmw[isou][ix][iz+1]-rmw[isou][ix][iz-1])+\
            wp[isou][ix][iz]*(rpw[isou][ix][iz+1]-rpw[isou][ix][iz-1])
  Parallel.loop(ns,Loop())

def taperAdjointSources(ru):
  sigma = 1.0/dx # half width for taper
  #sigma = 0.5/dx # half width for taper
  gt = GaussianTaper(sigma)
  points(gt.getWeightsForLength(nx))
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      gt.apply1(ru[isou],ru[isou])
  Parallel.loop(ns,Loop())

def computeGradient(offset,stride,ru,wave,src,rco,u4,a4,b4):
  class Reduce(Parallel.ReduceInt):
    def compute(self,isou):
      ksou = (isou-offset)/stride
      rui,ui,ai,bi = ru[isou],u4.get(ksou),a4.get(ksou),b4.get(ksou)
      wave.applyForward(src[isou],ui)
      #gp = wave.collapse(ui,ui,nabsorb) # gradient preconditioner
      #mul(1.0/max(gp),gp,gp)
      wave.applyAdjoint(Source.ReceiverSource(rco[isou]),ai)
      # TODO: second time derivative?
      wave.applyAdjoint(Source.WavefieldSource(ai,rui),bi)
      ga = wave.collapse(ui,bi,nabsorb) # source side
      wave.applyForward(Source.WavefieldSource(ui,rui),bi)
      gb = wave.collapse(bi,ai,nabsorb) # receiver side
      add(ga,gb,gb)
      #taperImage(isou,gb,gb,imageIsTransposed=False)
      if isou==nx/2:
        pixels(gb,cmap=rwb,sperc=100.0,title='gb')
      return gb
    def combine(self,ga,gb):
      return add(ga,gb)
  timer.start('gradient')
  #g = PartialParallel(np).reduce(offset,ns-offset,stride,Reduce()) 
  g = PartialParallel(np).reduce(2*offset,ns-2*offset,stride,Reduce()) 
  #g = PartialParallel(np).reduce(2*ns/5,3*ns/5,stride,Reduce()) 
  #g = PartialParallel(np).reduce(ns/3,2*ns/3,stride,Reduce()) 
  #g = PartialParallel(np).reduce(offset,1+offset,stride,Reduce()) # one shot
  #g = PartialParallel(np).reduce(nx/2,1+nx/2,stride,Reduce()) # one shot
  timer.stop('gradient')
  return g

def newstep():
  for isou in range(nsou):
    ksou = (1+isou)*ns/(1+nsou)

def findStepLength(s,m,p,nmig,offset,reco,u4,a4,warping):
  """Line search using only one shot."""
  amin = -0.05 # minimum step length
  #amax = 0.00 # maximum step length FIXME
  #amax = 0.01 # maximum step length
  amax = 0.05 # maximum step length
  atol = 0.01 # tolerance
  #atol = 0.005 # tolerance
  #atol = 0.001 # tolerance
  #atol = 0.0001 # tolerance
  hx = nx/2 # shot location for line search
  xss = [hx-offset,hx,hx+offset]
  zss = [0,0,0]
  src = zeros(3,Source)
  rco = zeros(3,Receiver)
  rcp = zeros(3,Receiver)
  r = zerofloat(nx,nz,3)
  for isou in range(3):
    ri = reco[xss[isou]]
    rco[isou] = Receiver(ri)
    rcp[isou] = Receiver(ri.getXIndices(),ri.getZIndices(),nt)
    src[isou] = Source.RickerSource(xss[isou],zss[isou],dt,fpeak)
  born = BornOperatorS(s,dx,dt,nabsorb,u4,a4)
  ref = RecursiveExponentialFilter(1.0/(fpeak*dx*sqrt(2.0)))
  bsolver = BornSolver(born,src,rcp,rco,ref,None)
  class Misfit(BrentMinFinder.Function):
    def evaluate(self,a):
      e = mul(a/max(abs(p)),p)
      add(s,e,e)
      born.setSlowness(e)
      bsolver.solve(nmig,r)
      normalizeImages(transpose(m),r,r)
      rm = transpose12(r[0])
      rr = transpose12(r[1])
      rp = transpose12(r[2])
      if m is not None:
        mul(m,rm,rm)
        mul(m,rr,rr)
        mul(m,rp,rp)
      taperImage(isou,rr,rr,imageIsTransposed=True)
      taperImage(isou,rm,rm,imageIsTransposed=True)
      taperImage(isou,rp,rp,imageIsTransposed=True)
      um = warping.findShifts(rr,rm)
      up = warping.findShifts(rr,rp)
      add(abs(um),abs(up),up)
      #add(um,up,up)
      if m is not None:
        mul(m,up,up)
        pass
      #pixels(rr,title='_rr_%f'%a)
      #pixels(up,cmap=rwb,sperc=100.0,title='_up_%f'%a)
      return sum(mul(up,up))
  aopt = BrentMinFinder(Misfit()).findMin(amin,amax,atol)
  print 'aopt=%f'%aopt
  return aopt/max(abs(p))

def xfindStepLength(s,p,m,nmig,offset,reco,u4,a4,warping):
  nsou = 21 # number of shots to migrate
  amin = -0.01 # minimum step length
  amax = 0.04 # maximum step length
  atol = 0.20*(amax-amin) # tolerance
  msou = nsou+2*offset
  r = zerofloat(nx,nz,msou)
  rr = zerofloat(nz,nx,nsou)
  rm = zerofloat(nz,nx,nsou)
  rp = zerofloat(nz,nx,nsou)
  ur = zerofloat(nz,nx,nsou)
  um = zerofloat(nz,nx,nsou)
  up = zerofloat(nz,nx,nsou)
  fxs = nx/2-nsou/2-offset
  xss = rampint(fxs,1,msou)
  zss = fillint(0,msou)
  src = zeros(msou,Source)
  rco = zeros(msou,Receiver)
  rcp = zeros(msou,Receiver)
  born = BornOperatorS(s,dx,dt,nabsorb,u4,a4)
  ref = RecursiveExponentialFilter(1.0/(fpeak*dx*sqrt(2.0)))
  bsolver = BornSolver(born,src,rcp,rco,ref,None) # solver
  for isou in range(msou):
    rco[isou] = Receiver(reco[fxs+isou])
    rcp[isou] = Receiver(xr,zr,nt)
    src[isou] = Source.RickerSource(xss[isou],zss[isou],dt,fpeak)
  class Misfit(BrentMinFinder.Function):
    def evaluate(self,a):
      e = mul(a/max(abs(p)),p)
      add(s,e,e)
      born.setSlowness(e)
      bsolver.solve(nmig,r)
      class Loop(Parallel.LoopInt):
        def compute(self,isou):
          transpose12(r[isou       ],rr[isou-offset])
          transpose12(r[isou-offset],rm[isou-offset])
          transpose12(r[isou+offset],rp[isou-offset])
      Parallel.loop(offset,msou-offset,Loop())
      class Loop(Parallel.LoopInt):
        def compute(self,isou):
          mul(m,rr[isou],rr[isou])
          mul(m,rm[isou],rm[isou])
          mul(m,rp[isou],rp[isou])
      Parallel.loop(nsou,Loop())
      warping.findShifts(rr,rm,um)
      warping.findShifts(rr,rp,up)
      add(um,up,up)
      class Loop(Parallel.LoopInt):
        def compute(self,isou):
          mul(m,up[isou],up[isou])
      Parallel.loop(nsou,Loop())
      return sum(mul(up,up))
  aopt = BrentMinFinder(Misfit()).findMin(amin,amax,atol)
  print 'aopt=%f'%aopt
  return aopt/max(abs(p))

def goMigration(niter=5):
  _,_,_,_,bs,src,rco,_,_,_,_ = getInputs()
  r = zerofloat(nx,nz,ns) # reflectivity image
  bs.solve(niter,r)
  if datDir:
    write(datDir+'r.dat',r)

def goGradient():
  doWarping = True
  doGradient = True
  ks = 20 # shot stride
  t,s,e,r,m = getModelAndMask()
  ms = nx # number of sources

  # Read images
  rt = zerofloat(nx,nz,ms) # reflectivity images
  #read('./dat/r100.dat',rt)
  read('./dat/r090.dat',rt)
  rr = zerofloat(nz,nx,ms) # reflectivity images
  rm = zerofloat(nz,nx,ms) # reflectivity images shifted by ks shots
  rp = zerofloat(nz,nx,ms) # reflectivity images shifted by ks shots
  um = zerofloat(nz,nx,ms) # shifts
  up = zerofloat(nz,nx,ms) # shifts
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      transpose12(rt[isou   ],rr[isou])
      transpose12(rt[isou-ks],rm[isou])
      transpose12(rt[isou+ks],rp[isou])
  Parallel.loop(ks,ms-ks,Loop())

  # Preconditioning
  mp = getModelPrecondition(rr)
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      mul(mp,rr[isou],rr[isou])
      mul(mp,rm[isou],rm[isou])
      mul(mp,rp[isou],rp[isou])
  Parallel.loop(ms,Loop())

  # Warping
  if doWarping:
    maxShift = 0.1 # max shift (km)
    strain1,strain2,strain3 = 0.5,0.5,1.0
    smooth1,smooth2,smooth3 = 4.0,4.0,4.0
    warping = ImageWarping(
      strain1,strain2,strain3,smooth1,smooth2,smooth3,maxShift,dz)
    timer.start('warping')
    warping.findShifts(rr,rm,um) # TODO: warp which image?
    warping.findShifts(rr,rp,up) # TODO: warp which image?
    #warping.filterAndFindShifts(rr,rm,um) # add noise
    #warping.filterAndFindShifts(rr,rp,up) # add noise
    timer.stop('warping')
  uu = add(um,up)
  ru = zerofloat(nx,nz,ms) # derivative of scaled images
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      for ix in range(nx):
        for iz in range(1,nz-1):
          ru[isou][iz][ix] =\
            0.5*uu[isou][ix][iz]*(rr[isou][ix][iz+1]-rr[isou][ix][iz-1])/dz
  Parallel.loop(ms,Loop())
  
  # Gradient
  if doGradient:
    timer.start('gradient')
    g = getGradientFromImages(ru)
    timer.start('gradient')
    pixels(g,cmap=rwb,sperc=100.0,title='g')
    if datDir:
      write(datDir+'g.dat',g)

  rw = warping.applyShifts(um,rm)
  pixels(um[ms/2],cmap=rwb,sperc=100.0,title='um')
  pixels(rm[ms/2],sperc=100.0,title='rm')
  pixels(rw[ms/2],sperc=100.0,title='rw')
  pixels(rr[ms/2],sperc=100.0,title='rr')
  pixels(ru[ms/2],sperc=100.0,title='ru')

#  display(rr,title='rr')
#  print max(rr)
#  print min(rr)
#  display(ru,title='ru')
#  display(uu,cmap=rwb,sperc=99.9,title='uu')

def getGradientFromImages(ru):
  _,_,wave,_,_,src,rco,u,a,_,_ = getInputs()
  b = SharedFloat4(nxp,nzp,nt,np)
  class Reduce(Parallel.ReduceInt):
    def compute(self,isou):
      ui,ai,bi = u.get(isou),a.get(isou),b.get(isou)
      rui = ru[Sampling(nx).indexOfNearest(xs[isou])]
      wave.applyForward(src[isou],ui)
      wave.applyAdjoint(Source.ReceiverSource(rco[isou]),ai)
      # TODO: second time derivative?

      # First half
      wave.applyAdjoint(Source.WavefieldSource(ai,rui),bi)
      ga = wave.collapse(ui,bi,nabsorb)

      # Second half
      wave.applyForward(Source.WavefieldSource(ui,rui),bi)
      gb = wave.collapse(bi,ai,nabsorb)

      #pixels(ga,cmap=rwb,sperc=100.0,title='ga')
      #pixels(gb,cmap=rwb,sperc=100.0,title='gb')
      return add(ga,gb)
    def combine(self,ga,gb):
      return add(ga,gb)
  return mul(1.0/ns,PartialParallel(np).reduce(ns,Reduce()))

def getModelPrecondition(rr):
  r = getStackedImage(rr)
  return getEnergy(r)

def getEnergy(r):
  e = abs(r)
  #e = mul(r,r)
  efilter(8.0,e,e)
  mul(1.0/max(e),e,e)
  return e

def efilter(sigma,f,g,niter=2):
  ref = RecursiveExponentialFilter(sigma/sqrt(niter))
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref.apply(f,g)
  for i in range(1,niter):
    ref.apply(g,g)

def getStackedImage(rr):
  n = len(rr)
  r = like2(rr[0])
  for i in range(n):
    add(rr[i],r,r)
  return r

def transpose12(f,g=None):
  n1,n2 = len(f[0]),len(f)
  if g is None:
    g = zerofloat(n2,n1)
  for i2 in range(n2):
    for i1 in range(n1):
      g[i1][i2] = f[i2][i1]
  return g

def getConjugateDirection(g,gm=None,pm=None):
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

#############################################################################

def makeBornModel(s):
  """
  sigma0: smoothing for background model
  sigma1: smoothing for perturbation
  """
  print 'sigma0=%f'%sigma0
  print 'sigma1=%f'%sigma1
  s0,s1 = like2(s),like2(s)
  ref0 = RecursiveExponentialFilter(sigma0/dx)
  ref0.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref0.apply(s,s0)
  t = copy(s)
  ref1 = RecursiveExponentialFilter(sigma1/dx)
  ref1.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref1.apply(s,t)
  sub(s,t,s1)
  r = sub(mul(s,s),mul(t,t))
  #if nx!=767:
  #  GaussianTaper.apply1(0.25,r,r)
  return s0,r

def like2(x):
  return zerofloat(len(x[0]),len(x))

def like3(x):
  return zerofloat(len(x[0][0]),len(x[0]),len(x))

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

def report(str,stopwatch):
  s = stopwatch.time()
  h = int(s/3600); s -= h*3600
  m = int(s/60); s -= m*60
  if h>0:
    print str+': %02d:%02d:%02d'%(h,m,s)
  else:
    print str+': %02d:%02d'%(m,s)

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
  if not plots:
    return
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if nx==189:
    sp.setSize(1000,650)
  elif nx==251:
    sp.setSize(1265,650)
  else:
    sp.setSize(1000,650)
  #sp.setFontSizeForSlide(1.0,1.0)
  cb = sp.addColorBar()
  cb.setWidthMinimum(90)
  #cb.setWidthMinimum(150)
  if title:
    sp.addTitle(title)
    pass
  if len(x)==nz:
    #pv = sp.addPixels(sz,sx,transpose(x))
    #sp.setHLabel("Distance (km)")
    #sp.setVLabel("Depth (km)")
    pv = sp.addPixels(transpose(x))
  elif len(x[0])==nt:
    #pv = sp.addPixels(st,sx,x)
    #sp.setHLabel("Distance (km)")
    #sp.setVLabel("Time (s)")
    pv = sp.addPixels(x)
  else:
    pv = sp.addPixels(x)
  pv.setColorModel(cmap)
  if perc<100.0:
    pv.setPercentiles(100.0-perc,perc)
  if sperc is not None: # symmetric percentile clip
    clips = Clips(100.0-sperc,sperc,x)
    clip = max(abs(clips.getClipMin()),abs(clips.getClipMax()))
    pv.setClips(-clip,clip)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if title and savDir:
    #write(savDir+title+'.dat',x)
    #sp.paintToPng(720,3.33,savDir+title+'.png')
    sp.paintToPng(360,3.33,savDir+title+'.png')

def points(x):
  if not plots:
    return
  SimplePlot.asPoints(x)

def display(x,cmap=gray,cbar=None,sperc=None,title=None):
  if not plots:
    return
  frame = SimpleFrame()
  ipg = frame.addImagePanels(x)
  ipg.setColorModel(cmap)
  colorBar = ColorBar(cbar) if cbar else ColorBar()
  colorBar.setWidthMinimum(70)
  ipg.addColorMapListener(colorBar)
  colorBar.setFont(colorBar.getFont().deriveFont(12.0))
  frame.add(colorBar,BorderLayout.EAST)
  if sperc is not None:
    clips = Clips(100-sperc,sperc,x)
    clip = max(abs(clips.getClipMin()),abs(clips.getClipMax()))
    ipg.setClips(-clip,clip)
  elif min(x)<max(x):
    ipg.setClips(min(x),max(x)) # fix colorbar
  if title:
    frame.setTitle(title)

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
