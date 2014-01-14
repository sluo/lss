#############################################################################
# Image-warping tomography

from imports import *

#############################################################################

savDir = None
#savDir = os.getenv('HOME')+'/Desktop/pngdat/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat2/'
savDir = '/Users/sluo/Dropbox/pngdat/'

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
  goInversion()

def setupForLayered():
  global sz,sx,st,nz,nx,nt,nxp,nzp,dz,dx,dt
  global zs,xs,zr,xr,ns,nr,fpeak,nabsorb,np
  global sigma0,sigma1
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
  #xs,zs = [0],[0]
  #xs,zs = [nx/2],[0]
  #xs,zs = [nx/2,1+nx/2],[0,0]
  #xs,zs = [nx/5,2*nx/5,3*nx/5,4*nx/5],[0,0,0,0]
  #xs,zs = rampint(2,5,38),fillint(0,38)
  xs,zs = rampint(0,1,nx),fillint(0,nx)
  xr,zr = rampint(0,1,nx),fillint(0,nx)
  ns,nr = len(xs),len(xr)
  fpeak = 20.0 # Ricker wavelet peak frequency
  nabsorb = 22 # absorbing boundary size
  nxp,nzp = nx+2*nabsorb,nz+2*nabsorb
  np = min(16,ns) # number of parallel sources
  sigma0 = 4.0/fpeak
  sigma1 = 1.0/fpeak
  print 'ns=%d'%ns
  #return getSingleLayeredModelAndMask()
  return getMultiLayeredModelAndMask()

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
  #GaussianTaper.apply1(t1,t1)
  if constantBackground:
    fill(tb,t0)
  m = fillfloat(1.0,nx,nz)
  for iz in range(nz/3-6):
    for ix in range(nx):
      m[iz][ix] = 0.0
  RecursiveExponentialFilter(1.0).apply2(m,m)
  #return t,t0,t1,m
  return t,t0,t1,None,None

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
  pixels(t,cmap=jet)
  return t,t0,t1,None,None

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

def addGaussianError(s,emax=0.010,m=None,sigma=None):
  g = zerofloat(nx,nz)
  g[nz/3][nx/2] = 1.0
  if sigma is None:
    RecursiveGaussianFilter(nz/8.0).apply00(g,g)
  else:
    RecursiveGaussianFilter(sigma).apply00(g,g)
  mul(emax/max(g),g,g)
  if m is not None:
    mul(m,g,g)
  return add(s,g)

def getModelAndMask():
  t,s,r,m,p = setupForLayered()
  e = copy(s)
  #s = mul(0.90,e) # erroneous background slowness
  #s = mul(1.10,e) # erroneous background slowness
  #s = addGaussianError(e,emax=0.05,m=m)
  #s = addGaussianError(e,emax=-0.05,m=m)
  s = addGaussianError(e,emax=0.05,m=m,sigma=nz/10.0)
  #s = addGaussianError(e,emax=-0.05,m=m,sigma=nz/10.0)
  #s = addGaussianError(e,emax=0.04,m=m)
  #s = addGaussianError(e,emax=-0.04,m=m)
  #s = addGaussianError(e,emax=-0.04,m=m)
  #pixels(t,cmap=jet,title='t')
  pixels(s,cmap=jet,title='s')
  pixels(e,cmap=jet,title='e')
  if m is not None:
    pixels(m,cmap=gray,title='m')
  pixels(sub(s,e),cmap=rwb,sperc=100.0,title='s-e')
  pixels(r,cmap=gray,sperc=100.0,title='r')
  return t,s,e,r,p

def xgetModelAndMask():
  t,s,r,m,p = setupForLayered()
  e = copy(s)
  e = mul(0.90,s) # erroneous background slowness
  #e = mul(1.10,s) # erroneous background slowness
  e = addGaussianError(s,emax=0.050,m=m)
  #pixels(t,cmap=jet,title='t')
  pixels(s,cmap=jet,title='s')
  pixels(e,cmap=jet,title='e')
  if m is not None:
    pixels(m,cmap=gray,title='m')
  pixels(sub(s,e),cmap=rwb,sperc=100.0,title='s-e')
  pixels(r,cmap=gray,sperc=100.0,title='r')
  return t,s,e,r,p

def getInputs():
  t,s,e,r,m = getModelAndMask()

  # Wavefields.
  timer.start('allocating')
  u4 = SharedFloat4(nxp,nzp,nt,np)
  a4 = SharedFloat4(nxp,nzp,nt,np)
  b4 = SharedFloat4(nxp,nzp,nt,np)
  timer.stop('allocating')

  # Sources and receivers.
  src = zeros(ns,Source)
  rco = zeros(ns,Receiver)
  rcp = zeros(ns,Receiver)
  for isou in range(ns):
    rco[isou] = Receiver(xr,zr,nt)
    rcp[isou] = Receiver(xr,zr,nt)
    src[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
  bornt = BornOperatorS(s,dx,dt,nabsorb,u4,a4) # true slowness
  timer.start('observed data')
  bornt.applyForward(src,r,rco) # observed data
  timer.stop('observed data')

  # Born modeling and solver.
  born = BornOperatorS(e,dx,dt,nabsorb,u4,a4)
  ref = RecursiveExponentialFilter(1.0/(fpeak*dx*sqrt(2.0)))
  bornsolver = BornSolver(born,src,rcp,rco,ref,m) # solver

  # Wave operator
  wave = WaveOperator(e,dx,dt,nabsorb)

  # Warping
  maxShift = 0.1 # max shift in km
  strain1,strain2,strain3 = 0.5,0.5,1.0
  smooth1,smooth2,smooth3 = 4.0,4.0,4.0
  warping = ImageWarping(
    strain1,strain2,strain3,smooth1,smooth2,smooth3,maxShift,dz)

  return e,wave,born,bornsolver,src,rco,u4,a4,b4,warping

def goInversion():
  s,wave,born,bs,src,rco,u4,a4,b4,warping = getInputs()
  doMigration = False
  doLineSearch = False
  ninv = 1 # inversion iterations
  nmig = 5 # migration iterations within each inversion iteration
  nlsr = 2 # migration iterations within each inversion iteration
  #offset,stride = 10,5 # shot offset for warping, stride for gradient
  offset,stride = 20,5 # shot offset for warping, stride for gradient
  #offset,stride = 30,5 # shot offset for warping, stride for gradient

  r = zerofloat(nx,nz,ns) # reflectivity images
  ru = zerofloat(nx,nz,ns) # adjoint source images
  rr = zerofloat(nz,nx,ns) # reflectivity images (transposed)
  rm = zerofloat(nz,nx,ns) # shifted images (transposed)
  rp = zerofloat(nz,nx,ns) # shifted images (transposed)
  um = zerofloat(nz,nx,ns) # shifts (transposed)
  up = zerofloat(nz,nx,ns) # shifts (transposed)
  dm = zerofloat(nz,nx,ns) # denominator for adjoint sources (transposed)
  dp = zerofloat(nz,nx,ns) # denominator for adjoint sources (transposed)
  gm,pm = None,None # gradient & conjugate-gradient directions

  for i in range(ninv):
    timer.start('ITERATION')

    # Migration
    timer.start('migration')
    if i>0 or doMigration:
      bs.solve(nmig,r)
      if i==0:
        write('./dat/r.dat',r)
    else:
      read('./dat/r.dat',r)
    timer.stop('migration')
    
    # Transpose and shift images, apply preconditioner
    class Loop(Parallel.LoopInt):
      def compute(self,isou):
        transpose12(r[isou       ],rr[isou])
        transpose12(r[isou-offset],rm[isou])
        transpose12(r[isou+offset],rp[isou])
    Parallel.loop(offset,ns-offset,Loop())
    mp = getModelPrecondition(rr) # preconditioner
    class Loop(Parallel.LoopInt):
      def compute(self,isou):
        mul(mp,rr[isou],rr[isou])
        mul(mp,rm[isou],rm[isou])
        mul(mp,rp[isou],rp[isou])
    Parallel.loop(offset,ns-offset,Loop())
    pixels(mp,title='mp_%d'%i)
    pixels(rm[ns/2],sperc=99.9,title='rm_%d'%i)
    pixels(rr[ns/2],sperc=99.9,title='rr_%d'%i)
    pixels(rp[ns/2],sperc=99.9,title='rp_%d'%i)

    # Warping
    timer.start('warping')
    warping.findShifts(rr,rm,um)
    warping.findShifts(rr,rp,up)
    timer.stop('warping')
    rm = warping.applyShifts(um,rm) # shift images
    rp = warping.applyShifts(up,rp) # shift images
    pixels(rm[ns/2],sperc=99.9,title='rmw_%d'%i)
    pixels(rp[ns/2],sperc=99.9,title='rpw_%d'%i)
    pixels(mul(mp,um[ns/2]),cmap=rwb,sperc=99.9,title='um_%d'%i)
    pixels(mul(mp,up[ns/2]),cmap=rwb,sperc=99.9,title='up_%d'%i)

    # Denominator for adjoint sources
    #rgfDeriv = RecursiveGaussianFilter(1.0)
    #rgfSmooth = RecursiveGaussianFilter(0.125/dx)
    #class Loop(Parallel.LoopInt):
    #  def compute(self,isou):
    #    dma,dmb = copy(rm[isou]),copy(rm[isou]) # shifted image
    #    dpa,dpb = copy(rp[isou]),copy(rp[isou]) # shifted image
    #    rgfDeriv.apply1X(dma,dma) # 1st derivative of shifted image
    #    rgfDeriv.apply1X(dpa,dpa) # 1st derivative of shifted image
    #    rgfDeriv.apply2X(dmb,dmb) # 2nd derivative of shifted image
    #    rgfDeriv.apply2X(dpb,dpb) # 2nd derivative of shifted image
    #    mul(dma,dma,dma) # 1st derivative squared, first term in denom
    #    mul(dpa,dpa,dpa) # 1st derivative squared, first term in denom
    #    mul(sub(rr[isou],rm[isou]),dmb,dmb) # second term in denom
    #    mul(sub(rr[isou],rp[isou]),dpb,dpb) # second term in denom
    #    sub(dma,dmb,dmb) # denominator
    #    sub(dpa,dpb,dpb) # denominator
    #    #copy(dma,dm[isou]) # denominator (ignoring second term)
    #    #copy(dpa,dp[isou]) # denominator (ignoring second term)
    #    rgfSmooth.apply00(dmb,dm[isou]) # smooth
    #    rgfSmooth.apply00(dpb,dp[isou]) # smooth
    #Parallel.loop(offset,ns-offset,Loop())
    #mul(1.0/max(abs(dm)),dm,dm) # normalize
    #mul(1.0/max(abs(dp)),dp,dp) # normalize
    #add(0.01,dm,dm) # stabilize for division
    #add(0.01,dp,dp) # stabilize for division
    ##mul(1.0/max(abs(dm)),dm,dm) # normalize
    ##mul(1.0/max(abs(dp)),dp,dp) # normalize
    #pixels(dm[ns/2],cmap=rwb,sperc=99.9,title='dm_%d'%i)
    #pixels(dp[ns/2],cmap=rwb,sperc=99.9,title='dp_%d'%i)
  
    # Denominator assuming center image is always used for adjoint source
    rgfDeriv = RecursiveGaussianFilter(1.0)
    rgfSmooth = RecursiveGaussianFilter(0.125/dx)
    class Loop(Parallel.LoopInt):
      def compute(self,isou):
        #dma,dmb = copy(rr[isou]),copy(rm[isou]) # shifted image
        #dpa,dpb = copy(rr[isou]),copy(rp[isou]) # shifted image
        dma,dmb = copy(rr[isou]),copy(rr[isou]) # shifted image
        dpa,dpb = copy(rr[isou]),copy(rr[isou]) # shifted image
        rgfDeriv.apply1X(dma,dma) # 1st derivative of shifted image
        rgfDeriv.apply1X(dpa,dpa) # 1st derivative of shifted image
        rgfDeriv.apply2X(dmb,dmb) # 2nd derivative of shifted image
        rgfDeriv.apply2X(dpb,dpb) # 2nd derivative of shifted image
        mul(dma,dma,dma) # 1st derivative squared, first term in denom
        mul(dpa,dpa,dpa) # 1st derivative squared, first term in denom
        mul(sub(rr[isou],rm[isou]),dmb,dmb) # second term in denom
        mul(sub(rr[isou],rp[isou]),dpb,dpb) # second term in denom
        #mul(sub(rm[isou],rr[isou]),dmb,dmb) # second term in denom
        #mul(sub(rp[isou],rr[isou]),dpb,dpb) # second term in denom
        sub(dma,dmb,dmb) # denominator
        sub(dpa,dpb,dpb) # denominator
        #copy(dma,dm[isou]) # denominator (ignoring second term)
        #copy(dpa,dp[isou]) # denominator (ignoring second term)
        rgfSmooth.apply00(dmb,dm[isou]) # smooth
        rgfSmooth.apply00(dpb,dp[isou]) # smooth
    Parallel.loop(offset,ns-offset,Loop())
    mul(1.0/max(abs(dm)),dm,dm) # normalize
    mul(1.0/max(abs(dp)),dp,dp) # normalize
    add(0.01,dm,dm) # stabilize for division
    add(0.01,dp,dp) # stabilize for division
    #mul(1.0/max(abs(dm)),dm,dm) # normalize
    #mul(1.0/max(abs(dp)),dp,dp) # normalize
    pixels(dm[ns/2],cmap=rwb,sperc=99.9,title='dm_%d'%i)
    pixels(dp[ns/2],cmap=rwb,sperc=99.9,title='dp_%d'%i)


    # Adjoint sources
    uu = add(div(um,dm),div(up,dp))
    #uu = add(um,up) # XXX
    class Loop(Parallel.LoopInt):
      def compute(self,isou):

        # Extra preconditioner
        #pp = copy(rr[isou])
        #abs(pp,pp)
        #efilter(8.0,pp,pp)
        #mul(1.0/max(pp),pp,pp)
        #add(0.01,pp,pp)
        #div(1.0,pp,pp)

        for ix in range(nx):
          for iz in range(1,nz-1):
            ru[isou][iz][ix] =\
              mp[ix][iz]*uu[isou][ix][iz]*\
              (rr[isou][ix][iz+1]-rr[isou][ix][iz-1])

    Parallel.loop(offset,ns-offset,Loop())
    pixels(mul(mp,uu[ns/2]),cmap=rwb,sperc=99.9,title='uu_%d'%i)
    pixels(ru[ns/2],sperc=99.9,title='ru_%d'%i)

    # Gradient and conjugate gradient
    class Reduce(Parallel.ReduceInt):
      def compute(self,isou):
        ksou = (isou-offset)/stride
        rui,ui,ai,bi = ru[isou],u4.get(ksou),a4.get(ksou),b4.get(ksou)
        wave.applyForward(src[isou],ui)
        wave.applyAdjoint(Source.ReceiverSource(rco[isou]),ai)
        # TODO: second time derivative?
        wave.applyAdjoint(Source.WavefieldSource(ai,rui),bi)
        ga = wave.collapse(ui,bi,nabsorb) # source side
        wave.applyForward(Source.WavefieldSource(ui,rui),bi)
        gb = wave.collapse(bi,ai,nabsorb) # receiver side
        return add(ga,gb)
      def combine(self,ga,gb):
        return add(ga,gb)
    timer.start('gradient')
    g = PartialParallel(np).reduce(offset,ns-offset,stride,Reduce()) 
    #g = PartialParallel(np).reduce(nx/2,1+nx/2,stride,Reduce()) # one shot
    print 'sum(g)=%f'%sum(g)
    timer.stop('gradient')
    RecursiveExponentialFilter(1.0/(fpeak*dx*sqrt(2.0))).apply(g,g) # XXX
    p = getConjugateDirection(g,gm,pm)
    gm,pm = g,p # rotate
    pixels(g,cmap=rwb,sperc=99.9,title='g_%d'%i)
    pixels(p,cmap=rwb,sperc=99.9,title='p_%d'%i)

    # Line search
    if doLineSearch:
      timer.start('line search')
      aopt = findStepLength(s,p,mp,nlsr,offset,rco,u4,a4,warping)
      timer.stop('line search')
    else:
      aopt = 0.0

    # Update slowness model
    add(s,mul(aopt,p),s)
    born.setSlowness(s)
    wave.setSlowness(s)
    pixels(s,cmap=jet,title='s_%d'%i)

    timer.stop('ITERATION')

def findStepLength(s,p,m,nmig,offset,reco,u4,a4,warping):
  """Line search using only one shot."""
  amin = -0.04 # minimum step length
  amax = 0.00 # maximum step length FIXME
  #atol = 0.01 # tolerance
  atol = 0.005 # tolerance
  #atol = 0.001 # tolerance
  #atol = 0.0001 # tolerance
  hx = nx/2
  xss = [hx-offset,hx,hx+offset]
  zss = [0,0,0]
  src = zeros(3,Source)
  rco = zeros(3,Receiver)
  rcp = zeros(3,Receiver)
  r = zerofloat(nx,nz,3)
  for isou in range(3):
    src[isou] = Source.RickerSource(xss[isou],zss[isou],dt,fpeak)
    rco[isou] = Receiver(reco[xss[isou]])
    rcp[isou] = Receiver(xr,zr,nt)
  born = BornOperatorS(s,dx,dt,nabsorb,u4,a4)
  ref = RecursiveExponentialFilter(1.0/(fpeak*dx*sqrt(2.0)))
  bsolver = BornSolver(born,src,rcp,rco,ref,None)
  class Misfit(BrentMinFinder.Function):
    def evaluate(self,a):
      e = mul(a/max(abs(p)),p)
      add(s,e,e)
      born.setSlowness(e)
      bsolver.solve(nmig,r)
      rm = transpose12(r[0]); mul(m,rm,rm)
      rr = transpose12(r[1]); mul(m,rr,rr)
      rp = transpose12(r[2]); mul(m,rp,rp)
      um = warping.findShifts(rr,rm)
      up = warping.findShifts(rr,rp)
      #add(um,up,up)
      add(abs(um),abs(up),up)
      mul(m,up,up)
      pixels(rr,title='_rr_%f'%a)
      pixels(up,cmap=rwb,sperc=100.0,title='_up_%f'%a)
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
  _,_,_,bs,src,rco,_,_,_,_ = getInputs()
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
  _,wave,_,_,src,rco,u,a,_,_ = getInputs()
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
  #mul(r,r,r)
  abs(r,r)
  efilter(8.0,r,r)
  mul(1.0/max(r),r,r)
  return r

def efilter(sigma,f,g,niter=2):
  ref = RecursiveExponentialFilter(sigma/sqrt(niter))
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref.apply(f,g)
  for i in range(1,niter):
    ref.apply(g,g)

def getStackedImage(rr):
  n = len(rr)
  r = like(rr[0])
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
  #if nx!=767:
  #  GaussianTaper.apply1(0.25,r,r)
  return s0,r

def like(x):
  return zerofloat(len(x[0]),len(x))

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
  #sp.setSize(1010,740)
  sp.setSize(1000,650)
  #sp.setFontSizeForSlide(1.0,1.0)
  cb = sp.addColorBar()
  cb.setWidthMinimum(100)
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
