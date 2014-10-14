#############################################################################
# Image-warping tomography

from imports import *

#############################################################################

savDir = None
#savDir = os.getenv('HOME')+'/Desktop/pngdat/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat2/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat3/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat4/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat5/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat6/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat7/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat8/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat9/'
#savDir = '/Users/sluo/Dropbox/pngdat/'
#savDir = '/Users/sluo/Dropbox/pngdat2/'
savDir = os.getenv('HOME')+'/Desktop/png/'

#datDir = './dat/multiLayer/nz401/plusMediumGaussian/'
#datDir = './dat/marmousi/' # ds=2,stride=4
#datDir = './dat2/marmousi/' # ds=4,stride=4
datDir = './dat3/marmousi/' # ds=8,stride=8
#datDir = './dat3/marmousiTrue/' # ds=8,stride=8, true background slowness

plots = True
#plots = False

#############################################################################

# TODO: change Ricker to Gaussian derivative
# TODO: second time derivative in adjoint source

timer = Timer()
def main(args):
  #getModelAndMask()
  #getInputs(modelChanged=True) # simulate data
  #goInversion()
  #plotSlowness()
  plotImages()

def plotSlowness():
  t,s,e,r,m = getModelAndMask()
  s0 = zerofloat(nx,nz)
  s1 = zerofloat(nx,nz)
  s2 = zerofloat(nx,nz)
  s3 = zerofloat(nx,nz)
  s4 = zerofloat(nx,nz)
  #read(os.getenv('HOME')+'/Desktop/iwt/marmousi/s_0.dat',s0)
  #read(os.getenv('HOME')+'/Desktop/iwt/marmousi/s_1.dat',s1)
  #read(os.getenv('HOME')+'/Desktop/iwt/marmousi/s_2.dat',s2)
  #read(os.getenv('HOME')+'/Desktop/iwt/marmousi/s_3.dat',s3)
  read(os.getenv('HOME')+'/Desktop/iwt/layered/s_0.dat',s0)
  read(os.getenv('HOME')+'/Desktop/iwt/layered/s_1.dat',s1)
  read(os.getenv('HOME')+'/Desktop/iwt/layered/s_2.dat',s2)
  read(os.getenv('HOME')+'/Desktop/iwt/layered/s_3.dat',s3)
  read(os.getenv('HOME')+'/Desktop/iwt/layered/s_4.dat',s4)
  smin,smax = 0.0,0.0
  #smin,smax = min(s),max(s)
  pixels(s,cmap=jet,cmin=smin,cmax=smax,cbar='Slowness (s/km)',title='s')
  pixels(s0,cmap=jet,cmin=smin,cmax=smax,cbar='Slowness (s/km)',title='s0')
  pixels(s1,cmap=jet,cmin=smin,cmax=smax,cbar='Slowness (s/km)',title='s1')
  pixels(s2,cmap=jet,cmin=smin,cmax=smax,cbar='Slowness (s/km)',title='s2')
  pixels(s3,cmap=jet,cmin=smin,cmax=smax,cbar='Slowness (s/km)',title='s3')

def plotImages():
  t,s,e,r,m = getModelAndMask()
  r0 = zerofloat(nx,nz)
  r1 = zerofloat(nx,nz)
  r2 = zerofloat(nx,nz)
  r3 = zerofloat(nx,nz)
  r4 = zerofloat(nx,nz)
  read(os.getenv('HOME')+'/Desktop/iwt/marmousi/g_0.dat',r0)
  read(os.getenv('HOME')+'/Desktop/iwt/marmousi/s_1.dat',r1)
  read(os.getenv('HOME')+'/Desktop/iwt/marmousi/s_2.dat',r2)
  read(os.getenv('HOME')+'/Desktop/iwt/marmousi/s_3.dat',r3)
  #read(os.getenv('HOME')+'/Desktop/iwt/layered/g_0.dat',r0)
  #read(os.getenv('HOME')+'/Desktop/iwt/layered/g_1.dat',r1)
  #read(os.getenv('HOME')+'/Desktop/iwt/layered/g_2.dat',r2)
  #read(os.getenv('HOME')+'/Desktop/iwt/layered/g_3.dat',r3)
  #read(os.getenv('HOME')+'/Desktop/iwt/layered/g_4.dat',r4)
  r0 = transpose(r0)
  r1 = transpose(r1)
  r2 = transpose(r2)
  r3 = transpose(r3)
  r4 = transpose(r4)
  #mul(dx*1000.0,r0,r0)
  #mul(dx*1000.0,r1,r1)
  #mul(dx*1000.0,r2,r2)
  #mul(dx*1000.0,r3,r3)
  #mul(dx*1000.0,r4,r4)
  #smin,smax,sp = 0.0,0.0,99.9
  smin,smax,sp = 0.0,0.0,100.0
  #smin,smax,sp = 0.0,0.0,None
  #smin,smax,sp = min(s),max(s),None
  #cbar,cmap = None,gray
  cbar,cmap = None,rwb
  #cbar,cmap = 'Slowness (s/km)',jet
  #cbar,cmap = 'Reflectivity',gray
  #cbar,cmap = 'Shift (m)',rwb
  pixels(r0,cmap,cmin=smin,cmax=smax,sperc=sp,cbar=cbar,title='r0')
  pixels(r1,cmap,cmin=smin,cmax=smax,sperc=sp,cbar=cbar,title='r1')
  pixels(r2,cmap,cmin=smin,cmax=smax,sperc=sp,cbar=cbar,title='r2')
  pixels(r3,cmap,cmin=smin,cmax=smax,sperc=sp,cbar=cbar,title='r3')
  pixels(r4,cmap,cmin=smin,cmax=smax,sperc=sp,cbar=cbar,title='r4')

def setupForLayered():
  global zs,xs,ns,ds
  global sz,sx,st,nz,nx,nt,nxp,nzp,dz,dx,dt
  global fpeak,nabsorb,np,sigma0,sigma1

  #sx = Sampling(501,0.016,0.0)
  ##sz = Sampling(201,0.016,0.0)
  ##st = Sampling(1803,0.0015,0.0)
  #sz = Sampling(401,0.016,0.0)
  #st = Sampling(2503,0.0015,0.0)
  #sigma0,sigma1 = 0.2,0.05 # background/reflectivity models
  #fpeak = 20.0 # Ricker wavelet peak frequency
  #np = 16 # number of parallel shots
  #ds = 2 # shot stride

  sz = Sampling(265,0.012,0.0)
  sx = Sampling(767,0.012,0.0)
  #st = Sampling(5001,0.0012,0.0)
  st = Sampling(6001,0.001,0.0)
  sigma0,sigma1 = 0.2,0.05 # background/reflectivity models
  fpeak = 15.0 # Ricker wavelet peak frequency
  np = 8 # number of parallel shots
  #ds = 2 # shot stride
  #ds = 4 # shot stride
  ds = 8 # shot stride

  nz,nx,nt = sz.count,sx.count,st.count
  dz,dx,dt = sz.delta,sx.delta,st.delta
  fz,fx,ft = sz.first,sx.first,st.first
  ns = (nx-1+ds)/ds
  xs,zs = rampint(0,ds,ns),fillint(0,ns)
  nabsorb = 22 # absorbing boundary size
  nxp,nzp = nx+2*nabsorb,nz+2*nabsorb
  np = min(np,ns) # number of parallel shots
  #sigma0 = 4.0/fpeak
  #sigma1 = 1.0/fpeak
  print 'ns=%d'%ns
  #return getSingleLayeredModelAndMask()
  #return getMultiLayeredModelAndMask()
  return getMarmousi()

def getSourcesAndReceivers():
  #maxOffset = 80
  maxOffset = nx # full offsets, fixed spread
  src = zeros(ns,Source)
  rco = zeros(ns,Receiver)
  rcp = zeros(ns,Receiver)
  for isou in range(ns):
    fxr,lxr = max(xs[isou]-maxOffset,0),min(xs[isou]+maxOffset,nx-1)
    nxr = 1+lxr-fxr
    xr,zr = rampint(fxr,1,nxr),fillint(0,nxr)
    rco[isou] = Receiver(xr,zr,nt)
    rcp[isou] = Receiver(xr,zr,nt)
    src[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
  return src,rcp,rco

def getSingleLayeredModelAndMask():
  t,t0,t1,m = getMultiLayeredModelAndMask()
  nlayer = 1+int(nz*dz)/0.25
  klayer = nlayer-2 # depth of single layer
  for iz in range(int((klayer-0.5)*nz/nlayer)):
    zero(t1[iz])
  for iz in range(int((klayer+0.5)*nz/nlayer),nz):
    zero(t1[iz])
  sigma = 3
  m = fillfloat(1,nx,nz)
  for iz in range(klayer*nz/nlayer-3*sigma):
    zero(m[iz])
  ref = RecursiveExponentialFilter(sigma)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref.apply2(m,m)
  return t,t0,t1,m

def getMultiLayeredModelAndMask():
  tb = 0.25 # background slowness
  nlayer = 1+int(nz*dz)/0.25
  vstart,vstep = 100.0,-1.0
  t = fillfloat(vstart,nx,nz)
  #for ilayer in range(2,nlayer):
  for ilayer in range(3,nlayer):
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

import socket
def getMarmousi():
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
  t = fillfloat(2.0/3.0,nz,nx)
  copy(248,767,0,0,3,3,p,17,0,1,1,t)
  t = transpose(t)
  t0,t1 = makeBornModel(t)
  m = makeTaperedMask(16)
  mul(m,t1,t1)
  return t,t0,t1,m

def makeTaperedMask(kz,sigma=2.0):
  m = fillfloat(1.0,nx,nz)
  for iz in range(kz):
    for ix in range(nx):
      m[iz][ix] = 0.0
  efilter2(sigma,m,m)
  return m

def makeLinearModel(t,fz=0):
  nx,nz = len(t[0]),len(t)
  ts = zerofloat(nz)
  for iz in range(nz):
    ts[iz] = sum(t[iz])
  mul(1.0/nx,ts,ts)
  ls = linearRegression(rampint(0,1,nz),ts)
  s = copy(t)
  for iz in range(fz,nz):
    fill(ls[iz],s[iz])
  efilter(8.0,s,s)
  return s

def linearRegression(x,y):
  """Fits a line to (x,y).
  Parameters:
    x - x-coordinates.
    y - y-coordinates.
  Returns:
    The sampled best fit line in an array the same size as x and y.
  """
  n = len(x)
  sumx = 0.0 # sum_n(x_i)
  sumy = 0.0 # sum_n(y_i)
  sumxx = 0.0 # sum_n(x_i*x_i)
  sumxy = 0.0 # sum_n(x_i*y_i)
  for i in range(n):
    xi = x[i]
    yi = y[i]
    sumx += xi
    sumy += yi
    sumxx += xi*xi
    sumxy += xi*yi
  beta = (sumxy-sumx*sumy/n)/(sumxx-sumx*sumx/n)
  alpha = (sumy-beta*sumx)/n
  z = zerofloat(n)
  for i in range(n):
    z[i] = alpha+beta*x[i]
  return z

#############################################################################

def addGaussianError(s,emax=0.010,sigma=None):
  g = zerofloat(nx,nz)
  #g[nz/3][nx/2] = 1.0
  #g[2*nz/5][nx/2] = 1.0
  g[nz/2][nx/2] = 1.0
  #g[100][nx/2] = 1.0
  if sigma is None:
    RecursiveGaussianFilter(nz/8.0).apply00(g,g)
  else:
    RecursiveGaussianFilter(sigma).apply00(g,g)
  mul(emax/max(g),g,g)
  return add(s,g)

def getModelAndMask():
  t,s,r,m = setupForLayered()
  e = copy(s)
  subDir = datDir.split('/')[-2]
  #s = mul(0.90,e) # erroneous background slowness
  #s = mul(1.10,e) # erroneous background slowness #
  #s = addGaussianError(e,emax=0.05)
  #s = addGaussianError(e,emax=-0.05)
  #s = addGaussianError(e,emax=0.05,sigma=nz/10.0) # plusGaussian
  #s = addGaussianError(e,emax=-0.05,sigma=nz/10.0) # minusGaussian
  #s = addGaussianError(e,emax=-0.05,sigma=nz/10.0)
  #s = addGaussianError(e,emax=0.04)
  #s = addGaussianError(e,emax=-0.04)
  #s = addGaussianError(e,emax=-0.04)
  if subDir=='plusGaussian':
    #s = addGaussianError(e,emax=0.05,sigma=nz/10.0) # plusGaussian
    s = addGaussianError(e,emax=0.05,sigma=20.0) # plusGaussian
  elif subDir=='plusLargeGaussian':
    #s = addGaussianError(e,emax=0.05,sigma=nz/3.0) # plusGaussian
    s = addGaussianError(e,emax=0.05,sigma=67.0) # plusGaussian
  elif subDir=='plusMediumGaussian':
    #s = addGaussianError(e,emax=0.025,sigma=nz/4.0) # plusGaussian
    s = addGaussianError(e,emax=0.025,sigma=50.0) # plusGaussian
  elif subDir=='mul110':
    s = mul(1.10,e) # erroneous background slowness #
  elif subDir=='marmousi':
    e = makeLinearModel(s,fz=20)
  elif subDir=='marmousiTrue':
    e = copy(s)
  pixels(t,cmap=jet,cbar='Slowness (s/km)',title='t')
  pixels(s,cmap=jet,cbar='Slowness (s/km)',title='s')
  pixels(e,cmap=jet,cbar='Slowness (s/km)',title='e')
  if m is not None:
    pixels(m,cmap=gray,title='m')
  pixels(sub(s,e),cmap=rwb,sperc=100.0,cbar='Slowness (s/km)',title='s-e')
  pixels(r,cmap=gray,sperc=100.0,cbar='Reflectivity',title='r')
  print 'sum(sub(s,e))=%f'%(sum(sub(s,e)))

  global GTRUE
  GTRUE = sub(e,s)

  return t,s,e,r,m

def getInputs(modelChanged=True):
  t,s,e,r,m = getModelAndMask()

  # Warping
  maxShift = 0.20 # max shift in km
  #bstrain1,bstrain2,bstrain3 = 16,8,8
  bstrain1,bstrain2,bstrain3 = 8,4,4
  #bstrain1,bstrain2,bstrain3 = 4,4,4 ###
  #bstrain1,bstrain2,bstrain3 = 4,2,2
  smooth1,smooth2,smooth3 = 4.0,4.0,4.0
  bstrain3 = max(bstrain3/ds,1)
  warping = ImageWarping(
    bstrain1,bstrain2,bstrain3,smooth1,smooth2,smooth3,maxShift,dz)
  warping.setWindowSizeAndOverlap(nx,nx,0.5,0.5)
  print 'bstrain1=%d'%bstrain1
  print 'bstrain2=%d'%bstrain2
  print 'bstrain3=%d'%bstrain3

  # Wavefields.
  timer.start('allocating')
  u4 = SharedFloat4(nxp,nzp,nt,np)
  a4 = SharedFloat4(nxp,nzp,nt,np)
  b4 = SharedFloat4(nxp,nzp,nt,np)
  timer.stop('allocating')

  # Sources and receivers.
  src,rcp,rco = getSourcesAndReceivers()
  if modelChanged: # if model has changed, resimulated data...
  #if False:
    bornt = BornOperatorS(s,dx,dt,nabsorb,u4,a4) # true slowness
    timer.start('observed data')
    bornt.applyForward(src,r,rco) # observed data
    timer.stop('observed data')
    timer.start('writing observed data')
    for isou in range(ns):
      write(datDir+'d%d.dat'%isou,rco[isou].getData())
    timer.stop('writing observed data')
  else: # ...else read from disk.
    timer.start('reading observed data')
    for isou in range(ns):
      read(datDir+'d%d.dat'%isou,rco[isou].getData())
    timer.stop('reading observed data')
  pixels(rco[ns/2].getData(),sperc=99.9,title='rco')

  # Born modeling and solver.
  born = BornOperatorS(e,dx,dt,nabsorb,u4,a4)
  ref = RecursiveExponentialFilter(1.0/(fpeak*dx*sqrt(2.0)))
  bornsolver = BornSolver(born,src,rcp,rco,ref,m) # solver

  # Wave operator
  wave = WaveOperator(e,dx,dt,nabsorb)

  return e,m,wave,born,bornsolver,src,rco,u4,a4,b4,warping

#############################################################################

def goInversion():
  modelChanged = False # if True, resimulate and remigrated observed data
  doLineSearch = False
  ninv = 5 # inversion iterations
  nmig = 2 # migration iterations within each inversion iteration
  nlsr = 2 # migration iterations within each inversion iteration
  #offset,stride = 10,4 # shot offset for warping, stride for gradient
  #offset,stride = 20,4 # shot offset for warping, stride for gradient
  #offset,stride = 16,8 # shot offset for warping, stride for gradient
  #offset,stride = 24,8 # shot offset for warping, stride for gradient
  offset,stride = 40,4 # shot offset for warping, stride for gradient
  s,m,wave,born,bs,src,rco,u4,a4,b4,warping = getInputs(modelChanged)
  Check.argument(offset%ds==0,'offset%ds==0')
  Check.argument(stride%ds==0,'stride%ds==0')
  offset /= ds; stride /= ds # divide by shot stride

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
  pscale = None # scale conjugate-gradient direction before line search

  for i in range(ninv):
    print ''
    timer.start('ITERATION %d'%i)

    # Migration
    if modelChanged or i>0:
      timer.start('migration')
      bs.solve(nmig,r)
      timer.stop('migration')
      if i==0:
        write(datDir+'r.dat',r)
    else:
      read(datDir+'r.dat',r)
    
    # Transpose images
    transposeImages(r,rr)

    # Model preconditioner TODO: here or below? (shouldn't matter)
    #mp = getModelPrecondition(rr) # preconditioner
    #applyModelPrecondition(mp,rr,rr)
    #pixels(mp,title='mp')

    # Find warping shifts
    findShifts(offset,warping,m,rr,rm,rp,um,up,i)
  
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
      #applyModelPrecondition(mp,rr,rr)
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
    #taperImages(offset,rr,rr)
    #taperImages(offset,wm,wm)
    #taperImages(offset,wp,wp)

    # Adjoint sources
    makeAdjointSources(wm,wp,rr,ru)

    class Loop(Parallel.LoopInt):
      def compute(self,isou):
        for iz in range(nz):
          mul(float(iz),ru[isou][iz],ru[isou][iz]) # XXX
          #mul(pow(float(iz),2.0),ru[isou][iz],ru[isou][iz]) # XXX
          #mul(pow(float(iz),4.0),ru[isou][iz],ru[isou][iz]) # XXX
          pass
    Parallel.loop(ns,Loop())

    #taperAdjointSources(ru)

    # Gradient and conjugate gradient
    g = computeGradient(offset,stride,ru,wave,src,rco,u4,a4,b4)
    efilter(1.0/(fpeak*dx),g,g) # XXX
    #efilter(5.0/(fpeak*dx),g,g) # XXX
    #mul(transpose(m),g,g) # mask XXX
    p = getConjugateDirection(g,gm,pm)
    gm,pm = g,p # rotate
    print 'sum(g)=%f'%sum(g)
    print 'sum(p)=%f'%sum(p)

#    # Line search
#    if doLineSearch or ninv>1:
#      timer.start('line search')
#      #aopt = findStepLength(s,m,p,nlsr,offset,rco,u4,a4,warping)
#      aopt = findStepLengthX(s,m,p,nlsr,offset,rco,u4,a4,warping) # test
#      timer.stop('line search')
#    else:
#      aopt = 0.0
#
#    # Update slowness model
#    add(s,mul(aopt,p),s)
#    born.setSlowness(s)
#    wave.setSlowness(s)

    if doLineSearch or ninv>1:
      timer.start('line search')
      updateSlownessModel(s,m,p,nlsr,offset,src,rco,u4,a4,warping)
      #if i==0: # if first iteration, set scale
      #  pscale = 1.0/max(abs(p))
      #updateSlownessModelX(s,p,pscale) # fixed step length
      timer.stop('line search')
      born.setSlowness(s)
      wave.setSlowness(s)

    #if i==(ninv-1) and (doLineSearch or ninv>1):
    if i==(ninv-1) and ninv>1:
      timer.start('migration')
      bs.solve(nmig,r)
      timer.stop('migration')
      pixels(getStackedImage(r),sperc=100.0,title='r_final')

    # Plot
    plots(i,rr,rn,rm,rp,rmw,rpw,um,up,sm,sp,wm,wp,ru,g,p,s)

    # Zero arrays
    zero(r)
    zero(rr)
    zero(rn)
    zero(rm)
    zero(rp)
    zero(rmw)
    zero(rpw)
    zero(um)
    zero(up)
    zero(sm)
    zero(sp)
    zero(wm)
    zero(wp)
    zero(ru)

    timer.stop('ITERATION %d'%i)

def plots(iiter,rr,rn,rm,rp,rmw,rpw,um,up,sm,sp,wm,wp,ru,g,p,s):
  #ks = [ns/2] # shots to plot
  ks = [ns/3,ns/2] # shots to plot
  for k in ks:
    pixels(getStackedImage(rr),title='r_%d'%iiter)
    pixels(rr[k],sperc=100.0,title='rr_%d_isou=%d'%(iiter,k))
    pixels(rn[k],sperc=100.0,title='rn_%d_isou=%d'%(iiter,k))
    #pixels(rm[k],sperc=100.0,title='rm_%d_isou=%d'%(iiter,k))
    pixels(rp[k],sperc=100.0,title='rp_%d_isou=%d'%(iiter,k))
    #pixels(rmw[k],sperc=100.0,title='rmw_%d_isou=%d'%(iiter,k))
    pixels(rpw[k],sperc=100.0,title='rpw_%d_isou=%d'%(iiter,k))
    #pixels(um[k],cmap=rwb,sperc=100.0,title='um_%d_isou=%d'%(iiter,k))
    pixels(up[k],cmap=rwb,sperc=100.0,title='up_%d_isou=%d'%(iiter,k))
    pixels(add(um[k],up[k]),cmap=rwb,sperc=100.0,
      title='um+up_%i_isou=%d'%(iiter,k))
    #pixels(sm[k],cmap=rwb,sperc=100.0,title='sm_%d_isou=%d'%(iiter,k))
    pixels(sp[k],cmap=rwb,sperc=100.0,title='sp_%d_isou=%d'%(iiter,k))
    #pixels(wm[k],cmap=rwb,sperc=100.0,title='wm_%d_isou=%d'%(iiter,k))
    #pixels(wp[k],cmap=rwb,sperc=100.0,title='wp_%d_isou=%d'%(iiter,k))
    #pixels(add(wm[k],wp[k]),cmap=rwb,sperc=100.0,
    #  title='wm+wp_%d_isou=%d'%(iiter,k))
    pixels(ru[k],sperc=100.0,title='ru_%d_isou=%d'%(iiter,k))
  pixels(g,cmap=rwb,sperc=100.0,title='g_%d'%iiter)
  pixels(p,cmap=rwb,sperc=100.0,title='p_%d'%iiter)
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

def applyModelPrecondition(mp,rr,ry):
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      mul(mp,rr[isou],ry[isou])
  Parallel.loop(ns,Loop())

addRandomNoise = False
def findShifts(offset,warping,m,rr,rm,rp,um,up,iii):
  rn = like3(rr)
  normalizeImages(m,rr,rn) # normalize images before warping
  offsetImages(offset,rn,rm,rp) # offset normalized images
  taperImages(offset,rn,rn) # taper images
  taperImages(offset,rm,rm) # taper images
  taperImages(offset,rp,rp) # taper images
  timer.start('warping')
  if addRandomNoise:
    warping.filterAndFindShifts(copy(rn),copy(rm),um)
    warping.filterAndFindShifts(copy(rn),copy(rp),up)
  else:
    warping.findShifts(rn,rm,um)
    warping.findShifts(rn,rp,up)
    #warping.findShifts(rm,rn,um); mul(-1.0,um,um)
    #warping.findShifts(rp,rn,up); mul(-1.0,up,up)
  timer.stop('warping')
  taperImages(offset,um,um) # taper shifts
  taperImages(offset,up,up) # taper shifts
  applyModelPrecondition(m,um,um)
  applyModelPrecondition(m,up,up)

def applyShifts(warping,um,up,rm,rp,rmw,rpw):
  warping.applyShifts(um,rm,rmw) # shift image
  warping.applyShifts(up,rp,rpw) # shift image

def makeScale(offset,m,rr,rmw,rpw,dm,dp):
  niter = 8 # number of times to apply ref
  #sigma = 0.064 # sigma in km
  #sigma = 0.125 # sigma in km
  sigma = 0.250 # sigma in km #
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
      #mul(dma,dma,dma) # raised to 4th power XXX
      #mul(dpa,dpa,dpa) # raised to 4th power XXX
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

def xmakeScaleN(m,rr,sm,sp=None,useAbs=True):
  print "not using scale!"
  nsou = len(rr)
  for isou in range(nsou):
    copy(m,sm[isou])
    if sp is not None:
      copy(m,sp[isou])

def makeScaleN(m,rr,sm,sp=None,useAbs=True):
  #sigma = 0.125 # sigma in km
  sigma = 0.250 # sigma in km
  nsou = len(rr)
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      smi = sm[isou]
      if useAbs:
        abs(rr[isou],smi)
      else:
        mul(rr[isou],rr[isou],smi)
      efilter(sigma/dx,smi,smi)
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
      taperImage(xs[isou],rr[isou],ry[isou])
  Parallel.loop(ns,Loop())

def taperImage(xsou,rx,ry,imageIsTransposed=True):
  #hx = 10 # half taper width
  #hx = 20 # half taper width
  #hx = 40 # half taper width
  #hx = 60 # half taper width
  #hx = 80 # half taper width
  hx = 100 # half taper width
  #hx = 3*offset # half taper width
  #hx = 4*offset # half taper width
  #hx = 6*offset # half taper width
  #hx = nx # no taper
  #sigma = 0.125/dx # half width for taper transition
  sigma = 8 # half width for taper transition
  w = zerofloat(nx)
  fx,lx = max(xsou-hx,0),min(xsou+hx,nx-1)
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
  #fsou = offset # first source
  fsou = 2*offset # first source
  #fsou = 4*offset # first source
  #fsou = 5*offset # first source
  class Reduce(Parallel.ReduceInt):
    def compute(self,isou):
      ksou = (isou-fsou)/stride
      rui,ui,ai,bi = ru[isou],u4.get(ksou),a4.get(ksou),b4.get(ksou)
      wave.applyForward(src[isou],ui) # source wavefield
      #gp = wave.collapse(ui,ui,nabsorb) # gradient preconditioner
      #mul(gp,gp,gp) # squared
      ##mul(1.0/max(gp),gp,gp) # normalized
      wave.applyAdjoint(Source.ReceiverSource(rco[isou]),ai) # rec wavefield
      # TODO: second time derivative?
      wave.applyAdjoint(Source.WavefieldSource(ai,rui),bi)
      ga = wave.collapse(ui,bi,nabsorb) # source side

      sourceSideOnly = False
      if sourceSideOnly:
        return ga
      else:
        wave.applyForward(Source.WavefieldSource(ui,rui),bi)
        gb = wave.collapse(bi,ai,nabsorb) # receiver side
        add(ga,gb,gb)
        #if isou==ns/2 or isou==2*offset:
        #  pixels(gb,cmap=rwb,sperc=100.0,title='g_isou=%d'%isou)

        ##for iz in range(nz):
        ##  mul(pow(float(iz),2.0),gb[iz],gb[iz])
        #div(gb,gp,gb)
        return gb
    def combine(self,ga,gb):
      return add(ga,gb)
  timer.start('gradient')
  g = PartialParallel(np).reduce(fsou,ns-fsou,stride,Reduce()) 
  #g = PartialParallel(np).reduce(ns/2,ns/2+1,stride,Reduce()) # one shot
  timer.stop('gradient')
  return g

def updateSlownessModel(s,m,p,nmig,offset,source,receiver,u4,a4,warping):
  hs = 4 # total number shots = 1+2*hs
  amin = -0.005 # minimum step length
  amax = 0.005 # maximum step length
  atol = 0.001 # tolerance
  nsou = 1+2*hs
  ks = rampint(ns/2-hs*offset,offset,nsou) # source indices
  zss = zeroint(nsou)
  src = zeros(nsou,Source)
  rco = zeros(nsou,Receiver)
  rcp = zeros(nsou,Receiver)
  r = zerofloat(nx,nz,nsou)
  for isou in range(nsou):
    srci = source[ks[isou]]
    reci = receiver[ks[isou]]
    src[isou] = srci
    rco[isou] = reci
    rcp[isou] = Receiver(reci.getXIndices(),reci.getZIndices(),nt)
  ss = mul(s,s) # slowness squared
  ref = RecursiveExponentialFilter(1.0/(fpeak*dx*sqrt(2.0)))
  born = BornOperatorS(s,dx,dt,nabsorb,u4,a4)
  bsolver = BornSolver(born,src,rcp,rco,ref,transpose(m))
  class Misfit(BrentMinFinder.Function):
    def evaluate(self,a):
      e = mul(a/max(abs(p)),p)
      add(ss,e,e) # update slowness squared
      sqrt(e,e)
      born.setSlowness(e)
      bsolver.solve(nmig,r)
      normalizeImages(transpose(m),r,r)
      class Reduce(Parallel.ReduceInt):
        def compute(self,isou):
          rm = transpose12(r[isou-1])
          rr = transpose12(r[isou  ])
          rp = transpose12(r[isou+1])
          taperImage(xs[ks[isou]],rr,rr,imageIsTransposed=True)
          taperImage(xs[ks[isou]],rm,rm,imageIsTransposed=True)
          taperImage(xs[ks[isou]],rp,rp,imageIsTransposed=True)
          um = warping.findShifts(rr,rm)
          up = warping.findShifts(rr,rp)
          taperImage(xs[ks[isou]],um,um,imageIsTransposed=True)
          taperImage(xs[ks[isou]],up,up,imageIsTransposed=True)
          mul(m,um,um)
          mul(m,up,up)
          if isou==1:
          #if isou==nsou/2:
            ##pixels(e,cmap=jet,title='_s_%f'%a)
            #pixels(rr,title='_rr_%f'%a)
            #pixels(rp,title='_rp_%f'%a)
            ##pixels(up,cmap=rwb,sperc=100.0,title='_up_%f'%a)
            pass
          mul(um,um,um)
          mul(up,up,up)
          return add(um,up)
        def combine(self,ua,ub):
          return add(ua,ub)
      return sum(Parallel.reduce(1,nsou-1,Reduce()))
  aopt = BrentMinFinder(Misfit()).findMin(amin,amax,atol)
  print 'aopt=%f'%aopt
  e = mul(aopt/max(abs(p)),p)
  add(ss,e,e)
  sqrt(e,s)
  #return aopt/max(abs(p))

def updateSlownessModelX(s,p,pscale=None):
  """Fixed step length."""
  #aopt = 0.002
  #aopt = 0.005
  #aopt = 0.010
  aopt = 0.050
  #aopt = -0.005
  print 'aopt=%f'%aopt
  #e = mul(aopt/max(abs(p)),p) if
  e = mul(pscale*aopt,p) if pscale else mul(aopt/max(abs(p)),p)
  add(mul(s,s),e,e)
  sqrt(e,s)

def getModelPrecondition(rr):
  r = getStackedImage(rr)
  return getEnergy(r)

def getEnergy(r):
  e = abs(r)
  #e = mul(r,r)
  efilter(8.0,e,e)
  mul(1.0/max(e),e,e)
  return e

def efilter(sigma,f,g,niter=4):
  ref = getRef(sigma/sqrt(niter))
  ref.apply(f,g)
  for i in range(1,niter):
    ref.apply(g,g)
def efilter1(sigma,f,g,niter=4):
  ref = getRef(sigma/sqrt(niter))
  ref.apply1(f,g)
  for i in range(1,niter):
    ref.apply1(g,g)
def efilter2(sigma,f,g,niter=4):
  ref = getRef(sigma/sqrt(niter))
  ref.apply2(f,g)
  for i in range(1,niter):
    ref.apply2(g,g)
def getRef(sigma):
  ref = RecursiveExponentialFilter(sigma)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  return ref

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

def getConjugateDirection(gi,gm=None,pm=None):
  """ Polak-Ribiere nonlinear conjugate gradient method.
  Parameters:
    gi - gradient ascent direction for current iteration.
    gm - gradient ascent direction for previous iteration.
    pm - conjugate ascent direction for previous iteration.
  Returns:
    conjugate ascent direction for the current iteration.
  """
  if gm is None or pm is None:
    return gi
  else:
    beta = sum(mul(gi,sub(gi,gm)))/sum(mul(gm,gm))
    orthogonality = abs(sum(mul(gi,gm)))/sum(mul(gi,gi))
    print 'beta=%f'%beta
    print 'orthogonality=%f'%orthogonality
    if beta<0.0:
      print "restarting CG (beta<0.0)"
      beta = 0.0
    elif orthogonality<0.1: # 0.1 threshold
      print "restarting CG (orthogonality<0.1)"
      beta = 0.0
    return add(gi,mul(beta,pm))
    #b = sum(mul(g,sub(g,gm)))/sum(mul(gm,gm))
    #if b<0.0:
    #  b = 0.0
    #  print "  CG DIRECTION RESET"
    #return add(g,mul(b,pm))

#############################################################################

def makeBornModel(s):
  """
  sigma0: smoothing for background model
  sigma1: smoothing for perturbation
  """
  print 'sigma0=%f'%sigma0
  print 'sigma1=%f'%sigma1
  s0,s1 = like2(s),like2(s)
  #ref0 = RecursiveExponentialFilter(sigma0/dx)
  #ref0.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  #ref0.apply(s,s0)
  efilter(sigma0/dx,s,s0)
  t = copy(s)
  #ref1 = RecursiveExponentialFilter(sigma1/dx)
  #ref1.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  #ref1.apply(s,t)
  efilter(sigma1/dx,s,t)
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
def pixels(x,cmap=gray,perc=100.0,sperc=None,
  cbar=None,cmin=0.0,cmax=0.0,title=None):
  if not plots:
    return
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if nx==189:
    sp.setSize(1000,650)
  elif nx==251:
    sp.setSize(1265,650)
  elif nx==501 and nz==201:
    sp.setSize(1530,640)
  elif nx==501 and nz==401:
    #sp.setSize(1210,900)
    sp.setSize(1270,900)
  else:
    sp.setSize(1000,650)
  if cbar is not None:
    cb = sp.addColorBar(cbar)
  else:
    cb = sp.addColorBar()
  #cb.setWidthMinimum(120)
  sp.setFontSizeForSlide(1.2,1.2)
  #cb.setWidthMinimum(190)
  cb.setWidthMinimum(150)
  if title:
    #sp.addTitle(title)
    pass
  if len(x)==nz:
    sp.setHLabel("Distance (km)")
    sp.setVLabel("Depth (km)")
    pv = sp.addPixels(sz,sx,transpose(x))
    #pv = sp.addPixels(transpose(x))
  elif len(x)==nx:
    sp.setHLabel("Distance (km)")
    sp.setVLabel("Depth (km)")
    pv = sp.addPixels(sz,sx,x)
    #pv = sp.addPixels(transpose(x))
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
    sp.paintToPng(360,3.33,savDir+title+'.png')
    #sp.paintToPng(720,3.33,savDir+title+'.png')

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
