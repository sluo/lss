#############################################################################
# Image-warping tomography

from imports import *

#############################################################################

savDir = None
#savDir = os.getenv('HOME')+'/Desktop/pngdat/'
#savDir = os.getenv('HOME')+'/Desktop/pngdat2/'

datDir = None
#datDir = './dat/'

plots = True
#plots = False

#############################################################################

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
  return getLayeredModelAndMask()

def getLayeredModelAndMask():
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

def addGaussianError(s,emax=0.010,m=None):
  g = zerofloat(nx,nz)
  g[nz/2][nx/2] = 1.0
  RecursiveGaussianFilter(nz/5.0).apply00(g,g)
  mul(emax/max(g),g,g)
  if m is not None:
    mul(m,g,g)
  return add(s,g)

def getModelAndMask():
  t,s,r,m,p = setupForLayered()
  e = copy(s)
  e = mul(0.90,s) # erroneous background slowness
  #e = addGaussianError(s,emax=0.025,m=m)
  pixels(t,cmap=jet,title='t')
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
  nvel = 1 # inversion iterations
  nmig = 1 # migration iterations within each inversion iteration
  offset = 20 # shot offset for warping
  #stride = 10 # shot stride for gradient computation
  stride = 1+nx/2 # shot stride for gradient computation
  r = zerofloat(nx,nz,ns) # reflectivity image
  rr = zerofloat(nz,nx,ns) # transposed image
  rm = zerofloat(nz,nx,ns) # transposed shifted imaged
  rp = zerofloat(nz,nx,ns) # transposed shifted imaged
  um = zerofloat(nz,nx,ns) # shifts
  up = zerofloat(nz,nx,ns) # shifts
  gm,pm = None,None # gradient & conjugate-gradient directions
  for i in range(nvel):
    timer.start('ITERATION')

    # Migration
    timer.start('migration')
    bs.solve(nmig,r)
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
    Parallel.loop(ns,Loop())

    # Warping
    timer.start('warping')
    warping.findShifts(rr,rm,um)
    warping.findShifts(rr,rp,up)
    timer.stop('warping')
    uu = add(um,up)
    ru = zerofloat(nx,nz,ns) # derivative of scaled images
    class Loop(Parallel.LoopInt):
      def compute(self,isou):
        for ix in range(nx):
          for iz in range(1,nz-1):
            ru[isou][iz][ix] =\
              0.5*uu[isou][ix][iz]*(rr[isou][ix][iz+1]-rr[isou][ix][iz-1])/dz
    Parallel.loop(ns,Loop())

    rw = warping.applyShifts(um,rm)
    pixels(um[ns/2],cmap=rwb,sperc=100.0,title='um')
    pixels(rm[ns/2],sperc=100.0,title='rm')
    pixels(rw[ns/2],sperc=100.0,title='rw')
    pixels(rr[ns/2],sperc=100.0,title='rr')
    pixels(ru[ns/2],sperc=100.0,title='ru')

    # Gradient and conjugate gradient
    class Reduce(Parallel.ReduceInt):
      def compute(self,isou):
        ui,ai,bi = u4.get(isou),a4.get(isou),b4.get(isou)
        rui = ru[Sampling(nx).indexOfNearest(xs[isou])]
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
    g = mul(1.0/ns,PartialParallel(np).reduce(stride,ns,stride,Reduce()))
    timer.stop('gradient')
    p = getConjugateDirection(g,gm,pm)
    gm,pm = g,p # rotate
    pixels(g,cmap=rwb,sperc=100.0,title='g%d'%i)
    pixels(p,cmap=rwb,sperc=100.0,title='p%d'%i)

    # Update slowness model
    aopt = findStepLength(s,p,mp,offset,rco,u4,a4,warping)
    add(s,mul(aopt,p),s)
    born.setSlowness(s)
    wave.setSlowness(s)
    pixels(s,cmap=jet,title='s%d'%i)

    timer.stop('ITERATION')

def findStepLength(s,p,mp,nmig,offset,rco,u4,a4,warping):
  ms = 1+int(ns/offset)
  amin,amax = -0.01,0.05; atol = 0.2*(amax-amin)
  q = mul(1.0/max(abs(p)),p) # normalized conjugate direction
  r = zerofloat(nx,nz,ms) # reflectivity image
  rr = zerofloat(nz,nx,ms) # transposed image
  rm = zerofloat(nz,nx,ms) # transposed shifted imaged
  rp = zerofloat(nz,nx,ms) # transposed shifted imaged
  um = zerofloat(nz,nx,ms) # shifts
  up = zerofloat(nz,nx,ms) # shifts
  class Misfit(BrentMinFinder.Function):
    def __init__(self):
      self.xs = rampint(0,ms,offset)
      self.zs = fillint(0,ms)
      self.src = zeros(ms,Source)
      self.rco = zeros(ms,Receiver)
      self.rcp = zeros(ms,Receiver)
      for isou in range(ms):
        self.src[isou] = Source.RickerSource(
          self.xs[isou],self.zs[isou],dt,fpeak)
        self.rco[isou] = Receiver(rco[Sampling(nx).indexOfNearest(xs[isou])])
        self.rcp[isou] = Receiver(xr,zr,nt)

    def evaluate(self,a):
      # Born modeling and solver.
      e = add(s,mul(a,q))
      born = BornOperatorS(e,dx,dt,nabsorb,u4,a4)
      ref = RecursiveExponentialFilter(1.0/(fpeak*dx*sqrt(2.0)))
      m = None
      bs = BornSolver(born,self.src,self.rcp,self.rco,ref,m)
      bs.solve(nmig,r)
      class Loop(Parallel.LoopInt):
        def compute(self,isou):
          transpose12(r[isou  ],rr[isou])
          transpose12(r[isou-1],rm[isou])
          transpose12(r[isou+1],rp[isou])
      Parallel.loop(1,ms-1,Loop())
      class Loop(Parallel.LoopInt):
        def compute(self,isou):
          mul(mp,rr[isou],rr[isou])
          mul(mp,rm[isou],rm[isou])
          mul(mp,rp[isou],rp[isou])
      Parallel.loop(ms,Loop())

      # Warping
      warping.findShifts(rr,rm,um)
      warping.findShifts(rr,rp,up)
      uu = add(um,up)
      class Loop(Parallel.LoopInt):
        def compute(self,isou):
          mul(mp,uu[isou])
      Parallel.loop(ms,Loop())

      return sum(mul(uu,uu))

  return BrentMinFinder(Misfit()).findMin(amin,amax,atol)

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
  #pixels(r)
  mul(r,r,r)
  #pixels(r)
  efilter(8.0,r,r)
  #pixels(r)
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

def transpose12(f,g):
  n1,n2 = len(f[0]),len(f)
  for i2 in range(n2):
    for i1 in range(n1):
      g[i1][i2] = f[i2][i1]

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
  if sperc is not None: # symmetric percentile clip
    clips = Clips(100-sperc,sperc,x)
    clip = max(abs(clips.getClipMin()),abs(clips.getClipMax()))
    pv.setClips(-clip,clip)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if title and savDir:
    #sp.paintToPng(360,3.33,savDir+title+'.png')
    sp.paintToPng(180,3.33,savDir+title+'.png')
    write(savDir+title+'.dat',x)

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
