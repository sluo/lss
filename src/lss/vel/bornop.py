#############################################################################
# Wavefield modeling using AcousticWaveOperator

from imports import *
from dnp import *

#############################################################################

#sz = Sampling(11,0.016,0.0)
#sx = Sampling(12,0.016,0.0)
#st = Sampling(13,0.0012,0.0)
#sz = Sampling(201,0.016,0.0)
#sx = Sampling(202,0.016,0.0)
#st = Sampling(2003,0.0012,0.0)
#sz = Sampling(401,0.0025,0.0) # for 40 Hz
#sx = Sampling(402,0.0025,0.0)
#st = Sampling(5003,0.0001,0.0)
sz = Sampling(201,0.005,0.0) # for 30 Hz
sx = Sampling(202,0.005,0.0)
st = Sampling(3003,0.0004,0.0)
#sz = Sampling(265,0.012,0.0); stride = 3
#sx = Sampling(767,0.012,0.0)
#st = Sampling(5001,0.0012,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
#xs,zs = [0],[0]
xs,zs = [nx/2],[0]
#xs,zs = [nx/2],[nz/2]
#xs,zs = [nx/4,nx/2,3*nx/4],[0,0,0]
#xs,zs = [nx/5,2*nx/5,3*nx/5,4*nx/5],[0,0,0,0]
#xs,zs = rampint(1,15,52),fillint(0,52)
xr,zr = rampint(0,1,nx),fillint(0,nx)
ns,nr = len(xs),len(xr)
fpeak = 30.0 # Ricker wavelet peak frequency
nabsorb = 22 # absorbing boundary size
nxp,nzp = nx+2*nabsorb,nz+2*nabsorb
np = min(16,ns) # number of parallel sources

# TODO
# Test adjoint with TransformQuadratic

def main(args):
  #goAcousticData()
  goBornData()
  #adjointTest()
  #adjointTestMultiSource()
  #adjointTestMultiSourceParallel()
  #goInversion()

def goInversion():
  #s = getLayeredModel()
  s = getMarmousi(stride)
  s0,rx = makeBornModel(s)

  # Allocate wavefields.
  print "allocating..."
  u = SharedFloat4(nxp,nzp,nt,np)
  a = SharedFloat4(nxp,nzp,nt,np)

  # Sources and receivers.
  source = BornWaveOperatorS.getSourceArray(ns)
  receiver = BornWaveOperatorS.getReceiverArray(ns)
  for isou in range(ns):
    source[isou] = Source.RickerSource(xs[isou],zs[isou],dt,fpeak)
    receiver[isou] = Receiver(xr,zr,nt)

  # Born operator.
  born = BornWaveOperatorS(s0,dx,dt,nabsorb,u,a)

  # RHS vector.
  print "computing RHS..."
  rr = like(rx)
  born.applyHessian(source,receiver,rx,rr)
  vb = VecArrayFloat2(rr)

  # Solution vector.
  ry = zerofloat(nx,nz)
  vx = VecArrayFloat2(ry)

  # CG solver.
  print "solving..."
  niter = 20
  ma = CgOperator(born,source,receiver)
  cg = CgSolver(0.0,niter);
  cg.solve(ma,vb,vx);

  plot(rx,sperc=99.5)
  plot(ry,sperc=99.5)

class CgOperator(CgSolver.A):
  def __init__(self,born,source,receiver):
    self.born = born
    self.source = source
    self.receiver = receiver
  def apply(self,vx,vy):
    x = vx.getArray()
    y = vy.getArray()
    self.born.applyHessian(self.source,self.receiver,x,y)

class QsOperator(LinearTransform):
  def __init__(self):
    pass

def goAcousticData():
  s = getLayeredModel()
  #s = getMarmousi(stride)
  awo = AcousticWaveOperator(s,dx,dt,nabsorb)
  source = Source.RickerSource(xs[0],zs[0],dt,fpeak)
  receiver = Receiver(xr,zr,nt)
  u = zerofloat(nxp,nzp,nt)
  sw = Stopwatch(); sw.start()
  awo.applyForward(source,receiver,u)
  print sw.time()
  d = receiver.getData()
  for it in range(900):
    for ir in range(nr):
      d[ir][it] = 0.0
  plot(d,cmap=gray,cmin=-0.35,cmax=0.35,title="data")
  #plot(d,cmap=gray,sperc=100.0,title="data")
#  plot(s,cmap=jet,title="slowness (s/km)")
#  #display(u,perc=99.0,title="wavefield")

def goBornData():
  s = getLayeredModel()
  #s = getMarmousi(stride)
  s0,s1 = makeBornModel(s)
  bwo = BornWaveOperator(
    s,dx,dt,nabsorb)
  receiver = Receiver(xr,zr,nt)
  u = zerofloat(nxp,nzp,nt)
  source = Source.Gaussian4Source(xs[0],zs[0],dt,fpeak)
  sw = Stopwatch(); sw.start()
  bwo.applyForward(source,u,s1,receiver)
  print sw.time()
  d = receiver.getData()
  #plot(d,cmap=gray,cmin=-0.35,cmax=0.35,title="data")
  plot(d,cmap=gray,sperc=99.5,title="data")
  plot(s0,cmap=jet,title="background slowness (s/km)")
  plot(s1,sperc=100.0,title="reflectivity")

def makeBornModel(s,sigma0=0.050,sigma1=None):
  """
  sigma0: smoothing for background model
  sigma1: smoothing for perturbation
  """
  s0,s1 = like(s),like(s)
  ref0 = RecursiveExponentialFilter(sigma0/dx)
  ref0.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref0.apply(s,s0)
  t = copy(s)
  if sigma1 is None:
    sigma1 = sigma0
  ref1 = RecursiveExponentialFilter(sigma1/dx)
  ref1.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref1.apply(s,t)
  sub(s,t,s1)
  r = sub(mul(s,s),mul(t,t))
  GaussianTaper.apply1(0.25,r,r)
  return s0,r

def like(x):
  return zerofloat(len(x[0]),len(x))

def adjointTest():
  s = fillfloat(0.25,nx,nz)
  add(mul(sub(randfloat(Random(0),nx,nz),0.5),0.05),s,s)
  random = Random(01234)
  #random = Random()
  b = randfloat(random,nxp,nzp,nt); sub(b,0.5,b) # background wavefield
  bwo = BornWaveOperator(s,dx,dt,nabsorb)
  ra = randfloat(random,nx,nz) # random reflectivity
  rb = zerofloat(nx,nz) 
  db = randfloat(random,nt,nr) # random data
  reca = Receiver(xr,zr,nt)
  recb = Receiver(xr,zr,db)
  bwo.applyForward(b,ra,reca)
  bwo.applyAdjoint(b,zerofloat(nxp,nzp,nt),recb,rb)
  sum1 = dot(ra,rb)
  sum2 = dot(reca.getData(),recb.getData())
  print "adjoint test:",AcousticWaveOperator.compareDigits(sum1,sum2)
  print sum1
  print sum2

def adjointTestMultiSource():
  nsou = 3
  #random = Random()
  random = Random(01234)
  s = fillfloat(0.25,nx,nz)
  add(mul(sub(randfloat(Random(0),nx,nz),0.5),0.05),s,s)
  b = AcousticWaveOperator.randfloat(random,nxp,nzp,nt,nsou)
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
  print "adjoint test:",AcousticWaveOperator.compareDigits(sum1,sum2)
  print sum1
  print sum2

def adjointTestMultiSourceParallel():
  nsou = 4 # number of sources
  psou = 2 # number of sources to compute in parallel
  #random = Random()
  random = Random(01234)
  s = fillfloat(0.25,nx,nz)
  add(mul(sub(randfloat(Random(0),nx,nz),0.5),0.05),s,s)
  u = SharedFloat4(AcousticWaveOperator.randfloat(random,nxp,nzp,nt,psou))
  a = SharedFloat4(nxp,nzp,nt,psou)
  born = BornWaveOperatorS(s,dx,dt,nabsorb,u,a)

  # Forward
  ra = randfloat(random,nx,nz) # random reflectivity
  da = zerofloat(nt,nr,nsou)
  reca = BornWaveOperatorS.getReceiverArray(nsou)
  for isou in range(nsou):
    reca[isou] = Receiver(xr,zr,da[isou])
  born.applyForward(ra,reca)

  # Adjoint
  rb = zerofloat(nx,nz)
  db = randfloat(random,nt,nr,nsou) # random data
  recb = BornWaveOperatorS.getReceiverArray(nsou)
  for isou in range(nsou):
    recb[isou] = Receiver(xr,zr,db[isou])
  born.applyAdjoint(recb,rb)

  # Dot product
  sum1 = dot(ra,rb)
  sum2 = dot(da,db)
  print "adjoint test:",AcousticWaveOperator.compareDigits(sum1,sum2)
  print sum1
  print sum2

def dot(u,a):
  return AcousticWaveOperator.dot(u,a)

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

#############################################################################
# plotting

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,cmap=gray,perc=100.0,sperc=None,cmin=0.0,cmax=0.0,title=None):
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

def display(image,cmap=gray,cmin=0,cmax=0,perc=100,title=None):
  world = World()
  addImageToWorld(world,image,cmap,cmin,cmax,perc)
  makeFrame(world,title)

def addImageToWorld(world,image,cmap=gray,cmin=0,cmax=0,perc=100):
  ipg = ImagePanelGroup(image)
  ipg.setColorModel(cmap)
  #ipg.setSlices(k1,k2,k3)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
  if perc<100:
    ipg.setPercentiles(100-perc,perc)
  world.addChild(ipg)
  return ipg

def addColorBar(frame,label):
  cbar = ColorBar(label)
  cbar.setFont(cbar.getFont().deriveFont(64.0))
  frame.add(cbar,BorderLayout.EAST)
  #frame.viewCanvas.setBackground(frame.getBackground())
  return cbar

def makeFrame(world,title=None):
  frame = SimpleFrame(world)
  if title:
    frame.setTitle(title)
  view = frame.getOrbitView()
  view.setAxesScale(1.0,10.0,10.0)
  view.setScale(1.2)
  #view.setAzimuth(250)
  #view.setElevation(50)
  #frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(800,600)
  frame.setVisible(True)
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
