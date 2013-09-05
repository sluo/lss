#############################################################################
# Wavefield modeling using AcousticWaveOperator

from imports import *

#############################################################################

#sz = Sampling(11,0.016,0.0)
#sx = Sampling(12,0.016,0.0)
#st = Sampling(13,0.0012,0.0)
sz = Sampling(201,0.016,0.0)
sx = Sampling(202,0.016,0.0)
st = Sampling(2003,0.0012,0.0)
#sz = Sampling(265,0.012,0.0); stride = 3
#sx = Sampling(767,0.012,0.0)
#st = Sampling(5001,0.0012,0.0)
#sz = Sampling(199,0.016,0.0); stride = 4
#sx = Sampling(576,0.016,0.0)
#st = Sampling(5001,0.0012,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
xs,zs = [nx/2],[0]
#xs,zs = [nx/2],[nz/2]
#xs,zs = [nx/4,nx/2,3*nx/4],[0,0,0]
#xs,zs = rampint(25,25,22),fillint(0,22)
#xs,zs = rampint(2,11,53),fillint(0,53)
xr,zr = rampint(0,1,nx),fillint(0,nx)
ns,nr = len(xs),len(xr)
fpeak = 12.0 # Ricker wavelet peak frequency
nabsorb = 12 # absorbing boundary size
nxp,nzp = nx+2*nabsorb,nz+2*nabsorb

def main(args):
  #goAcousticData()
  goBornData()
  #adjointTest()
  #adjointTestMulti()

def goAcousticData():
  s = getLayeredModel()
  #s = getMarmousi(stride)
  awo = AcousticWaveOperator(s,dx,dt,nabsorb)
  source = AcousticWaveOperator.RickerSource(xs[0],zs[0],dt,fpeak)
  receiver = AcousticWaveOperator.Receiver(xr,zr,nt)
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
    #AcousticWaveOperator.RickerSource(xs[0],zs[0],dt,fpeak),
    AcousticWaveOperator.Gaussian4Source(xs[0],zs[0],dt,fpeak),
    s,dx,dt,nabsorb)
  receiver = AcousticWaveOperator.Receiver(xr,zr,nt)
  sw = Stopwatch(); sw.start()
  bwo.applyForward(s1,receiver)
  print sw.time()
  d = receiver.getData()
  plot(d,cmap=gray,cmin=-0.35,cmax=0.35,title="data")
  #plot(d,cmap=gray,sperc=100.0,title="data")
  #plot(s0,cmap=jet,title="background slowness (s/km)")
  #plot(s1,sperc=100.0,title="reflectivity")

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
  bwo = BornWaveOperator(b,s,dx,dt,nabsorb)
  ra = randfloat(random,nx,nz) # random reflectivity
  rb = zerofloat(nx,nz) 
  da = zerofloat(nt,nr)
  db = randfloat(random,nt,nr) # random data
  reca = AcousticWaveOperator.Receiver(xr,zr,nt)
  recb = AcousticWaveOperator.Receiver(xr,zr,nt)
  copy(da,reca.getData())
  copy(db,recb.getData())
  bwo.applyForward(ra,reca)
  bwo.applyAdjoint(zerofloat(nxp,nzp,nt),recb,rb)
  sum1 = dot(ra,rb)
  sum2 = dot(da,db)
  sum2 = dot(reca.getData(),recb.getData())
  print "adjoint test:",AcousticWaveOperator.compareDigits(sum1,sum2)
  print sum1
  print sum2

def adjointTestMulti():
  s = fillfloat(0.25,nx,nz)
  add(mul(sub(randfloat(Random(0),nx,nz),0.5),0.05),s,s)
  random = Random(01234)
  #random = Random()
  ra = randfloat(random,nx,nz) # random reflectivity
  rb = zerofloat(nx,nz) 

  # TODO
  #for isou in range(1):
  #b = randfloat(random,nxp,nzp,nt); sub(b,0.5,b) # background wavefield
  #bwo = BornWaveOperator(b,s,dx,dt,nabsorb)
  #da = zerofloat(nt,nr)
  #db = randfloat(random,nt,nr) # random data
  #reca = AcousticWaveOperator.Receiver(xr,zr,nt)
  #recb = AcousticWaveOperator.Receiver(xr,zr,nt)
  #copy(da,reca.getData())
  #copy(db,recb.getData())
  #bwo.applyForward(ra,reca)
  #bwo.applyAdjoint(zerofloat(nxp,nzp,nt),recb,rb)
  #sum1 = dot(ra,rb)
  #sum2 = dot(da,db)
  #sum2 = dot(reca.getData(),recb.getData())
  #print "adjoint test:",AcousticWaveOperator.compareDigits(sum1,sum2)
  #print sum1
  #print sum2

def dot(u,a):
  return AcousticWaveOperator.dot(u,a)

def getLayeredModel():
  """Make slowness (s/km) model."""
  s = fillfloat(0.5,nx,nz)
  for iz in range(nz/3,2*nz/3):
    for ix in range(nx):
      s[iz][ix] = 0.35
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
