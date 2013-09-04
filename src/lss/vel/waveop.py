#############################################################################
# Wavefield modeling using AcousticWaveOperator

from imports import *

#############################################################################

#sz = Sampling(11,0.016,0.0)
#sx = Sampling(12,0.016,0.0)
#st = Sampling(13,0.0012,0.0)
#sz = Sampling(201,0.016,0.0)
#sx = Sampling(202,0.016,0.0)
#st = Sampling(2003,0.0012,0.0)
#sz = Sampling(101,0.032,0.0)
#sx = Sampling(102,0.032,0.0)
#st = Sampling(1003,0.0024,0.0)
#sz = Sampling(265,0.012,0.0); stride = 3
#sx = Sampling(767,0.012,0.0)
#st = Sampling(5001,0.0012,0.0)
sz = Sampling(199,0.016,0.0); stride = 4
sx = Sampling(576,0.016,0.0)
st = Sampling(5001,0.0012,0.0)
#sz = Sampling(159,0.020,0.0); stride = 5
#sx = Sampling(461,0.020,0.0)
#st = Sampling(2501,0.0024,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
#xs,zs = [nx/2],[0]
#xs,zs = [nx/2],[nz/2]
#xs,zs = [nx/4,nx/2,3*nx/4],[0,0,0]
#xs,zs = rampint(25,25,22),fillint(0,22)
xs,zs = rampint(2,11,53),fillint(0,53)
xr,zr = rampint(0,1,nx),fillint(0,nx)
ns,nr = len(xs),len(xr)
fpeak = 12.0 # Ricker wavelet peak frequency
nabsorb = 10 # absorbing boundary size
nxp,nzp = nx+2*nabsorb,nz+2*nabsorb

def main(args):
  #plot(getMarmousi(5))
  #goForward()
  goMigration()
  #ajointTest()
  #absorbAdjointTest()

def goForward():
  s = fillfloat(0.25,nx,nz)
  #s = getLayeredModel() # slowness model
  #s = getMarmousi(stride)
  u = zerofloat(nxp,nzp,nt) # wavefield
  awo = AcousticWaveOperator(s,dx,dt,nabsorb)
  source = AcousticWaveOperator.RickerSource(xs[0],zs[0],dt,fpeak)
  receiver = AcousticWaveOperator.Receiver(xr,zr,nt)
  sw = Stopwatch(); sw.start()
  awo.applyForward(source,receiver,u)
  print sw.time()
  d = receiver.getData()
  plot(u[500])
#  plot(s,cmap=jet,title="slowness (s/km)")
#  plot(d,cmap=gray,perc=99.0,title="data")
#  #display(u,perc=99.0,title="wavefield")

def goMigration():
  #s = getLayeredModel() # true slowness
  s = getMarmousi(stride) # true slowness
  b = refSmooth(0.1,s) # smooth slowness
  #b = copy(s) # smooth slowness
  u = zerofloat(nxp,nzp,nt) # forward wavefield
  a = zerofloat(nxp,nzp,nt) # adjoint wavefield
  r = zerofloat(nx,nz) # image
  awoS = AcousticWaveOperator(s,dx,dt,nabsorb) # awo with true slowness
  awoB = AcousticWaveOperator(b,dx,dt,nabsorb) # awo with smooth slowness
  receiver = AcousticWaveOperator.Receiver(xr,zr,nt)
  sw = Stopwatch(); sw.start()
  for isou in range(ns):
    print 'isou=%d'%isou
    zero(u); zero(a);
    awoS.applyForward(
      AcousticWaveOperator.RickerSource(xs[isou],zs[isou],dt,fpeak),receiver)
    awoB.applyForward(
      AcousticWaveOperator.RickerSource(xs[isou],zs[isou],dt,fpeak),u)
    awoB.applyAdjoint(AcousticWaveOperator.ReceiverSource(receiver),a)
    add(AcousticWaveOperator.collapse(u,a,nabsorb),r,r)
  sw.stop(); print 'time:',sw.time()
  #print sum(u); print sum(a); print sum(r)
  for iz in range(int(0.2/dx)):
    for ix in range(nx):
      r[iz][ix] = 0.0
      pass
  plot(s)
  plot(b)
  plot(r,sperc=98.0,cmin=0,cmax=0)

def refSmooth(sigma,x):
  y = like(x)
  ref = RecursiveExponentialFilter(sigma/dx)
  #ref = RecursiveExponentialFilter(1.0)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref.apply(x,y)
  return y

def adjointTest():
  s = fillfloat(0.25,nx,nz)
  add(mul(add(randfloat(Random(0),nx,nz),-0.5),0.05),s,s)
  random = Random(01234)
  #random = Random()
  sa = add(randfloat(random,nxp,nzp,nt),0.0)
  sb = add(randfloat(random,nxp,nzp,nt),0.0)
  va = copy(sa)
  vb = copy(sb)
  ua = zerofloat(nxp,nzp,nt)
  ub = zerofloat(nxp,nzp,nt)
  awo = AcousticWaveOperator(s,dx,dt,nabsorb)
  sw = Stopwatch(); sw.start()
  awo.applyForward(AcousticWaveOperator.WavefieldSource(sa),ua)
  print 'time:',sw.time(); sw.restart()
  for i in range(1):
    awo.applyAdjoint(AcousticWaveOperator.WavefieldSource(sb),ub)
  sw.stop(); print 'time:',sw.time()
  print 'sum(ua)=%f'%sum(ua)
  print 'sum(ub)=%f'%sum(ub)
  sum1 = dot(ua,vb)
  sum2 = dot(ub,va)
  print "adjoint test:",AcousticWaveOperator.compareDigits(sum1,sum2)
  print sum1
  print sum2

def absorbAdjointTest():
  s = fillfloat(0.25,nx,nz)
  add(mul(add(randfloat(Random(0),nx,nz),-0.5),0.05),s,s)
  random = Random(01234)
  #random = Random()
  ua = add(randfloat(random,nxp,nzp,3),0.0)
  ub = add(randfloat(random,nxp,nzp,3),0.0)
  va = copy(ua)
  vb = copy(ub)
  awo = AcousticWaveOperator(s,dx,dt,nabsorb)
  sw = Stopwatch(); sw.start()
  awo.forwardAbsorb(ua[0],ua[1],ua[2])
  awo.adjointAbsorb(ub[0],ub[1],ub[2])
  sw.stop(); print 'time:',sw.time()
  print 'sum(ua)=%f'%sum(ua)
  print 'sum(ub)=%f'%sum(ub)
  sum1 = dot(ua,vb)
  sum2 = dot(ub,va)
  print "adjoint test:",AcousticWaveOperator.compareDigits(sum1,sum2)
  print sum1
  print sum2

def dot(u,a):
  #return pdot(u,a)
  return jdot(u,a)

def jdot(u,a):
  return AcousticWaveOperator.dot(u,a)

def pdot(u,a):
  n1 = len(u[0][0])
  n2 = len(u[0])
  n3 = len(u)
  sum = 0.0
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        sum += u[i3][i2][i1]*a[i3][i2][i1]
  return sum

def dotTest():
  n1,n2,n3 = 1,2,3
  random = Random(012345)
  u = randfloat(random,n1,n2,n3)
  a = randfloat(random,n1,n2,n3)
  pd = pdot(u,a)
  jd = jdot(u,a)
  print pd
  print jd
  print abs(pd-jd)

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
