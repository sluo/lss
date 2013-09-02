#############################################################################
# Wavefield modeling using AcousticWaveOperator

from imports import *

#############################################################################

#sz = Sampling(11,0.016,0.0)
#sx = Sampling(12,0.016,0.0)
#st = Sampling(13,0.0012,0.0)
sz = Sampling(201,0.016,0.0)
sx = Sampling(203,0.016,0.0)
st = Sampling(2001,0.0012,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
xs,zs = dx*(nx/2),0.0 # source location
xr,zr = rampfloat(0.0,dx,nx),fillfloat(0.0,nx)
fpeak = 10.0 # Ricker wavelet peak frequency
nabsorb = 12 # absorbing boundary size
nxp,nzp = nx+2*nabsorb,nz+2*nabsorb

def main(args):
  #goForward()
  goMigration()
  #adjointTest()

def goForward():
  s = getLayeredModel() # slowness model
  u = zerofloat(nzp,nxp,nt) # wavefield
  awo = AcousticWaveOperator(s,dx,dt,nabsorb)
  source = AcousticWaveOperator.RickerSource(xs,zs,dt,fpeak)
  receiver = AcousticWaveOperator.Receiver(xr,zr,nt)
  sw = Stopwatch(); sw.start()
  awo.applyForward(source,receiver,u)
  print sw.time()
  d = receiver.getData()
  plot(s,cmap=jet,title="slowness (s/km)")
  plot(d,cmap=gray,perc=99.0,title="data")
  #display(u,perc=99.0,title="wavefield")

def goMigration():
  s = getLayeredModel() # slowness model
  u = zerofloat(nzp,nxp,nt) # forward wavefield
  a = zerofloat(nzp,nxp,nt) # adjoint wavefield
  receiver = AcousticWaveOperator.Receiver(xr,zr,nt)
  awo = AcousticWaveOperator(s,dx,dt,nabsorb)
  awo.applyForward(AcousticWaveOperator.RickerSource(xs,zs,dt,fpeak),receiver,u)
  awo.applyAdjoint(AcousticWaveOperator.ReceiverSource(receiver),a)
  r = AcousticWaveOperator.collapse(u,a,nabsorb)
  for ix in range(nx):
    for iz in range(20):
      r[ix][iz] = 0.0
  print sum(r)
  plot(r)

def adjointTest():
  s = fillfloat(0.5,nz,nx)
  #add(mul(add(randfloat(Random(01),nz,nx),-0.5),0.05),s,s)
  random = Random(0123)
  #random = Random()
  ua = randfloat(random,nzp,nxp,nt)
  ub = randfloat(random,nzp,nxp,nt)
  va = copy(ua)
  vb = copy(ub)
  awo = AcousticWaveOperator(s,dx,dt,nabsorb)
#  awo.applyForward(AcousticWaveOperator.WavefieldSource(ua),ua)
#  awo.applyAdjoint(AcousticWaveOperator.WavefieldSource(ub),ub)
#  print sum(ua)
#  print sum(ub)
#  sum1 = dot(ua,vb)
#  sum2 = dot(ub,va)
#  print "adjoint test:",AcousticWaveOperator.compareDigits(sum1,sum2)
#  print sum1
#  print sum2

def dot(u,a):
  n1 = len(u[0][0])
  n2 = len(u[0])
  n3 = len(u)
  sum = 0.0
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        sum += u[i3][i2][i1]*a[i3][i2][i1]
  return sum

def getLayeredModel():
  """Make slowness (s/km) model."""
  s = fillfloat(0.5,nz,nx)
  for ix in range(nx):
    for iz in range(nz/3,2*nz/3):
      s[ix][iz] = 0.35
    for iz in range(2*nz/3,nz):
      s[ix][iz] = 0.25
  return s


#############################################################################
# plotting

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,cmap=gray,perc=100.0,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  sp.setSize(600,600)
  if title:
    sp.addTitle(title)
  pv = sp.addPixels(x)
  pv.setColorModel(cmap)
  if perc<100.0:
    pv.setPercentiles(100.0-perc,perc)

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
