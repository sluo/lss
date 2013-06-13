#############################################################################
# Wavefield modeling demo

from imports import *

#############################################################################

sz = Sampling(201,0.016,0.0)
sx = Sampling(201,0.016,0.0)
st = Sampling(2001,0.0012,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
kzs,kxs = [0],[nx/2] # source location
kxr,kzr = rampint(0,1,nx),fillint(0,nx) # receiver locations
ns,nr = len(kxs),len(kxr)
fpeak = 10.0 # Ricker wavelet peak frequency

def main(args):
  s = getLayeredModel() # slowness model
  u = zerofloat(nz,nx,nt) # wavefield
  d = zerofloat(nt,nr) # data
  wave = Wavefield(sz,sx,st)
  wave.modelAcousticDataAndWavefield(
    Wavefield.RickerSource(fpeak,kzs[0],kxs[0]),
    Wavefield.Receiver(kzr,kxr),
    s,d,u)
  display(u,perc=99.0,title="wavefield")
  plot(d,cmap=gray,perc=99.0,title="data")
  plot(s,cmap=jet,title="slowness (s/km)")

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
