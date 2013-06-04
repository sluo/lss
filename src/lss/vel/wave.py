"""
Acoustic wavefield modeling
"""
from imports import *

#############################################################################
# Parameters

#sz = Sampling(201,16.0,0.0)
#sx = Sampling(201,16.0,0.0)
sz = Sampling(201,0.016,0.0)
sx = Sampling(201,0.016,0.0)
st = Sampling(2001,0.0012,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
v = fillfloat(2000.0,nz,nx) # velocity
kzs,kxs = nz/2,nx/2 # source location
fpeak = 5.0 # Ricker wavelet peak frequency

#############################################################################

def main(args):
  #wavefield = AcousticWavefield(sz,sx,st)
  #wavelet = AcousticWavefield.SourceWavelet.RICKER;
  #wavefield.forwardPropagate(wavelet,fpeak,kzs,kxs,v)
  #u = wavefield.getWavefield()
  #display(u)

  #wavefield = AcousticWavefieldY(sz,sx,st)
  #source = AcousticWavefieldY.RickerSource(fpeak,kzs,kxs)
  #wavefield.propagate(source,v)
  #u = wavefield.getWavefield()

  s = div(1000.0,v) # slowness (s/km)
  u = zerofloat(nz,nx,nt)
  wave = Wavefield(sz,sx,st)
  wave.modelAcousticWavefield(
    Wavefield.RickerSource(fpeak,kzs,kxs),
    s,u)

  display(u)
  plot(u[900]);
  plot(u[1000]);
  plot(u[1800]);
  print sum(abs(u[1800]));

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,cmap=gray,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  sp.setSize(600,600)
  if title:
    sp.addTitle(title)
  pv = sp.addPixels(x)
  pv.setColorModel(cmap)

def display(image,cmap=gray,cmin=0,cmax=0,perc=100,name=None):
  world = World()
  addImageToWorld(world,image,cmap,cmin,cmax,perc)
  makeFrame(world,name)

def display2(image1,image2,cmap1=gray,cmap2=gray,name=None):
  world = World()
  addImageToWorld(world,image1,cmap1)
  addImageToWorld(world,image2,cmap2)
  makeFrame(world,name)

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

def makeFrame(world,name=None):
  frame = SimpleFrame(world)
  if name:
    frame.setTitle(name)
  view = frame.getOrbitView()
  view.setAxesScale(1.0,1.0,1.0)
  view.setScale(2.5)
  #view.setAzimuth(250)
  #view.setElevation(50)
  #frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(800,600)
  frame.setVisible(True)
  return frame

#############################################################################
# Do everything on Swing thread.
import sys,time
class RunMain(Runnable):
  def run(self):
    start = time.time()
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02ds'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
