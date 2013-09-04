#############################################################################
# Wavefield modeling using AcousticWaveOperator

from imports import *

#############################################################################

sz = Sampling(11,0.016,0.0)
sx = Sampling(12,0.016,0.0)
st = Sampling(13,0.0012,0.0)
#sz = Sampling(201,0.016,0.0)
#sx = Sampling(202,0.016,0.0)
#st = Sampling(2003,0.0012,0.0)
#sz = Sampling(265,0.012,0.0); stride = 3
#sx = Sampling(767,0.012,0.0)
#st = Sampling(5001,0.0012,0.0)
sz = Sampling(199,0.016,0.0); stride = 4
sx = Sampling(576,0.016,0.0)
st = Sampling(5001,0.0012,0.0)
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
  ajointTest()

def adjointTest():
  # TODO
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
  awo.applyForward(ua,AcousticWaveOperator.WavefieldSource(sa))
  print 'time:',sw.time(); sw.restart()
  for i in range(1):
    awo.applyAdjoint(ub,AcousticWaveOperator.WavefieldSource(sb))
  sw.stop(); print 'time:',sw.time()
  print 'sum(ua)=%f'%sum(ua)
  print 'sum(ub)=%f'%sum(ub)
  sum1 = dot(ua,vb)
  sum2 = dot(ub,va)
  print "adjoint test:",AcousticWaveOperator.compareDigits(sum1,sum2)
  print sum1
  print sum2

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
