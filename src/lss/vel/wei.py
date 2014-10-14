"""
Traveltime-based wave equation inversion (Zhang, 2009)
"""
from imports import *

#############################################################################
# Parameters

sz = Sampling(126,16.0,0.0)
sx = Sampling(376,16.0,0.0)
st = Sampling(3001,0.0012,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
fpeak = 10.0 # Ricker wavelet peak frequency

#kzs,kxs = nz/2,nx/2 # source location
kzs,kxs = nz/10,nx/2 # source location
#kzr,kxr = fillint(nz-1,nx),rampint(0,1,nx)
kzr,kxr = fillint(nz/10,nx),rampint(0,1,nx)
#kzr,kxr = fillint(9*nz/10,nx),rampint(0,1,nx)

ft = 1.0 # plane wave delay

pngDir = None
#pngDir = '/Users/sluo/Desktop/'
#pngDir = '/Users/sluo/Desktop/png/'
#pngDir = '/Users/sluo/Desktop/png'+str(nt)+'/'
#############################################################################

def main(args):
  #v,c = getModel(); plot(v,cmap=jet,title='v')
  #r = getReflectorModel(v); plot(r,cmap=jet,title='r')
  goForward()
  #goMovie()
  #goForwardBorn()
  #goBornImaging()

def goForward():

  # TODO smooth model for Born modeling

  # Focus time
  tf = ft
  kt = tf/dt

  alpha = 0.0 # plane wave angle
  bornModeling = True

  # Recorded data
  v,c = getModel()
  #v = fillfloat(2500.0,nz,nx); c = copy(v)
  #source = AcousticWavefield.RickerSource(fpeak,kzs,kxs)
  source = AcousticWavefield.PlaneWaveSource(alpha,ft,fpeak,dx,kzs,kxs,v)
  aw = AcousticWavefield(sz,sx,st)
  aw.forwardPropagate(source,v)
  d = aw.getWavefield(kzr,kxr)
  aw.forwardPropagate(source,c)
  e = aw.getWavefield(kzr,kxr) # direct arrival
  sub(d,e,d) # remove direct arrival
  if bornModeling:
    r = getReflectorModel(v)
    bw = BornWavefield(sz,sx,st,r)

  plot(d,perc=100,title='data',png='d')
  plot(v,cmap=jet,title='true velocity',png='v')

  # Taper
  taper = fillfloat(1.0,nx)
  RecursiveGaussianFilter(nx/4).apply0(taper,taper)

  # Taper data before back-propagating
  for ix in range(nx):
    #mul(taper[ix],d[ix],d[ix])
    pass

  # Initial model
  #mul(0.9,v,v)
  #mul(0.9,c,c)
  for ix in range(nx/2):
    for iz in range(9*nz/10):
      v[ix][iz] *= 1.1
      c[ix][iz] *= 1.1
      #v[ix][iz] *= 0.9
      #c[ix][iz] *= 0.9
  for ix in range(nx/2,nx):
    for iz in range(9*nz/10):
      v[ix][iz] *= 0.9
      c[ix][iz] *= 0.9
  """
  r = zerofloat(nz,nx)
  r[nx/2][nz/2] = 1.0
  RecursiveGaussianFilter(16).apply00(r,r)
  div(r,max(r),r)
  mul(r,-400.0,r)
  add(r,v,v)
  """
  plot(v,cmap=jet,title='initial velocity',png='c')

  # Model without reflector for born modeling
  b = zerofloat(nz,nx)
  for ix in range(nx):
    fill(v[ix][0],b[ix])
  plot(b,cmap=jet,title='b')

  # Back-propagate delayed shot-gather
  source = AcousticWavefield.AdjointSource(dt,kzr,kxr,d)
  if bornModeling:
  #if False:
    print 'Born back-propagating'
    bw.backPropagate(source,b)
    d = bw.getWavefield(kzr,kxr)
    u = bw.getWavefield()
    mul(-1.0,d,d) # FIXME
  else:
    aw.backPropagate(source,v)
    d = aw.getWavefield(kzr,kxr)
    u = aw.getWavefield()
    aw.backPropagate(source,c)
    e = aw.getWavefield(kzr,kxr); sub(d,e,d) # remove direct arrival
    #y = aw.getWavefield(); sub(u,y,u) # remove downgoing wave
  plot(d,title='backpropagated data',png='b')

  # Penalty function
  """
  p = 300 # penalty half-width
  w = zerofloat(nt,nx)
  #for ix in range(nx):
  #  w[ix][int(kt)] = 1.0
  #RecursiveGaussianFilter(p).apply0X(w,w)
  for ix in range(nx):
    for it in range(kt-p,kt+p):
    #for it in range(0,int(4.0*kt)):
      w[ix][it] = it-kt
  RecursiveGaussianFilter(0.25*p).apply0X(w,w)
  #pow(w,3,w)
  div(w,max(abs(w)),w)
  #SimplePlot.asPoints(w[0])
  """
  w = makePenalty(alpha,v[kxs][kzs])
  RecursiveGaussianFilter(50).apply0X(w,w)
  plot(w,title='penalty',png='w')

  # Penalized back-propagated data (pseudosource)
  mul(w,d,d)
  plot(d,title='penalized backpropagated',png='p')

  # Taper pseudosource (adjoint source)
  for ix in range(nx):
    #mul(taper[ix],d[ix],d[ix])
    pass

  """
  cmin = 0.5*min(u)
  cmax = 0.5*max(u)
  #for i in range(0,nt,100):
  for i in range(nt-td-20,nt-td+20,2):
  #for i in range(td-40,td+10,2):
    plot(u[i],cmap=jet,cmin=cmin,cmax=cmax,title=str(i))
  """

  # Forward-propagate pseudosource
  jzs,jxs = fillint(kzs,nx),rampint(0,1,nx)
  source = AcousticWavefield.DefaultSource(dt,jzs,jxs,d)
  if False:
    print 'Born forward-propagating'
    bw.forwardPropagate(source,b)
    a = bw.getWavefield()
    e = bw.getWavefield(kzr,kxr)
    plot(e)
  else:
    aw.forwardPropagate(source,v)
    a = aw.getWavefield()
    #aw.forwardPropagate(source,c)
    #e = aw.getWavefield(); sub(a,e,a) # remove downgoing wave

  # Gradient
  g = correlate(u,a)
  div(g,max(abs(g)),g)
  plot(g,cmap=rwb,cmin=-0.9,cmax=0.9,title='gradient',png='g')

def makePenalty(alpha,v0):
  """Penalty function plane-wave angle alpha."""
  s = 0.36 # penalty half-width in seconds
  p = sin(DBL_PI*alpha/180.0)/v0
  def f(t):
    if -s<=t and t<=s:
      return t/s
    else:
      return 0.0
  w = zerofloat(nt,nx)
  for ix in range(nx):
    px = p*(kxs-ix)*dx
    for it in range(nt):
      t = it*dt-ft
      w[ix][it] = f(t-px)
  #plot(w)
  return w

def goMovie():
  tf = ft
  kt = tf/dt

  alpha = 0.0 # plane wave angle

  v,c = getModel()
  wavefield = AcousticWavefield(sz,sx,st)
  source = AcousticWavefield.PlaneWaveSource(alpha,tf,fpeak,dx,kzs,kxs,v)
  wavefield.forwardPropagate(source,v)
  d = wavefield.getWavefield(kzr,kxr)
  wavefield.forwardPropagate(source,c)
  e = wavefield.getWavefield(kzr,kxr) # direct arrival
  sub(d,e,d) # remove direct arrival

  plot(d,perc=100,title='data',png='d')
  plot(v,cmap=jet,title='true',png='v')

  # Taper
  taper = fillfloat(1.0,nx)
  RecursiveGaussianFilter(nx/4).apply0(taper,taper)

  # Taper data before back-propagating
  for ix in range(nx):
    mul(taper[ix],d[ix],d[ix])
    pass

  # Initial model
  for ix in range(nx/2):
    for iz in range(9*nz/10):
      v[ix][iz] *= 1.1
  for ix in range(nx/2,nx):
    for iz in range(9*nz/10):
      v[ix][iz] *= 0.9
  """
  c = zerofloat(nz,nx)
  c[nx/2][nz/2] = 1.0
  RecursiveGaussianFilter(16).apply00(c,c)
  div(c,max(c),c)
  mul(c,-400.0,c)
  add(c,v,v)
  """
  plot(v,cmap=jet,title='initial',png='c')

  source = AcousticWavefield.AdjointSource(dt,kzr,kxr,d)
  wavefield.backPropagate(source,v)
  d = wavefield.getWavefield(kzr,kxr)
  plot(d,title='backpropagated',png='b')
  u = wavefield.getWavefield()

  cmax = 0.5*max(u)
  cmin = 0.5*min(u)
  for it in range((int)(kt),3000,100):
    plot(u[it],cmap=gyr,cmax=cmax,cmin=cmin,
      title='u['+str(it)+']',png='u'+str(it))
    pass

def goBornImaging():
  v,c = getModel();
  #r = getReflectorModel(v)
  s,r = getBackgroundAndPerturbationModel(v)
  aw = AcousticWavefield(sz,sx,st)
  bw = BornWavefield(sz,sx,st,r)
  rickerSource = AcousticWavefield.RickerSource(fpeak,kzs,kxs)
  plot(v,cmap=jet,title='v',png='v')
  plot(c,cmap=jet,cmin=min(v),cmax=max(v),title='c',png='c')
  plot(r,cmap=jet,title='r',png='r')

  # Recorded data
  aw.forwardPropagate(rickerSource,v)
  d = aw.getWavefield(kzr,kxr)
  aw.forwardPropagate(rickerSource,c)
  e = aw.getWavefield(kzr,kxr) # direct arrival
  sub(d,e,d) # remove direct arrival
  plot(d,title='data')

  # Adjoint source
  taper = fillfloat(1.0,nx)
  RecursiveGaussianFilter(nx/8).apply0(taper,taper)
  for ix in range(nx):
    mul(taper[ix],d[ix],d[ix])
    pass
  adjointSource = AcousticWavefield.AdjointSource(dt,kzr,kxr,d)

  # RTM
  """
  """
  #aw.forwardPropagate(rickerSource,v); u = aw.getWavefield()
  #aw.backPropagate(adjointSource,v); a = aw.getWavefield()
  aw.forwardPropagate(rickerSource,c); u = aw.getWavefield()
  aw.backPropagate(adjointSource,c); a = aw.getWavefield()

  # First half
  """
  aw.forwardPropagate(rickerSource,c); u = aw.getWavefield()

  aw.backPropagate(adjointSource,c); a = aw.getWavefield()
  #bw.backPropagate(adjointSource,c); a = bw.getWavefield()

  #aw.backPropagate(adjointSource,c); a = aw.getWavefield()
  #aw.backPropagate(adjointSource,v); b = aw.getWavefield()
  #sub(b,a,a)
  """

  # Second half
  """
  bw.forwardPropagate(rickerSource,c); u = bw.getWavefield()
  aw.backPropagate(adjointSource,c); a = aw.getWavefield()
  """

  g = correlate(u,a)
  plot(g,perc=99,title='migrated image',png='g')

def correlate(u,a):
  """Zero-lag correlation."""
  g = zerofloat(nz,nx)
  for it in range(1,nt-1):
    #t = add(add(u[it-1],u[it+1]),mul(-2.0,u[it])) # 2nd time derivative
    #add(mul(a[it],t),g,g)
    add(mul(a[it],u[it]),g,g)
  return g

def goForwardBorn():
  v,c = getModel()
  r = getReflectorModel(v)
  wavefield = BornWavefield(sz,sx,st,r)
  wavelet = BornWavefield.SourceWavelet.RICKER;
  wavefield.forwardPropagate(wavelet,fpeak,kzs,kxs,v)
  d = wavefield.getWavefield(kzr,kxr)
  plot(d,title='d',perc=100)
  plot(v,cmap=jet,title='v')
  plot(r,cmap=jet,title='r')

def getReflectorModel(v):
  b,p = getBackgroundAndPerturbationModel(v)
  return p

def getBackgroundAndPerturbationModel(v):
  """Make a reflector model for use by BornWavefield."""
  a = sum(v)/nx/nz # average velocity
  w = a/fpeak/dx # dominant wavelength in samples
  #w = w*0.1
  #w = w*0.5
  s = copy(v)
  #RecursiveGaussianFilter(w).apply00(v,s) # smoothed velocity
  RecursiveExponentialFilter(w).apply(v,s) # smoothed velocity
  """
  print 'average =',a
  print 'wavelength =',w
  plot(s,cmap=jet,title='s')
  """
  return s,sub(v,s)

def getModel():
  return layeredModel()
  #return gaussianModel()

def gaussianModel():
  v = fillfloat(2000.0,nz,nx)
  c = copy(v)
  t = zerofloat(nz,nx)
  t[nx/2][nz/2] = 1.0
  RecursiveGaussianFilter(20.0).apply00(t,t)
  t = normalize(t)
  t = mul(t,500.0)
  v = add(t,v)
  #plot(v)
  return v,c

def layeredModel():
  v = fillfloat(2000.0,nz,nx) # velocity
  c = copy(v)
  for ix in range(nx):
    for iz in range(nz/10,9*nz/10):
      v[ix][iz] = 2000.0
      c[ix][iz] = 2000.0
    for iz in range(9*nz/10,nz):
      v[ix][iz] = 6000.0
  return v,c

def normalize(x):
  return div(x,max(abs(x)))

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
gyr = ColorMap.GRAY_YELLOW_RED
def plot(x,cmap=gray,cmin=0,cmax=0,perc=100,title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  cb = sp.addColorBar()
  cb.setWidthMinimum(80)
  if len(x[0])>nz:
    sp.setSize(1000,700)
  else:
    sp.setSize(1000,500)
  if title:
    sp.addTitle(title)
  pv = sp.addPixels(x)
  pv.setColorModel(cmap)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if perc<100:
    pv.setPercentiles(100-perc,perc)
  if png and pngDir:
    sp.paintToPng(360.0,3.0,pngDir+png+'.png')


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
  frame.setSize(1000,400)
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
