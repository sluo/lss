"""
Unfaulting F3D with Dave's fault throws
"""
from imports import *

subset = 'a'
#subset = 'b'

pngDir = None
#pngDir = "/Users/sluo/Desktop/"
#pngDir = "/Users/sluo/Desktop/png/"
dataDir = '/data/seis/f3/sub'+subset+'/'
n1,n2,n3 = 90,221,220
"""
d1,d2,d3 = 0.004,0.025,0.025
f1,f2,f3 = 1.204 if subset=='a' else 3.044,1.25,2.5
s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
"""
s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
k1,k2,k3 = 174,400,500
azimuth=300; elevation=50 # for 3D views
h1 = 26 if subset=="a" else 53 # horizons
h1s = [20,40,60]

#############################################################################

def main(args):
  show()
  #goInterpolateThrows()
  #goCheckInterpolation()
  #goUnfault()
  #goFlatten()
  #goMapping()

def show():
  """
  seismic: input, unfaulted, and flattened images
  nearest: nearest neighbor interpolation of fault throws
  blended: blended neighbor interpolation of fault throws
  shiftsr: flattening shifts for unfaulted image
  shiftss: composite shifts
  mapping: mapping x(u) from composite shifts
  horizon: horizons
  """
  seismic = False
  nearest = False
  blended = False
  shiftsr = False
  shiftss = False
  mapping = False
  horizon = True
  if seismic:
    g = read("g")
    h = read("h")
    f = read("f")
    #display(g,name="input") # input
    #display(h,name="unfaulted") # unfaulted
    #display(f,name="flat") # unfaulted and flattened
    display2(g,h,name="input and unfaulted")
    display2(h,f,name="unfaulted and flat")
  if nearest or blended:
    t1 = read('t1')
    t2 = read('t2')
    t3 = read('t3')
    display(t1,cmap=jet,name='t1')
    display(t2,cmap=jet,name='t2')
    display(t3,cmap=jet,name='t3')
  if nearest:
    p1 = read('p1')
    p2 = read('p2')
    p3 = read('p3')
    display(mul(-1,p1),cmap=jet,name='p1')
    display(mul(-1,p2),cmap=jet,name='p2')
    display(mul(-1,p3),cmap=jet,name='p3')
  if blended:
    q1 = read('q1')
    q2 = read('q2')
    q3 = read('q3')
    display(mul(-1,q1),cmap=jet,name='q1')
    display(mul(-1,q2),cmap=jet,name='q2')
    display(mul(-1,q3),cmap=jet,name='q3')
    #display(mul(-1,q1),cmap=prism,perc=99.5,name='q1')
    #display(mul(-1,q2),cmap=prism,perc=99.5,name='q2')
    #display(mul(-1,q3),cmap=prism,perc=99.5,name='q3')
  if shiftsr or horizon:
    r1,r2,r3 = read('r1'),read('r2'),read('r3')
  if mapping or horizon:
    x1,x2,x3 = read('x1'),read('x2'),read('x3')
  if shiftsr:
    display(r1,cmap=jet,name='r1')
    display(r2,cmap=jet,name='r2')
    display(r3,cmap=jet,name='r3')
  if shiftss:
    s1,s2,s3 = read('s1'),read('s2'),read('s3')
    display(s1,cmap=jet,name='s1')
    display(s2,cmap=jet,name='s2')
    display(s3,cmap=jet,name='s3')
  if mapping:
    #x1,x2,x3 = sin(x1),sin(x2),sin(x3)
    display(x1,cmap=jet,name='x1')
    display(x2,cmap=jet,name='x2')
    display(x3,cmap=jet,name='x3')
  if horizon:
    g = read('g')
    r = array(r1,r2,r3)
    x = array(x1,x2,x3)
    world1 = World()
    ipg = addImageToWorld(world1,g)
    #ipg.setSlices(n1-1,0,n3-1)
    ipg.setSlices(n1-1,7,n3-8)

    """
    xyz = getHorizonVertices(h1,r,x)
    tg = QuadGroup(True,xyz)
    tg.setColor(Color.ORANGE)
    """
    xyz,rgb = getHorizonVertices(h1,r,x)
    tg = QuadGroup(True,xyz,rgb)

    world1.addChild(tg)
    frame = makeFrame(world1,name='horizon')
    colorbar = addColorBar(frame,'Amplitude')
    colorbar.setWidthMinimum(250)
    ipg.addColorMapListener(colorbar)
    """
    world2 = World()
    addImageToWorld(world2,g).setSlices(n1-1,7,n2-8)
    colors = [Color.BLUE,Color.RED,Color.CYAN,Color.YELLOW]
    for i in range(len(h1s)):
      xyz = getHorizonVertices(h1s[i],r,x)
      tg = QuadGroup(True,xyz)
      tg.setColor(colors[i])
      world2.addChild(tg)
    makeFrame(world2,name='horizons')
    """

def goMapping():
  """Computes the composite mapping x(u) = u-s(u)."""
  s1,s2,s3 = read('s1'),read('s2'),read('s3')
  x1,x2,x3 = zeros(n1,n2,n3),zeros(n1,n2,n3),zeros(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        x1[i3][i2][i1] = i1-s1[i3][i2][i1]
        x2[i3][i2][i1] = i2-s2[i3][i2][i1]
        x3[i3][i2][i1] = i3-s3[i3][i2][i1]
  #display(x1,cmap=jet,name='x1')
  #display(x2,cmap=jet,name='x2')
  #display(x3,cmap=jet,name='x3')
  write('x1',x1)
  write('x2',x2)
  write('x3',x3)
  
def getHorizonVertices(k1,r=None,x=None):
  """Computes horizon vertex coordinates for use in QuadGroup."""
  # Coordinates: x = x(x)
  if r:
    r1,r2,r3 = r[0],r[1],r[2]
  else:
    r1,r2,r3 = read('r1'),read('r2'),read('r3')
    r = array(r1,r2,r3)
  if x:
    x1,x2,x3 = x[0],x[1],x[2]
  else:
    x1,x2,x3 = read('x1'),read('x2'),read('x3')
    x = array(x1,x2,x3)
  # Fault location: w = w(x)
  t1,t2,t3,z1,z2,z3 = getGriddedThrows()
  e = Faulter.getGriddedThrows(0.0,t1,t2,t3,z1,z2,z3)
  z1,z2,z3 = e[3],e[4],e[5] # live value coords
  w = zerofloat(n1,n2,n3)
  for i in range(len(z1)):
    i1,i2,i3 = int(z1[i]),int(z2[i]),int(z3[i])
    w[i3][i2][i1] = 1.0
  w = FlattenerUtil.applyShiftsR(w,r) # w(u) = w(u-r(u))
  display(w,cmap=jet,cmin=0,cmax=1)
  # Slice
  x1 = slice23(k1,x1) # x1(u=k1)
  x2 = slice23(k1,x2) # x2(u=k1)
  x3 = slice23(k1,x3) # x3(u=k1)
  w = slice23(k1,w)   # w(u=k1)
  #return FlattenerUtil.makeQuadVertices(x1,x2,x3,w)
  xyzrgb = FlattenerUtil.makeQuadVerticesAndRgbFloats(x1,x2,x3,w)
  return xyzrgb[0],xyzrgb[1]

def goFlatten():
  """Flattens the unfaulted image."""
  smoothImage = True # sos before lof
  findSlopes = True
  findShifts = True
  timer = Stopwatch()
  h = read('h')
  #h = read('f0')
  u1,u2,u3,ep = like(h),like(h),like(h),like(h)
  if findSlopes:
    lof = LocalOrientFilter(8.0,2.0)
    timer.restart(); print 'slopes...'
    if smoothImage: # Structure-oriented smoothing
      hs = like(h)
      et = lof.applyForTensors(h)
      p0 = 1.0 # amplitude
      p1 = 1.0 # linearity
      p2 = 10.0 # planarity
      et.invertStructure(p0,p1,p2)
      c = 2.0
      LocalSmoothingFilter().apply(et,c,h,hs)
    lof.applyForNormalPlanar(hs,u1,u2,u3,ep)
    timer.stop(); print 'slopes in %.2fs'%timer.time()
    write("u1",u1)
    write("u2",u2)
    write("u3",u3)
    write("ep",ep)
  flat = FlattenerRT(6.0,6.0)
  if findShifts:
    read("u1",u1)
    read("u2",u2)
    read("u3",u3)
    read("ep",ep)
    pow(ep,8.0,ep)
    p = array(u1,u2,u3,ep)
    timer.restart(); print 'shifts...'
    r = flat.findShifts(p)
    timer.stop(); print 'shifts in %.2fs'%timer.time()
    q = array(read("q1"),read("q2"),read("q3"))
    s = shiftShifts(r,q) # composite shifts
    #f = applyShifts(g,s)
    f = applyShifts(h,r) # faults look cleaner
    write("r1",r[0])
    write("r2",r[1])
    write("r3",r[2])
    write("s1",s[0])
    write("s2",s[1])
    write("s3",s[2])
    write("f",f)
    #write("f0",f)
  goMapping() # composite mapping x(u) = u-s(u) 
  f = read('f')
  display(f,name="f")

def shiftShifts(r,q):
  """Combine interpolated fault throws and flattening shifts
     r = flattening shifts, q = unfaulting shifts
  """
  s1 = add(FlattenerUtil.applyShiftsRLinear(q[0],r),r[0])
  s2 = add(FlattenerUtil.applyShiftsRLinear(q[1],r),r[1])
  s3 = add(FlattenerUtil.applyShiftsRLinear(q[2],r),r[2])
  return array(s1,s2,s3)

def goUnfault():
  """Unfault the image."""
  smoothFaults = True
  g,q1,q2,q3 = read("g"),read("q1"),read("q2"),read("q3")
  #h = FlattenerUtil.applyShiftsR(g,[q1,q2,q3])
  h = applyShifts(g,[q1,q2,q3])
  if smoothFaults:
    h0 = copy(h)
    t1,t2,t3,x1,x2,x3 = getGriddedThrows()
    s = zerofloat(n1,n2,n3)
    for i in range(len(x1)):
      i1,i2,i3 = int(x1[i]),int(x2[i]),int(x3[i])
      r1,r2,r3 = t1[i],t2[i],t3[i]
      r = sqrt(r1*r1+r2*r2+r3*r3)
      #s[i3][i2][i1] = 1.0*r
      s[i3][i2][i1] = 1.0
    RecursiveGaussianFilter(2.0).apply000(s,s)
    t3 = LocalOrientFilter(8.0).applyForTensors(h)
    t3.invertStructure(1.0,1.0,1.0)
    LocalSmoothingFilter().apply(t3,2.0,s,h,h)
    display(s,name="scale")
    display2(h0,h,name="h0 and h")
  else:
    display(h,name="h")
  display(g,name="g")
  #display(q1,cmap=jet,name="q1")
  #display(q2,cmap=jet,name="q2")
  #display(q3,cmap=jet,name="q3")
  write("h",h)
  return h

def applyShifts(x,r):
  """Undo gain before applying shifts."""
  def slog(x):
    return mul(sgn(x),sub(exp(abs(x)),1.0))
  def sexp(x):
    return mul(sgn(x),log(add(abs(x),1.0)))
  return sexp(FlattenerUtil.applyShiftsR(slog(x),r))

def goInterpolateThrows():
  """Interpolate throw vectors."""
  t1,t2,t3,x1,x2,x3 = getGriddedThrows()
  sw = Stopwatch()
  print "gridding..."
  sw.restart()
  q1,p1 = gridBlended(t1,x1,x2,x3)
  q2,p2 = gridBlended(t2,x1,x2,x3)
  q3,p3 = gridBlended(t3,x1,x2,x3)
  sw.stop()
  print "done:",sw.time()
  display(q1,cmap=jet,name="q1")
  display(q2,cmap=jet,name="q2")
  display(q3,cmap=jet,name="q3")
  write("p1",p1)
  write("p2",p2)
  write("p3",p3)
  write("q1",q1)
  write("q2",q2)
  write("q3",q3)
  return q1,q2,q3

def gridBlended(f1,x1,x2,x3):
  """Blended neighbor interpolation."""
  grid = Grid3(f1,x1,x2,x3)
  d,q = zerofloat(n1,n2,n3),zerofloat(n1,n2,n3)
  p = grid.gridNearest(s1,s2,s3,d)
  display(d,cmap=jet,name="d")
  Faulter.adjustDistances(1.0,d)
  blend = BlendedGridder3()
  blend.setBlendingKernel(
    LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D21))
  #blend.setTimeMax(20.0)
  blend.gridBlended(d,p,q)
  return q,p

def adjustDistances(d):
  """Flood distances to correct for finite difference stencil."""
  #sqrt3 = sqrt(3.0)
  sqrt3 = 1.0 # if using D21 stencil
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        d[i3][i2][i1] = max(0.0,d[i3][i2][i1]-sqrt3)
  return d

def getGriddedThrows():
  null = -0.12345
  t1,t2,t3 = read("t1"),read("t2"),read("t3")
  f1 = SimpleGridder3.getGriddedSamples(null,s1,s2,s3,t1)
  f2 = SimpleGridder3.getGriddedSamples(null,s1,s2,s3,t2)
  f3 = SimpleGridder3.getGriddedSamples(null,s1,s2,s3,t3)
  return mul(-1.0,f1[0]),mul(-1.0,f2[0]),mul(-1.0,f3[0]),f1[1],f1[2],f1[3]

def goCheckInterpolation():
  """Checks interpolation condition."""
  q1,q2,q3 = read("q1"),read("q2"),read("q3")
  #q1,q2,q3 = interpolateThrows()
  t1,t2,t3,x1,x2,x3 = getGriddedThrows()
  e1,e2,e3 = 0.0,0.0,0.0
  for i in range(ns):
    i1 = int(x1[i])
    i2 = int(x2[i])
    i3 = int(x3[i])
    e1 = e1+abs(t1[i]-q1[i3][i2][i1])
    e2 = e2+abs(t2[i]-q2[i3][i2][i1])
    e3 = e3+abs(t3[i]-q3[i3][i2][i1])
  print "error1 =",e1
  print "error2 =",e2
  print "error3 =",e3

def array(x1,x2,x3=None,x4=None):
  if x3 and x4:
    return jarray.array([x1,x2,x3,x4],Class.forName('[[[F'))
  elif x3:
    return jarray.array([x1,x2,x3],Class.forName('[[[F'))
  else:
    return jarray.array([x1,x2],Class.forName('[[[F'))

def zeros(n1,n2=None,n3=None):
  if n2 and n3:
    return zerofloat(n1,n2,n3)
  elif n2:
    return zerofloat(n1,n2)
  else:
    return zerofloat(n1)

def normal(x,y):
  div(x,max(abs(x)),y)

def like(x):
  return zerofloat(len(x[0][0]),len(x[0]),len(x))

def slog(x):
  return mul(sgn(x),sub(exp(abs(x)),1.0))
def sexp(x):
  return mul(sgn(x),log(add(abs(x),1.0)))

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
prism = ColorMap.PRISM
def plot(x,cmap=gray,cmin=0.0,cmax=0.0,perc=100.0,
         cbar=None,nearest=False,name=None):
  pan = panel(cbar)
  pix = pan.addPixels(x)
  pix.setColorModel(cmap)
  if cmin<cmax:
    pix.setClips(cmin,cmax)
  if perc<100:
    pix.setPercentiles(100-perc,perc)
  if nearest:
    pix.setInterpolation(PixelsView.Interpolation.NEAREST)
  frame(pan,name)

def panel(cbar=None):
  p = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  cb = p.addColorBar()
  if cbar:
    cb.setLabel(cbar)
  return p

def frame(panel,name=None):
  frame = PlotFrame(panel)
  frame.setBackground(Color(204,204,204,255))
  #frame.setFontSizeForSlide(1.0,1.0)
  #frame.setSize(1200,600)
  #frame.setSize(1200,1000)
  frame.setSize(1500,800)
  if name:
    frame.setTitle(name)
  frame.setVisible(True)
  if name and pngDir:
    frame.paintToPng(360,3.0,pngDir+name+'.png')
    #frame.paintToPng(1080,12.0,pngDir+name+'.png')

def read(name,image=None):
  if not image:
    image = zerofloat(n1,n2,n3)
  fileName = dataDir+name+".dat"
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def write(name,image):
  fileName = dataDir+name+'.dat'
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()

#############################################################################
# graphics

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
  ipg = ImagePanelGroup(s1,s2,s3,image)
  ipg.setColorModel(cmap)
  ipg.setSlices(k1,k2,k3)
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
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  frame = SimpleFrame(world)
  #frame.setBackground(Color(204,204,204,255))
  frame.setBackground(Color.WHITE)
  if name:
    frame.setTitle(name)
  view = frame.getOrbitView()
  zscale = 0.6*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.3)
  view.setAzimuth(azimuth)
  view.setElevation(elevation)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  frame.viewCanvas.setBackground(frame.getBackground())
  #frame.setSize(1020,750)
  frame.setSize(1800,1200)
  frame.setVisible(True)
  if pngDir and name:
    frame.paintToFile(pngDir+name+'.png')
  return frame

def slice12(k3,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n2)
  SimpleFloat3(f).get12(n1,n2,0,0,k3,s)
  return s

def slice13(k2,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n3)
  SimpleFloat3(f).get13(n1,n3,0,k2,0,s)
  return s

def slice23(k1,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n2,n3)
  SimpleFloat3(f).get23(n2,n3,k1,0,0,s)
  return s

#############################################################################
# Run the function main on the Swing thread
import time
class RunMain(Runnable):
  def run(self):
    start = time.time()
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02d'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
