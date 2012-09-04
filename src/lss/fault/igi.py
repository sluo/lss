"""
Image-guided interpolation
"""
from imports import *
from java.awt.image import IndexColorModel

subset = "a"
#subset = "b"

pngDir = None
#pngDir = "/Users/sluo/Desktop/"

#  s1 = Sampling(90,0.004,1.484)
#  s2 = Sampling(221,0.025,1.25)
#  s3 = Sampling(220,0.025,2.5)

if subset=='a':
  dataDir = '/data/seis/f3/suba/'
  s1 = Sampling(90)
  s2 = Sampling(221)
  s3 = Sampling(220)
  n1,n2,n3 = s1.count,s2.count,s3.count
  k1,k2,k3 = 22,48,151
  h1s = [32]
elif subset=='b':
  dataDir = '/data/seis/f3/subb/'
  s1 = Sampling(90)
  s2 = Sampling(221)
  s3 = Sampling(220)
  n1,n2,n3 = s1.count,s2.count,s3.count
  k1,k2,k3 = 56,48,150
  h1s = [55]
h1s = [32] if subset=="a" else [75] # horizons

#############################################################################

def main(args):
  igi()

def igi():
  semblance = True # replace eigenvalues with local semblance
  k3 = 126
  h = slice12(k3,read("h"))
  #h = slice12(k3,read("g"))
  g = slice12(k3,read("g"))
  f,x1,x2 = fakeSamples()
  et = LocalOrientFilter(8.0).applyForTensors(h)
  if semblance:
    eps = 0.001
    lsf1 = LocalSemblanceFilter(2,2)
    lsf2 = LocalSemblanceFilter(2,8)
    m1 = lsf1.semblance(LocalSemblanceFilter.Direction2.V,et,g)
    m2 = lsf2.semblance(LocalSemblanceFilter.Direction2.UV,et,g)
    pow(m2,4.0,m2)
    fill(eps,m2)
    m1 = clip(eps,1.0,m1)
    m2 = clip(eps,1.0,m2)
    et.setEigenvalues(m2,m1)
  else:
    et.invertStructure(1.0,2.0)
  bg2 = BlendedGridder2(et,f,x1,x2)
  y = bg2.grid(s1,s2)
  def getShiftsS():
    #t1,t2,t3,x1,x2,x3 = getGriddedThrows()
    #b1 = gridBlended(t1,x1,x2,x3)
    #b2 = gridBlended(t2,x1,x2,x3)
    #b3 = gridBlended(t3,x1,x2,x3)
    b1 = read('q1')
    b2 = read('q2')
    b3 = read('q3')
    return b1,b2,b3
  b1,b2,b3 = getShiftsS()
  s = [slice12(k3,b1),slice12(k3,b2),slice12(k3,b3)]
  plot2(h,y,name='h')
  plot2(g,FlattenerUtil.applyShiftsSLinear(y,s),name='g')

def fakeSamples():
  k2 = n2/3
  f = rampfloat(2.0,1.0/n1,n1)
  #f = zerofloat(n1)
  for i in range(n1):
    f[i] += (i*8/n1)*0.5
      
  #random = Random()
  #for i in range(n1):
  #  f[i] += 0.05*random.nextGaussian()
  x1 = rampfloat(0.0,1.0,n1)
  x2 = fillfloat(k2,n1)
  return f,x1,x2

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,cmap=gray,cmin=0.0,cmax=0.0,cbar=None,name=None):
  pan = panel()
  s1 = Sampling(90,0.004,1.484)
  s2 = Sampling(221,0.025,1.25)
  s3 = Sampling(220,0.025,2.5)
  pix = pan.addPixels(s1,s2,x)
  pix.setColorModel(cmap)
  if cmin<cmax:
    pix.setClips(cmin,cmax)
  pix.setInterpolation(PixelsView.Interpolation.LINEAR)
  frame(pan,name)

def plot2(f1,f2,cmap1=gray,cmap2=jet,name=None):
  pan = panel()
  s1 = Sampling(90,0.004,1.484)
  s2 = Sampling(221,0.025,1.25)
  s3 = Sampling(220,0.025,2.5)
  pix1 = pan.addPixels(s1,s2,f1)
  pix1.setColorModel(cmap1)
  def trans(alpha,cmap):
    r,g,b = zerobyte(256),zerobyte(256),zerobyte(256)
    a = fillbyte(int(alpha*255),256)
    cmap.getReds(r)
    cmap.getGreens(g)
    cmap.getBlues(b)
    return IndexColorModel(8,256,r,g,b,a)
  pix2 = pan.addPixels(s1,s2,f2)
  pix2.setColorModel(trans(0.0,cmap2))
  frame(pan,name)

def panel():
  p = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p.setVLabel("Time (s)")
  p.setVInterval(0.1);
  p.setHLabel("Inline (km)")
  cb = p.addColorBar("Velocity (km/s)")
  cb.setWidthMinimum(120)
  #if cbar:
  #  cb.setLabel(cbar)
  return p

def frame(panel,name=None):
  frame = PlotFrame(panel)
  #frame.setBackground(Color(204,204,204,255))
  frame.setFontSizeForSlide(1.0,1.0)
  #frame.setSize(1200,600)
  frame.setSize(1024,700)
  if name:
    frame.setTitle(name)
  frame.setVisible(True)
  if name and pngDir:
    frame.paintToPng(360,3.0,pngDir+name+'.png')

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

def display(image,cmap=gray,name=None):
  world = World()
  addImageToWorld(world,image,cmap)
  makeFrame(world,name)

def display2(image1,image2,cmap1=gray,cmap2=gray,name=None):
  world = World()
  addImageToWorld(world,image1,cmap1)
  addImageToWorld(world,image2,cmap2)
  makeFrame(world,name)

def addImageToWorld(world,image,cmap=gray,cmin=0,cmax=0):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  ipg.setColorModel(cmap)
  ipg.setSlices(k1,k2,k3)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
  world.addChild(ipg)
  return ipg

def makeFrame(world,name=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  frame = SimpleFrame(world)
  frame.setBackground(Color(204,204,204,255))
  if name:
    frame.setTitle(name)
  view = frame.getOrbitView()
  zscale = 0.6*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.4)
  view.setAzimuth(azimuth)
  view.setElevation(elevation)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(1100,800)
  frame.setVisible(True)
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
import sys,time
class RunMain(Runnable):
  def run(self):
    start = time.time()
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02d'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
