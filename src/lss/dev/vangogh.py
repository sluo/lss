"""
Testing linear/nonlinear isotropic/anisotropic smoothing filters
"""
from imports import *
from java.awt.image import *
from java.io import *
from javax.imageio import *
from ldf import BilateralFilter

pngDir = None
dataDir = '/Users/sluo/Desktop/vangogh/'
n1,n2 = 1765,2228

#############################################################################

def main(args):
  #makeImage()
  f = read('f'); plot(f,name='input')
  goLinearIsotropic(f)
  goLinearAnisotropic(f)
  goNonlinearIsotropic(f)
  goNonlinearAnisotropic(f)
  goBilateral(f)

def goLinearIsotropic(f):
  g = like(f)
  RecursiveGaussianFilter(12.0).apply00(f,g)
  plot(g,cmin=-1,cmax=1,name='linear isotropic')

def goLinearAnisotropic(f):
  g = like(f)
  RecursiveGaussianFilter(12.0).apply0X(f,g)
  RecursiveGaussianFilter(1.0).applyX0(g,g)
  plot(g,cmin=-1,cmax=1,name='linear anisotropic')

def goNonlinearIsotropic(f):
  g = like(f)
  s = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      if f[i2][i1]>0.5:
        s[i2][i1] = 1.0
  lsf = LocalSmoothingFilter()
  c = 1000.0
  lsf.apply(c,s,f,g)
  plot(g,cmin=-1,cmax=1,name='nonlinear isotropic')

def goNonlinearAnisotropic(f):
  g = like(f)
  lof = LocalOrientFilter(24.0)
  t = lof.applyForTensors(f)
  t.invertStructure(1.0,10.0)
  lsf = LocalSmoothingFilter()
  c = 100.0
  lsf.apply(t,c,f,g)
  plot(g,cmin=-1,cmax=1,name='nonlinear anisotropic')

def goBilateral(f):
  g = like(f)
  sigmaS = 24.0 # spatial filter sigma
  sigmaR = computeSigmaR(f)
  blf = BilateralFilter(sigmaS,sigmaR)
  blf.apply(f,g)
  plot(g,cmin=-1,cmax=1,name='bilateral filter')

def computeSigmaR(g):
  sigmaR = 0.5*(Quantiler.estimate(0.75,g)-Quantiler.estimate(0.25,g))
  print 'sigmaR=%f'%sigmaR
  return sigmaR

def like(f):
  return zerofloat(len(f[0]),len(f))

def makeImage():
  ffile = File('/Users/sluo/Desktop/vangogh/small.jpg')
  img = ImageIO.read(ffile)
  n1,n2 = img.getHeight(),img.getWidth()
  f = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      f[i2][i1] = img.getRGB(i2,i1)
  f = sub(f,min(f))
  f = div(f,max(abs(f)))
  f = mul(2.0,sub(f,0.5))
  plot(f)
  write('f',f)
  return f

def read(name,image=None,dir=None):
  if not image:
    image = zerofloat(n1,n2)
  if not dir:
    dir = dataDir
  fileName = dir+name+".dat"
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

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
bwr = ColorMap.BLUE_WHITE_RED
gyr = ColorMap.GRAY_YELLOW_RED
def plot(x,cmap=gyr,cmin=0,cmax=0,perc=100,cbar=None,name=None):
  pan = panel(cbar)
  pix = pan.addPixels(x)
  pix.setColorModel(cmap)
  if cmin<cmax:
    pix.setClips(cmin,cmax)
  if perc<100:
    pix.setPercentiles(100-perc,perc)
  pix.setInterpolation(PixelsView.Interpolation.LINEAR)
  frame(pan,name)

def panel(cbar=None):
  p = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  #p.setHLabel("index i2")
  #p.setVLabel("index i1")
  cb = p.addColorBar()
  if cbar:
    cb.setLabel(cbar)
    cb.setWidthMinimum(120)
  else:
    cb.setWidthMinimum(80)
  return p

def frame(panel,name=None):
  frame = PlotFrame(panel)
  #frame.setBackground(Color(204,204,204,255))
  #frame.setFontSizeForSlide(1.0,1.0)
  frame.setSize(800,600)
  if name:
    frame.setTitle(name)
  frame.setVisible(True)
  if name and pngDir:
    frame.paintToPng(360,3.0,pngDir+name+'.png')

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
    print '%02d:%02d:%02ds'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
