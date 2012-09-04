'''
Interpolation of missing data by POCS
'''
from imports import *

pngDir=None
#pngDir='/Users/sluo/Desktop/'
#seisDir = './bp2011/'
seisDir = '/Users/sluo/Desktop/others/tony/bp2011/'
s1,s2 = Sampling(707),Sampling(1027)
n1,n2 = s1.count,s2.count

#############################################################################

def main(args):
  interpf()

def interpf():
  h,f,m = read('h'),read('f'),read('m')
  fi = PocsInterpolator(n1,n2)
  g = fi.interpolate(0.0,m,f)
  c = 0.1*max(abs(f))
  plot(h,cmin=-c,cmax=c,name='h')
  plot(g,cmin=-c,cmax=c,name='g')
  plot(f,cmin=-c,cmax=c,name='f')
  plot(m,cmap=jet,name='m')

def like(x):
  return zerofloat(len(x[0]),len(x))

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,et=None,cmap=gray,cmin=0,cmax=0,perc=100,cbar=None,name=None):
  pan = panel(cbar)
  pix = pan.addPixels(x)
  pix.setColorModel(cmap)
  pix.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pix.setClips(cmin,cmax)
  if perc<100:
    pix.setPercentiles(100-perc,perc)
  if et:
    tv = TensorsView(et)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(2)
    tv.setEllipsesDisplayed(20)
    tv.setScale(0.9)
    pan.getTile(0,0).addTiledView(tv)
  frame(pan,name)

def panel(cbar=None):
  p = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  cb = p.addColorBar()
  cb.setWidthMinimum(80)
  if cbar:
    cb.setLabel(cbar)
  return p

def frame(panel,name=None):
  frame = PlotFrame(panel)
  frame.setBackground(Color(204,204,204,255))
  frame.setFontSizeForSlide(1.5,1.5)
  frame.setSize(1024,700)
  if name:
    frame.setTitle(name)
  frame.setVisible(True)
  if name and pngDir:
    frame.paintToPng(360,3.0,pngDir+name+'.png')

def read(name,image=None):
  if not image:
    image = zerofloat(n1,n2)
  fileName = seisDir+name+'.dat'
  ais = ArrayInputStream(fileName)
  #ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image

def write(name,image,directory=seisDir):
  fileName = directory+name+'.dat'
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()

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
