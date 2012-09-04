from imports import *

pngDir = None
seismicDir = '/data/seis/f3/'
s1 = Sampling(462,0.004,0.004)
s2 = Sampling(651,0.025,0.000)
s3 = Sampling(951,0.025,0.000)
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
f1,f2,f3 = s1.first,s2.first,s3.first
l1,l2,l3 = s1.last,s2.last,s3.last

############################################################################

def main(args):
  f = read('f3')
  display(f,perc=99.5)

def read(name,image=None):
  if image==None:
    image = zerofloat(n1,n2,n3)
  fileName = seismicDir+name+".dat"
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def write(name,image):
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def display(f,cmap=ColorMap.GRAY,cmin=0,cmax=0,perc=100):
  world = World()
  ipg = addImageToWorld(world,f,cmap,cmin,cmax,perc)
  makeFrame(world)

def addImageToWorld(world,image,cmap=gray,cmin=0,cmax=0,perc=100):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  ipg.setColorModel(cmap)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
  if perc<100:
    ipg.setPercentiles(100-perc,perc)
  world.addChild(ipg)
  return ipg

def addImage2ToWorld(world,image1,image2):
  ipg = ImagePanelGroup2(s1,s2,s3,image1,image2)
  ipg.setColorModel1(ColorMap.getGray())
  #ipg.setColorModel2(ColorMap.getJet(0.3))
  ipg.setColorModel2(ColorMap.getHue(0.0,20.0,0.3))
  world.addChild(ipg)
  return ipg

def makeFrame(world):
  frame = SimpleFrame(world)
  #frame.setBackground(Color(255,254,255,255))
  view = frame.getOrbitView()
  zscale = 0.6*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.2)
  #view.setAzimuth(270)
  #view.setElevation(50)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
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
    print '%02d:%02d:%02d'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
