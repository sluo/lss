'''
Non-vertical flattening using transformed rotation equations
'''
from imports import *

pngDir = None
#pngDir = '/Users/sluo/Desktop/'
seisDir = '/data/scratch/'
s1,s2,s3 = Sampling(201),Sampling(201),Sampling(201)
n1,n2,n3 = s1.count,s2.count,s3.count
k1,k2,k3 = 190,190,190; azimuth=240; elevation=20 # for 3D views
timer = Stopwatch()

#############################################################################

def main(args):
  show()
  #makeFake3d()
  #goFlatten()

def show():
  f = read('f'); display(f,name='f')
  g = read('g'); display(g,name='g')

def goFlatten():
  slopes = False
  shifts = False
  f = read('f')
  if slopes:
    u1,u2,u3,ep = like(f),like(f),like(f),like(f)
    timer.restart(); print 'slopes...'
    LocalOrientFilter(1.0,1.0).applyForNormalPlanar(f,u1,u2,u3,ep)
    timer.stop(); print 'slopes in %.2fs'%timer.time()
    fill(u1[n3/2][n2/2][n1/2],u1)
    fill(u2[n3/2][n2/2][n1/2],u2)
    fill(u3[n3/2][n2/2][n1/2],u3)
    fill(1.0,ep)
    #display(u1,cmap=jet,name='u1')
    #display(u2,cmap=jet,name='u2')
    #display(u3,cmap=jet,name='u3')
    #display(ep,cmap=jet,name='ep')
    write('u1',u1)
    write('u2',u2)
    write('u3',u3)
    write('ep',ep)
  if shifts:
    u1,u2,u3,ep = read('u1'),read('u2'),read('u3'),read('ep')
    timer.restart(); print 'shifts...'
    r = FlattenerRT(6.0,6.0).findShifts([u1,u2,u3,pow(ep,8.0)])
    timer.stop(); print 'shifts in %.2fs'%timer.time()
    r1,r2,r3 = r[0],r[1],r[2]
    #display(r1,cmap=jet,name='r1')
    #display(r2,cmap=jet,name='r2')
    #display(r3,cmap=jet,name='r3')
    g = FlattenerUtil.applyShiftsR(f,r)
    write('g',g)
    write('r1',r1)
    write('r2',r2)
    write('r3',r3)
  g,r1,r2,r3 = read('g'),read('r1'),read('r2'),read('r3')
  display(f,name='f')
  display(g,name='g')
  display(r1,cmap=jet,name='r1')
  display(r2,cmap=jet,name='r2')
  display(r3,cmap=jet,name='r3')

def like(x):
  return zerofloat(len(x[0][0]),len(x[0]),len(x))

def makeFake3d():
  n1p = int(1.8*n1)
  g = FakeData.seismic3d2010A(n1p,n2,n3,0.0,0.0,0.0,1.0,0.0)
  s = like(g)
  for i3 in range(n3):
    for i2 in range(n2):
      r = abs(i3+i2)*0.3
      s[i3][i2] = fillfloat(r,n1p)
  zero = zerofloat(n1p,n2,n3)
  t = FlattenerUtil.applyShiftsS(g,[s,zero,zero])
  f = copy(n1,n2,n3,int(0.5*(n1p-n1)),0,0,t)
  display(f)
  write('f',f)

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,cmap=gray,cmin=0.0,cmax=0.0,cbar=None,name=None):
  pan = panel(cbar)
  pix = pan.addPixels(x)
  pix.setColorModel(cmap)
  if cmin<cmax:
    pix.setClips(cmin,cmax)
  pix.setInterpolation(PixelsView.Interpolation.LINEAR)
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
  frame.setSize(1200,600)
  if name:
    frame.setTitle(name)
  frame.setVisible(True)
  if name and pngDir:
    frame.paintToPng(360,3.0,pngDir+name+'.png')

def read(name,image=None):
  if not image:
    image = zerofloat(n1,n2,n3)
  fileName = seisDir+name+'.dat'
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def write(name,image,directory=seisDir):
  fileName = directory+name+'.dat'
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()

from org.python.util import PythonObjectInputStream
def readTensors(name):
  fis = FileInputStream(seisDir+name+".dat")
  ois = PythonObjectInputStream(fis)
  tensors = ois.readObject()
  fis.close()
  return tensors

def writeTensors(name,tensors):
  fos = FileOutputStream(seisDir+name+".dat")
  oos = ObjectOutputStream(fos)
  oos.writeObject(tensors)
  fos.close()

#############################################################################
# graphics

def display(image,tensors=None,cmap=gray,cbar=None,
            cmin=0,cmax=0,perc=100,name=None):
  #return
  world = World()
  ipg = addImageToWorld(world,image,cmap,cmin,cmax,perc)
  if tensors:
    addTensorsToIpg(ipg,tensors)
  frame = makeFrame(world,name)
  if cbar:
    colorbar = addColorBar(frame,cbar)
    ipg.addColorMapListener(colorbar)

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

def addTensorsToIpg(ipg,mt):
  def add(ip,mt,esize=20):
    #tp = TensorsPanel(s1,s2,s3,mt)
    tp = TensorsPanel(mt)
    tp.setEllipsoidSize(esize)
    ip.getFrame().addChild(tp)
    return tp
  add(ipg.getImagePanel(Axis.X),mt)
  add(ipg.getImagePanel(Axis.Y),mt)
  add(ipg.getImagePanel(Axis.Z),mt)

def addColorBar(frame,label):
  cbar = ColorBar(label)
  cbar.setFont(cbar.getFont().deriveFont(24.0))
  frame.add(cbar,BorderLayout.EAST)
  #frame.viewCanvas.setBackground(frame.getBackground())
  return cbar

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
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.3)
  view.setAzimuth(azimuth)
  view.setElevation(elevation)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  frame.viewCanvas.setBackground(frame.getBackground())
  #frame.setSize(1460,980)
  frame.setSize(1020,750)
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

def normal(x,y):
  div(x,max(abs(x)),y)

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
