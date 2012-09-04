'''
Metric tensors computed from shifts r(u)
'''
from imports import *

pngDir = None
#pngDir = '/Users/sluo/Desktop/'
s1,s2 = Sampling(301),Sampling(501)
n1,n2 = s1.count,s2.count
timer = Stopwatch()

#############################################################################

def main(args):
  #flatten()
  goTensors()

def goTensors():
  epow = False # raise eigenvalues to some power
  f = fakeImage()
  g,r,u1,u2,_ = flatten()
  u = [u1,u2]
  st = FlattenerUtil.getMetricTensors(r)
  et = st.asEigenTensors()
  if epow:
    #et.invertStructure(1.0,1.0)
    a = zerofloat(2)
    for i2 in range(n2):
      for i1 in range(n1):
        et.getEigenvalues(i1,i2,a)
        pow(a,1.0,a)
        et.setEigenvalues(i1,i2,a)
  a11,a12,a22 = like(f),like(f),like(f)
  st.getTensors(a11,a12,a22)
  plot(a11,cmap=rwb,cmin=0.5,cmax=1.5,name="a11")
  plot(a12,cmap=rwb,cmin=-0.5,cmax=0.5,name="a12")
  plot(a22,cmap=rwb,cmin=0.5,cmax=1.5,name="a22")
  #z = mul(a12,mul(sub(a11,1.0),sub(a22,1.0))); cz = max(abs(z))
  #z = mul(sub(a11,1.0),sub(a22,1.0)); cz = max(abs(z))
  z = sub(a22,a11); cz = max(abs(z))
  plot(z,cmap=rwb,cmin=-cz,cmax=cz,name="z")
  la22 = laplacian(a22); cl = 0.02*max(abs(la22))
  plot(la22,cmap=jet,perc=99.9,name="laplacian(a22)")
  plot(g,name="g") # metric tensors
  plot(g,et,name="tensors") # metric tensors
  d = sub(1.0,determinant(st))
  plot(d,cmap=rwb,cmin=-0.5,cmax=0.5,name="determinant")
  x1,x2 = zerofloat(n1,n2,2),zerofloat(n1,n2,2)
  FlattenerUtil.getFrame(r,x1,x2)
  FlattenerUtil.normalize(x1,x1)
  angle = dot(u,x1); ca = max(abs(angle))
  plot(angle,cmap=rwb,name="angle")
  plot(laplacian(angle),cmap=rwb,name="laplacian(angle)")

def dot(x,y):
  z = like(y[0])
  for i2 in range(n2):
    for i1 in range(n1):
      z[i2][i1] = x[0][i2][i1]*y[0][i2][i1]+x[1][i2][i1]*y[1][i2][i1]
  return z

def determinant(st):
  d = zerofloat(n1,n2)
  a = zerofloat(3)
  for i2 in range(n2):
    for i1 in range(n1):
      st.getTensor(i1,i2,a)
      a11,a12,a22 = a[0],a[1],a[2]
      d[i2][i1] = a11*a22-a12*a12
  return d

def laplacian(a):
  y = like(a)
  for i2 in range(1,n2-1):
    for i1 in range(1,n1-1):
      y[i2][i1] = a[i2-1][i1]+a[i2][i1-1]+\
                  a[i2+1][i1]+a[i2][i1+1]-4.0*a[i2][i1]
  RecursiveGaussianFilter(6.0).apply00(y,y)
  return y

def xlaplacian(a):
  a11,a22 = like(a),like(a)
  rgf = RecursiveGaussianFilter(8.0)
  rgf.apply20(a,a11)
  rgf.apply02(a,a22)
  return add(a11,a22)

def displayTensors():
  et = readTensors('et')
  g,m,el = read('g'),read('m'),read('el')
  mask = ZeroMask(m)
  mask.apply(0.0001,el)
  pow(el,8.0,el)
  et.scale(el)
  display(g,tensors=et,name='et')

def flatten():
  f = fakeImage()
  u1,u2,el = like(f),like(f),like(f)
  LocalOrientFilter(1.0,1.0).applyForNormalLinear(f,u1,u2,el)
  #fill(u1[n2/2][n1/2],u1),fillfloat(u2[n2/2][n1/2],u2),fill(1.0,el)
  a = fillfloat(1.0,n1,n2) # rotation
  p = [u1,u2,pow(el,8.0),a]
  #r = FlattenerR(6.0,6.0).findShifts(p)
  r = FlattenerRT(6.0,6.0).findShifts(p)
  g = FlattenerUtil.applyShiftsR(f,r)
  plot(f,name='f')
  plot(g,name='g')
  return g,r,u1,u2,el

def fakeImage():
  #return constant()
  return fakedata()

def constant():
  maxshift = 200.0
  ds2 = maxshift/n2
  m1,m2 = 2*n1,2*n2
  f = FakeData.seismic2d2011A(m1,m2,0.0)
  s = rampfloat(0.0,0.0,ds2,m1,m2)
  sub(s,sum(s)/m1/m2,s)
  t = FlattenerCg().applyShifts(f,s)
  g = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      g[i2][i1] = t[i2+n2/2][i1+n1/2]
  return g

def like(x):
  return zerofloat(len(x[0]),len(x))

def fakedata():
  maxdip = 40.0
  f = FakeData.seismic2d2011A(2*n1,n2,maxdip)
  g = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      g[i2][i1] = f[i2][i1+n1/2]
  return g

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,et=None,cmap=gray,cmin=0,cmax=0,perc=100,cbar=None,name=None):
  pan = panel(cbar)
  pix = pan.addPixels(x)
  pix.setColorModel(cmap)
  pix.setInterpolation(PixelsView.Interpolation.LINEAR)
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
  if cbar:
    cb.setLabel(cbar)
  return p

def frame(panel,name=None):
  frame = PlotFrame(panel)
  frame.setBackground(Color(204,204,204,255))
  #frame.setFontSizeForSlide(1.0,1.0)
  frame.setSize(1200,690)
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

def write(name,image):
  fileName = seisDir+name+'.dat'
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
