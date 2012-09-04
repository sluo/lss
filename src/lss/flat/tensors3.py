'''
Metric and curvature tensors computed from shifts r(u)
'''
from imports import *

pngDir = None
#pngDir = '/Users/sluo/Desktop/'
seisDir = '/data/seis/tp/csm/seismict/subt_251_4_500/flatr/'
#seisDir = '/data/seis/fake/t1/' # TP = 1.0
#seisDir = '/data/seis/fake/t0/' # TP = 0.0
s1,s2,s3 = Sampling(251),Sampling(357),Sampling(161)
n1,n2,n3 = s1.count,s2.count,s3.count
k1,k2,k3 = 197,310,110; azimuth=240; elevation=20 # for 3D views
timer = Stopwatch()

#############################################################################

def main(args):
  show()
  #makeFake3d()
  #flatten()
  #metricTensors()
  #curvatureTensors()
  #gaussianCurvature()
  #compareNormals()

def show():
  f = read('f')
  g = read('g')
  k = read('k')
  a13 = read('a13')
  a33 = read('a33')
  t = readTensors('mt');
  display(f,cbar='Amplitude',name='f')
  display(g,cbar='Amplitude',name='g')
  display(k,cmap=rwb,cmin=-0.1,cmax=0.1, cbar='K',name='k')
  display(a13,cmap=rwb,cmin=-0.5,cmax=0.5, cbar='(x_u).(x_t)',name='a13')
  display(a33,cmap=rwb,cmin=0.8,cmax=1.2, cbar='||x_u||^2',name='a33')
  displayTensors(t)

def compareNormals():
  """Compares normal vectors computed with LocalOrientFilter 
  to those computed from metric tensors.
  """
  computeNormals = False
  if computeNormals:
    r1,r2,r3 = read('r1'),read('r2'),read('r3')
    r = [r1,r2,r3]
    x2 = [like(r1),like(r1),like(r1)]
    x3 = [like(r1),like(r1),like(r1)]
    v =  [like(r1),like(r1),like(r1)]
    FlattenerUtil.getFrame(r,None,x2,x3)
    FlattenerUtil.cross(x3,x2,v)
    FlattenerUtil.normalize(v,v)
    write('v1',v[0])
    write('v2',v[1])
    write('v3',v[2])
  v1,v2,v3 = read('v1'),read('v2'),read('v3')
  u1,u2,u3 = read('u1'),read('u2'),read('u3')
  display(sub(v1,u1),cmap=rwb,cmin=-0.2,cmax=0.2,name='v1-u1')
  display(sub(v2,u2),cmap=rwb,cmin=-0.2,cmax=0.2,name='v2-u2')
  display(sub(v3,u3),cmap=rwb,cmin=-0.2,cmax=0.2,name='v3-u3')

def gaussianCurvature():
  computeCurvature = False
  f = read('f')
  k = like(f)
  if computeCurvature:
    mt,ct = readTensors('mt'),readTensors('ct')
    mti,cti = zerofloat(6),zerofloat(6)
    for i3 in range(n3):
      for i2 in range(n2):
        for i1 in range(n1):
          mt.getTensor(i1,i2,i3,mti)
          ct.getTensor(i1,i2,i3,cti)
          m22,m23,m33 = mti[3],mti[4],mti[5]
          c22,c23,c33 = cti[3],cti[4],cti[5]
          num = c22*c33-c23
          den = m22*m33-m23
          k[i3][i2][i1] = num/den
    #mask(k)
    write('k',k)
  read('k',k)
  print 'min =',min(k),' max =',max(k)
  display(k,cmap=rwb,cmin=-0.1,cmax=0.1,cbar="Gaussian curvature",name="k")

def metricTensors():
  makeTensors = False
  if makeTensors:
    r = (read('r1'),read('r2'),read('r3'))
    timer.restart(); print 'tensors...'
    t = FlattenerUtil.getMetricTensors(r)
    a11 = zerofloat(n1,n2,n3)
    a12 = zerofloat(n1,n2,n3)
    a13 = zerofloat(n1,n2,n3)
    a22 = zerofloat(n1,n2,n3)
    a23 = zerofloat(n1,n2,n3)
    a33 = zerofloat(n1,n2,n3)
    t.getTensors(a11,a12,a13,a22,a23,a33)
    mt = t.asEigenTensors()
    timer.stop(); print 'tensors in %.2fs'%timer.time()
    #writeTensors('mt',mt)
    write('a11',a11)
    write('a12',a12)
    write('a13',a13)
    write('a22',a22)
    write('a23',a23)
    write('a33',a33)
  mt = readTensors('mt')
  displayTensors(mt)

def curvatureTensors():
  makeTensors = False
  if makeTensors:
    r = [read('r1'),read('r2'),read('r3')]
    u1 = FlattenerUtil.applyShiftsRLinear(read('u1'),r)
    u2 = FlattenerUtil.applyShiftsRLinear(read('u2'),r)
    u3 = FlattenerUtil.applyShiftsRLinear(read('u3'),r)
    u = [u1,u2,u3]
    timer.restart(); print 'tensors...'
    ct = FlattenerUtil.getCurvatureTensors(r,u).asEigenTensors()
    timer.stop(); print 'tensors in %.2fs'%timer.time()
    writeTensors('ct',ct)
  ct = readTensors('ct')
  displayTensors(ct)

def displayTensors(t):
  g,ep = read('g'),read('ep')
  mask(ep)
  r1,r2,r3 = read('r1'),read('r2'),read('r3')
  ep = FlattenerUtil.applyShiftsRLinear(ep,[r1,r2,r3])
  pow(ep,8.0,ep)
  t.scale(ep)
  display(g,tensors=t,name='tensors')

def flatten():
  slopes = False
  shifts = False
  f = read('f')
  if slopes:
    u1,u2,u3,ep = like(f),like(f),like(f),like(f)
    lof = LocalOrientFilter(8.0,2.0)
    timer.restart(); print 'slopes...'
    lof.applyForNormalPlanar(f,u1,u2,u3,ep)
    timer.stop(); print 'slopes in %.2fs'%timer.time()
    #display(u1,cmap=jet,name='u1')
    #display(u2,cmap=jet,name='u2')
    #display(u3,cmap=jet,name='u3')
    #display(ep,cmap=jet,name='ep')
    write('u1',u1)
    write('u2',u2)
    write('u3',u3)
    write('ep',ep)
  flat = FlattenerR(8.0,8.0)
  if shifts:
    u1,u2,u3,ep = read('u1'),read('u2'),read('u3'),read('ep')
    a = fillfloat(1.0,n1,n2,n3) # flattening by rotation
    p = [u1,u2,u3,pow(ep,8.0),a]
    timer.restart(); print 'shifts...'
    r = flat.findShifts(p)
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
  g,r = read('g'),[read('r1'),read('r2'),read('r3')]
  display(f,name='f')
  display(g,name='g')

def mask(x,value=0.0001):
  try:
    m = read('m')
  except:
    pass
  else:
    ZeroMask(m).apply(value,x)

def like(x):
  return zerofloat(len(x[0][0]),len(x[0]),len(x))

def makeFake3d():
  n1p = int(1.8*n1)
  g = FakeData.seismic3d2010A(n1p,n2,n3,0.0,0.0,0.0,1.0,0.0)
  s = like(g)
  for i3 in range(n3):
    for i2 in range(n2):
      r = abs(i3+i2)*0.2
      s[i3][i2] = fillfloat(r,n1p)
  zero = zerofloat(n1p,n2,n3)
  t = FlattenerUtil.applyShiftsS(g,[s,zero,zero])
  f = copy(n1,n2,n3,int(0.5*(n1p-n1)),0,0,t)
  display(f)
  write('f',f,directory='/data/seis/fake/')

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
