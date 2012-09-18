"""
Nonlinear structure tensors (See Brox, 2006)
"""
from imports import *

pngDir = None
#pngDir = '/Users/sluo/Desktop/'
#seisImage = 'fakeA'
#seisImage = 'fakeB'
seisImage = 'fakeC'
#seisImage = 'f3d'

#############################################################################

def main(args):
  #f = getImage(); plot(f)
  #goLof()
  #goIsotropic()
  #goAnisotropic()
  #goFlattenLof()
  #goFlattenIsotropic()
  goFlattenAnisotropic()

def goFlattenLof():
  f = getImage()
  sigma1,sigma2 = 8.0,1.0
  u1,u2 = goLof(f,sigma1,sigma2)
  flattenS(f,u1,u2)

def goFlattenIsotropic():
  f = getImage()
  u1,u2 = goIsotropic(f)
  flattenS(f,u1,u2)

def goFlattenAnisotropic():
  f = getImage()
  u1,u2 = goAnisotropic(f)
  flattenS(f,u1,u2)

def flattenS(f,u1,u2):
  p2 = getSlopesFromNormals(u1,u2)
  fl = FlattenerS(6.0,6.0)
  s = fl.findShifts(p2)
  g = fl.applyShifts(f,s)
  plot(f,name='f')
  plot(g,name='g')

def goIsotropic(f=None):
  if f==None:
    f = getImage(); plot(f,name='f')
  g11,g12,g22 = getGradientOuterProducts(f)
  d = getDiffusionScalars(f)
  plot(d,cmap=jet,name='diffusivity')
  t11,t12,t22 = like(g11),like(g12),like(g22)
  lsf = LocalSmoothingFilter()
  c = 1.0
  lsf.apply(c,d,g11,t11)
  lsf.apply(c,d,g12,t12)
  lsf.apply(c,d,g22,t22)
  u1,u2,_,_,_,_ = getEigenFromTensors(t11,t12,t22)
  plot(u1,cmap=jet,name='u1')
  plot(u2,cmap=jet,name='u2')
  return u1,u2
    
def goAnisotropic(f=None):
  if f==None:
    f = getImage(); plot(f,name='f')
  f = getImage(); plot(f,name='f')
  g11,g12,g22 = getGradientOuterProducts(f)
  t11,t12,t22 = like(g11),like(g12),like(g22)
  t = getDiffusionTensors(f)
  lsf = LocalSmoothingFilter()
  c = 1.0
  lsf.apply(t,c,g11,t11)
  lsf.apply(t,c,g12,t12)
  lsf.apply(t,c,g22,t22)
  u1,u2,_,_,_,_ = getEigenFromTensors(t11,t12,t22)
  plot(u1,cmap=jet,name='u1')
  plot(u2,cmap=jet,name='u2')
  return u1,u2

def getDiffusionScalars(f):
  """Scalar diffusivity."""
  g11,g12,g22 = getGradientOuterProducts(f)
  h11 = getGradientInnerProducts(g11)
  h12 = getGradientInnerProducts(g12)
  h22 = getGradientInnerProducts(g22)
  h = getGradientInnerProducts(f) # extra channel
  plot(h11,cmap=jet,name='h11')
  plot(h12,cmap=jet,name='h12')
  plot(h22,cmap=jet,name='h22')
  return diffusivity(add(add(h11,h12),h22))
  #return diffusivity(add(add(h11,h12),add(h22,h)))

def getDiffusionTensors(f):
  """Anisotropic diffusion tensor."""
  g11,g12,g22 = getGradientOuterProducts(f)
  g1111,g1112,g1122 = getGradientOuterProducts(g11)
  g1211,g1212,g1222 = getGradientOuterProducts(g12)
  g2211,g2212,g2222 = getGradientOuterProducts(g22)
  t11 = add(add(g1111,g1211),g2211)
  t12 = add(add(g1112,g1212),g2212)
  t22 = add(add(g1122,g1222),g2222)
  u1,u2,v1,v2,eu,ev = getEigenFromTensors(t11,t12,t22)
  du = diffusivity(eu)
  dv = diffusivity(ev)
  plot(eu,cmap=jet,name='eu')
  plot(ev,cmap=jet,name='ev')
  plot(du,cmap=jet,name='du')
  plot(dv,cmap=jet,name='dv')
  #t11,t12,t22 = getTensorsFromEigen(u1,u2,v1,v2,du,dv)
  #return t11,t12,t22
  return EigenTensors2(u1,u2,du,dv)

def diffusivity(g,eps=0.1):
  """Diffusivity function."""
  d = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      gi = g[i2][i1]
      d[i2][i1] = 1.0/sqrt(eps*eps+gi)
  return d

def getSlopesFromNormals(u1,u2):
  p2min,p2max = -2.0,2.0
  n1,n2 = len(u1[0]),len(u1)
  p2 = like(u1)
  for i2 in range(n2):
    for i1 in range(n1):
      u1i = u1[i2][i1]
      u2i = u2[i2][i1]
      if -u2i<p2min*u1i:
        u2i = -p2min*u1i
      if -u2i>p2max*u1i:
        u2i = -p2max*u1i
      if u1i==0.0:
        p2[i2][i1] = p2max if u2i<0.0 else p2min
      else:
        p2[i2][i1] = -u2i/u1i
  return p2

def getEigenFromTensors(t11,t12,t22):
  """Eigendecomposition of symmetric 2x2 matrices."""
  n1,n2 = len(t11[0]),len(t11)
  u1,u2 = zerofloat(n1,n2),zerofloat(n1,n2)
  v1,v2 = zerofloat(n1,n2),zerofloat(n1,n2)
  eu,ev = zerofloat(n1,n2),zerofloat(n1,n2)
  a = zerodouble(2,2)
  v = zerodouble(2,2) # eigenvectors
  d = zerodouble(2)   # eigenvalues
  for i2 in range(n2):
    for i1 in range(n1):
      a[0][0] = t11[i2][i1]
      a[0][1] = a[1][0] = t12[i2][i1]
      a[1][1] = t22[i2][i1]
      Eigen.solveSymmetric22(a,v,d)
      u1[i2][i1] = v[0][0]
      u2[i2][i1] = v[0][1]
      v1[i2][i1] = v[1][0]
      v2[i2][i1] = v[1][1]
      eu[i2][i1] = d[0]
      ev[i2][i1] = d[1]
  return u1,u2,v1,v2,eu,ev

def getTensorsFromEigen(u1,u2,v1,v2,eu,ev):
  """ A = QDQ' """
  a11 = zerofloat(n1,n2)
  a12 = zerofloat(n1,n2)
  a22 = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      u1i = u1[i2][i1]
      u2i = u2[i2][i1]
      v1i = v1[i2][i1]
      v2i = v2[i2][i1]
      eui = eu[i2][i1]
      evi = ev[i2][i1]
      a11[i2][i1] = eui*u1i*u1i+evi*v1i*v1i
      a12[i2][i1] = eui*u1i*u2i+evi*v1i*v2i
      a22[i2][i1] = eui*u2i*u2i+evi*v2i*v2i
  return a11,a12,a22

def goLof(f=None,sigma1=12.0,sigma2=12.0):
  if f==None:
    f = getImage()
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof = LocalOrientFilter(sigma1,sigma2)
  lof.applyForNormalLinear(f,u1,u2,el)
  el = pow(el,8)
  #plot(f)
  plot(u1,cmap=jet,name='u1 lof')
  plot(u2,cmap=jet,name='u2 lof')
  #plot(el,cmap=jet,name='el')
  #plot(div(u2,u1),cmap=jet,perc=100)
  return u1,u2

def getGradientOuterProducts(f):
  """Unsmoothed outer product of image gradients."""
  g1,g2 = getGradients(f)
  return mul(g1,g1),mul(g1,g2),mul(g2,g2)

def getGradientInnerProducts(f):
  """Inner product of image gradients."""
  g1,g2 = getGradients(f)
  return add(mul(g1,g1),mul(g2,g2))

def getGradients(f):
  g1 = zerofloat(n1,n2)
  g2 = zerofloat(n1,n2)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply10(f,g1)
  rgf.apply01(f,g2)
  return g1,g2

def addNoise(nrms,f):
  n1,n2 = len(f[0]),len(f)
  r = Random(31415)
  nrms *= max(abs(f))
  g = mul(2.0,sub(randfloat(r,n1,n2),0.5))
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply10(g,g)
  frms = sqrt(sum(mul(f,f))/n1/n2)
  grms = sqrt(sum(mul(g,g))/n1/n2)
  g = mul(g,nrms*frms/grms)
  return add(f,g)

def getImage():
  if seisImage=='fakeA':
    return fakeImageA()
  elif seisImage=='fakeB':
    return fakeImageB()
  elif seisImage=='fakeC':
    return fakeImageC()
  elif seisImage=='f3d':
    return read('f3d400')

def fakeImageA():
  f = FakeData.seismic2d2012A(n1,n2,0.0,0.0,0.0,1.0,0.0)
  r = like(f)
  t = copy(f[0])
  si = SincInterpolator()
  si.setUniform(n1,1.0,0.0,t)
  slope = 0.5
  for i2 in range(0,2*n2/3):
    for i1 in range(n1/3,n1):
      f[i2][i1] = si.interpolate(slope*i2+i1-n1/3)
      r[i2][i1] = 1.0
  #return f
  #return f,r
  #return addNoise(0.1,f),r
  return addNoise(0.1,f)

def fakeImageB():
  str = 0.0 # maximum vertical shift when adding structure.
  fau = 0.0 # maximum displacement when adding a fault.
  ero = 0.0 # sample at which to create an erosional unconformity.
  amp = 0.8 # minimum scale factor when scaling amplitudes.
  noi = 0.0 # rms of added noise, relative to rms signal.
  f = FakeData.seismic2d2012A(n1,n2,str,fau,ero,amp,noi)
  r = like(f)
  t = copy(f[0])
  si = SincInterpolator()
  si.setUniform(n1,1.0,0.0,t)
  slope = 0.1
  for i2 in range(0,n2):
    for i1 in range(n1/3,n1):
      f[i2][i1] = si.interpolate(slope*i2+i1-n1/3)
  return addNoise(0.1,f)

def fakeImageC():
  str = 24.0 # maximum vertical shift when adding structure.
  fau = 8.0 # maximum displacement when adding a fault.
  ero = 50 # sample at which to create an erosional unconformity.
  amp = 1.0 # minimum scale factor when scaling amplitudes.
  noi = 0.0 # rms of added noise, relative to rms signal.
  return FakeData.seismic2d2012A(n1,n2,str,fau,ero,amp,noi)

def like(f):
  return zerofloat(len(f[0]),len(f))

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

def setup():
  global s1,s2,n1,n2,dataDir
  if seisImage=='fakeA':
    s1,s2 = Sampling(201),Sampling(201)
    dataDir = None
  if seisImage=='fakeB':
    s1,s2 = Sampling(401),Sampling(801)
    dataDir = None
  if seisImage=='fakeC':
    s1,s2 = Sampling(301),Sampling(501)
    dataDir = None
  elif seisImage=='f3d':
    s1,s2 = Sampling(462),Sampling(951)
    dataDir = '/data/seis/f3/'
  n1,n2 = s1.count,s2.count

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,cmap=gray,cmin=0,cmax=0,perc=100,cbar=None,name=None):
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
    setup()
    start = time.time()
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02ds'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
