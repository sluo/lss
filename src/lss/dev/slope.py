"""
Testing different methods for slope estimation
"""
from imports import *
from ldf import BilateralFilter

pngDir = None
#pngDir = '/Users/sluo/Desktop/'
#seisImage = 'fakeA'
#seisImage = 'fakeB'
#seisImage = 'fakeC'
#seisImage = 'art'
seisImage = 'f3d'
#seisImage = 'tp'

#############################################################################

def main(args):
  #goBlf()
  #goGradientOuterProducts()
  #goNormals()
  goFlatten()

def goNormals():
  f = getImage(); plot(f,perc=99.5,name='f')
  goLof(f,24.0,24.0)
  goLofX(f) # lof using lsf
  goNormalABC(f,niter=1) # van de Weijer
  goIsotropicNst(f) # Brox
  goAnisotropicNst(f) # Brox

def goFlatten():
  f = getImage(); plot(f,perc=99.5,name='f')
  #p2 = goLof(f,8.0,2.0)
  #p2 = goLofX(f)
  #p2 = goNormalABC(f,niter=1)
  #p2 = goIsotropicNst(f)
  #p2 = goAnisotropicNst(f)
  p2 = goDsf(f) 
  fl = FlattenerS(6.0,6.0)
  s = fl.findShifts(p2)
  g = fl.applyShifts(f,s)
  plot(g,perc=99.5,name='g')
  plot(p2,cmap=jet,name='p2')

#############################################################################

def goGradientOuterProducts():
  f = getImage(); plot(f,perc=99.5,name='f')
  g11,g12,g22 = getGradientOuterProducts(f)
  plot(g11,cmap=jet,name='g11')
  plot(g12,cmap=jet,name='g12')
  plot(g22,cmap=jet,name='g22')

def goBlf():
  """Bilateral filter example."""
  f = fillfloat(-1.0,n1,n2)
  for i2 in range(n2/2,n2):
    fill(1.0,f[i2])
  f = addNoise(0.25,f)
  sigma = 24.0
  sigmaS = sigma
  sigmaR = computeSigmaR(f)
  gblf,grgf = like(f),like(f)
  BilateralFilter(sigmaS,sigmaR).apply(f,gblf)
  RecursiveGaussianFilter(sigma).apply00(f,grgf)
  cmin,cmax = -1.5,1.5
  plot(f,cmin=cmin,cmax=cmax,name='f')
  plot(gblf,cmin=cmin,cmax=cmax,name='blf')
  plot(grgf,cmin=cmin,cmax=cmax,name='rgf')

  # Tukey's biweight function
  dt = 0.05*sigmaR
  nt = 1+int((4.0*sigmaR)/dt)
  ft = -2.0*sigmaR
  st = Sampling(nt,dt,ft)
  f = zerofloat(nt)
  for it in range(nt):
    t = ft+it*dt
    if abs(t)<sigmaR:
      s = t/sigmaR
      r = 1.0-s*s
      f[it] = r*r
  SimplePlot.asPoints(st,f)

#############################################################################

def goIsotropicNst(f):
  """Isotropic nonlinear structure tensors."""
  g11,g12,g22 = getGradientOuterProducts(f)
  d = getDiffusionScalars(f)
  #plot(d,cmap=jet,name='diffusivity')
  t11,t12,t22 = like(g11),like(g12),like(g22)
  lsf = LocalSmoothingFilter()
  c = 2.0
  lsf.apply(c,d,g11,t11)
  lsf.apply(c,d,g12,t12)
  lsf.apply(c,d,g22,t22)
  u1,u2,_,_,_,_ = getEigenFromTensors(t11,t12,t22)
  #plot(u1,cmap=jet,name='u1 (isotropic nst)')
  plot(u2,cmap=jet,cmin=-1,cmax=1,name='u2 (isotropic nst)')
  return getSlopesFromNormals(u1,u2)

def getDiffusionScalars(f):
  """Scalar diffusivity."""
  g11,g12,g22 = getGradientOuterProducts(f)
  h11 = getGradientInnerProducts(g11)
  h12 = getGradientInnerProducts(g12)
  h22 = getGradientInnerProducts(g22)
  h = getGradientInnerProducts(f) # extra channel
  #plot(h11,cmap=jet,name='h11')
  #plot(h12,cmap=jet,name='h12')
  #plot(h22,cmap=jet,name='h22')
  return diffusivity(add(add(h11,h12),h22))
  #return diffusivity(add(add(h11,h12),add(h22,h)))
    
def goAnisotropicNst(f):
  """Anisotropic nonlinear structure tensors."""
  g11,g12,g22 = getGradientOuterProducts(f)
  t11,t12,t22 = like(g11),like(g12),like(g22)
  t = getDiffusionTensors(f)
  lsf = LocalSmoothingFilter()
  c = 2.0
  lsf.apply(t,c,g11,t11)
  lsf.apply(t,c,g12,t12)
  lsf.apply(t,c,g22,t22)
  u1,u2,_,_,_,_ = getEigenFromTensors(t11,t12,t22)
  #plot(u1,cmap=jet,name='u1 (anisotropic nst)')
  plot(u2,cmap=jet,cmin=-1,cmax=1,name='u2 (anisotropic nst)')
  return getSlopesFromNormals(u1,u2)

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
  #plot(eu,cmap=jet,name='eu')
  #plot(ev,cmap=jet,name='ev')
  #plot(du,cmap=jet,name='du')
  #plot(dv,cmap=jet,name='dv')
  #t11,t12,t22 = getTensorsFromEigen(u1,u2,v1,v2,du,dv)
  #return t11,t12,t22
  #return EigenTensors2(u1,u2,du,dv)
  return EigenTensors2(v1,v2,du,dv) # TODO: check this

def diffusivity(g,eps=0.01):
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
      u1i = v[0][0]
      u2i = v[0][1]
      if u1i<0.0:
        u1i = -u1i
        u2i = -u2i
      u1[i2][i1] = u1i
      u2[i2][i1] = u2i
      v1[i2][i1] = -u2i
      v2[i2][i1] = u1i
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

def getAnglesFromGradients(g1,g2):
  """Angle between gradient vector and vertical."""
  n1,n2 = len(g1[0]),len(g1)
  a = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      g1i = g1[i2][i1]
      g2i = g2[i2][i1]
      den = sqrt(g1i*g1i+g2i*g2i)
      a[i2][i1] = g1i/den # cos of angle
      #a[i2][i1] = acos(g1i/den) # angle
  a = mul(a,a)
  #plot(a,cmap=jet,name='a')
  #plot(g1,cmap=jet,name='g1')
  #plot(g2,cmap=jet,name='g2')
  return a

def goLof(f,sigma1=12.0,sigma2=12.0):
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof = LocalOrientFilter(sigma1,sigma2)
  lof.applyForNormalLinear(f,u1,u2,el)
  #plot(u1,cmap=jet,name='u1 (lof)')
  plot(u2,cmap=jet,cmin=-1,cmax=1,name='u2 (lof)')
  #plot(el,cmap=jet,name='el')
  return getSlopesFromNormals(u1,u2)

def goLofX(f):
  """Estimates normal vectors using a modified lof."""
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof = LocalOrientFilter(12.0,12.0)
  lof.applyForNormalLinear(f,u1,u2,el)
  #el = pow(el,8)
  c = 100.0
  lofx = LocalOrientFilterX(c,el)
  lofx.applyForNormal(f,u1,u2)
  #plot(u1,cmap=jet,name='u1 (lofx)')
  plot(u2,cmap=jet,cmin=-1,cmax=1,name='u2 (lofx)')
  return getSlopesFromNormals(u1,u2)

def goNormalABC(f,niter=1):
  """Estimate normal vectors using a bilateral filter."""
  g1,g2 = getGradients(f)
  g11 = mul(g1,g1)
  g12 = mul(g1,g2)
  g22 = mul(g2,g2)
  g = getAnglesFromGradients(g1,g2)
  def goOnce(v1,v2):
    t11,t12,t22 = like(g11),like(g12),like(g22)
    v = getAnglesFromGradients(v1,v2)
    sigmaS = 24.0 # spatial filter sigma
    sigmaR = computeSigmaR(v1)
    blf = BilateralFilter(sigmaS,sigmaR)
    blf.setType(BilateralFilter.Type.TUKEY_ANGLE)
    blf.applyABC(g,v,g11,t11)
    blf.applyABC(g,v,g12,t12)
    blf.applyABC(g,v,g22,t22)
    u1,u2,_,_,_,_ = getEigenFromTensors(t11,t12,t22)
    return u1,u2
  for ii in range(niter):
    if ii==0:
      u1,u2 = like(f),like(f)
      #LocalOrientFilter(8.0,8.0).applyForNormal(f,u1,u2)
      #LocalOrientFilter(12.0,12.0).applyForNormal(f,u1,u2)
      LocalOrientFilter(24.0,24.0).applyForNormal(f,u1,u2)
      #LocalOrientFilter(32.0,8.0).applyForNormal(f,u1,u2)
      #LocalOrientFilter(64.0,32.0).applyForNormal(f,u1,u2)
      cmin1,cmax1 = getClips(u1)
      cmin2,cmax2 = getClips(u2)
      #plot(u1,cmap=jet,cmin=cmin1,cmax=cmax1,name='u1 (lof)') # lof
      #plot(u2,cmap=jet,cmin=cmin2,cmax=cmax2,name='u2 (lof)') # lof
    u1,u2 = goOnce(u1,u2)
  u1,u2 = getNormalsFromEigenvectors(u1,u2)
  #plot(f,name='f')
  #plot(u1,cmap=jet,cmin=cmin1,cmax=cmax1,name='u1 (blf)') # blf
  plot(u2,cmap=jet,cmin=-1,cmax=1,name='u2 (blf)') # blf
  return getSlopesFromNormals(u1,u2)

def computeSigmaR(g):
  sigmaR = 0.25*(Quantiler.estimate(0.75,g)-Quantiler.estimate(0.25,g))
  print 'sigmaR=%f'%sigmaR
  return sigmaR

def getNormalsFromEigenvectors(u1,u2):
  """Postive u1 only."""
  n1,n2 = len(u1[0]),len(u1)
  for i2 in range(n2):
    for i1 in range(n1):
      u1i = u1[i2][i1]
      u2i = u2[i2][i1]
      if u1i<0.0:
        u1i = -u1i
        u2i = -u2i
      u1[i2][i1] = u1i
      u2[i2][i1] = u2i
  return u1,u2

def goDsf(f):
  """Dynamic slope finder."""
  dsf = DynamicSlopeFinder(12.0)
  p2 = like(f)
  dsf.findSlopes(f,p2)
  plot(p2,cmap=jet,name='p2 (dsf)')
  return p2

def getClips(f):
  minf,maxf = min(f),max(f)
  rangef = 0.1*(maxf-minf)
  cmin = minf-rangef
  cmax = maxf+rangef
  return cmin,cmax

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

def getImage():
  if seisImage=='fakeA':
    return fakeImageA()
  elif seisImage=='fakeB':
    return fakeImageB()
  elif seisImage=='fakeC':
    return fakeImageC()
  elif seisImage=='f3d':
    f = read('f3d400')
    return div(f,max(abs(f)))
  elif seisImage=='tp':
    f = read('tp73')
    return div(f,max(abs(f)))
  elif seisImage=='art':
    f = read('art')
    return addNoise(0.5,f)

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
  elif seisImage=='art':
    s1,s2 = Sampling(600),Sampling(600)
    dataDir = '/Users/sluo/Desktop/'
  elif seisImage=='tp':
    s1,s2 = Sampling(251),Sampling(357)
    dataDir = '/data/seis/tp/'
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
  #pix.setInterpolation(PixelsView.Interpolation.LINEAR)
  pix.setInterpolation(PixelsView.Interpolation.NEAREST)
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
