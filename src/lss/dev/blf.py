"""
Bilateral filtering for slope estimation
"""
from imports import *
from ldf import *
#from lss.dsp import *

pngDir = None
#pngDir = '/Users/sluo/Desktop/'
#seisImage = 'fakeA'
#seisImage = 'fakeB'
seisImage = 'f3d'

#############################################################################

def main(args):
  #f = getImage(); plot(f)
  #goBlf()
  #goLof()
  #goLofX()
  #goGradientProduct()
  #goError()
  #goNormalAB()
  goNormalABC(1)
  #goFlatten()

def goFlatten():
  bilateralFilter = True
  f = getImage()
  if bilateralFilter:
    u1,u2 = goNormalABC(1)
  else:
    u1,u2 = goLof(8.0,2.0)
    plot(u1,cmap=jet,name='u1')
    plot(u2,cmap=jet,name='u2')
  p2 = getSlopesFromNormals(u1,u2)
  fl = FlattenerS(6.0,6.0)
  s = fl.findShifts(p2)
  g = fl.applyShifts(f,s)
  plot(f,perc=99.5,name='f')
  plot(g,perc=99.5,name='g')

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

def goNormalABC(niter=1):
  """Estimate normal vectors using a bilateral filter."""
  f = getImage()
  g1,g2 = goGradient()
  g11 = mul(g1,g1)
  g12 = mul(g1,g2)
  g22 = mul(g2,g2)
  g = getAnglesFromGradients(g1,g2)
  """
  #sigmaS = 6.0 # spatial filter sigma
  sigmaS = 24.0 # spatial filter sigma
  sigmaR = computeSigmaR(g)
  blf = BilateralFilter(sigmaS,sigmaR)
  blf.setType(BilateralFilter.Type.TUKEY_ANGLE)
  """
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
    return getEigenvectors(t11,t12,t22)
  for ii in range(niter):
    if ii==0:
      u1,u2 = goLof(8.0,8.0)
      #u1,u2 = goLof(12.0,12.0)
      #u1,u2 = goLof(16.0,4.0)
      #u1,u2 = goLof(32.0,8.0)
      #u1,u2 = goLof(64.0,32.0)
      cmin1,cmax1 = getClips(u1)
      cmin2,cmax2 = getClips(u2)
      plot(u1,cmap=jet,cmin=cmin1,cmax=cmax1,name='u1 (lof)') # lof
      plot(u2,cmap=jet,cmin=cmin2,cmax=cmax2,name='u2 (lof)') # lof
    u1,u2 = goOnce(u1,u2)
  u1,u2 = getNormalsFromEigenvectors(u1,u2)
  plot(f,name='f')
  plot(u1,cmap=jet,cmin=cmin1,cmax=cmax1,name='u1') # blf
  plot(u2,cmap=jet,cmin=cmin2,cmax=cmax2,name='u2') # blf
  return u1,u2

def goNormalAB():
  """Estimate normal vectors using a bilateral filter."""
  doCholesky = True
  f,r = getImage()
  g1,g2 = goGradient()
  g11 = mul(g1,g1)
  g12 = mul(g1,g2)
  g22 = mul(g2,g2)
  a = getAnglesFromGradients(g1,g2) # for use in range function
  plot(f,name='f')
  #plot(r,name='r')
  plot(a,name='a')
  #plot(g1,cmap=jet,name='g1')
  #plot(g2,cmap=jet,name='g2')
  #plot(g11,cmap=jet,name='g11')
  #plot(g12,cmap=jet,name='g12')
  #plot(g22,cmap=jet,name='g22')
  sigmaS = 24.0 # spatial filter
  #sigmaR = 0.5 # range function
  sigmaR = computeSigmaR(a)
  blf = BilateralFilter(sigmaS,sigmaR)
  blf.setType(BilateralFilter.Type.TUKEY)
  #blf.setType(BilateralFilter.Type.GAUSS)
  #blf.setType(BilateralFilter.Type.HUBER)
  if doCholesky:
    """Smooth Cholesky decomposition."""
    l11,l12,l22 = cholesky(g11,g12,g22) # cholesky decomposition
    t11,t12,t22 = like(l11),like(l12),like(l22)
    #blf.apply(l11,t11)
    #blf.apply(l12,t12)
    #blf.apply(l22,t22)
    blf.applyAB(a,l11,t11)
    blf.applyAB(a,l12,t12)
    blf.applyAB(a,l22,t22)
    """
    rgf = RecursiveGaussianFilter(sigmaS)
    rgf.apply00(l11,t11)
    rgf.apply00(l12,t12)
    rgf.apply00(l22,t22)
    """
    for i2 in range(n2):
      for i1 in range(n1):
        """ LL' """
        t11i = t11[i2][i1]
        t12i = t12[i2][i1]
        t22i = t22[i2][i1]
        t11[i2][i1] = t11i*t11i
        t12[i2][i1] = t11i*t12i
        t22[i2][i1] = t12i*t12i+t22i*t22i
    #plot(l11,cmap=jet,name='l11')
    #plot(l12,cmap=jet,name='l12')
    #plot(l22,cmap=jet,name='l22')
  else:
    """Smooth gradient outer product."""
    t11,t12,t22 = like(g11),like(g12),like(g22)
    #blf.apply(g11,t11)
    #blf.apply(g12,t12)
    #blf.apply(g22,t22)
    blf.applyAB(a,g11,t11)
    blf.applyAB(a,g12,t12)
    blf.applyAB(a,g22,t22)
  u1,u2 = getEigenvectors(t11,t12,t22)
  v1,v2 = goLof()
  cmin1,cmax1 = 0.88,1.01
  cmin2,cmax2 = -0.05,0.5
  #cmin1,cmax1 = getClips(v1)
  #cmin2,cmax2 = getClips(v2)
  #cmin1,cmax1 = min(v1),max(v1);
  #cmin2,cmax2 = min(v2),max(v2);
  plot(v1,cmap=jet,cmin=cmin1,cmax=cmax1,name='v1') # lof
  plot(v2,cmap=jet,cmin=cmin2,cmax=cmax2,name='v2') # lof
  plot(u1,cmap=jet,cmin=cmin1,cmax=cmax1,name='u1') # blf
  plot(u2,cmap=jet,cmin=cmin2,cmax=cmax2,name='u2') # blf

def getClips(f):
  minf,maxf = min(f),max(f)
  rangef = 0.1*(maxf-minf)
  cmin = minf-rangef
  cmax = maxf+rangef
  return cmin,cmax

def getEigenvectors(t11,t12,t22):
  """Eigendecomposition of structure tensors."""
  n1,n2 = len(t11[0]),len(t11)
  u1,u2 = zerofloat(n1,n2),zerofloat(n1,n2)
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
  return u1,u2

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

def computeSigmaR(g):
  sigmaR = 0.25*(Quantiler.estimate(0.75,g)-Quantiler.estimate(0.25,g))
  print 'sigmaR=%f'%sigmaR
  return sigmaR

def cholesky(g11,g12,g22):
  """Cholesky decomposition of structure tensors."""
  n1,n2 = len(g11[0]),len(g11)
  l11 = zerofloat(n1,n2)
  l12 = zerofloat(n1,n2)
  l22 = zerofloat(n1,n2)
  a = zerodouble(2,2)
  l = zerodouble(2,2)
  for i2 in range(n2):
    for i1 in range(n1):
      a[0][0] = g11[i2][i1]
      a[1][0] = a[0][1] = g12[i2][i1]
      a[1][1] = g22[i2][i1]
      d = DMatrix(a)
      DMatrixChd(d).getL().get(l)
      l11[i2][i1] = l[0][0]
      l12[i2][i1] = l[1][0]
      l22[i2][i1] = l[1][1]
  return l11,l12,l22

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

def like(f):
  return zerofloat(len(f[0]),len(f))

def goError():
  """Compare error in gradient vector to structure tensor eigenvectors.
     See van de Weijer, 2005, equation 49."""
  g1,g2 = goGradient()
  u1,u2 = goLof()
  plot(g1,jet,name='g1')
  plot(g2,jet,name='g2')
  plot(u1,jet,name='u1')
  plot(u2,jet,name='u2')
  """
  p = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      p[i2][i1] = -u2[i2][i1]/u1[i2][i1]
  plot(p,jet,name='p')
  """
  t = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      g1i = g1[i2][i1]
      g2i = g2[i2][i1]
      u1i = u1[i2][i1]
      u2i = u2[i2][i1]
      normg = sqrt(g1i*g1i+g2i*g2i)
      normu = sqrt(u1i*u1i+u2i*u2i)
      g11i = g1i*g1i
      g12i = g1i*g2i
      g22i = g2i*g2i
      u11i = u1i*u1i
      u12i = u1i*u2i
      u22i = u2i*u2i
      #t[i2][i1] = g11i+g22i-g11i*u11i-2.0*g12i*u12i-g22i*u22i
      t[i2][i1] = acos(abs((g1i*u1i+g2i*u2i)/(normg*normu)))
      #t[i2][i1] = acos((g1i*u1i+g2i*u2i))
      #t[i2][i1] = abs(atan2(g1i,g2i))
  plot(t,cmap=jet,name='error',perc=99)

def goBlf():
  #f = getImage()
  #f = zerofloat(n1,n2); f[n2/2][n1/2] = 1.0
  f = sub(randfloat(n1,n2),0.5)
  RecursiveGaussianFilter(4.0).apply00(f,f)
  gb = zerofloat(n1,n2)
  gr = zerofloat(n1,n2)
  sigmaSpace = 1.0 # half-width of spatial filter
  sigmaRange = 1.0e6 # half-width of range function
  blf = BilateralFilter(sigmaSpace,sigmaRange)

  sigma1 = 10.0
  sigma2 = 1.0
  et = makeTensors(sigma1,sigma2)
  blf.apply(et,f,gb)
  rgf1 = RecursiveGaussianFilter(sigma1)
  rgf2 = RecursiveGaussianFilter(sigma2)
  rgf1.apply0X(f,gr)
  rgf2.applyX0(gr,gr)

  cmin = min(min(gb),min(gr))
  cmax = max(max(gb),max(gr))
  plot(f)
  plot(gb,cmin=cmin,cmax=cmax)
  plot(gr,cmin=cmin,cmax=cmax)
  print sum(sub(gb,gr))

def makeTensors(sigma1,sigma2):
  u1 = fillfloat(1.0,n1,n2) # 1st component of eigenvector u
  u2 = fillfloat(0.0,n1,n2) # 2nd component of eigenvector u
  au = fillfloat(500.0,n1,n2) # eigenvalues for eigenvector u
  av = fillfloat(50.0,n1,n2) # eigenvalues for eigenvector v
  return EigenTensors2(u1,u2,au,av)

def goLof(sigma1=24.0,sigma2=None):
  f = getImage()
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  if sigma2==None:
    lof = LocalOrientFilter(sigma1)
  else:
    lof = LocalOrientFilter(sigma1,sigma2)
  lof.applyForNormalLinear(f,u1,u2,el)
  el = pow(el,8)
  #plot(f)
  #plot(u1,cmap=jet,name='u1 lof')
  #plot(u2,cmap=jet,name='u2 lof')
  #plot(el,cmap=jet,name='el')
  #plot(div(u2,u1),cmap=jet,perc=100)
  return u1,u2

def goLofX():
  f = getImage()
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof = LocalOrientFilter(24.0,24.0)
  lof.applyForNormalLinear(f,u1,u2,el)
  el = pow(el,64)

  plot(f)
  #plot(u1,cmap=jet,cmin=0.86,cmax=1.00,name='u1')
  #plot(u2,cmap=jet,cmin=0.00,cmax=0.50,name='u2')
  #plot(el,cmap=jet,name='el')
  plot(div(u2,u1),cmap=jet,name='old')

  #el = fillfloat(1.0,n1,n2)
  lofx = LocalOrientFilterX(800.0,el)
  lofx.applyForNormal(f,u1,u2)

  plot(el,cmap=jet,name='el')
  #plot(u1,cmap=jet,cmin=0.86,cmax=1.00,name='u1')
  #plot(u2,cmap=jet,cmin=0.00,cmax=0.50,name='u2')
  plot(div(u2,u1),cmap=jet,name='new')

  return u1,u2

def goGradientProduct():
  """Unsmoothed outer product of image gradients."""
  f = getImage()
  g1,g2 = goGradient()
  """
  for i2 in range(n2):
    for i1 in range(n1):
      g1i = g1[i2][i1]
      g2i = g2[i2][i1]
      g11[i2][i1] = g1i*g1i
      g12[i2][i1] = g1i*g2i
      g22[i2][i1] = g2i*g2i
  """
  g11 = mul(g1,g1)
  g12 = mul(g1,g2)
  g22 = mul(g2,g2)
  #plot(f,name='f')
  #plot(g11,cmap=jet,name='g11')
  #plot(g12,cmap=jet,name='g12')
  #plot(g22,cmap=jet,name='g22')
  return g11,g12,g22

def goGradient():
  f = getImage()
  g1 = zerofloat(n1,n2)
  g2 = zerofloat(n1,n2)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply10(f,g1)
  rgf.apply01(f,g2)
  # TODO: ok to flip orientation?
  """
  for i2 in range(n2):
    for i1 in range(n1):
      g1i = g1[i2][i1]
      g2i = g2[i2][i1]
      if g1i<0.0:
        g1[i2][i1] = -g1i
        g2[i2][i1] = -g2i
  """
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

def getClipsForPlot3(f,r=1.0):
  s12 = slice12(k3,f)
  s13 = slice13(k2,f)
  s23 = slice23(k1,f)
  #cmin = 0.0
  cmin = r*min(min(s12),min(min(s13),min(s23)))
  cmax = r*max(max(s12),max(max(s13),max(s23)))
  return cmin,cmax

def getLinearAlphaColorMap(cmap):
  alpha = rampfloat(0.0,1.0/255.0,256)
  return ColorMap.setAlpha(cmap,alpha)

def getNotchAlphaColorMap(cmap,knotch):
  alpha = fillfloat(1.0,256)
  alpha[knotch] = 0.0
  if knotch>0:
    alpha[knotch-1] = 0.0
  if knotch<255:
    alpha[knotch+1] = 0.0
  return ColorMap.setAlpha(cmap,alpha)

def plot3(f,cmap=gray,cmin=0.0,cmax=0.0,perc=100.0,
          t=None,cmap1=jet,cmin1=0.0,cmax1=0.0,perc1=0.0,
          nearest=False,cbar=None,name=None):
  o = PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT
  a = PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM
  panel = PlotPanelPixels3(o,a,s1,s2,s3,f)
  if t:
    s12 = slice12(k3,t)
    s13 = slice13(k2,t)
    s23 = slice23(k1,t)
    p12 = panel.addPixels(1,0,s1,s2,s12)
    p13 = panel.addPixels(1,1,s1,s3,s13)
    p23 = panel.addPixels(0,0,s3,s2,transpose(s23))
    p12.setInterpolation(PixelsView.Interpolation.NEAREST)
    p13.setInterpolation(PixelsView.Interpolation.NEAREST)
    p23.setInterpolation(PixelsView.Interpolation.NEAREST)
    if cmin1<cmax1:
      p12.setClips(cmin1,cmax1)
      p13.setClips(cmin1,cmax1)
      p23.setClips(cmin1,cmax1)
    p12.setColorModel(cmap1)
    p13.setColorModel(cmap1)
    p23.setColorModel(cmap1)
  if nearest:
    panel.setInterpolation(PixelsView.Interpolation.NEAREST)
  panel.setColorModel(cmap)
  panel.setSlices(k1,k2,k3)
  panel.setLabel1("Time (s)")
  panel.setLabel2("Inline (km)")
  panel.setLabel3("Crossline (km)")
  panel.setLineColor(Color.YELLOW)
  if (cbar):
    panel.addColorBar(cbar)
    panel.setColorBarWidthMinimum(90)
  panel.setBackground(Color(255,254,255,255))
  panel.setInterval1(0.1)
  panel.setInterval2(0.5)
  panel.setInterval3(0.5)
  if cmin<cmax:
    panel.setClips(cmin,cmax)
  if perc<100.0:
    panel.setPercentiles(100-perc,perc)
  #panel.mosaic.setWidthElastic(1,100)
  #panel.mosaic.setHeightElastic(0,205)
  panel.mosaic.setHeightElastic(0,230)
  frame = PlotFrame(panel)
  frame.setFontSizeForSlide(1.0,1.0)
  #frame.setSize(1150,800)
  frame.setSize(1200,795)
  frame.setVisible(True)
  if name and pngDir:
    frame.paintToPng(720,7.0,pngDir+name+".png")
  return frame

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

def array(x1,x2,x3=None,x4=None):
  if x3 and x4:
    return jarray.array([x1,x2,x3,x4],Class.forName('[[[F'))
  elif x3:
    return jarray.array([x1,x2,x3],Class.forName('[[[F'))
  else:
    return jarray.array([x1,x2],Class.forName('[[[F'))

def zeros(n1,n2,n3):
  return zerofloat(n1,n2,n3)

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

def addImage2ToWorld(world,image1,image2):
  ipg = ImagePanelGroup2(s1,s2,s3,image1,image2)
  ipg.setColorModel1(ColorMap.getGray())
  #ipg.setColorModel2(ColorMap.getJet(0.3))
  ipg.setColorModel2(ColorMap.getHue(0.0,20.0,0.3))
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
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.55)
  view.setAzimuth(azimuth)
  view.setElevation(elevation)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3-1.0,f2,f1,l3,l2,l1)))
  frame.viewCanvas.setBackground(frame.getBackground())
  #frame.setSize(1020,750)
  frame.setSize(1800,1200)
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

def setup():
  global s1,s2,n1,n2,dataDir
  if seisImage=='fakeA':
    s1,s2 = Sampling(201),Sampling(201)
    dataDir = None
  if seisImage=='fakeB':
    s1,s2 = Sampling(401),Sampling(801)
    dataDir = None
  elif seisImage=='f3d':
    s1,s2 = Sampling(462),Sampling(951)
    dataDir = '/data/seis/f3/'
  n1,n2 = s1.count,s2.count

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
