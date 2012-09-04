"""
Reflection FWI
"""
from imports import *
from fault import DynamicWarping

pngDir=None
#pngDir='/Users/sluo/Desktop/'
#pngDir='/Users/sluo/Dropbox/png/'
datDir='./dat/'
sz = Sampling(201,0.015,0.0)
sx = Sampling(801,0.015,0.0)
st = Sampling(1504,0.005,0.0)
#st = Sampling(804,0.005,0.0)
nt,nx,nz = st.count,sx.count,sz.count
dt,dx,dz = st.delta,sx.delta,sz.delta
ft,fx,fz = st.first,sx.first,sz.first
fpeak = 10.0

#kxs = [nx/2]
#kxs = rampint(5,5,(nx-5)/5)
#kxs = rampint(10,10,(nx-10)/10)
#kxs = rampint(25,25,(nx-25)/25)
kxs = rampint(50,50,(nx-50)/50)
#kxs = rampint(100,100,(nx-100)/100)
kzs = fillint(0,len(kxs))
kzr,kxr = fillint(0,nx),rampint(0,1,nx)
ns = len(kxs)

#############################################################################

def main(args):
  v,_,_ = makeModel(); plot(v,jet,name='v')
  #goIterations()

def goIterations():
  niter = 10
  global a,b
  v,c,r = makeModel()
  a = AcousticWavefield(sz,sx,st)
  b = BornWavefield(sz,sx,st,r)

  #def findUpdate(c,g,cmin,cmax):
  #  #global misfit
  #  """Line search."""
  #  ks = ns/2 # which source to compare misfit
  #  pmax = 0.1 # maximum percent slowness change
  #  def misfitDecreased(cp):
  #    global misfit
  #    ds = makeData(kzs[ks],kxs[ks],cp)
  #    m = findMisfit(dos[ks],ds)
  #    print m
  #    if m<misfit:
  #      misfit = m
  #      print 'decreased'
  #      return True
  #    else:
  #      print 'increased'
  #      return False
  #  s = div(1.0,c)
  #  h = div(1.0,g)
  #  savg = sum(s)/nx/nz
  #  hrms = rms(h)
  #  #maxUpdate =  pmax/(savg*(1.0-pmax))
  #  #minUpdate = -pmax/(savg*(1.0+pmax))
  #  for i in range(3):
  #    p = pmax*pow(0.5,i)
  #    a = p*savg/hrms
  #    sp = add(s,mul(a,h))
  #    cp = div(1.0,sp)
  #    plot(cp,jet)
  #    if misfitDecreased(cp):
  #      return clip(cp,min(v),max(v))
  #  return clip(cp,min(v),max(v))

  def findUpdate(v,g,vmin,vmax):
    """Line search."""
    sl = 0.050 # initial step length (km/s)
    ks = ns/2 # which source to compare misfit
    pmax = 0.1 # maximum percent slowness change
    g = mul(pmax*vmin,g) # scaled gradient
    def misfitDecreased(vp):
      global misfit
      ds = makeData(kzs[ks],kxs[ks],vp)
      m = findMisfit(dos[ks],ds)
      print m
      if m<misfit:
        misfit = m
        print 'decreased'
        return True
      else:
        print 'increased'
        return False
    ntry = 5 # max number of times step length is halved
    for i in range(ntry):
      print sl
      vp = add(v,mul(sl,g))
      if i==ntry-1 or misfitDecreased(vp):
        #return clip(vp,vmin,vmax)
        return vp
      else:
        sl *= 0.5

  dos = makeObservedData(v)
  for iiter in range(niter):
    print 'ITER %d/%d'%(iiter+1,niter)
    g = makeGradient(dos,c)
    c = findUpdate(c,g,min(v),max(v))
    gn = 'g'+str(iiter+1)
    cn = 'c'+str(iiter+1)
    write(gn,g)
    write(cn,c)
    plot(g,jet,name=gn)
    plot(c,jet,name=cn)

def clip(c,cmin,cmax):
  n1,n2 = len(c[0]),len(c)
  for i2 in range(n2):
    for i1 in range(n1):
      if c[i2][i1]>cmax:
        c[i2][i1] = cmax
      if c[i2][i1]<cmin:
        c[i2][i1] = cmin
  return c

#############################################################################
# Gradient

def makeGradient(dos,c):
  global misfit
  g = zeros(nz,nx)
  for i in range(ns):
    print "\n%d/%d\n---"%(i+1,ns)
    kzsrc,kxsrc = kzs[i],kxs[i]
    ds,us = makeDataAndWavefield(kzsrc,kxsrc,c)
    ad = makeAdjointSource(dos[i],ds)
    if i==ns/2:
      misfit = findMisfit(dos[i],ds)
      print misfit
    a.backPropagate(ad,kzr,kxr,c)
    ur = a.getWavefield()
    add(imagingCondition(us,ur),g,g)
  g = div(g,mul(c,mul(c,c)))
  RecursiveGaussianFilter(1.0).apply00(g,g)
  r = sum(abs(g))/nx/nz
  return div(g,r)

def imagingCondition(us,ur):
  r = zeros(nz,nx)
  for it in range(1,nt-1):
    rs = us[it]
    rr = add(add(ur[it-1],ur[it+1]),mul(-2.0,ur[it]))
    #add(r,mul(rs,rr),r)
    sub(r,mul(rs,rr),r)
  #return laplacian(r)
  return r

def makeAdjointSource(do,ds):
  ad = zeros(nt,nx)
  s = findShifts(do,ds)
  #h = applyShifts(s,do)
  h = ds # XXX
  w = fillfloat(1.0,nx) # Taper
  #RecursiveGaussianFilter(8.0).apply0(w,w)
  #RecursiveGaussianFilter(nx/2).apply0(w,w)
  RecursiveGaussianFilter(nx/4).apply0(w,w)
  sub(w,w[0],w)
  for ix in range(nx):
    for it in range(1,nt-1):
      hm = h[ix][it-1]
      hp = h[ix][it+1]
      ad[ix][it] = w[ix]*(s[ix][it]*(hp-hm))
  return ad

def makeObservedData(c):
  dos = zeros(nt,nx,ns)
  for i in range(ns):
    dos[i] = makeData(kzs[i],kxs[i],c)
  #write('dos',dos)
  return dos

def makeData(kzsrc,kxsrc,c):
  b.forwardPropagate(fpeak,kzsrc,kxsrc,c)
  return b.getWavefield(kzr,kxr)

def makeWavefield(kzsrc,kxsrc,c):
  b.forwardPropagate(fpeak,kzsrc,kxsrc,c)
  return b.getWavefield()

def makeDataAndWavefield(kzsrc,kxsrc,c):
  b.forwardPropagate(fpeak,kzsrc,kxsrc,c)
  return b.getWavefield(kzr,kxr),b.getWavefield()

def findShifts(f,g):
  s = findWarping(g,f)
  return s

def findMisfit(f,g):
  s = findShifts(f,g)
  return sum(mul(s,s))/nt/nx

#############################################################################
# Dynamic warping

def findWarping(f,g):
  global smoothShifts; smoothShifts = True
  global smoothSigma; smoothSigma = 2.0
  shiftMax = 16
  stretchMax = 0.25
  dw = DynamicWarping(-shiftMax,shiftMax)
  dw.setStretchMax(stretchMax)
  dw = DynamicWarping(-shiftMax,shiftMax)
  dw.setStretchMax(stretchMax)
  e = dw.computeErrors(f,g)
  u = shifts121(dw,e)
  #h = applyShifts(g,u)
  return u

def shifts121(dw,e):
  e = dw.accumulate1(e)
  e = normalize(e)
  e = dw.accumulate2(e)
  e = normalize(e)
  d = dw.accumulateForward1(e)
  u = dw.findShiftsReverse1(d,e)
  return smooth(u)

def smooth(u):
  v = copy(u)
  if smoothShifts:
    RecursiveGaussianFilter(smoothSigma).apply00(u,v)
  return v

def normalize(x):
  xmin = min(x)
  xmax = max(x)
  return mul(sub(x,xmin),1.0/(xmax-xmin))

def applyShifts(g,u):
  n1,n2 = len(u[0]),len(u)
  si = SincInterpolator()
  si.setUniformSampling(n1,1.0,0.0)
  h = copy(g)
  r = rampfloat(0.0,1.0,n1)
  for i2 in range(n2):
    t = add(r,u[i2])
    si.setUniformSamples(g[i2])
    si.interpolate(n1,t,h[i2])
  return h
  
#############################################################################
# Model

def makeModel(sigma=8.0):
  def checkers():
    w = 50 # checker width in samples
    ps = 0.05 # percent slowness anomaly
    vp = 2.0/(1.0-ps)
    vm = 2.0/(1.0+ps)
    v = fill(2.0,nz,nx)
    c = fill(2.0,nz,nx)
    r = zeros(nz,nx)
    def setv(ix,v1,v2):
      z0 = nz/3
      z1,z2 = z0,z0+w
      count = 0
      while z1<nz:
        for iz in range(z1,z2):
          v[ix][iz] = v1 if count%2==0 else v2
        r[ix][z1-1] = -0.5
        r[ix][z1  ] =  1.0
        r[ix][z1+1] = -0.5
        """
        num = v[ix][z1+1]-v[ix][z1-1]
        den = v[ix][z1+1]+v[ix][z1-1]
        r[ix][z1  ] = num/den
        """
        z1 = z2
        z2 = min(z2+w,nz)
        count += 1
    x1,x2 = 0,w/2
    count = 0
    while x1<nx:
      for ix in range(x1,x2):
        if count%2==0:
          setv(ix,vp,vm)
        else:
          setv(ix,vm,vp)
      x1 = x2
      x2 = min(x2+w,nx)
      count += 1
    return v,c,r
  return checkers()

#############################################################################
# Utilities

def smoothSlowness(sigma,v):
  v = div(1.0,v)
  ref = RecursiveExponentialFilter(sigma)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref.apply(v,v)
  return div(1.0,v)

def laplacian(x):
  n1,n2 = len(x[0]),len(x)
  x1 = zeros(n1,n2)
  x2 = zeros(n1,n2)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply2X(x,x1)
  rgf.applyX2(x,x2)
  return add(x1,x2)

def zeros(n1,n2=None,n3=None):
  if n3 and n2:
    return zerofloat(n1,n2,n3)
  elif n2:
    return zerofloat(n1,n2)
  else:
    return zerofloat(n1)

def fill(v,n1,n2,n3=None):
  if n3:
    return fillfloat(v,n1,n2,n3)
  else:
    return fillfloat(v,n1,n2)

def rms(x):
  n1,n2 = len(x[0]),len(x)
  r = 0.0
  for i2 in range(n2):
    for i1 in range(n1):
      xi = x[i2][i1]
      r += xi*xi
  return sqrt(r/n1/n2)

def read(name,image=None):
  if not image:
    image = zerofloat(nz,nx)
  fileName = datDir+name+'.dat'
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def write(name,image,directory=datDir):
  fileName = directory+name+'.dat'
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()

#############################################################################
# Graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,cmap=gray,sz=None,sx=None,cmin=0,cmax=0,perc=100,
         cbar=None,name=None):
  pan = panel(cbar)
  pix = pan.addPixels(x)
  if sz and sx:
    pix.set(sz,sx,x)
  pix.setColorModel(cmap)
  pix.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pix.setClips(cmin,cmax)
  if perc<100:
    pix.setPercentiles(100-perc,perc)
  size = [1400,min(int(1.05*1400*len(x[0])/len(x)),800)]
  frame(pan,size,name)

def panel(cbar=None):
  p = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  cb = p.addColorBar()
  cb.setWidthMinimum(100)
  if cbar:
    cb.setLabel(cbar)
  return p

def frame(panel,size=None,name=None):
  frame = PlotFrame(panel)
  #frame.setFontSizeForSlide(1.0,1.0)
  if size:
    frame.setSize(size[0],size[1])
  else:
    frame.setSize(1200,1200)
  if name:
    frame.setTitle(name)
  frame.setVisible(True)
  if name and pngDir:
    frame.paintToPng(360,3.0,pngDir+name+'.png')

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
  ipg = ImagePanelGroup(sz,sx,st,image)
  ipg.setColorModel(cmap)
  ipg.setSlices(nz/2,nx/2,nt/2)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
  if perc<100:
    ipg.setPercentiles(100-perc,perc)
  world.addChild(ipg)
  return ipg

def addColorBar(frame,label):
  cbar = ColorBar(label)
  cbar.setFont(cbar.getFont().deriveFont(64.0))
  frame.add(cbar,BorderLayout.EAST)
  #frame.viewCanvas.setBackground(frame.getBackground())
  return cbar

def makeFrame(world,name=None):
  n1,n2,n3 = sz.count,sx.count,st.count
  d1,d2,d3 = sz.delta,sx.delta,st.delta
  f1,f2,f3 = sz.first,sx.first,st.first
  l1,l2,l3 = sz.last,sx.last,st.last
  frame = SimpleFrame(world)
  #frame.setBackground(Color(204,204,204,255))
  frame.setBackground(Color.WHITE)
  if name:
    frame.setTitle(name)
  view = frame.getOrbitView()
  #zscale = 0.6*max(n2*d2,n3*d3)/(n1*d1)
  #view.setAxesScale(1.0,1.0,zscale)
  view.setAxesScale(1.0,1.0,1.0)
  view.setScale(1.3)
  view.setAzimuth(290)
  view.setElevation(50)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(1200,800)
  frame.setVisible(True)
  if pngDir and name:
    frame.paintToFile(pngDir+name+'.png')
  return frame

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
