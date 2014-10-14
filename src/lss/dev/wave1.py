from imports import *
from fd import *

True = 1
False = 0

#############################################################################
# parameters

fontSize = 24
width = 800
height = 600
#pngDir = "./png"
pngDir = None

#############################################################################
# functions

def main(args):
  tmp()
  #test1()
  return

def tmp():

  findShifts = True

  d,v = makeModel()
  nx = len(d)
  vs = [v,fillfloat(1.0,nx)]
  #SimplePlot.asPoints(d)
  SimplePlot.asPoints(v)

  dt = 0.5
  dx = 1.0
  nt = 2001
  f = zerofloat(nx,nt,2)
  for i in range(2):
    v = vs[i]
    wave1 = Wave1(Wave1.Method.SYMMETRIC,dt,dx,d,v)
    for it in range(nt):
      fi = wave1.step(1) # step once
      copy(fi,f[i][it])
    #SimplePlot.asPixels(f[i])
    
  uo = copy(1+nx/2,nt,f[0]) # observed
  us = copy(1+nx/2,nt,f[1]) # simulated
  uo = addRandomNoise(1.0e1,uo)
  us = addRandomNoise(1.0e1,us)
  plot(uo)
  plot(us)

  if findShifts:
    shiftMax = 50
    strainMax = 0.25
    dw = DynamicWarping(-shiftMax,shiftMax)
    dw.setStrainMax(strainMax)
    dw.setShiftSmoothing(1.0,1.0); # shift smoothing
    dw.setErrorSmoothing(2) # number of error smoothings
    s = dw.findShifts(uo,us)
    plot(s,jet)

gray = ColorMap.GRAY
jet = ColorMap.JET
def plot(x,cmap=gray):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  pv = sp.addPixels(x)
  pv.setColorModel(cmap)

def xplot(w,png):
  nw = len(w)
  nx = len(w[0])
  colors = [Color.RED,Color.GREEN,Color.BLUE]
  pp = PlotPanel(1,1,
    PlotPanel.Orientation.X1RIGHT_X2UP,
    PlotPanel.AxesPlacement.LEFT_BOTTOM)
  pp.setVLimits(-30.0,30.0)
  dx = 1.0
  fx = -0.5*(nx-1)*dx
  sx = Sampling(nx,dx,fx)
  for iw in range(nw):
    pv = pp.addPoints(sx,w[iw])
    pv.setLineColor(colors[iw%3])
    pv.setLineWidth(3)
  pf = PlotFrame(pp)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  pf.setFontSizeForSlide(1.0,0.9)
  pf.setSize(width,height)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(720,3.3,pngDir+"/"+png+".png")


def addRandomNoise(r,x):
  n1,n2 = len(x[0]),len(x)
  # xmax = max(abs(p)) # peak signal
  xrms = sqrt(sum(mul(x,x))/n1/n2) # rms of signal
  s = sub(randfloat(Random(3),n1,n2),0.5)
  #s = sub(randfloat(n1,n2),0.5)
  sigma = 1.0
  RecursiveGaussianFilter(sigma).apply00(s,s); # bandlimited noise
  srms = sqrt(sum(mul(s,s))/n1/n2) # rms of noise
  return add(mul(xrms/(srms*r),s),x); # r = rms-signal / rms-noise

def makeModel():
  p = 0.1 # max slowness change
  nx = 2001
  d = fillfloat(1.0,nx)
  s = zerofloat(nx) # slowness
  s[  nx/6] =  1.0
  s[2*nx/6] = -1.0
  #RecursiveGaussianFilter(nx/32).apply0(v,v)
  RecursiveGaussianFilter(32.0).apply0(s,s)
  #div(v,4.0*max(v),v)
  mul(s,p/max(s),s)
  add(s,1.0,s)
  v = zerofloat(nx)
  div(1.0,s,v)
  #add(v,1.0,v)
  return d,v

def xmakeModel(dl,dr):
  nx = 801
  d = zerofloat(nx)
  v = zerofloat(nx)
  for ix in range(nx):
    v[ix] = 1.0
    if ix<3*nx/4:
      d[ix] = dl
    else:
      d[ix] = dr
  return d,v

def test1():
  nmodel = 4
  lmodel = ["05","03","02","01"]
  dllist = [1.0,1.0,1.0,1.0] # densities left
  drlist = [0.5,0.3,0.2,0.1] # densities right
  for imodel in range(nmodel):
    dl = dllist[imodel]
    dr = drlist[imodel]
    d,v = makeModel(dl,dr)
    nx = len(d)
    rc = (dr-dl)/(dr+dl)
    tc = 1.0+rc
    print "rc =",rc,"  tc =",tc
    dt = 0.5
    dx = 1.0
    nt = 5
    mt = 125
    methods = [
      Wave1.Method.SYMMETRIC,
      Wave1.Method.PRODUCT1,
      Wave1.Method.PRODUCT2
    ]
    nmethod = len(methods)
    f = zerofloat(nx,nmethod,nt)
    for imethod in range(nmethod):
      wave1 = Wave1(methods[imethod],dt,dx,d,v)
      for it in range(nt):
        fi = wave1.step(mt)
        ft = copy(nx/2,fi)
        copy(fi,f[it][imethod])
        fmax = max(ft)
        #print "fmax =",fmax,"  fr =",fmax*rc,"  ft =",fmax*tc
    prefix = "wave1"+lmodel[imodel]
    if imodel==0:
      for it in range(nt):
        png = prefix+str(it)
        plot(f[it],png)
    else:
      it = nt-1
      png = prefix+str(it)
      plot(f[it],png)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
