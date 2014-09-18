"""
Correlated leakage method for 2D images
"""
from imports import *
import ipf

noise = 0.00 # noise-to-signal ratio
#noise = 0.05 # noise-to-signal ratio
n1,n2 = 800,200
shiftMax = 2.0 # maximum shift in samples (true shift)
shiftTrial = 0.1 # user-specified trial shift (see abstract)
#window1,window2 = 10,1 # window half widths
window1,window2 = 12,2 # window half widths
#window1,window2 = 24,4 # window half widths
#window1,window2 = 20,5 # window half widths
#window1,window2 = 20,1 # window half widths
#window1,window2 = 50,5 # window half widths
#window1,window2 = 50,10 # window half widths
#window1,window2 = n1,1 # window half widths

#############################################################################

timer = Timer()
def main(args):
  f,g,s = getImagesAndShifts()
  pixels(f,sperc=99.5,title='f')
  pixels(g,sperc=99.5,title='g')
  smin,smax = -1.1*shiftMax,1.1*shiftMax
  pixels(s,cmap=rwb,cmin=smin,cmax=smax,title='true')
  cls = findClShifts(f,g); pixels(cls,cmap=rwb,cmin=smin,cmax=smax,title='cl')
  ccs = findCcShifts(f,g); pixels(ccs,cmap=rwb,cmin=smin,cmax=smax,title='cc')

def findCcShifts(f,g):
  """Local crosscorrelation"""
  sigma1 = max(sqrt(window1*(window1+1.0)/3.0),1.0)
  sigma2 = max(sqrt(window2*(window2+1.0)/3.0),1.0)
  lsf = LocalShiftFinder(sigma1,sigma2)
  u = zerofloat(n1,n2)
  min1,max1 = -int(2.0*shiftMax),int(2.0*shiftMax)
  lsf.find1(min1,max1,f,g,u)
  return mul(-1.0,u)

def findClShifts(f,g):
  """Correlated leakage method"""
  null = -99.99
  s = fillfloat(shiftTrial,n1,n2)
  fs = applyShifts(f,s)
  gs = applyShifts(g,s)

  # Points for cross plot
  x = zerofloat(n1,n2)
  y = zerofloat(n1,n2)
  class Loop(Parallel.LoopInt):
    def compute(self,i2):
      for i1 in range(n1):
        x[i2][i1] = 0.5*(fs[i2][i1]+gs[i2][i1])-0.5*(f[i2][i1]+g[i2][i1])
        y[i2][i1] = g[i2][i1]-f[i2][i1]
  Parallel.loop(n2,Loop())

  # Linear regressions
  r = zerofloat(n1,n2)
  class Loop(Parallel.LoopInt):
    def compute(self,i2):
      #print 'i2 =',i2
      for i1 in range(n1):
        f1,l1 = max(i1-window1,0),min(i1+window1,n1-1)
        f2,l2 = max(i2-window2,0),min(i2+window2,n2-1)
        r[i2][i1] = shiftTrial*linearRegression(x,y,f1,l1,f2,l2)
  timer.start('find shifts')
  Parallel.loop(n2,Loop())
  timer.stop('find shifts')
  return r

def getImagesAndShifts():
  random = Random(12345)
  f = getImage()
  s = getShifts()
  g = applyShifts(f,s)
  if noise>0.0:
    return addNoise(noise,random,f),addNoise(noise,random,g),s
  else:
    return f,g,s

def addNoise(nrms,random,f):
  nrms = nrms*max(abs(f))
  g = mul(2.0,sub(randfloat(random,n1,n2),0.5))
  RecursiveGaussianFilter(1.0).apply10(g,g)
  frms = sqrt(sum(mul(f,f))/n1/n2)
  grms = sqrt(sum(mul(g,g))/n1/n2)
  g = mul(g,nrms*frms/grms)
  return add(f,g)

def applyShifts(f,s):
  g = zerofloat(n1,n2)
  si = SincInterpolator()
  rr = rampfloat(0.0,1.0,n1)
  for i2 in range(n2):
    si.interpolate(n1,1.0,0.0,f[i2],n1,add(rr,s[i2]),g[i2])
  return g

def getImage():  
  g = zerofloat(n1,n2)
  f = ipf.FakeData.seismicAndSlopes2d2014A(0.0)[0]
  rr = rampfloat(0.0,501.0/(n1-1),n1)
  si = SincInterpolator()
  for i2 in range(n2):
    si.interpolate(501,1.0,0.0,f[i2],n1,rr,g[i2])
  return g

def xgetImage():  
  structure = 20.0 # maximum vertical shift when adding structure
  fault = 10.0 # maximum displacement when adding a fault.
  erosion = 0.25*n1/4 # sample at which to create an erosional unconformity.
  amplitude = 1.0 # minimum scale factor when scaling amplitudes.
  f = FakeData.seismic2d2012A(
    n1/2,n2,structure,fault,erosion,amplitude,0.0)
  g = zerofloat(n1,n2)
  rr = rampfloat(0.0,0.25,n1)
  si = SincInterpolator()
  for i2 in range(n2):
    si.interpolate(n1/4,1.0,0.0,f[i2],n1,rr,g[i2])
  return g

def getShifts():
  #return getShiftsConstant()
  #return getShiftsGaussian()
  return getShiftsRandom()
def getShiftsConstant():
  return fillfloat(shiftMax,n1,n2)
def getShiftsGaussian():
  sigma1 = 0.1*n1
  sigma2 = 0.1*n2
  s = zerofloat(n1,n2)
  s[n2/2][n1/2] = 1.0
  RecursiveGaussianFilter(sigma1).apply0X(s,s)
  RecursiveGaussianFilter(sigma2).applyX0(s,s)
  mul(shiftMax/max(s),s,s)
  return s
def getShiftsRandom():
  sigma1 = 0.02*n1
  sigma2 = 0.02*n2
  s = sub(randfloat(Random(12345),n1,n2),0.5)
  RecursiveGaussianFilter(sigma1).apply0X(s,s)
  RecursiveGaussianFilter(sigma2).applyX0(s,s)
  mul(shiftMax/max(abs(s)),s,s)
  return s

def linearRegression(x,y,f1,l1,f2,l2):
  """Fits a line to (x,y).
  Parameters:
    x - x-coordinates.
    y - y-coordinates.
  Returns:
    The sampled best fit line in an array the same size as x and y.
  """
  n = len(x)
  sumx = 0.0 # sum_n(x_i)
  sumy = 0.0 # sum_n(y_i)
  sumxx = 0.0 # sum_n(x_i*x_i)
  sumxy = 0.0 # sum_n(x_i*y_i)
  count = 0
  for i2 in range(f2,l2+1):
    for i1 in range(f1,l1+1):
      xi = x[i2][i1]
      yi = y[i2][i1]
      sumx += xi
      sumy += yi
      sumxx += xi*xi
      sumxy += xi*yi
      count += 1
  beta = (sumxy-sumx*sumy/count)/(sumxx-sumx*sumx/count)
  alpha = (sumy-beta*sumx)/count
  #z = zerofloat(n)
  #for i in range(n):
  #  if null is None or x[i]!=null:
  #    z[i] = alpha+beta*x[i]
  #print 'slope =',beta
  #return z
  return beta

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def pixels(f,cmap=gray,cmin=0,cmax=0,perc=100,sperc=None,cbar=None,title=None):
  panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  cb = panel.addColorBar()
  if cbar:
    cb.setLabel(cbar)
  cb.setWidthMinimum(100)
  pixel = panel.addPixels(f)
  pixel.setColorModel(cmap)
  if cmin<cmax:
    pixel.setClips(cmin,cmax)
  if perc<100:
    pixel.setPercentiles(100-perc,perc)
  if sperc is not None:
    clips = Clips(100.0-sperc,sperc,f)
    clip = max(abs(clips.getClipMin()),abs(clips.getClipMax()))
    pixel.setClips(-clip,clip)
  pixel.setInterpolation(PixelsView.Interpolation.LINEAR)
  frame = PlotFrame(panel)
  frame.setFontSizeForSlide(1.5,1.5)
  if title:
    frame.setTitle(title)
  frame.setSize(1100,1000)
  frame.setVisible(True)

#############################################################################
import sys,time
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
