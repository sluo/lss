"""
Correlated leakage method
"""
from imports import *

nt = 501
dt = 1.0
freq = 0.10
shiftMax = 2.0
shiftDrv = 1.0 

#nt = 1001
#dt = 0.001
#freq = 50
#shiftMax = 0.003
#shiftDrv = 0.001

#############################################################################

def main(args):
  f,g = getSequences()
  points(f.array(),g.array())
  makeCrossPlot(f,g)

def makeCrossPlot(f,g):
  ti = 1+int(shiftDrv/dt)
  x = zerofloat(nt-ti)
  y = zerofloat(nt-ti)
  for it in range(nt-ti):
    x[it] = 0.5*(
      f.evaluate(it+shiftDrv/dt)-f.evaluate(it)+
      g.evaluate(it+shiftDrv/dt)-g.evaluate(it))
    #x[it] = f.evaluate(it+shiftDrv/dt)-f.evaluate(it)
    #x[it] = g.evaluate(it+shiftDrv/dt)-g.evaluate(it)
    y[it] = g.evaluate(it)-f.evaluate(it)

  # Best-fit line
  z = linearRegression(x,y)

  # Average
  #avg = 0.0
  #for it in range(nt-ti):
  #  avg += y[it]/x[it]
  #print 'avg =',avg/(nt-ti)

  sp = SimplePlot()
  sp.setHLabel("f'(t)")
  sp.setVLabel("a f'(t)")
  pxy = sp.addPoints(x,y)
  pxy.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
  pxy.setLineStyle(PointsView.Line.NONE)
  pxz = sp.addPoints(x,z)
  pxz.setMarkStyle(PointsView.Mark.NONE)
  pxz.setLineWidth(2.0)
  pxz.setLineColor(Color.RED)
  sp.setSize(800,600)

def getSequences():
  f = zerofloat(nt)
  g = zerofloat(nt)
  for it in range(nt):
    f[it] = cos(freq*it*dt)
    g[it] = cos(freq*(it*dt+shiftMax))
  return Sequence(f),Sequence(g)

class Sequence():
  def __init__(self,y):
    self.y = y
    self.n = len(y) 
    self.si = SincInterpolator()
  def evaluate(self,x):
    return self.si.interpolate(self.n,1.0,0.0,self.y,x)
  def array(self):
    return self.y

def linearRegression(x,y):
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
  for i in range(n):
    xi = x[i]
    yi = y[i]
    sumx += xi
    sumy += yi
    sumxx += xi*xi
    sumxy += xi*yi
  beta = (sumxy-sumx*sumy/n)/(sumxx-sumx*sumx/n)
  alpha = (sumy-beta*sumx)/n
  z = zerofloat(n)
  for i in range(n):
    z[i] = alpha+beta*x[i]
  print 'slope =',beta
  print 'intercept =',alpha
  return z

def points(f,g):
  sp = SimplePlot()
  pf = sp.addPoints(f)
  pf.setLineColor(Color.RED)
  pg = sp.addPoints(g)
  pg.setLineColor(Color.BLUE)
  sp.setSize(1000,400)

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
  cb.setWidthMinimum(80)
  return p

def frame(panel,name=None):
  frame = PlotFrame(panel)
  #frame.setBackground(Color(204,204,204,255))
  #frame.setFontSizeForSlide(1.0,1.0)
  #frame.setSize(800,600)
  frame.setSize(1000,700)
  if name:
    frame.setTitle(name)
  frame.setVisible(True)
  if name and pngDir:
    frame.paintToPng(360,3.0,pngDir+name+'.png')

#############################################################################
import sys,time
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
