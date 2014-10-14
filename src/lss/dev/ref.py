"""
Testing different methods for slope estimation
"""
from imports import *

#############################################################################

def main(args):
  a = 0.9
  n = 101
  random = Random(314)
  xa = randfloat(random,n)
  #xa = zerofloat(n); xa[n/2] = 1.0
  xb = copy(xa)
  #ya = zerofloat(n)
  #yb = zerofloat(n)
  ya = xa
  yb = xb
  xsmoothDave(a,xa,ya)
  #smooth(a,xb,yb)
  #smoothSimon(a,xb,yb)
  smoothDave(a,xb,yb)
  SimplePlot.asPoints(ya)
  SimplePlot.asPoints(yb)

  #for i in range(n):
  #  print ya[i]-yb[i]
  print sum(xa)
  print sum(xb)
  print sum(ya)
  print sum(yb)

def xsmoothDave(a,x,y):
  n = len(x)
  b = 1.0-a
  sx = 1.0
  sy = a
  yi = 0.0
  y[0] = yi = sy*yi+sx*x[0]
  for i in range(1,n-1):
    y[i] = yi = a*yi+b*x[i]
  sx /= 1.0+a
  sy /= 1.0+a
  y[n-1] = yi = sy*yi+sx*x[n-1]
  for i in range(n-2,-1,-1):
    y[i] = yi = a*yi+b*y[i]

def smoothDave(a,x,y):
  n = len(x)
  b = 1.0-a
  sx = 1.0
  sy = a
  yi = 0.0
  y[0] = sx*x[0]
  for i in range(1,n-1):
    y[i] = a*y[i-1]+b*x[i]
  sx /= 1.0+a
  sy /= 1.0+a
  y[n-1] = sy*y[n-2]+sx*x[n-1]
  for i in range(n-2,-1,-1):
    y[i] = a*y[i+1]+b*y[i]

def smoothSimon(a,x,y):
  """Zero-value boundary conditions"""
  n = len(x)
  s = (1.0-a)/(1.0+a)
  r = (1.0-a)*(1.0-a)
  y[0] = x[0]
  for i in range(1,n):
    y[i] = a*y[i-1]+x[i]
  y[n-1] *= s
  for i in range(n-2,-1,-1):
    y[i] = a*y[i+1]+r*y[i]

def smooth(a,x,y):
  n = len(x)
  t = zerofloat(n)

  y[0] = (1.0-a)*x[0]
  t[0] = y[0]
  for i in range(1,n-1):
    y[i] = a*y[i-1]+(1.0-a)*x[i]
    t[i] = y[i]

  y[n-1] = a/(1.0+a)*y[n-2]+(1.0-a)/(1.0+a)*x[n-1]

#  for i in range(n-2,-1,-1):
#    #print i,y[i]-(a*y[i-1]+(1.0-a)*x[i])
#    y[i] = a*y[i+1]+(1.0-a)*y[i]
#    #if i>0:
#    #  yi = a*t[i-1]+(1.0-a)*x[i]
#    #  y[i] = a*y[i+1]+(1.0-a)*yi
#    #  #y[i] = a*y[i+1]+(1.0-a)*(a*t[i-1]+(1.0-a)*x[i])
#    #else:
#    #  yi = (1.0-a)*x[i]
#    #  y[i] = a*y[i+1]+(1.0-a)*yi
#    #  #y[i] = a*y[i+1]+(1.0-a)*((1.0-a)*x[i])

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
    print '%02d:%02d:%02ds'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
