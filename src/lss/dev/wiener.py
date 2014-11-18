##############################################################################
# Wiener filtering and deconvolution

from imports import *

##############################################################################

def main(args):
  #goAdjointTest()
  goDecon()

def goAdjointTest():
  n1 = 11
  random = Random(12345)
  aa = randfloat(random,n1)
  fa,ga = randfloat(random,n1),zerofloat(n1)
  fb,gb = randfloat(random,n1),zerofloat(n1)
  Conv.conv(n1,0,aa,
            n1,-n1/2,fa,
            n1,0,ga)
  Conv.xcor(n1,0,aa,
            n1,0,fb,
            n1,-n1/2,gb)
  #applyLhs(aa,fa,ga)
  #applyLhs(aa,fb,gb)
  check(fa,fb,ga,gb)

def goDecon():
  nt = 1001
  random = Random(12345)
  do = wavelet(nt,nt/2)
  dp = wavelet(nt,nt/4)
  #dc = deconFft(do,dp)
  dc = deconWiener(do,dp)
  SimplePlot.asPoints(dc)

  # Weights
  wt = zerofloat(nt)
  for it in range(nt):
    wt[it] = float(it-nt/2)
  #wd = mul(wt,dc)
  #dnorm = sum(mul(dc,dc))
  #wnorm = sum(mul(wd,wd))
  #for it in range(nt):
  #  wt[it] = (wt[it]-wnorm/dnorm)/dnorm
  SimplePlot.asPoints(wt)
  mul(wt,dc,dc)

  t = zerofloat(nt)
  Conv.conv(nt,0,do,
            nt,-nt/2,dc,
            nt,0,t)
  SimplePlot.asPoints(t)

  # Adjoint source
  x = zerofloat(nt)
  vx = VecArrayFloat1(x)
  vt = VecArrayFloat1(t)
  a1 = A1(do,epsilon=1.0e6)
  cg = CgSolver(0.001,100)
  cg.solve(a1,vt,vx)
  SimplePlot.asPoints(vx.getArray())

def deconWiener(do,dp):
  n = len(do)
  x = zerofloat(n); vx = VecArrayFloat1(x)
  b = zerofloat(n); vb = VecArrayFloat1(b)
  makeRhs(do,dp,b)
  a1 = A1(do)
  cg = CgSolver(0.001,100)
  cg.solve(a1,vb,vx)
  return vx.getArray()
  #return b

class A1(CgSolver.A):
  def __init__(self,do,epsilon=None):
    self.do = do
    # Stabilization
    if epsilon:
      self.epsilon = epsilon
    else:
      self.epsilon = 0.0
      #self.epsilon = 1.0E0
      #self.epsilon = 1.0E6
  def apply(self,vx,vy):
    x = vx.getArray()
    y = vy.getArray()
    applyLhs(self.do,x,y)
    if self.epsilon>0.0:
      vy.add(1.0,vx,self.epsilon)

def applyLhs(a,x,y):
  n = len(y)
  h = -n/2
  t = zerofloat(n)
  Conv.conv(n,0,a,
            n,h,x,
            n,0,t)
  Conv.xcor(n,0,a,
            n,0,t,
            n,h,y)

def makeRhs(a,x,y):
  n = len(y)
  h = -n/2
  Conv.xcor(n,0,a,
            n,0,x,
            n,h,y)

def deconFft(do,dp):
  epsilon = 1.0E-1
  #epsilon = 1.0E6
  fft = Fft(do)
  fp = fft.applyForward(dp)
  fo = fft.applyForward(do)
  co = cconj(fo)
  den = cadd(cmul(co,fo),Cfloat(epsilon))
  num = cmul(co,fp)
  fc = cdiv(num,den)
  dc = fft.applyInverse(fc)
  #BandPassFilter(0.00,0.03,0.01,0.01).apply(dc,dc)
  return dc

def wavelet(nt,kt=None):
  w = zerofloat(nt)
  if kt:
    w[kt] = -1.0
  else:
    w[nt/2] = -1.0
  RecursiveGaussianFilter(16.0).apply2(w,w)
  mul(1.0/max(w),w,w)
  return w

def check(fa,fb,ga,gb):
  f = sum(mul(fa,gb))
  g = sum(mul(ga,fb))
  print "adjoint test:",compare(f,g)
  print f
  print g

def compare(a,b):
  digits = 20
  equal = False
  while not equal and digits>0:
    almost = Almost(digits)
    equal = almost.equal(a,b)
    if not equal:
      digits -= 1
  return digits

##############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

