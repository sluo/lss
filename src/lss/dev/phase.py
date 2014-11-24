##############################################################################
# Hilbert transform and phase shifts

from imports import *

##############################################################################

nt = 1001
nfft = FftComplex.nfftFast(nt)
fft = FftComplex(nfft)
print 'nt =',nt
print 'nfft =',nfft

def main(args):
  f = wavelet()
  g = rotate(f)
  #h = hilbert(f)
  #g = rotate(rotate(f))
  #h = hilbert(hilbert(f))
  pointsReal(f,title='input')
  #points(mul(-1.0,f),title='negative')
  points(g,title='rotate')
  pointsReal(g,title='rotate')
  #points(h,title='hilbert')

def rotate(f):
  #sf = fft.getFrequencySampling1()
  #nf = sf.count
  g = zerofloat(2*nfft)
  fft.complexToComplex(-1,f,g)
  r = zerofloat(2*nfft)
  for i in range(0,nfft/2):
    #print 'f =',sf.getValue(i)
    r[2*i  ] = 0.0
    r[2*i+1] = 1.0
  for i in range(nfft/2,nfft):
    r[2*i  ] = 0.0
    r[2*i+1] = -1.0
  cmul(r,g,g)
  h = zerofloat(2*nfft)
  fft.complexToComplex(1,g,h) # inverse FFT
  fft.scale(nt,h)
  return h

def hilbert(f):
  n = len(f)
  g = zerofloat(n)
  HilbertTransformFilter().apply(n,f,g)
  return g

def wavelet():
  sigma = 50
  w = zerofloat(nt)
  w[nt/2] = 1.0
  RecursiveGaussianFilter(sigma).apply2(w,w)
  mul(-1.0/max(abs(w)),w,w)
  c = zerofloat(2*nfft) # complex
  for it in range(nt):
    c[2*it] = w[it]
  return c

def points(f,title=None):
  sp = SimplePlot.asPoints(f)
  if title:
    sp.addTitle(title)

def pointsReal(f,title=None):
  g = zerofloat(nt)
  for it in range(nt):
    g[it] = f[2*it]
  points(g,title)

##############################################################################
import sys
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

