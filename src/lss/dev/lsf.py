##############################################################################
# Adjoint tests for LocalSmoothingFilter

from imports import *

##############################################################################

def main(args):
  n1,n2 = 101,102
  random = Random(12345)
  fa,ga = randfloat(random,n1,n2),zerofloat(n1,n2)
  fb,gb = randfloat(random,n1,n2),zerofloat(n1,n2)
  """
  lsf = LocalSmoothingFilter()
  lsf.apply(fa,ga)
  lsf.apply(fb,gb)
  """
  lsfa = LocalSmoothingFilter(0.0,1)
  lsfb = LocalSmoothingFilter(0.0,2)
  lsfa.apply(fa,ga)
  lsfb.apply(fb,gb)
  check(fa,fb,ga,gb)

def check(fa,fb,ga,gb):
  f = dot(fa,gb)
  g = dot(ga,fb)
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

def dot(a,b):
  return sum(mul(a,b))

  

##############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

