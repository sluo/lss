#############################################################################
# Finite-difference coefficients

from imports import *
from edu.mines.jtk.lapack import *

#############################################################################

def main(args):
  derivativeOrder = 1
  halfStencilOrder = 1 # half the stencil order
  computeCoefficients(derivativeOrder,halfStencilOrder)

def computeCoefficients(derivativeOrder,halfStencilOrder):
  stencilOrder = 2*halfStencilOrder
  n = stencilOrder+1
  mata = DMatrix(n,n)
  r = halfStencilOrder
  for j in range(n):
    for i in range(n):
      mata.set(i,j,pow(r,i)/factorial(i))
    r -= 1
  mate = DMatrix(n,1)
  if derivativeOrder==1:
    mate.set(1,0,1.0)
  elif derivativeOrder==2:
    mate.set(2,0,1.0)
  matx = mata.solve(mate)
  x = matx.getArray()
  print x

def factorial(k):
  return 1 if k==0 else k*factorial(k-1)


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
