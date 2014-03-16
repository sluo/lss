##############################################################################
# Adjoint tests

from java.util import *
from edu.mines.jtk.util.ArrayMath import *
from lss.mod import *

##############################################################################

nabsorb = 20
nx,nz,nt = 10,11,402
dx,dz,dt = 0.01,0.01,0.001
nxp,nzp = nx+2*nabsorb,nz+2*nabsorb
random = Random(12345)

##############################################################################

def main(args):
  #testWaveOperator()
  #testBornOperator()
  testBornOperatorWithShifts()

def testWaveOperator():
  """Adjoint test for forward/adjoint wave propagation"""
  wave = WaveOperator(getSlowness(),dx,dt,nabsorb)
  fa,fb = rfloat(nxp,nzp,nt),zfloat(nxp,nzp,nt)
  ga,gb = rfloat(nxp,nzp,nt),zfloat(nxp,nzp,nt)
  wave.applyForward(Source.WavefieldSource(fa),gb)
  wave.applyAdjoint(Source.WavefieldSource(ga),fb)
  checkDotProducts(fa,fb,ga,gb)

def testBornOperator(useTimeShifts=False):
  """Adjoint test for Born modeling and migration"""
  xr = rampint(0,1,min(nx,nz)) # receiver x-coordinates
  zr = rampint(0,1,min(nx,nz)) # receiver z-coordinates
  nr = len(xr)
  born = BornOperator(getSlowness(),dx,dt,nabsorb)
  fa,fb = rfloat(nx,nz),zfloat(nx,nz) # reflectivity
  ga,gb = rfloat(nt,nr),zfloat(nt,nr) # data...
  ra,rb = Receiver(xr,zr,ga),Receiver(xr,zr,gb) # ...for receivers
  bwf = rfloat(nxp,nzp,nt) # random background wavefield
  awf = zfloat(nxp,nzp,nt) # array for adjoint wavefield
  if useTimeShifts:
    ts = rfloat(nt,nr) # random time shifts
    born.applyForward(bwf,fa,ts,rb)
    born.applyAdjoint(bwf,awf,ra,ts,fb)
  else:
    born.applyForward(bwf,fa,rb)
    born.applyAdjoint(bwf,awf,ra,fb)
  checkDotProducts(fa,fb,ga,gb)

def testBornOperatorWithShifts():
  testBornOperator(True)

##############################################################################

def getSlowness():
  r = randfloat(random,nx,nz)
  mul(0.1,r,r)
  s = fillfloat(0.2,nx,nz)
  return add(s,r)

def checkDotProducts(fa,fb,ga,gb):
  f = dot(fa,fb)
  g = dot(ga,gb)
  print "adjoint test:",WaveOperator.compareDigits(f,g)
  print f
  print g

def dot(u,v):
  return sum(mul(u,v))

def zfloat(n1,n2,n3=None):
  return zerofloat(n1,n2) if n3 is None else zerofloat(n1,n2,n3)

def rfloat(n1,n2,n3=None):
  r = randfloat(random,n1,n2) if n3 is None else randfloat(random,n1,n2,n3)
  sub(r,0.5,r)
  mul(2.0,r,r)
  return r

##############################################################################
# Do everything on Swing thread.

import sys
from java.lang import Runnable
from javax.swing import SwingUtilities
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
