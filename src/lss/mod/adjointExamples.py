from java.util import Random
from edu.mines.jtk.util import Almost
from edu.mines.jtk.util.ArrayMath import *

##############################################################################

nx,nt = 10,11
random = Random(12345)

##############################################################################

def main(args):
  fa,ga = rfloat(nx,nt),zfloat(nx,nt)
  fb,gb = rfloat(nx,nt),zfloat(nx,nt)
  applyForward(fa,ga)
  applyAdjoint(fb,gb)
  #applyForward1A(fa,ga)
  #applyAdjoint1A(fb,gb)
  #applyForward1B(fa,ga)
  #applyAdjoint1B(fb,gb)
  #applyForward2A(fa,ga)
  #applyAdjoint2A(fb,gb)
  #applyForward2B(fa,ga)
  #applyAdjoint2B(fb,gb)
  #applyForward3A(fa,ga)
  #applyAdjoint3A(fb,gb)
  #applyForward3B(fa,ga)
  #applyAdjoint3B(fb,gb)
  check(fa,fb,ga,gb)

##############################################################################
# Source injection, plus different time indexing for forward and adjoint

def applyForward(f,g):
  up = zerofloat(nt)
  ui = zerofloat(nt)
  um = zerofloat(nt)
  for ix in range(nx):
    um[ix] = f[0][ix]
    ui[ix] = f[1][ix]
    g[0][ix] = f[0][ix]
    g[1][ix] = f[1][ix]
  for it in range(2,nt):
    for ix in range(nx):
      up[ix] += f[it][ix]
    for ix in range(nx):
      up[ix] -= um[ix]
      up[ix] += ui[ix]*2.0
    for ix in range(nx):
      g[it][ix] = up[ix]
    copy(ui,um)
    copy(up,ui)
    zero(up)

def applyAdjoint(f,g):
  up = zerofloat(nt)
  ui = zerofloat(nt)
  um = zerofloat(nt)
  for ix in range(nx):
    um[ix] = f[nt-1][ix]
    ui[ix] = f[nt-2][ix]
    g[nt-1][ix] = f[nt-1][ix]
    g[nt-2][ix] = f[nt-2][ix]
  for it in range(nt-3,-1,-1):
    for ix in range(nx):
      up[ix] += f[it][ix]
    for ix in range(nx):
      up[ix] -= um[ix]
      ui[ix] += um[ix]*2.0
    for ix in range(nx):
      g[it  ][ix] = up[ix]
      g[it+1][ix] = ui[ix]
    copy(ui,um)
    copy(up,ui)
    zero(up)

##############################################################################
# Same time indexing for forward and adjoint
  
def applyForward1A(f,g):
  copy(f,g)
  for it in range(1,nt-1):
    for ix in range(nx):
      g[it+1][ix] -= g[it-1][ix]
      g[it+1][ix] += g[it  ][ix]*2.0

def applyAdjoint1A(f,g):
  copy(f,g)
  for it in range(nt-2,0,-1):
    for ix in range(nx):
      g[it-1][ix] -= g[it+1][ix]
      g[it  ][ix] += g[it+1][ix]*2.0

def applyForward1B(f,g):
  copy(f,g)
  for it in range(1,nt-1):
    up = g[it+1]
    ui = g[it  ]
    um = g[it-1]
    for ix in range(nx):
      up[ix] -= um[ix]
      up[ix] += ui[ix]*2.0

def applyAdjoint1B(f,g):
  copy(f,g)
  for it in range(nt-2,0,-1):
    up = g[it+1]
    ui = g[it  ]
    um = g[it-1]
    for ix in range(nx):
      um[ix] -= up[ix]
      ui[ix] += up[ix]*2.0

##############################################################################
# Same time indexing for forward and adjoint
  
def applyForward2A(f,g):
  copy(f,g)
  for it in range(2,nt):
    for ix in range(nx):
      g[it][ix] -= g[it-2][ix]
      g[it][ix] += g[it-1][ix]*2.0

def applyAdjoint2A(f,g):
  copy(f,g)
  for it in range(nt-1,1,-1):
    for ix in range(nx):
      g[it-2][ix] -= g[it][ix]
      g[it-1][ix] += g[it][ix]*2.0

def applyForward2B(f,g):
  copy(f,g)
  for it in range(2,nt):
    up = g[it  ]
    ui = g[it-1]
    um = g[it-2]
    for ix in range(nx):
      up[ix] -= um[ix]
      up[ix] += ui[ix]*2.0

def applyAdjoint2B(f,g):
  copy(f,g)
  for it in range(nt-1,1,-1):
    up = g[it  ]
    ui = g[it-1]
    um = g[it-2]
    for ix in range(nx):
      um[ix] -= up[ix]
      ui[ix] += up[ix]*2.0

##############################################################################
# Different time indexing for forward and adjoint
  
def applyForward3A(f,g):
  copy(f,g)
  for it in range(2,nt):
    for ix in range(nx):
      g[it][ix] -= g[it-2][ix]
      g[it][ix] += g[it-1][ix]*2.0

def applyAdjoint3A(f,g):
  copy(f,g)
  for it in range(nt-3,-1,-1):
    print it
    for ix in range(nx):
      g[it  ][ix] -= g[it+2][ix]
      g[it+1][ix] += g[it+2][ix]*2.0

def applyForward3B(f,g):
  copy(f,g)
  for it in range(2,nt):
    up = g[it]
    ui = g[it-1]
    um = g[it-2]
    for ix in range(nx):
      up[ix] -= um[ix]
      up[ix] += ui[ix]*2.0

def applyAdjoint3B(f,g):
  copy(f,g)
  for it in range(nt-3,-1,-1):
    up = g[it]
    ui = g[it+1]
    um = g[it+2]
    for ix in range(nx):
      up[ix] -= um[ix]
      ui[ix] += um[ix]*2.0

##############################################################################

def check(fa,fb,ga,gb):
  f = dot(fa,gb)
  g = dot(ga,fb)
  #print "adjoint test:",WaveOperator.compareDigits(f,g)
  print "adjoint test:",compareDigits(f,g)
  print f
  print g

def compareDigits(xa,xb):
  digits = 20
  equal = False
  while not equal and digits>0:
    almost = Almost(digits)
    equal = almost.equal(xa,xb)
    if not equal:
      digits -= 1
  return digits

def dot(u,v):
  return sum(mul(u,v))

def zfloat(n1,n2):
  return zerofloat(n1,n2)

def rfloat(n1,n2):
  r = randfloat(random,n1,n2)
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
