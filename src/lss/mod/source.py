#############################################################################
# Source tests

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util.ArrayMath import *
from lss.mod import *

#############################################################################

def main(args):
  si = SincInterpolator()
  table = si.getTable()
  nsinc = len(table)
  lsinc = len(table[0])
  print 'nsinc =',nsinc
  print 'lsinc =',lsinc
  pixels(table,cmap=jet)
  #for isinc in range(0,nsinc,100):
  #  SimplePlot.asPoints(table[isinc])
  SimplePlot.asPoints(table[nsinc/2])





#############################################################################
# plotting

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def pixels(x,cmap=gray,perc=100.0,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  sp.setSize(600,600)
  if title:
    sp.addTitle(title)
  pv = sp.addPixels(x)
  pv.setColorModel(cmap)
  if perc<100.0:
    pv.setPercentiles(100.0-perc,perc)

#############################################################################
# Do everything on Swing thread.

import sys
from java.lang import Runnable
from javax.swing import SwingUtilities
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
