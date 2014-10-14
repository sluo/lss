from imports import *

pngDir = None
n1,n2 = 201,301

#############################################################################

def main(args):
  f = getImage(); plot(f)
  itp = Interpret()
  e = itp.computeErrors(f[10],f[20])
  plot(e,jet)

def getImage():
  slope = 0.5
  return FakeData.seismic2d2012B(n1,n2,slope)

jet = ColorMap.JET
gray = ColorMap.GRAY
def plot(x,cmap=gray,cmin=0,cmax=0,perc=100,cbar=None,name=None):
  panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pixel = panel.addPixels(x)
  pixel.setColorModel(cmap)
  #pixel.setInterpolation(PixelsView.Interpolation.LINEAR)
  pixel.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pixel.setClips(cmin,cmax)
  if perc<100:
    pixel.setPercentiles(100-perc,perc)
  colorbar = panel.addColorBar()
  if cbar:
    colorbar.setLabel(cbar)
  frame = PlotFrame(panel)
  #frame.setFontSizeForSlide(1.0,1.0)
  frame.setSize(800,600)
  if name:
    frame.setTitle(name)
  frame.setVisible(True)
  if name and pngDir:
    frame.paintToPng(360,3.0,pngDir+name+'.png')

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
