import sys
from math import *
from java.awt import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from sem import *

#############################################################################
# parameters

#dir = '/Users/sluo/Home/doc/Research/Semblance/temp/'
dir = '/Users/sluo/Home/Desktop/'
datadir = '/Users/sluo/Home/box/idh/trunk/bench/src/sim/cdps/'

st = Sampling(1001,0.004,0.0) # sampling for sythetic data
sx = Sampling(60,0.050,0.050)
sv = Sampling(101,0.020,1.5)
#sv = Sampling(1001,0.002,1.5) # for cCurve

sst = Sampling(1500,0.004,0.0) # sampling for Vikin Graben data
ssx = Sampling(60,0.050,-3.212)
ssv = Sampling(101,0.020,1.5)

vp = (2.00,3.00)
vm = (1.98,2.70)

fpeak = 25.0
tsigma = 8.0
snr = 1.0e6

contour = 0.20
contour1 = 0.40
contourwidth = 8.0

perc = -0.008
width = 600
height = 1000
hint = 1.0
vint = 1.0 

tperc = 0.2 # for teaser plots
#tperc1 = -0.052
tperc1 = 0.052
hlo = 2.0
hhi = 3.5
vlo = 2.25
vhi = 3.75

n0 = 1062.546630859375
n1 = 94.75059509277344
h0 = 63752.8125
h1 = 62174.5
d0 = 168.18316650390625
d1 = 99.5364990234375

#n0 = 1740.4544677734375
#n1 = 294.77069091796875
#h0 = 104427.265625
#h1 = 101308.09375
#d0 = 177.6163330078125
#d1 = 116.84996795654297

#############################################################################
# functions

def main(args):
	makeSynthetic()
  #makeReal()
  #test()

def test():
  v1 = Velan.makeLinearVelocity(2.0,3.0,st)
  v2 = Velan.makeLinearVelocity(2.0,3.0,st)
  p1 = Velan.makeRickerGather(fpeak,v1,st,sx)
  p2 = Velan.makeRickerGather(fpeak,v2,st,sx)
  sc1 = Velan.velocitySpectrum(st,sx,p1,sv,tsigma,2)
  sc2 = Velan.velocitySpectrum(st,sx,p2,sv,tsigma,2)
  sc3 = Velan.velocitySpectrum(st,sx,add(p1,p2),sv,tsigma,2)
  #cPlot(sc1,st,sv)
  #cPlot(sc2,st,sv)
  cPlot(add(sc1,sc2),st,sv)
  cPlot(sc3,st,sv)

def makeSynthetic():
  vps = Velan.makeLinearVelocity(vp[0],vp[1],st)
  vms = Velan.makeLinearVelocity(vm[0],vm[1],st)
  p = Velan.makeRickerGather(fpeak,vps,st,sx)
  p = add(p,Velan.makeRickerGather(fpeak,vms,st,sx)) # multiples

  p = Velan.addRandomNoise(snr,p) # noise
  cmpPlot(p,st,sx)

  s = Velan.velocitySpectrum(st,sx,p,sv,tsigma,3)
  bPlot(s,st,sv)

  sc = Velan.velocitySpectrum(st,sx,p,sv,tsigma,2)
  cPlot(sc,st,sv)

  #sw = Velan.velocitySpectrum(st,sx,p,sv,tsigma,1)
  sw = Velan.newSpectrum(st,sx,p,sv,tsigma)
  wPlot(sw,st,sv)

  #cCurve(sc,sw,3.2,st,sv)
  #bCurve(n0,n1,h0,h1,d0,d1)
  #teaserPlots(p,sc,sw,st,sv,sx)

def makeReal():
  #p = BinaryIO.readFloats(1500,60,datadir+'cdp_noheader.dat');
  #p = BinaryIO.readFloats(1500,60,datadir+'cdp1200.noheader.su');

  p = BinaryIO.readFloats(1500,60,datadir+'cdp600.dat');
  #p = BinaryIO.readFloats(1500,60,datadir+'cdp500.radon.noheader.su');

  cmpPlot(p,sst,ssx)
  s = Velan.velocitySpectrum(sst,ssx,p,ssv,tsigma,3)
  bPlot(s,sst,ssv)
  sc = Velan.velocitySpectrum(sst,ssx,p,ssv,tsigma,2)
  cPlot(sc,sst,ssv)
  #sw = Velan.velocitySpectrum(sst,ssx,p,ssv,tsigma,1)
  sw = Velan.newSpectrum(sst,ssx,p,ssv,tsigma)
  wPlot(sw,sst,ssv)

#############################################################################
# plots

def cmpPlot(array,samplingv,samplingh):
	pp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
	pp.setHLabel('Offset (km)')
	pp.setVLabel('Time (s)')
	pp.setHLimits(samplingh.getFirst(),samplingh.getLast())
	pp.setVLimits(samplingv.getFirst(),samplingv.getLast())
	pp.setHInterval(hint)
	pp.setVInterval(vint)
	pv = pp.addPixels(samplingv,samplingh,array)
	pv.setColorModel(ColorMap.GRAY)
	pv.setInterpolation(PixelsView.Interpolation.LINEAR)
	pv.setClips(-6.0,6.0)
	pf = PlotFrame(pp)
	pf.setSize((int)(width-width*perc),height)
	pf.setFontSizeForPrint(8,115-115*perc)
	pf.setVisible(True)
	pf.paintToPng(720,1.67,dir+'1.png')

def cPlot(array,samplingv,samplingh):
	pp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
	pp.setHLabel('Velocity (km/s)')
	pp.setVLabel('Time (s)')
	pp.setHLimits(samplingh.getFirst(),samplingh.getLast())
	pp.setVLimits(samplingv.getFirst(),samplingv.getLast())
	pp.setHInterval(hint)
	pp.setVInterval(vint)
	pv = pp.addPixels(samplingv,samplingh,array)
	pv.setColorModel(ColorMap.JET)
	pv.setInterpolation(PixelsView.Interpolation.LINEAR)
	pv.setClips(0.0,1.0)
	cv = pp.addContours(samplingv,samplingh,array)
	cv.setLineColor(Color.MAGENTA)
	cv.setLineWidth(contourwidth)
	cv.setContours(Sampling(1,0.1,contour))
	pf = PlotFrame(pp)
	pf.setSize((int)(width-width*perc),height)
	pf.setFontSizeForPrint(8,115-115*perc)
	pf.setVisible(True)
	pf.paintToPng(720,1.67,dir+'3.png')

def wPlot(array,samplingv,samplingh):
  pp = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,PlotPanel.AxesPlacement.TOP)
  pp.setHLabel('Velocity (km/s)')
  pp.setHLimits(samplingh.getFirst(),samplingh.getLast())
  pp.setVLimits(samplingv.getFirst(),samplingv.getLast())
  pp.setHInterval(hint)
  cb = pp.addColorBar('Semblance')
  cb.setInterval(1.0)
  cb.setWidthMinimum(75)
  pv = pp.addPixels(samplingv,samplingh,array)
  pv.setColorModel(ColorMap.JET)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setClips(0.0,1.0)
  #pv.setClips(0.0,2.0)
  cv = pp.addContours(samplingv,samplingh,array)
  cv.setLineColor(Color.MAGENTA)
  cv.setLineWidth(contourwidth)
  cv.setContours(Sampling(1,0.1,contour))
  pf = PlotFrame(pp)
  pf.setSize((int)(width+width*perc),height)
  pf.setFontSizeForPrint(8,115+115*perc)
  pf.setVisible(True)
  pf.paintToPng(720,1.67,dir+'4.png')

def bPlot(array,samplingv,samplingh):
  pp = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,PlotPanel.AxesPlacement.TOP)
  pp.setHLabel('Velocity (km/s)')
  pp.setHLimits(samplingh.getFirst(),samplingh.getLast())
  pp.setVLimits(samplingv.getFirst(),samplingv.getLast())
  pp.setHInterval(hint)
  cb = pp.addColorBar('b value')
  cb.setInterval(1.0)
  cb.setWidthMinimum(75)
  pv = pp.addPixels(samplingv,samplingh,array)
  pv.setColorModel(ColorMap.JET)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  #pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #pv.setClips(0.80,1.0) #TODO: check b plot for real data
  pv.setClips(0.0,1.0)
  #pv.setClips(-10.0,10.0)
  #cv = pp.addContours(samplingv,samplingh,array)
  #cv.setLineColor(Color.BLACK)
  #cv.setLineWidth(2.0)
  #cv.setContours(Sampling(1,0.1,contour))
  pf = PlotFrame(pp)
  pf.setSize((int)(width+width*perc),height)
  pf.setFontSizeForPrint(8,115+115*perc)
  pf.setVisible(True)
  pf.paintToPng(720,1.67,dir+'2.png')

def cCurve(semc,semw,t,st,sv):
  ccc = Velan.coherenceCurve(t,semc,sv,st)
  ccw = Velan.coherenceCurve(t,semw,sv,st)
  pp = PlotPanel()
  pp.setHLabel('Velocity (km/s)')
  pp.setVLabel('Semblance')
  pp.setHLimits(2.0,3.5)
  pp.setVLimits(0.0,1.0)
  pp.setHInterval(0.5)
  pp.setVInterval(0.2)
  pvc = pp.addPoints(sv,ccc)
  pvw = pp.addPoints(sv,ccw)
  pvw.setLineStyle(PointsView.Line.DASH)
  pvc.setLineWidth(3.0)
  pvw.setLineWidth(3.0)
  pf = PlotFrame(pp)
  pf.setSize(height,4*height/5)
  pf.setFontSizeForPrint(8,240)
  pf.setVisible(True)
  pf.paintToPng(720,3.33,dir+'5.png')

def bCurve(n0,n1,h0,h1,d0,d1):
  bb = Velan.bCurve(n0,n1,h0,h1,d0,d1)
  pp = PlotPanel()
  pp.setHLabel('b value')
  pp.setVLabel('Semblance')
  pp.setHLimits(0.0,4.0)
  pp.setVLimits(-3.0,3.0)
  pp.setHInterval(1.0)
  pp.setVInterval(2.0)
  gv = pp.addGrid()
  gv.setHorizontal(GridView.Horizontal.ZERO);
  gv.setVertical(GridView.Vertical.ZERO);
  
  sb = Sampling(650,0.01,-1.5)
  
  pv = pp.addPoints(sb,bb)
  pv.setLineStyle(PointsView.Line.DASH)
  pv.setLineWidth(3.0);
  #pv.setLineColor(Color.red)
  pf = PlotFrame(pp)
  pf.setSize(1500,1200)
  pf.setFontSizeForPrint(8,240)
  pf.setVisible(True)
  pf.paintToPng(720,3.33,dir+'bcurve.png')

def teaserPlots(cmp,semc,semw,st,sv,sx):
  pp1 = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pp1.setHLabel('Offset (km)')
  pp1.setVLabel('Time (s)')
  pp1.setHLimits(sx.getFirst(),sx.getLast())
  pp1.setVLimits(vlo,vhi)
  pp1.setHInterval(1.0)
  pp1.setVInterval(0.5)
  pv1 = pp1.addPixels(st,sx,cmp)
  pv1.setColorModel(ColorMap.GRAY)
  pv1.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv1.setClips(-6.0,6.0)
  pf1 = PlotFrame(pp1)
  pf1.setSize((int)(800-800*tperc),800)
  pf1.setFontSizeForPrint(8,(int)(150-150*tperc))
  pf1.setVisible(True)
  pf1.paintToPng(300,2.0,dir+'teaser1.png')
  
  width1 = (int)(800+800*tperc)
  fsize = (int)(150+150*tperc)

  pp2 = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,PlotPanel.AxesPlacement.TOP)
  pp2.setHLabel('Velocity (km/s)')
  pp2.setHLimits(hlo,hhi)
  pp2.setVLimits(vlo,vhi)
  pp2.setHInterval(0.5)
  pv2 = pp2.addPixels(st,sv,semc)
  pv2.setColorModel(ColorMap.JET)
  pv2.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv2.setClips(0.0,1.0)
  cv2 = pp2.addContours(st,sv,semc)
  cv2.setLineColor(Color.BLACK)
  cv2.setLineWidth(2.0)
  cv2.setContours(Sampling(1,0.1,contour))
  cv21 = pp2.addContours(st,sv,semc)
  cv21.setLineColor(Color.BLACK)
  cv21.setLineWidth(2.0)
  cv21.setContours(Sampling(1,0.1,contour1))
  cb2 = pp2.addColorBar()
  cb2.setInterval(10.0)
  cb2.setWidthMinimum(10)
  pf2 = PlotFrame(pp2)
  pf2.setSize((int)(width1-width1*tperc1),800)
  pf2.setFontSizeForPrint(8,(int)(fsize-fsize*tperc1))
  pf2.setVisible(True)
  pf2.paintToPng(300,3.0,dir+'teaser2.png')

  pp3 = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,PlotPanel.AxesPlacement.TOP)
  pp3.setHLabel('Velocity (km/s)')
  pp3.setHLimits(hlo,hhi)
  pp3.setVLimits(vlo,vhi)
  pp3.setHInterval(0.5)
  pv3 = pp3.addPixels(st,sv,semw)
  pv3.setColorModel(ColorMap.JET)
  pv3.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv3.setClips(0.0,1.0)
  cv3 = pp3.addContours(st,sv,semw)
  cv3.setLineColor(Color.BLACK)
  cv3.setLineWidth(2.0)
  cv3.setContours(Sampling(1,0.1,contour))
  cv31 = pp3.addContours(st,sv,semw)
  cv31.setLineColor(Color.BLACK)
  cv31.setLineWidth(2.0)
  cv31.setContours(Sampling(1,0.1,contour1))
  cb3 = pp3.addColorBar('Semblance')
  cb3.setInterval(1.0)
  cb3.setWidthMinimum(110)
  pf3 = PlotFrame(pp3)
  pf3.setSize((int)(width1+width1*tperc1),800)
  pf3.setFontSizeForPrint(8,(int)(fsize+fsize*tperc1))
  pf3.setVisible(True)
  pf3.paintToPng(300,3.0,dir+'teaser3.png')


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
