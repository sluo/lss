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

pngDir = '/Users/sluo/Home/doc/Research/Semblance/temp/'
dataDir = '/Users/sluo/Home/box/idh/trunk/bench/src/sim/cdps/'

st = Sampling(1001,0.004,0.0) # sampling for sythetic data
sx = Sampling(60,0.050,0.050)
#sv = Sampling(101,0.020,1.5)
sv = Sampling(401,0.005,1.5)

sst = Sampling(1500,0.004,0.0) # sampling for Vikin Graben data
ssx = Sampling(60,0.050,-3.212)
ssv = Sampling(101,0.020,1.5)

vp = (2.00,3.00)
vm = (1.98,2.70)

fpeak = 25.0
tsigma = 8.0
snr = 1.0e6

contours = [0.20]
contour = 0.10
contour1 = 0.40

contourwidth = 7.0

perc = 0.03 
#width = 600
width = 1000
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
	spectrum()

def spectrum():
	makeSynthetic()
  #makeReal()

  #test(100)
  #n = 10
  #r = zerofloat(n)
  #for i in range(n):
  #  real = 100
  #  for ir in range(real):
  #    r[i] += test(i+1,ir)/float(real)
  #  print('r=%f'%r[i]) 

  #pp = PlotPanel()
  #pp.setHLabel('Number of traces')
  #pp.setVLabel('Semblance expectation')
  ##pp.setHLimits(2.0,3.5)
  #pp.setVLimits(0.9,1.5)
  ##pp.setHInterval(0.5)
  ##pp.setVInterval(0.2)
  #point = pp.addPoints(Sampling(n),r)
  #point.setLineWidth(3.0)
  ##pvw.setLineStyle(PointsView.Line.DASH)
  ##pvw.setLineWidth(3.0)
  #pf = PlotFrame(pp)
  #pf.setSize(height+200,height)
  #pf.setFontSizeForSlide(0.8,1.0)
  #pf.setVisible(True)
  #pf.paintToPng(600,5,pngDir+'expec.png')
    
  #n = 20
  #e1 = zerofloat(n)
  #e2 = zerofloat(n)
  #for i in range(n):
  #  e1[i] = test(i+1,0,1)
  #  e2[i] = test(i+1,0,2)
  #pp = PlotPanel()
  #pp.setHLabel('Number of traces')
  #pp.setVLabel('Semblance expectation')
  ##pp.setHLimits(2.0,3.5)
  #pp.setVLimits(0.0,1.1)
  ##pp.setHInterval(0.5)
  ##pp.setVInterval(0.2)
  #point = pp.addPoints(Sampling(n),e1)
  #point.setLineWidth(3.0)
  #point2 = pp.addPoints(Sampling(n),e2)
  #point2.setLineStyle(PointsView.Line.DASH)
  #point2.setLineWidth(3.0)
  #pf = PlotFrame(pp)
  #pf.setSize(height+200,height)
  #pf.setFontSizeForSlide(0.8,1.0)
  #pf.setVisible(True)
  #pf.paintToPng(600,5,pngDir+'expec1.png')
    

# random noise analysis (Fomel)
def test(n,seed=None,mode=None):
  nt = 10000
  q = zerofloat(nt,n)
  rann = RandomFloat(seed)

  for it in range(nt):
    #ran = RandomFloat(it)
    ran = RandomFloat(int(nt*rann.uniform()))
    for i in range(n):
      q[i][it] = ran.normal()

  cs = Velan.semblance(0.0,q)
  e1 = sum(cs)/len(cs)
  #print('e=%f'%e1)

  #vnmo = 100.0
  #st = Sampling(nt,1.004,0.004)
  #sx = Sampling(n,1.000,0.050)
  #ws = Velan.semblance(st,sx,vnmo,0.0,q,"",True,False)
  #e2 = sum(ws)/len(ws)
  ##print('e=%f'%e)
  
  e2 = 0.0
  for it in range(nt):
    qq = zerofloat(1,n)
    for i in range(n):
      qq[i][0] = q[i][it]

    vnmo = 2.0
    st = Sampling(1,1.0,1.0)
    sx = Sampling(n,0.050,0.050)
    ws = Velan.semblance(st,sx,vnmo,0.0,qq,"",True,False)
    e2 += ws[0]/float(nt)
  #print('e=%f'%e2)

  if mode==1:
    return e1
  elif mode==2:
    return e2
  else:
    return e1/e2
  

def makeSynthetic():
  vps = Velan.makeLinearVelocity(vp[0],vp[1],st)
  vms = Velan.makeLinearVelocity(vm[0],vm[1],st)
  p = Velan.makeRickerGather(fpeak,vps,st,sx)
  p = add(p,Velan.makeRickerGather(fpeak,vms,st,sx)) # multiples
  p = Velan.addRandomNoise(snr,p) # noise

  cmpPlot(p,st,sx)
  #cmpPlot(Velan.nmo(3.0,st,sx,p),st,sx)
  b = Velan.velocitySpectrum(st,sx,p,sv,tsigma,3)
  bPlot(b,st,sv,'2')

  cs = Velan.velocitySpectrum(st,sx,p,sv,tsigma,2)
  plot(cs,st,sv,'3')

  ws = Velan.newSpectrum(st,sx,p,sv,tsigma)
  plot(ws,st,sv,'4')

  #cCurve(cs,ws,3.2,st,sv)
  #bCurve(n0,n1,h0,h1,d0,d1)
  #teaserPlots(p,cs,ws,st,sv,sx)

def makeReal():
  p = BinaryIO.readFloats(1500,60,dataDir+'cdp_noheader.dat')
  #p = BinaryIO.readFloats(1500,60,dataDir+'cdp1200.noheader.su')
  #p = BinaryIO.readFloats(1500,60,dataDir+'cdp1200_a.su')
  #p = BinaryIO.readFloats(1500,60,dataDir+'cdp900.dat')
  cmpPlot(p,sst,ssx)
  b = Velan.velocitySpectrum(sst,ssx,p,ssv,tsigma,3)
  bPlot(b,sst,ssv,'2')
  cs = Velan.velocitySpectrum(sst,ssx,p,ssv,tsigma,2)
  plot(cs,sst,ssv,'3')

  #ws = Velan.velocitySpectrum(sst,ssx,p,ssv,tsigma,1)
  ws = Velan.newSpectrum(sst,ssx,p,ssv,tsigma)
  plot(ws,sst,ssv,'4')

#############################################################################
# plots

def cmpPlot(f,s1,s2):
  pp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pp.setHLabel('Offset (km)')
  pp.setVLabel('Time (s)')
  pp.setHLimits(s2.getFirst(),s2.getLast())
  pp.setVLimits(s1.getFirst(),s1.getLast())
  pp.setHInterval(hint)
  pp.setVInterval(vint)
  pv = pp.addPixels(s1,s2,f) 
  pv.setColorModel(ColorMap.GRAY)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setClips(-6.0,6.0)
  pf = PlotFrame(pp)
  #pf.setSize((int)(width-width*perc),height)
  pf.setSize((int)(600-600*perc),height)
  pf.setFontSizeForSlide(0.8,1.0)
  pf.setVisible(True)
  pf.paintToPng(300,5,pngDir+'1.png')

def cmpPlot1(f,s1,s2):
  pp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pp.setHLabel('Offset (km)')
  pp.setVLabel('Time (s)')
  pp.setHLimits(sx.getFirst(),sx.getLast())
  pp.setVLimits(s1.getFirst(),s1.getLast())
  pp.setHInterval(hint)
  pp.setVInterval(vint)
  # reverses x-axis for viking graben cmp 
  nx = 60 
  nt = 1500
  temp = [[0.0 for x in xrange(nt)] for y in xrange(nx)] 
  ixx = nx
  for ix in range(nx):
    ixx = ixx-1;
    for it in range(nt):
      temp[ix][it] = f[ixx][it]
  cb = pp.addColorBar('Amplitude')
  cb.setInterval(6.0)
  cb.setWidthMinimum(120)
  pv = pp.addPixels(s1,sx,temp) # for viking graben cmp
  pv.setColorModel(ColorMap.GRAY)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setClips(-6.0,6.0)
  pf = PlotFrame(pp)
  pf.setSize((int)(width-width*perc),height)
  pf.setFontSizeForSlide(0.8,1.0)
  pf.setVisible(True)
  pf.paintToPng(300,5,pngDir+'1.png')

def plot(f,s1,s2,png):
  p = semblancePanel()
  p.setHLimits(s2.getFirst(),s2.getLast())
  p.setVLimits(s1.getFirst(),s1.getLast())
  pixel = p.addPixels(s1,s2,f)
  pixel.setColorModel(ColorMap.JET)
  pixel.setInterpolation(PixelsView.Interpolation.LINEAR)
  pixel.setClips(0.0,1.0)
  contour(p,f,s1,s2)
  frame(p,png)

def semblancePanel():
  p = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  #p.setHLimits(1.5,3.5)
  #p.setVLimits(0.0,4.0)
  p.setHLabel('Velocity (km/s)')
  p.setVLabel('Time (s)')
  p.setHInterval(hint)
  p.setVInterval(vint)
  cb = p.addColorBar('Semblance')
  cb.setInterval(1.0)
  cb.setWidthMinimum(120)
  return p

def contour(panel,f,s1,s2):
  contour = panel.addContours(s1,s2,f)
  contour.setLineColor(Color.MAGENTA)
  contour.setLineWidth(contourwidth)
  contour.setContours(contours)

def frame(panel,png):
  frame = PlotFrame(panel)
  frame.setSize((int)(width-width*perc),height)
  frame.setFontSizeForSlide(0.8,1.0)
  frame.setVisible(True)
  frame.paintToPng(300,5,pngDir+png+'.png')

def bPlot(f,s1,s2,png):
  p = semblancePanel()
  p.setHLimits(s2.getFirst(),s2.getLast())
  p.setVLimits(s1.getFirst(),s1.getLast())
  cb = p.addColorBar('b')
  cb.setInterval(1.0)
  cb.setWidthMinimum(120)
  pv = p.addPixels(s1,s2,f)
  pv.setColorModel(ColorMap.JET)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setClips(0.0,1.0)
  frame(p,png)

def cCurve(semc,semw,t,st,sv):
  ccc = Velan.coherenceCurve(t,semc,sv,st)
  ccw = Velan.coherenceCurve(t,semw,sv,st)
  pp = PlotPanel()
  pp.setHLabel('Velocity (km/s)')
  pp.setVLabel('Semblance')
  pp.setHLimits(2.0,3.5)
  pp.setHInterval(0.5)
  pp.setVInterval(0.2)
  pvc = pp.addPoints(sv,ccc)
  pvw = pp.addPoints(sv,ccw)
  pvc.setLineWidth(3.0)
  pvw.setLineStyle(PointsView.Line.DASH)
  pvw.setLineWidth(3.0)
  pf = PlotFrame(pp)
  pf.setSize(height+200,height)
  pf.setFontSizeForSlide(0.8,1.0)
  pf.setVisible(True)
  pf.paintToPng(600,5,pngDir+'5.png')

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
  pv.setLineWidth(2);
  pv.setLineColor(Color.red)
  pf = PlotFrame(pp)
  pf.setSize(1500,900)
  pf.setFontSizeForPrint(8,225)
  pf.setVisible(True)
  pf.paintToPng(1200,3.33,pngDir+'bcurve.png')

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
  pf1.paintToPng(300,2.0,pngDir+'teaser1.png')
  
  width1 = (int)(800+800*tperc)
  fsize = (int)(150+150*tperc)

  pp2 = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,
                      PlotPanel.AxesPlacement.TOP)
  pp2.setHLabel('Velocity (km/s)')
  pp2.setHLimits(hlo,hhi)
  pp2.setVLimits(vlo,vhi)
  pp2.setHInterval(0.5)
  pv2 = pp2.addPixels(st,sv,semc)
  pv2.setColorModel(ColorMap.JET)
  pv2.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv2.setClips(0.0,1.0)
  contour(pp2,semc,st,sv)
  cb2 = pp2.addColorBar()
  cb2.setInterval(10.0)
  cb2.setWidthMinimum(10)
  pf2 = PlotFrame(pp2)
  pf2.setSize((int)(width1-width1*tperc1),800)
  pf2.setFontSizeForPrint(8,(int)(fsize-fsize*tperc1))
  pf2.setVisible(True)
  pf2.paintToPng(300,3.0,pngDir+'teaser2.png')

  pp3 = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.TOP)
  pp3.setHLabel('Velocity (km/s)')
  pp3.setHLimits(hlo,hhi)
  pp3.setVLimits(vlo,vhi)
  pp3.setHInterval(0.5)
  pv3 = pp3.addPixels(st,sv,semw)
  pv3.setColorModel(ColorMap.JET)
  pv3.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv3.setClips(0.0,1.0)
  contour(pp3,semw,st,sv)
  cb3 = pp3.addColorBar('Semblance')
  cb3.setInterval(1.0)
  cb3.setWidthMinimum(110)
  pf3 = PlotFrame(pp3)
  pf3.setSize((int)(width1+width1*tperc1),800)
  pf3.setFontSizeForPrint(8,(int)(fsize+fsize*tperc1))
  pf3.setVisible(True)
  pf3.paintToPng(300,3.0,pngDir+'teaser3.png')


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
