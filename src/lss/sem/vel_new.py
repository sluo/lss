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
from fwi import *

#############################################################################
# parameters

#dir = '/Users/sluo/Home/doc/Research/Semblance/temp/'
dir = '/Users/sluo/Desktop/'
datadir = '/Users/sluo/Home/box/idh/trunk/bench/src/sim/cdps/'

st = Sampling(1001,0.004,0.0) # sampling for sythetic data
sx = Sampling(60,0.050,0.050)
#sv = Sampling(101,0.020,1.5)
sv = Sampling(501,0.004,1.5) # for smooth semblance curve
nt,nx,nv = st.count,sx.count,sv.count
dt,dx,dv = st.delta,sx.delta,sv.delta
ft,fx,fv = st.first,sx.first,sv.first

sst = Sampling(1500,0.004,0.0) # sampling for Vikin Graben data
ssx = Sampling(60,0.050,-3.212)
ssv = Sampling(101,0.020,1.5)

vp = (2.00,3.00)
#vm = (1.98,2.70)
vm = (1.98,2.50)
#vm = None

#vp = (2.00,2.91) # XXX
#vm = (1.98,2.70) # XXX
#vm = None

fpeak = 25.0
tsigma = 8.0 # smoothing window
#tsigma = 12.0
snr = 1.0e8
#snr = 1.0e0

width = 600
height = 1000

#############################################################################
# functions

def main(args):
	makeSynthetic()
  #makeReal()
  #test()

def test():
  old = zerofloat(nt,nv)
  new = zerofloat(nt,nv)
  read("/Users/sluo/Desktop/swOld.dat",old)
  read("/Users/sluo/Desktop/swNew.dat",new)
  dif = sub(old,new)
  print max(dif)
  plot(dif,st,sv,0.0,0.0,True,jet,'Velocity (km/s)','Difference','dif')

def makeSynthetic():
  vps = Velan.makeLinearVelocity(vp[0],vp[1],st)
  p = Velan.makeRickerGather(fpeak,vps,st,sx)
  if vm:
    vms = Velan.makeLinearVelocity(vm[0],vm[1],st)
    p = add(p,Velan.makeRickerGather(fpeak,vms,st,sx)) # multiples
  p = Velan.addRandomNoise(snr,p) # noise

  b,sc,sw = zerofloat(nt,nv),zerofloat(nt,nv),zerofloat(nt,nv)
  sc = Velan.velocitySpectrum(st,sx,p,sv,tsigma,False)
  sw = Velan.velocitySpectrum(st,sx,p,sv,tsigma,True,b)

  # plots
  r = 0.5
  cmin = r*min(p)
  cmax = r*max(p)
  plot(p,st,sx,cmin,cmax,False,gray,'Offset (km)','Amplitude','1') # cmp
  #plot(b,st,sv,0.0,1.0,False,jet,'Velocity (km/s)','b value','2') # b
  #plot(sc,st,sv,0.0,1.0,True,jet,'Velocity (km/s)','Semblance','3') # sc
  #plot(sw,st,sv,0.0,1.0,True,jet,'Velocity (km/s)','Semblance','4') # sw
  plot(b, st,sv,0.0,1.0,False,gray,'Velocity (km/s)','b value','2') # b
  plot(sc,st,sv,0.0,1.0,False,gray,'Velocity (km/s)','Semblance','3') # sc
  plot(sw,st,sv,0.0,1.0,False,gray,'Velocity (km/s)','Semblance','4') # sw

  #semblanceCurve(sc,sw,1.5,"curve1")
  semblanceCurve(sc,sw,3.2,"curve1")
  #semblanceCurve(sc,sw,3.5,"curve2")

  write("sw",sw)

  #cCurve(sc,sw,3.2,st,sv)
#  #bCurve(n0,n1,h0,h1,d0,d1)
#  #teaserPlots(p,sc,sw,st,sv,sx)

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

pixels = 690
points = pixels*120.0/600.0
inches = pixels*1.667/600.0
fwidth,fheight = 716,1000
hint,vint = 1,1

gray = ColorMap.GRAY
jet = ColorMap.JET
def plot(f,s1,s2,cmin,cmax,contour,cmap,hlabel,cbarlabel,name):
  pp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pp.setHLabel(hlabel)
  pp.setVLabel('Time (s)')
  pp.setHLimits(s2.getFirst(),s2.getLast())
  pp.setVLimits(s1.getFirst(),s1.getLast())
  pp.setHInterval(hint)
  pp.setVInterval(vint)
  cb = pp.addColorBar(cbarlabel)
  cb.setInterval(1.0)
  cb.setWidthMinimum(80)
  pv = pp.addPixels(s1,s2,f)
  pv.setColorModel(cmap)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if (contour):
    cv = pp.addContours(s1,s2,f)
    cv.setLineColor(Color.MAGENTA)
    cv.setLineWidth(8.0)
    cv.setContours([0.2])

    #pvp = pp.addPoints([0.0,4.0],[2.00,3.0])
    #pvm = pp.addPoints([0.0,4.0],[1.98,2.7])
    #pvp.setLineWidth(4.0)
    #pvm.setLineWidth(4.0)

  pf = PlotFrame(pp)
  pf.setSize(fwidth,fheight)
  pf.setFontSizeForPrint(8,points)
  pf.setVisible(True)
  pf.paintToPng(720,inches,dir+name+'.png')

def semblanceCurve(sc,sw,tzero,name='curve'):
  c = zerofloat(nv)
  w = zerofloat(nv)
  kt = int(tzero/dt)
  print 'kt=',kt
  for iv in range(nv):
    c[iv] = sc[iv][kt]
    w[iv] = sw[iv][kt]
  pp = PlotPanel()
  pp.setHLabel('Velocity (km/s)')
  pp.setVLabel('Semblance')
  pp.setHLimits(2.0,3.5)
  pp.setVLimits(0.0,1.0)
  pp.setHInterval(0.5)
  pp.setVInterval(0.2)
  pvc = pp.addPoints(sv,c)
  pvw = pp.addPoints(sv,w)
  
  # true NMO velocity
  #m = zerofloat(nv)
  #v2 = vp[0]+(tzero/4.0)*(vp[1]-vp[0])
  #kv2 = int((v2-fv)/dv)+1
  #kv1 = 0
  #if vm:
  #  v1 = vm[0]+(tzero/4.0)*(vm[1]-vm[0])
  #  kv1 = int((v1-fv)/dv)
  #for iv in range(kv1):
  #  m[iv] = -10
  #for iv in range(kv1,kv2):
  #  m[iv] = 10
  #for iv in range(kv2,nv):
  #  m[iv] = -10
  #pvm = pp.addPoints(sv,m)
  ##pvm.setLineStyle(PointsView.Line.DOT)
  #pvm.setLineWidth(3.0)

  pvw.setLineStyle(PointsView.Line.DASH)
  pvc.setLineWidth(3.0)
  pvw.setLineWidth(3.0)
  pf = PlotFrame(pp)
  pf.setSize(1000,800)
  pf.setFontSizeForPrint(8,240)
  pf.setVisible(True)
  pf.paintToPng(720,3.33,dir+name+'.png')

def read(name,image):
  InversionUtil.readFloats(name,image)

def write(name,image):
  saveDir = '/Users/sluo/Desktop/'
  filename = saveDir+name+'.dat'
  InversionUtil.writeFloats(filename,image)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
