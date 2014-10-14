"""
Semblance for Oz Yilmaz's shot gathers
"""
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

from lss.sem import *
from lss.util import *

#############################################################################

tsigma = 8.0
datDir = '/data/seis/oz/'
pngDir = None
#pngDir = '/Users/sluo/Desktop/'

#############################################################################

def main(args):
  goSemblance()

def goSemblance():
  go('oz03',Sampling(501,0.010,1.0),0.4,6.0)
  go('oz08',Sampling(401,0.010,3.0),0.0,5.0)
  go('oz16',Sampling(201,0.010,1.0),0.0,5.0,reverse=True)
  go('oz30',Sampling(201,0.010,1.0),1.0,8.0,reverse=True)

def go(name,sv,vmin=0,vmax=0,reverse=False):
  print name
  st,sx = getParms(name)
  nt,nx,nv = st.count,sx.count,sv.count
  p = zerofloat(nt,nx)
  readImage(name+'.bin',p)
  if reverse==True:
    p = reverseGather(p)
  gain(p) # Gain

  b = zerofloat(nt,nv)
  sc = Velan.velocitySpectrum(st,sx,p,sv,tsigma,False)
  sw = Velan.velocitySpectrum(st,sx,p,sv,tsigma,True,b)

  r = 0.03
  cmin,cmax = r*min(p),r*max(p)
  png1,png2,png3,png4 = name+'_1',name+'_2',name+'_3',name+'_4',
  plot(p,st,sx,'Offset (km)','Amplitude',vmin,vmax,cmin,cmax,png=png1)
  #plot(b,st,sv,'Velocity (km/s)','b',vmin,vmax,0.0,1.0,jet,png2)
  #plot(sc,st,sv,'Velocity (km/s)','Semblance',vmin,vmax,0.0,0.5,jet,png3)
  #plot(sw,st,sv,'Velocity (km/s)','Semblance',vmin,vmax,0.0,0.5,jet,png4)
  plot(b,st,sv,'Velocity (km/s)','b value',vmin,vmax,0.0,1.0,png=png2)
  plot(sc,st,sv,'Velocity (km/s)','Semblance',vmin,vmax,0.0,0.5,png=png3)
  plot(sw,st,sv,'Velocity (km/s)','Semblance',vmin,vmax,0.0,0.5,png=png4)

def getParms(fname):
  infile = open(datDir+fname+'.H','r')
  text = infile.read()
  infile.close()
  def get(par):
    index = text.find(par)
    start = index+3
    end = start
    while text[end]!=' ':
      end += 1
    #print(text[start:end])
    return text[start:end]
  n1,d1,f1 = int(get('n1')),float(get('d1')),float(get('f1'))
  n2,d2,f2 = int(get('n2')),float(get('d2')),float(get('f2'))
  #print('n1=%d, d1=%f, f1=%f'%(n1,d1,f1))
  #print('n2=%d, d2=%f, f2=%f'%(n2,d2,f2))
  return Sampling(n1,d1,f1),Sampling(n2,d2,f2)

def reverseGather(f):
  g = copy(f)
  n2 = len(f)
  for i in range(n2):
    f[i] = g[n2-1-i]
  return f

def gain(f):
  n2 = len(f)
  n1 = len(f[0])
  for i in range(n2):
    for j in range(n1):
      r = j+1.0
      #f[i][j] = f[i][j]*r
      f[i][j] = f[i][j]*r*r

#############################################################################
# plots

pixels = 690
points = pixels*120.0/600.0
inches = pixels*1.667/600.0
fwidth,fheight = 716,1000
hint,vint = 1,1

gray = ColorMap.GRAY
jet = ColorMap.JET
def plot(f,s1,s2,hlabel,cbarlabel,vmin,vmax,cmin,cmax,cmap=gray,png=None):
  pp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pp.setHLabel(hlabel)
  pp.setVLabel('Time (s)')
  if vmin<vmax:
    pp.setVLimits(vmin,vmax)
  pp.setHInterval(hint)
  pp.setVInterval(vint)
  cb = pp.addColorBar(cbarlabel)
  cb.setInterval(cmax)
  cb.setWidthMinimum(80)
  pv = pp.addPixels(s1,s2,f)
  pv.setColorModel(cmap)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  pf = PlotFrame(pp)
  pf.setSize(fwidth,fheight)
  pf.setFontSizeForPrint(8,points)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(720,inches,pngDir+png+'.png')

def readImage(name,y):
  ais = ArrayInputStream(datDir+name,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(y)
  ais.close()

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
