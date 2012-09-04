"""
Add faults to images
"""
from imports import *
import jarray

pngDir = None
#pngDir = '/Users/sluo/Desktop/'
seismicDir = "/data/seis/tp/"
n1,n2 = 251,751

#############################################################################

def main(args):
  show()
  #fault1()
  figure()

def figure():
  f = fakeData()
  [g,r1,r2,s1,s2,[f1,f2,x1,x2]] = Faulter.addOneFault(500.0,f)
  plot(f,cmap=gray,title="f",png="f")
  g = FlattenerUtil.applyShiftsR(f,mul(-1.0,array(r1,r2)))
  plot(g,cmap=gray,title="g",png="g")
  #h = copy(80,230,145,300,g)
  h = copy(60,200,170,320,g)
  def smallplot(x):
    sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
    sp.setFontSizeForPrint(8.0,240.0)
    #cb = sp.addColorBar()
    #cb.setWidthMinimum(130)
    sp.setSize(1450,900)
    pv = sp.addPixels(x)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    #pv.setColorModel(gray)
    if pngDir:
      sp.paintToPng(720,3.33,pngDir+'fault.png')
  smallplot(h)

def fault1():
  f = fakeData()
  [g,r1,r2,s1,s2,[f1,f2,x1,x2]] = Faulter.addOneFault(500.0,f)
  t1 = fillfloat(0.5*(max(f1)+min(f1)),n1,n2)
  t2 = fillfloat(0.5*(max(f2)+min(f2)),n1,n2)
  for i in range(len(f1)):
    i1,i2 = int(x1[i]),int(x2[i])
    t1[i2][i1] = f1[i]
    t2[i2][i1] = f2[i]
  g1,g2 = Grid2(f1,x1,x2),Grid2(f2,x1,x2)
  print "sibson..."
  q1 = g1.gridSibson(Sampling(n1),Sampling(n2))
  q2 = g2.gridSibson(Sampling(n1),Sampling(n2))
  print "done"
  h = Faulter.applyShiftsR(g,[r1,r2])
  k = Faulter.applyShiftsR(g,[q1,q2])
  plot(f,cmap=gray,title="input",png="f")
  plot(g,cmap=gray,title="faulted",png="g")
  plot(h,cmap=gray,title="unfaulted (true shift)",png="h")
  plot(k,cmap=gray,title="unfaulted (interpolated throw)",png="k")
  plot(r1,cmap=jet,title="r1",png="r1")
  plot(r2,cmap=jet,title="r2",png="r2")
  plot(s1,cmap=jet,title="s1",png="s1")
  plot(s2,cmap=jet,title="s2",png="s2")
  plot(q1,cmap=jet,title="interpolated throw1",png="q1")
  plot(q2,cmap=jet,title="interpolated throw2",png="q2")
  plot(t1,cmap=rwb,cmin=min(f1),cmax=max(f1),title="throw1",png="t1")
  plot(t2,cmap=rwb,cmin=min(f2),cmax=max(f2),title="throw2",png="t2")

def show():
  f = fakeData()
  [g,r1,r2,s1,s2,[f1,f2,x1,x2]] = Faulter.addOneFault(500.0,f)
  plot(f,cmap=gray,title="f",png="f")
  plot(g,cmap=gray,title="g",png="g")

def fakeData():
  return FakeData.seismic2d2011A(n1,n2,10.0)

def array(x1,x2,x3=None,x4=None):
  if x3 and x4:
    return jarray.array([x1,x2,x3,x4],Class.forName('[[F'))
  elif x3:
    return jarray.array([x1,x2,x3],Class.forName('[[F'))
  else:
    return jarray.array([x1,x2],Class.forName('[[F'))

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,cmap=gray,perc=100,cmin=0.0,cmax=0.0,\
         colorbar=None,title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #sp.setBackground(Color(204,204,204,255))
  #sp.setFontSizeForSlide(1.0,1.0)
  cb = sp.addColorBar()
  #cb.setWidthMinimum(130)
  if colorbar:
    cb.setLabel(colorbar)
  sp.setSize(1400,550)
  if title:
    sp.addTitle(title)
  pv = sp.addPixels(x)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cmap)
  if perc<100:
    pv.setPercentiles(100-perc,perc)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if png and pngDir:
    sp.paintToPng(360,3.0,pngDir+png+'.png')

def readImage(name,image):
  fileName = seismicDir+name+".dat"
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()

def writeImage(name,image):
  fileName = saveDir+name+'.dat'
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()

#############################################################################
# Run the function main on the Swing thread
import sys
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
run(main)
