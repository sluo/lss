'''
Semblance
'''
from imports import *

pngDir=None
#pngDir='/Users/sluo/Desktop/'
seisDir = './bp2011/'
st = Sampling(707,0.0255,0.0)
#sx = Sampling(129,0.312,-20.007) 
sx = Sampling(151,0.312,-23.517) # for cmp1400
sv = Sampling(301,0.020,1.0)
nt,nx,nv = st.count,sx.count,sv.count
dt,dx,dv = st.delta,sx.delta,sv.delta
ft,fx,fv = st.first,sx.first,sv.first

#############################################################################

def main(args):
  #show()
  #test()
  spectrum()

def spectrum():
  #c = read('cmp1')
  #c = read('cmp2')
  #c = read('cmp3')
  c = read('cmp1400')
  ss = SemblanceSpectrum(st,sx,sv)
  s = ss.apply(c)
  plot(c,perc=99,name='cmp')
  plot(s,cmap=jet,perc=99.0,s1=st,s2=sv,name='semblance')

def test():
  vp = [2.00,3.00]
  fpeak = 5.0;
  snr = 1.0e6;
  vps = SemblanceSpectrum.makeLinearVelocity(vp[0],vp[1],st);
  p = SemblanceSpectrum.makeRickerGather(fpeak,snr,vps,st,sx);
  ss = SemblanceSpectrum(st,sx,sv);
  s = ss.apply(p);
  plot(p,name='cmp')
  plot(s,cmap=jet,s1=st,s2=sv,name='semblance')

def show():
  c1 = read('cmp1')
  c2 = read('cmp2')
  c3 = read('cmp3')
  plot(c1,perc=99,name='cmp1')
  plot(c2,perc=99,name='cmp2')
  plot(c3,perc=99,name='cmp3')

def like(x):
  return zerofloat(len(x[0]),len(x))

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,et=None,cmap=gray,s1=None,s2=None,
         cmin=0,cmax=0,perc=100,cbar=None,name=None):
  pan = panel(cbar)
  pix = pan.addPixels(x)
  if s1 and s2:
    pix.set(s1,s2,x)
  pix.setColorModel(cmap)
  pix.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pix.setClips(cmin,cmax)
  if perc<100:
    pix.setPercentiles(100-perc,perc)
  if et:
    tv = TensorsView(et)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(2)
    tv.setEllipsesDisplayed(20)
    tv.setScale(0.9)
    pan.getTile(0,0).addTiledView(tv)
  frame(pan,name)

def panel(cbar=None):
  p = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  cb = p.addColorBar()
  cb.setWidthMinimum(120)
  if cbar:
    cb.setLabel(cbar)
  return p

def frame(panel,name=None):
  frame = PlotFrame(panel)
  frame.setBackground(Color(204,204,204,255))
  #frame.setFontSizeForSlide(1.0,1.0)
  frame.setSize(1200,800)
  if name:
    frame.setTitle(name)
  frame.setVisible(True)
  if name and pngDir:
    frame.paintToPng(360,3.0,pngDir+name+'.png')

def read(name,image=None):
  if not image:
    image = zerofloat(nt,nx)
  fileName = seisDir+name+'.dat'
  ais = ArrayInputStream(fileName)
  #ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image

def write(name,image,directory=seisDir):
  fileName = directory+name+'.dat'
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()

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
    print '%02d:%02d:%02d'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
