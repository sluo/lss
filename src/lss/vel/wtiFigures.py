from imports import *
from wtiDiving import *

#############################################################################
sz = Sampling(265,0.012,0.0)
sx = Sampling(767,0.012,0.0)
st = Sampling(4000,0.0015,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta

pngDir = None
pngDir = '/Users/sluo/Desktop/png/'

gfile = None
pfile = None
sfile = None
mfile = None
#gfile = '/Users/sluo/Dropbox/save/fwi/marmousi/5hz/2000m_100p/fwi/g_0.dat'
#pfile = '/Users/sluo/Dropbox/save/fwi/marmousi/5hz/2000m_95p/alt/p_iter4.dat'
#sfile = '/Users/sluo/Dropbox/save/fwi/marmousi/5hz/2000m_100p/split4a/s_19.dat'
#sfile = '/Users/sluo/Dropbox/save/lsm/marmousi/95p/amp4/s0_true.dat'
mfile = '/Users/sluo/Dropbox/save/lsm/marmousi/100p/fwi2/s1_19.dat'
nfile = '/Users/sluo/Dropbox/save/lsm/marmousi/100p/fwi2/s1_true.dat'

#############################################################################

#widthPoints,widthInches = 240.0,3.33 # 1 column
widthPoints,widthInches = 260.0,3.61 # 1/2 page
#widthPoints,widthInches = 469.0,6.51 # 2 column (full page)

def main(args):
  plotFiles()
  #plotObjectiveFunction() # TODO

def plotFiles():
  setModel()
  if sfile is not None:
    s = read(sfile)
    plot(s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',title='s')
  if gfile is not None:
    g = read(gfile)
    div(g,max(abs(g))/1.0,g)
    plot(g,cmap=rwb,cmin=-1.0,cmax=1.0,title='g')
    #plot(g,cmap=rwb,sperc=99.0,title='g')
  if pfile is not None:
    p = read(pfile)
    div(p,max(abs(p)),p)
    plot(p,cmap=rwb,cmin=-1.0,cmax=1.0,title='p')
    #plot(p,cmap=rwb,sperc=99.0,title='p')
  if mfile is not None:
    m = read(mfile)
    #n = read(mfile.split('/')[:-1]+'s1_true.dat')
    n = read(nfile)
    clip = max(abs(n))
    #plot(m,cmap=rwb,sperc=100,cbar='Reflectivity s2/km2',title='m')
    plot(m,cmap=rwb,cmin=-clip,cmax=clip,cbar='Reflectivity',title='m')
    plot(n,cmap=rwb,cmin=-clip,cmax=clip,cbar='Reflectivity',title='n')

def plotObjectiveFunction():
  setModel()
  do = sub(modelData(_t),modelDirectArrival(_t))[ns/2] # observed data
  niter = 4
  r = zerofloat(1+niter,2)

  r[0][0] = r[0][1] = sum(mul(do,do))

  for i in range(1,niter+1):
    pre = '/Users/sluo/Dropbox/save/marmousi/5hz/2000m_100p/'

    s = read(pre+'fwi/dat/s_'+str(i)+'.dat')
    ds = sub(modelData(s),modelDirectArrival(s))[ns/2]  # simulated data
    t = sub(ds,do)
    r[0][i] = sum(mul(t,t))

    s = read(pre+'split4/dat/s_'+str(i)+'.dat')
    ds = sub(modelData(s),modelDirectArrival(s))[ns/2]  # simulated data
    t = sub(ds,do)
    r[1][i] = sum(mul(t,t))

  div(r,max(r),r)

  sp = SimplePlot()
  sp.addPoints(r[0]).setMarkStyle(PointsView.Mark.FILLED_SQUARE)
  sp.addPoints(r[1])

def setModel():
  global _t,_s
  #_t,_s = getGaussian()
  #_t,_s = getMarmousi(0.2,1.0)
  #_t,_s = getMarmousi(0.5,1.0)
  _t,_s = getMarmousi(2.0,1.0)
  #_t,_s = getMarmousi(2.0,0.95)
  #plot(_t,cmap=jet,cbar='Slowness (s/km)',title='s_true')
  #plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
  #     title='s_init')

import socket
def getMarmousi(sigmaC=0.5,smul=1.0):
  p = zerofloat(751,2301)
  if socket.gethostname()=='backus.Mines.EDU':
    read("/data/sluo/marmousi/marmousi.dat",p)
  else:
    read("/data/seis/marmousi/marmousi.dat",p)
  p = copy(743,2301,8,0,p)
  div(1000.0,p,p) # convert to slowness
  q = copy(p)
  refC = RecursiveExponentialFilter(sigmaC/0.004)
  refC.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  refC.apply(q,q)
  mul(smul,q,q) # scale slowness, e.g., so that simulated data arrives earlier
  v = fillfloat(2.0/3.0,nz,nx)
  c = copy(v)
  copy(248,767,0,0,3,3,p,17,0,1,1,v)
  copy(248,767,0,0,3,3,q,17,0,1,1,c)
  return v,c

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(f,cmap=gray,cmin=0,cmax=0,perc=100,sperc=None,cbar=None,title=None):
  panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  cb = panel.addColorBar()
  if cbar:
    cb.setLabel(cbar)
  if widthPoints==240.0:
    cb.setWidthMinimum(150)
  if widthPoints==260.0:
    cb.setWidthMinimum(140)
  elif widthPoints==469.0:
    cb.setWidthMinimum(80)
  if len(f[0])==nz and len(f)==nx:
    pixel = panel.addPixels(sz,sx,f)
    panel.setHLabel('Distance (km)')
    panel.setVLabel('Depth (km)')
  elif len(f[0])==nt and len(f)==nr:
    pixel = panel.addPixels(st,Sampling(len(f)),f)
    panel.setHLabel('Receiver')
    panel.setVLabel('Time (s)')
  else:
    pixel = panel.addPixels(f)
  pixel.setColorModel(cmap)
  if cmin<cmax:
    pixel.setClips(cmin,cmax)
  if perc<100:
    pixel.setPercentiles(100-perc,perc)
  if sperc is not None: # symmetric percentile clip (for plotting gradients)
    clips = Clips(100-sperc,sperc,f)
    clip = max(abs(clips.getClipMin()),abs(clips.getClipMax()))
    pixel.setClips(-clip,clip)
  pixel.setInterpolation(PixelsView.Interpolation.LINEAR)
  frame = PlotFrame(panel)
  #frame.setFontSizeForSlide(1.5,1.5)
  frame.setFontSizeForPrint(8.0,widthPoints)
  if (len(f[0])==nz):
    if widthPoints==240.0:
      frame.setSize(1200,468)
    if widthPoints==260.0:
      frame.setSize(1200,465)
    elif widthPoints==469.0:
      frame.setSize(1200,445)
  else:
    frame.setSize(1200,1000)
  if title:
    frame.setTitle(title)
  frame.setVisible(True)
  if title and pngDir:
    frame.paintToPng(720,widthInches,pngDir+title+'.png')
  return panel

def read(name,image=None):
  if not image:
    image = zerofloat(nz,nx)
  fileName = name
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def write(fname,image):
  aos = ArrayOutputStream(fname)
  aos.writeFloats(image)
  aos.close()

def cleanDir(dir):
  os.chdir(dir)
  for f in os.listdir(dir):
    os.remove(f)

#############################################################################
# Do everything on Swing thread.
import sys,time
class RunMain(Runnable):
  def run(self):
    start = time.time()
    if pngDir is not None:
      print 'cleaning pngDir'
      cleanDir(pngDir)
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02d'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
