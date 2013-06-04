from imports import *
import wtiDiving
import lsrtm

#############################################################################
#sz = Sampling(301,0.012,0.0)
#sx = Sampling(501,0.012,0.0)
#st = Sampling(2750,0.0015,0.0)
#nz,nx,nt = sz.count,sx.count,st.count
#dz,dx,dt = sz.delta,sx.delta,st.delta

sz = Sampling(265,0.012,0.0)
sx = Sampling(767,0.012,0.0)
st = Sampling(2700,0.0015,0.0)
#st = Sampling(4001,0.0015,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta

pngDir = None
pngDir = '/Users/sluo/Desktop/pngdat/'

widthPoints = None # slides
#widthPoints = 240.0 # 1 column
#widthPoints = 260.0 # 1/2 page
#widthPoints = 469.0 # 2 column (full page)

#############################################################################

def main(args):
  plotFwi()
  #plotLsm()
  #plotData()
  #plotBornData()
  #plotObjectiveFunction()
  #plotFwiObjectiveFunctions()
  #plotLsmObjectiveFunctions()
  #plotFiles()

def plotData():
  ddir,iiter = '/Users/sluo/Dropbox/save/fwi/marmousi/2000m_100p/fwi/',None
  t = read(ddir+'s_true.dat')
  s = read(ddir+'s_'+iiter+'.dat') if iiter else read(ddir+'s_init.dat')
  do = sub(modelData(t),modelDirectArrival(t))
  ds = sub(modelData(s),modelDirectArrival(s))
  cmin = -0.8*max(abs(do))
  cmax =  0.8*max(abs(do))
  plot(do,cmin=cmin,cmax=cmax,cbar='Amplitude',title='do')
  plot(ds,cmin=cmin,cmax=cmax,cbar='Amplitude',title='ds')
  plot(t,cmap=jet,cbar='Slowness (s/km)',title='t')
  plot(s,cmap=jet,cbar='Slowness (s/km)',title='s')
  ra,rt,rc,dw,v = like(do),like(do),like(do),like(do),like(do)
  wtiDiving.makeWarpedResidual(ds,do,ra,rt,rc,dw,v,reverseOrder=False)
  plotResiduals(do,ds,ra,rt,rc,dw,v)
def modelData(s):
  fpeak = 5.0
  kxs,kzs = nx/2,0
  kxr,kzr = rampint(0,1,nx),fillint(0,nx)
  d = zerofloat(nt,nx)
  wave = Wavefield(sz,sx,st)
  source = Wavefield.RickerSource(fpeak,kzs,kxs)
  receiver = Wavefield.Receiver(kzr,kxr)
  wave.modelAcousticData(source,receiver,s,d)
  return d
def modelDirectArrival(s):
  t = copy(s)
  # Find water bottom then flood below
  t00 = t[0][0]
  zbot = 0
  while t[0][zbot]==t00:
    zbot += 1
  for ix in range(nx):
    for iz in range(zbot,nz):
      t[ix][iz] = t[ix][zbot] # flood
  # Scale and smooth model
  t = transpose(t)
  sigma = 20.0
  for iz in range(zbot+3,nz):
    r = 1.0+0.1*(iz-zbot-3.0)/(nz-zbot-10.0)
    mul(r,t[iz],t[iz])
    ref = RecursiveExponentialFilter(sigma)
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
    ref.apply(t[iz],t[iz])
    sigma += 1.0
  t = transpose(t)
  return modelData(t)

def plotBornData():
  #ddir,iiter = '/Users/sluo/Dropbox/save/lsm/marmousi/95p/amp5/','19'
  ddir,iiter = '/Users/sluo/Dropbox/save/lsm/layered1/95p/amp/','9'
  t0 = read(ddir+'s0_true.dat')
  t1 = read(ddir+'s1_true.dat')
  s0 = read(ddir+'s0_init.dat')
  s1 = read(ddir+'s1_'+iiter+'.dat')
  if nx==501:
    mask(s1,int(nz/4))
  do = modelBornData(t0,t1)
  if nx==767:
    GaussianTaper.apply(0.5,s1,s1)
  ds = modelBornData(s0,s1)
  plot(t0,cmap=jet,cbar='Slowness (s/km)',title='t0')
  plot(t1,cmap=rwb,cmin=-max(abs(t1)),cmax=max(abs(t1)),
    cbar='Reflectivity',title='t1')
  plot(s0,cmap=jet,cmin=min(t0),cmax=max(t0),cbar='Slowness (s/km)',title='s0')
  plot(s1,cmap=rwb,cmin=-max(abs(t1)),cmax=max(abs(t1)),
    cbar='Reflectivity',title='s1')
  ra,rt,rc,dw,v = like(do),like(do),like(do),like(do),like(do)
  lsrtm.makeWarpedResidual(ds,do,ra,rt,rc,dw,v)
  plotResiduals(do,ds,ra,rt,rc,dw,v)
def modelBornData(s0,s1):
  fpeak = 10.0
  kxs,kzs = nx/2,0
  kxr,kzr = rampint(0,1,nx),fillint(0,nx)
  d = zerofloat(nt,nx)
  u = zerofloat(nz,nx,nt)
  wave = Wavefield(sz,sx,st)
  wave.modelAcousticWavefield(
    Wavefield.RickerSource(fpeak,kzs,kxs),s0,u)
  wave.modelAcousticData(
    Wavefield.WavefieldSource(dt,s1,u),
    Wavefield.Receiver(kzr,kxr),
    s0,d)
  return d

def plotResiduals(do,ds,ra,rt,rc,dw,v):
  cmin,cmax = -0.3*max(abs(do)),0.3*max(abs(do))
  #rmin,rmax = cmin,cmax
  rmin,rmax = -1.0*max(abs(ra)),1.0*max(abs(ra))
  plot(do,cmin=cmin,cmax=cmax,cbar='Amplitude',title='do')
  plot(ds,cmin=cmin,cmax=cmax,cbar='Amplitude',title='ds')
  plot(dw,cmin=cmin,cmax=cmax,cbar='Amplitude',title='dw')
  plot(mul(dt,v),cmap=rwb,sperc=100.0,cbar='Traveltime shift (s)',title='v')
  plot(sub(ds,do),cmin=rmin,cmax=rmax,cbar='Amplitude',title='rd')
  plot(ra,cmin=rmin,cmax=rmax,cbar='Amplitude',title='ra')
  plot(rt,cmin=rmin,cmax=rmax,cbar='Amplitude',title='rt')
  plot(rc,cmin=rmin,cmax=rmax,cbar='Amplitude',title='rc')

def plotLsm():
  ddir,iiter = '/Users/sluo/Dropbox/save/lsm/marmousi/95p/fwi3/','19'
  #ddir,iiter = '/Users/sluo/Dropbox/save/lsm/layered1/95p/fwi/','9'
  t0 = read(ddir+'s0_true.dat')
  t1 = read(ddir+'s1_true.dat')
  s0 = read(ddir+'s0_init.dat')
  s1 = read(ddir+'s1_'+iiter+'.dat')
  g = read(ddir+'g_'+iiter+'.dat'); mul(g,1.0/max(abs(g)),g)
  p = read(ddir+'p_'+iiter+'.dat'); mul(p,1.0/max(abs(p)),p)
  mul(s1,0.15/max(abs(s1)),s1)
  plot(t0,cmap=jet,cmin=min(t0),cmax=max(t0),cbar='Slowness (s/km)',title='t0')
  plot(t1,cmap=rwb,cmin=-max(abs(t1)),cmax=max(abs(t1)),
    cbar='Reflectivity',title='t1')
  plot(s0,cmap=jet,cmin=min(t0),cmax=max(t0),cbar='Slowness (s/km)',title='s0')
  plot(s1,cmap=rwb,cmin=-max(abs(t1)),cmax=max(abs(t1)),
    cbar='Reflectivity',title='s1')
  plot(g,cmap=rwb,cmin=-1.0,cmax=1.0,title='g')
  plot(p,cmap=rwb,cmin=-1.0,cmax=1.0,title='p')

def plotFwi():
  iiter = '19'
  ddir = '/Users/sluo/Dropbox/save/fwi/marmousi/2000m_100p/dres2/'
  t = read(ddir+'s_true.dat')
  s = read(ddir+'s_init.dat')
  r = read(ddir+'s_'+iiter+'.dat')
  g = read(ddir+'g_'+iiter+'.dat'); mul(g,1.0/max(abs(g)),g)
  p = read(ddir+'p_'+iiter+'.dat'); mul(p,1.0/max(abs(p)),p)
  plot(g,cmap=rwb,cmin=-1.0,cmax=1.0,title='g')
  plot(p,cmap=rwb,cmin=-1.0,cmax=1.0,title='p')
  plot(t,cmap=jet,cbar='Slowness (s/km)',title='t')
  plot(s,cmap=jet,cmin=min(t),cmax=max(t),cbar='Slowness (s/km)',title='s_init')
  plot(r,cmap=jet,cmin=min(t),cmax=max(t),cbar='Slowness (s/km)',
    title='s_'+iiter)

def plotFiles():
  setModel()
  gfile = None
  pfile = None
  sfile = None
  mfile = None
  #gfile = '/Users/sluo/Dropbox/save/fwi/marmousi/5hz/2000m_100p/fwi/g_0.dat'
  #pfile = '/Users/sluo/Dropbox/save/fwi/marmousi/5hz/2000m_95p/alt/p_iter4.dat'
  sfile = '/Users/sluo/Dropbox/save/lsm/marmousi/100p/fwi2/s0_true.dat'
  mfile = '/Users/sluo/Dropbox/save/lsm/marmousi/95p/fwi2/s1_19.dat'
  nfile = '/Users/sluo/Dropbox/save/lsm/marmousi/100p/fwi2/s1_true.dat'
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
    #mul(m,0.12/max(abs(m)),m)
    #n = read(mfile.split('/')[:-1]+'s1_true.dat')
    n = read(nfile)
    clip = max(abs(n))
    #plot(m,cmap=rwb,sperc=100,cbar='Reflectivity s2/km2',title='m')
    plot(m,cmap=rwb,cmin=-clip,cmax=clip,cbar='Reflectivity',title='m')
    plot(n,cmap=rwb,cmin=-clip,cmax=clip,cbar='Reflectivity',title='n')

def plotObjectiveFunction():
  res = zerofloat(21)
  read('/Users/sluo/Dropbox/save/lsm/marmousi/100p/fwi2/dres.dat',res)
  div(res,res[0],res)
  panel = PlotPanel()
  panel.setHLabel('Iteration')
  #panel.setVLabel('Error')
  panel.setVLimits(0.0,1.1)
  pa = panel.addPoints(res)
  pa.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pa.setLineWidth(1.5)
  frame = PlotFrame(panel)
  frame.setFontSizeForSlide(1.0,1.0)
  frame.setSize(1000,700)
  frame.setVisible(True)
  if pngDir:
    frame.paintToPng(1024,2.0,pngDir+'obj.png')

def plotFwiObjectiveFunctions():
  dres_dres = zerofloat(21) # data residual for waveform inversion
  cres_dres = zerofloat(21) # data residual for combined inversion
  cres_cres = zerofloat(21) # combined residual for combined inversion
  ddir = '/Users/sluo/Dropbox/save/fwi/marmousi/2000m_100p/'
  read(ddir+'/dres2/dres.dat',dres_dres)
  read(ddir+'/cres6/dres.dat',cres_dres)
  read(ddir+'/cres6/cres.dat',cres_cres)
  def normal(x):
    div(x,x[0],x)
  normal(dres_dres)
  normal(cres_dres)
  normal(cres_cres)

  panel = PlotPanel()
  panel.setHLabel('Iteration')
  #panel.setVLabel('Error')
  panel.setVLimits(0.0,1.4)

  pd = panel.addPoints(cres_cres)
  pd.setLineColor(Color.RED)
  pd.setMarkColor(Color.RED)
  pd.setLineWidth(1.5)
  pd.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)

  pc = panel.addPoints(cres_dres)
  pc.setLineColor(Color.RED)
  pc.setMarkColor(Color.RED)
  pc.setLineWidth(1.5)
  pc.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)

  pa = panel.addPoints(dres_dres)
  pa.setLineColor(Color.BLACK)
  pa.setMarkColor(Color.BLACK)
  pa.setLineWidth(1.5)
  pa.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)

  frame = PlotFrame(panel)
  frame.setFontSizeForSlide(1.0,1.0)
  frame.setSize(1000,500)
  frame.setVisible(True)
  if pngDir:
    frame.paintToPng(1024,2.0,pngDir+'obj.png')
  plotTemplatesForLegend()

def plotLsmObjectiveFunctions():
  alsm_ra95p = zerofloat(21)
  alsm_rd95p = zerofloat(21)
  lsm_rd100p = zerofloat(21)
  lsm_rd95p = zerofloat(21)
  ddir = '/Users/sluo/Dropbox/save/lsm/marmousi/'
  read(ddir+'95p/ares5/ares.dat',alsm_ra95p)
  read(ddir+'95p/ares5/dres.dat',alsm_rd95p)
  read(ddir+'100p/dres2/dres.dat',lsm_rd100p)
  #read(ddir+'95p/dres2/dres.dat',lsm_rd95p)
  read(ddir+'95p/dres3/dres.dat',lsm_rd95p)

  def normal(x):
    div(x,x[0],x)
  normal(alsm_ra95p)
  normal(alsm_rd95p)
  normal(lsm_rd100p)
  normal(lsm_rd95p)
  #points(alsm_ra95p)
  #points(alsm_rd95p)
  #points(lsm_rd100p)
  #points(lsm_rd95p)

  panel = PlotPanel()
  panel.setHLabel('Iteration')
  #panel.setVLabel('Error')
  panel.setVLimits(0.0,1.7)

  pd = panel.addPoints(alsm_ra95p)
  pd.setLineColor(Color.RED)
  pd.setMarkColor(Color.RED)
  pd.setLineWidth(1.5)
  pd.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)

  pc = panel.addPoints(alsm_rd95p)
  pc.setLineColor(Color.RED)
  pc.setMarkColor(Color.RED)
  pc.setLineWidth(1.5)
  pc.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)

  pb = panel.addPoints(lsm_rd95p)
  pb.setLineColor(Color.BLUE)
  pb.setMarkColor(Color.BLUE)
  pb.setLineWidth(1.5)
  pb.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)

  pa = panel.addPoints(lsm_rd100p)
  pa.setLineColor(Color.BLACK)
  pa.setMarkColor(Color.BLACK)
  pa.setLineWidth(1.5)
  pa.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)

  frame = PlotFrame(panel)
  frame.setFontSizeForSlide(1.0,1.0)
  frame.setSize(1000,500)
  frame.setVisible(True)
  if pngDir:
    frame.paintToPng(1024,2.0,pngDir+'obj.png')
  plotTemplatesForLegend()

def plotTemplatesForLegend():
  def template(color,markStyle,title):
    t = zerofloat(3)
    panel = PlotPanel()
    panel.setVLimits(-1.0,1.0)
    point = panel.addPoints(t)
    point.setLineColor(color)
    point.setMarkColor(color)
    #point.setLineWidth(1.5)
    point.setMarkStyle(markStyle)
    frame = PlotFrame(panel)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.setSize(1000,500)
    frame.setVisible(True)
    if pngDir:
      frame.paintToPng(1024,2.0,pngDir+title+'.png')
  template(Color.BLACK,PointsView.Mark.FILLED_CIRCLE,'black_filled')
  template(Color.BLUE,PointsView.Mark.FILLED_CIRCLE,'blue_filled')
  template(Color.RED,PointsView.Mark.FILLED_CIRCLE,'red_filled')
  template(Color.RED,PointsView.Mark.HOLLOW_CIRCLE,'red_hollow')

def points(f,title=None):
  panel = PlotPanel()
  panel.setHLabel('Iteration')
  panel.setVLabel('Error')
  panel.setVLimits(0.0,1.1*max(f))
  pv = panel.addPoints(f)
  pv.setLineColor(Color.BLACK)
  pv.setLineWidth(2.0)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  frame = PlotFrame(panel)
  frame.setFontSizeForSlide(1.5,1.5)
  frame.setSize(1000,600)
  if title:
    frame.setTitle(title)
  frame.setVisible(True)
  if title and pngDir:
    frame.paintToPng(720,3.08,pngDir+title+'.png')
  
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

def mask(g,l1=0):
  n1,n2 = len(g[0]),len(g)
  for i2 in range(n2):
    for i1 in range(l1):
      g[i2][i1] = 0.0

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
  if len(f[0])==nz:
    cb.setInterval(0.1)
  if cbar:
    cb.setLabel(cbar)
  if widthPoints is None:
    cb.setWidthMinimum(140)
  elif widthPoints==240.0:
    cb.setWidthMinimum(150)
  elif widthPoints==260.0:
    cb.setWidthMinimum(140)
  if len(f[0])==nz and len(f)==nx:
    pixel = panel.addPixels(sz,sx,f)
    panel.setHLabel('Distance (km)')
    panel.setVLabel('Depth (km)')
  elif len(f[0])==nt and len(f)==nx:
    pixel = panel.addPixels(st,sx,f)
    panel.setHLabel('Distance (km)')
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
  if widthPoints is None:
    frame.setFontSizeForSlide(1.0,1.0)
  else:
    frame.setFontSizeForPrint(8.0,widthPoints)
  if (len(f[0])==nz):
    if widthPoints is None:
      #frame.setSize(1200,770)
      frame.setSize(1024,670)
    elif widthPoints==240.0:
      frame.setSize(1200,468)
    elif widthPoints==260.0:
      frame.setSize(1200,465)
    elif widthPoints==469.0:
      frame.setSize(1200,445)
  else:
    if widthPoints is None:
      #frame.setSize(1200,770)
      frame.setSize(1024,670)
    else:
      frame.setSize(1200,1000)
  if title:
    frame.setTitle(title)
  frame.setVisible(True)
  if title and pngDir:
    if widthPoints is None:
      frame.paintToPng(1024,2.0,pngDir+title+'.png')
    else:
      frame.paintToPng(720.0,widthPoints/72.0,pngDir+title+'.png')
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
