from imports import *
import wtiDiving
import lsm

#############################################################################
#sz = Sampling(301,0.012,0.0)
#sx = Sampling(501,0.012,0.0)
#st = Sampling(2750,0.0015,0.0)
#nz,nx,nt = sz.count,sx.count,st.count
#dz,dx,dt = sz.delta,sx.delta,st.delta

sz = Sampling(265,0.012,0.0)
sx = Sampling(767,0.012,0.0)
st = Sampling(2650,0.0015,0.0)
#st = Sampling(4001,0.0015,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta

pngDir = None
#pngDir = '/Users/sluo/Desktop/png/'
pngDir = '/Users/sluo/Dropbox/png/'

plotVelocity = True # convert slowness to velocity

widthPoints = None # slides
#widthPoints = 175.0 # 1/3 column
#widthPoints = 240.0 # 1 column
#widthPoints = 260.0 # 1/2 page
#widthPoints = 504.0 # 2 column (full page)

#############################################################################

def main(args):
  #plotFwi()
  #plotLsm()
  #plotFwiData()
  #plotLsmData() # Born data
  #plotObjectiveFunction()
  #plotFwiObjectiveFunctions()
  plotLsmObjectiveFunctions()
  #plotModelMisfits()
  #plotFiles()

  #plotDataFromFile()

def plotModelMisfits():
  nm = 20 # number of iterations run
  slides = True if widthPoints is None else False

  lsmDir = '/Users/sluo/Dropbox/save/lsm/'
  ddir1 = lsmDir+'marmousi/100p/dres2/'
  ddir2 = lsmDir+'marmousi/95p/dres3/'
  ddir3 = lsmDir+'marmousi/95p/ares5/'
  #ddir1 = lsmDir+'marmousi/100p/dres2/'
  #ddir2 = lsmDir+'marmousi/random/10p/plus/dres/'
  #ddir3 = lsmDir+'marmousi/random/10p/plus/ares/'
  sname = 's1'

  #fwiDir = '/Users/sluo/Dropbox/save/fwi/'
  #ddir1 = fwiDir+'marmousi/2000m_100p/dres2/'
  #ddir2 = fwiDir+'marmousi/2000m_100p/cres4/'
  #ddir3 = None
  #sname = 's'

  def computeMisfit(ffile,t1):
    s1 = read(ffile)
    sub(t1,s1,s1)
    return sum(mul(s1,s1))
  mres1 = zerofloat(nm+1)
  mres2 = zerofloat(nm+1)
  mres3 = zerofloat(nm+1)
  t1 = read(ddir1+sname+'_true.dat')
  tt = read(ddir1+sname+'_init.dat') if sname=='s' else like(t1)
  mres1[0] = mres2[0] = mres3[0] = sum(mul(sub(t1,tt),sub(t1,tt)))
  for im in range(nm):
    mres1[im+1] = computeMisfit(ddir1+sname+'_'+str(im)+'.dat',t1)
    mres2[im+1] = computeMisfit(ddir2+sname+'_'+str(im)+'.dat',t1)
    if ddir3:
      mres3[im+1] = computeMisfit(ddir3+sname+'_'+str(im)+'.dat',t1)

  def normal(x):
    div(x,x[0],x)
  normal(mres1)
  normal(mres2)
  normal(mres3)

  panel = PlotPanel()
  panel.setHLabel('Iteration')
  if slides:
    panel.setVLimits(0.0,1.7)
  else:
    panel.setVLabel('Normalized misfit')
    panel.setVInterval(1.0)
    panel.setVLimits(0.0,1.3)
    #panel.setVLimits(0.0,2.0)

  # Reverse order for colors and styles
  arrays = [mres1,mres2,mres3]
  if slides:
    if ddir3:
      colors = [Color.BLACK,Color.BLUE,Color.RED]
    else:
      colors = [Color.BLACK,Color.RED,Color.RED]
    styles = [PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE,
              PointsView.Mark.FILLED_CIRCLE]
  else:
    colors = [Color.BLACK,Color.BLACK,Color.BLACK]
    styles = [PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.HOLLOW_CIRCLE,
              PointsView.Mark.FILLED_SQUARE]
  for i in range(3 if ddir3 else 2):
    p = panel.addPoints(arrays[i])
    p.setLineColor(colors[i])
    p.setMarkColor(colors[i])
    if slides:
      p.setLineWidth(1.5)
    else:
      p.setLineWidth(3.0)
    p.setMarkStyle(styles[i])

  frame = PlotFrame(panel)
  if slides:
    frame.setFontSizeForSlide(1.0,1.0,16.0/9.0)
    frame.setSize(1000,500)
  else:
    frame.setFontSizeForPrint(8.0,widthPoints)
    frame.setSize(1000,600)
  frame.setVisible(True)
  if pngDir:
    frame.paintToPng(1080,3.5,pngDir+'mres.png')
  #plotTemplatesForLegend(slides,colors,styles)

def plotDataFromFile():
  lsmDir = '/Users/sluo/Dropbox/save/lsm/'
  #ddir,iiter = lsmDir+'marmousi/100p/dres2/','19'
  #ddir,iiter = lsmDir+'marmousi/95p/dres3/','19'
  ddir,iiter = lsmDir+'marmousi/95p/ares5/','19'
  #ddir,iiter = lsmDir+'marmousi/random/pos/dres/','9'
  #ddir,iiter = lsmDir+'marmousi/random/pos/ares/','9'
  do = zerofloat(4001,nx)
  read(ddir+'do.dat',do)
  cmin,cmax = -0.3*max(abs(do)),0.3*max(abs(do))
  plot(do,cmin=cmin,cmax=cmax,cbar='Amplitude',grid=True,title='do')

def plotFwiData():
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

def plotLsmData():
  lsmDir = '/Users/sluo/Dropbox/save/lsm/'
  ddir,iiter = lsmDir+'marmousi/100p/dres2/','19'
  #ddir,iiter = lsmDir+'marmousi/95p/dres3/','19'
  #ddir,iiter = lsmDir+'marmousi/95p/ares5/','19'
  #ddir,iiter = lsmDir+'marmousi/random/25p/plus/dres/','9'
  #ddir,iiter = lsmDir+'marmousi/random/25p/plus/ares/','9'
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
  dw,v,rt = like(do),like(do),like(do)
  warp(ds,do,dw,v,rt) # right order
  ra = sub(ds,dw)
  rc = sub(ds,dw)
  add(rt,rc,rc)
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
  grid,hint,vint = True,2.0,1.0
  #grid,hint,vint = False,None,None
  rd = sub(ds,do)
  cmin,cmax = -0.3*max(abs(do)),0.3*max(abs(do))
  rmin,rmax = cmin,cmax
  #rmin,rmax = -1.0*max(abs(rd)),1.0*max(abs(rd))
  #rmin,rmax = -1.0*max(abs(ra)),1.0*max(abs(ra))
  #rmin,rmax = cmin,cmax
  plot(do,cmin=cmin,cmax=cmax,cbar='Amplitude',hint=hint,vint=vint,
    grid=grid,title='do')
  plot(ds,cmin=cmin,cmax=cmax,cbar='Amplitude',hint=hint,vint=vint,
    grid=grid,title='ds')
  plot(dw,cmin=cmin,cmax=cmax,cbar='Amplitude',hint=hint,vint=vint,
    grid=grid,title='dw')
  plot(mul(dt,v),cmap=rwb,sperc=100.0,cbar='Traveltime shift (s)',
    hint=hint,vint=vint,cint=0.1,grid=grid,title='v')
  plot(ra,cmin=rmin,cmax=rmax,cbar='Amplitude',title='ra')
  #plot(rt,cmin=rmin,cmax=rmax,cbar='Amplitude',title='rt')
  #plot(rc,cmin=rmin,cmax=rmax,cbar='Amplitude',title='rc')
  plot(sub(ds,do),cmin=rmin,cmax=rmax,cbar='Amplitude',title='rd')

def warp(ds,do,dw=None,v=None,rt=None):
  nt,nr = len(ds[0]),len(ds)
  doRmsFilter = True # filter by local rms amplitudes
  doEgain = False # exponential gain from first arrivals
  doAgc = False # agc
  td = 5 # time decimation
  rd = 1 # receiver decimation
  maxShift = 0.5 # max shift in seconds
  sw = Stopwatch(); sw.start()
  doc = copy(do)
  if doRmsFilter:
    ds,do = rmsFilter(ds,do,sigma=0.5*maxShift)
  if doEgain:
    ds = egain(ds)
    do = egain(do)
  if doAgc:
    ds,do = agc(ds,do)
  ds,do = addRandomNoise(10.0,ds,do,sigma=1.0)
  ds = copy(nt/td,nr/rd,0,0,td,rd,ds)
  do = copy(nt/td,nr/rd,0,0,td,rd,do)
  strainMax1,strainMax2 = 1.00,0.50
  #strainMax1,strainMax2 = 0.50,0.20
  #strainMax1,strainMax2 = 0.25,0.10
  #strainMax1,strainMax2 = 0.10,0.05
  shiftMax = int(maxShift/(td*dt))
  warp = DynamicWarping(-shiftMax,shiftMax)
  warp.setErrorExponent(1.0)
  warp.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.AVERAGE)
  warp.setStrainMax(strainMax1,strainMax2)
  warp.setShiftSmoothing(32.0/td,8.0/rd) # shift smoothing
  warp.setErrorSmoothing(2) # number of smoothings of alignment errors
  u = warp.findShifts(ds,do)
  mul(td,u,u) # scale shifts to compensate for decimation
  li = LinearInterpolator()
  li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT)
  li.setUniform(nt/td,td*dt,0.0,nr/rd,rd*dx,0.0,u)
  if v is None:
    v = zerofloat(nt,nr)
  for ir in range(nr):
    z = ir*dz
    for it in range(nt):
      t = it*dt
      v[ir][it] = li.interpolate(t,z) # interpolate shifts
  if dw is None:
    dw = zerofloat(nt,nr)
  warp.applyShifts(v,doc,dw)
  if rt is not None:
    warp.applyShifts(v,timeDerivative(doc),rt)
    mul(v,rt,rt)
  #report('warping',sw)
  #plot(ds,title='f') # used for dynamic warping
  #plot(do,title='g') # used for dynamic warping
  #plot(dw,title='h') # warped
  #plot(v,cmap=rwb,sperc=100.0,title='shifts')
  return dw,v
def rmsFilter(ds,do,sigma=0.1):
  x,y = copy(ds),copy(do)
  rmsx = rms(x)
  rmsy = rms(y)
  # equalize rms
  if rmsx>rmsy:
    mul(rms(y)/rms(x),x,x)
  else:
    mul(rms(x)/rms(y),y,y)
  xx = mul(x,x)
  yy = mul(y,y)
  rgf = RecursiveGaussianFilter(sigma/dt)
  rgf.apply00(xx,xx)
  rgf.apply00(yy,yy)
  num = mul(mul(2.0,xx),yy)
  den = add(mul(xx,xx),mul(yy,yy))
  add(1.0e-6,den,den)
  div(num,den,den)
  mul(den,x,x)
  mul(den,y,y)
  #plot(den,cmap=jet,title='rms_weights')
  return x,y
def rms(x):
  return sqrt(sum(mul(x,x))/len(x[0])/len(x))
def timeDerivative(f):
  """Derivative in time, assumed to be the 1st dimension."""
  n = len(f)
  #odt = 0.5/dt
  odt = 0.5
  g = like(f)
  for i in range(n):
    for it in range(1,nt-1):
      g[i][it] = odt*(f[i][it+1]-f[i][it-1])
  return g
def addRandomNoise(snr,f,g,sigma=1.0):
  n1,n2 = len(f[0]),len(f)
  frms = sqrt(sum(mul(f,f))/n1/n2)
  grms = sqrt(sum(mul(g,g))/n1/n2)
  xrms = 0.5*(frms+grms)  # (average) rms of signal
  random = Random(314159)
  s = sub(randfloat(random,n1,n2),0.5)
  RecursiveGaussianFilter(sigma).apply00(s,s) # bandlimited noise
  srms = sqrt(sum(mul(s,s))/n1/n2) # rms of noise
  mul(xrms/(srms*snr),s,s)
  return add(f,s),add(g,s)

def plotLsm():
  lsmDir = '/Users/sluo/Dropbox/save/lsm/'
  ddir,iiter = lsmDir+'marmousi/100p/dres2/','19'
  #ddir,iiter = lsmDir+'marmousi/100p/ares/','19'
  #ddir,iiter = lsmDir+'marmousi/95p/dres3/','19'
  #ddir,iiter = lsmDir+'marmousi/95p/ares5/','19'
  #ddir,iiter = lsmDir+'marmousi/random/10p/plus/dres/','19'
  #ddir,iiter = lsmDir+'marmousi/random/10p/plus/ares/','19'
  t = read(ddir+'s_true.dat')
  t0 = read(ddir+'s0_true.dat')
  t1 = read(ddir+'s1_true.dat')
  s0 = read(ddir+'s0_init.dat')
  s1 = read(ddir+'s1_'+iiter+'.dat')
  if plotVelocity:
    div(1.0,t,t)
    div(1.0,t0,t0)
    div(1.0,s0,s0)
  g = read(ddir+'g_'+iiter+'.dat'); mul(g,1.0/max(abs(g)),g)
  p = read(ddir+'p_'+iiter+'.dat'); mul(p,1.0/max(abs(p)),p)
  #mul(s1,0.10/max(abs(s1)),s1) # XXX
  #mul(s1,0.15/max(abs(s1)),s1) # XXX
  clipMin0,clipMax0 = min(t),max(t)
  #clipMin0,clipMax0 = min(t0),max(t0)
  #clipMin0,clipMax0 = 1.5,max(t0)
  cbar0 = 'Velocity (km/s)' if plotVelocity else 'Slowness (s/km)'
  #cmap1,cint,clip1 = rwb,0.1,max(abs(t1))
  cmap1,cint,clip1 = gray,None,0.05
  #cmap1,cint,clip1 = gray,0.05,0.05
  #cmapDiff,clipDiff = rwb,0.0
  #cmapDiff,clipDiff = rwb,0.04
  cmapDiff,clipDiff = jet,0.04
  plot(t,cmap=jet,cbar=cbar0,title='t')
  plot(t0,cmap=jet,cmin=clipMin0,cmax=clipMax0,cbar=cbar0,title='t0')
  plot(t1,cmap=cmap1,cmin=-clip1,cmax=clip1,
    cint=cint,cbar='Reflectivity',title='t1')
  #plot(sub(t0,s0),cmap=cmapDiff,sperc=100.0,
  plot(sub(s0,t0),cmap=cmapDiff,sperc=100.0,
    cbar='Velocity (km/s)' if plotVelocity else 'Slowness (s/km)',
    title='t0-s0')
  plot(s0,cmap=jet,cmin=clipMin0,cmax=clipMax0,cbar=cbar0,title='s0')
  plot(s1,cmap=cmap1,cmin=-clip1,cmax=clip1,cbar='Reflectivity',
    cint=cint,title='s1')
  #plot(g,cmap=rwb,cmin=-1.0,cmax=1.0,title='g')
  #plot(p,cmap=rwb,cmin=-1.0,cmax=1.0,title='p')

  zmin,zmax,xmin,xmax,width,height = 1.00,2.02,4.5,6.5,1024,494
  #zmin,zmax,xmin,xmax,width,height = 1.00,2.00,4.0,6.0,1024,494
  #zmin,zmax,xmin,xmax,width,height = 0.25,1.25,3.5,5.5,1024,494
  plotSubset(s1,zmin,zmax,xmin,xmax,width,height,-clip1,clip1,'s1_sub')
  plotSubset(t1,zmin,zmax,xmin,xmax,width,height,-clip1,clip1,'t1_sub')

def plotSubset(f,zmin,zmax,xmin,xmax,width,height,cmin,cmax,title):
  widthPoints = 240.0
  j1,n1 = int(zmin/dz),int((zmax-zmin)/dz)
  j2,n2 = int(xmin/dx),int((xmax-xmin)/dx)
  s1 = Sampling(n1,dz,zmin)
  s2 = Sampling(n2,dx,xmin)
  g = copy(n1,n2,j1,j2,f)
  panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  panel.setHLabel('Distance (km)')
  panel.setVLabel('Depth (km)')
  panel.setHInterval(0.5)
  panel.setVInterval(0.5)
  cb = panel.addColorBar('Reflectivity')
  cb.setWidthMinimum(150)
  cb.setInterval(0.05)
  pixel = panel.addPixels(s1,s2,g)
  pixel.setColorModel(gray)
  pixel.setClips(cmin,cmax)
  pixel.setInterpolation(PixelsView.Interpolation.LINEAR)
  frame = PlotFrame(panel)
  frame.setFontSizeForPrint(8.0,widthPoints)
  #frame.setSize(1024,670)
  frame.setSize(width,height)
  frame.setTitle(title)
  frame.setVisible(True)
  if title and pngDir:
    frame.paintToPng(720.0,widthPoints/72.0,pngDir+title+'.png')

def plotFwi():
  iiter = '4'
  ddir = '/Users/sluo/Dropbox/save/fwi/marmousi/2000m_100p/cres4/'
  #ddir = '/Users/sluo/Dropbox/save/fwi/marmousi/2000m_100p/dres2/'
  #ddir = '/Users/sluo/Dropbox/save/fwi/marmousi/200m_100p/dres/'
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
  frame.setFontSizeForSlide(1.0,1.0,16.0/9.0)
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
  frame.setFontSizeForSlide(1.0,1.0,16.0/9.0)
  frame.setSize(1000,500)
  frame.setVisible(True)
  if pngDir:
    frame.paintToPng(1024,2.0,pngDir+'obj.png')
  plotTemplatesForLegend()

def plotLsmObjectiveFunctions():

  slides = True if widthPoints is None else False

  alsm_ra95p = zerofloat(21)
  alsm_rd95p = zerofloat(21)
  lsm_rd95p = zerofloat(21)
  lsm_rd100p = zerofloat(21)
  ddir = '/Users/sluo/Dropbox/save/lsm/marmousi/'
  #read(ddir+'95p/ares5/ares.dat',alsm_ra95p)
  #read(ddir+'95p/ares5/dres.dat',alsm_rd95p)
  #read(ddir+'95p/dres3/dres.dat',lsm_rd95p)
  #read(ddir+'100p/dres2/dres.dat',lsm_rd100p)
  read(ddir+'random/10p/plus/ares/ares.dat',alsm_ra95p)
  read(ddir+'random/10p/plus/ares/dres.dat',alsm_rd95p)
  read(ddir+'random/10p/plus/dres/dres.dat',lsm_rd95p)
  read(ddir+'100p/dres2/dres.dat',lsm_rd100p)
  
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
  if not slides:
    panel.setVLabel('Normalized misfit')
    panel.setVInterval(1.0)
  if slides:
    #panel.setVLimits(0.0,1.7)
    panel.setVLimits(0.0,1.1)
  else:
    panel.setVLimits(0.0,1.3)
    #panel.setVLimits(0.0,2.0)

  # Reverse order for colors and styles
  #arrays = [alsm_ra95p,alsm_rd95p,lsm_rd95p,lsm_rd100p]
  arrays = [alsm_ra95p,alsm_rd95p,lsm_rd95p]
  #arrays = [lsm_rd95p,lsm_rd100p]
  if slides:
    #colors = [Color.RED,Color.RED,Color.BLUE,Color.BLACK]
    #styles = [PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.HOLLOW_CIRCLE,
    #          PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE]
    colors = [Color.RED,Color.RED,Color.BLUE]
    styles = [PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.HOLLOW_CIRCLE,
              PointsView.Mark.FILLED_CIRCLE]
    #colors = [Color.BLUE,Color.BLACK]
    #styles = [PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE]
  else:
    colors = [Color.BLACK,Color.BLACK,Color.BLACK,Color.BLACK]
    styles = [PointsView.Mark.FILLED_SQUARE,PointsView.Mark.CROSS,
              PointsView.Mark.HOLLOW_CIRCLE,PointsView.Mark.FILLED_CIRCLE]
  for i in range(len(colors)):
    p = panel.addPoints(arrays[i])
    p.setLineColor(colors[i])
    p.setMarkColor(colors[i])
    if slides:
      p.setLineWidth(1.5)
    else:
      p.setLineWidth(3.0)
    p.setMarkStyle(styles[i])

#  pd = panel.addPoints(alsm_ra95p)
#  pd.setLineColor(Color.RED)
#  pd.setMarkColor(Color.RED)
#  pd.setLineWidth(1.5)
#  pd.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
#
#  pc = panel.addPoints(alsm_rd95p)
#  pc.setLineColor(Color.RED)
#  pc.setMarkColor(Color.RED)
#  pc.setLineWidth(1.5)
#  pc.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
#
#  pb = panel.addPoints(lsm_rd95p)
#  pb.setLineColor(Color.BLUE)
#  pb.setMarkColor(Color.BLUE)
#  pb.setLineWidth(1.5)
#  pb.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
#
#  pa = panel.addPoints(lsm_rd100p)
#  pa.setLineColor(Color.BLACK)
#  pa.setMarkColor(Color.BLACK)
#  pa.setLineWidth(1.5)
#  pa.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)

  frame = PlotFrame(panel)
  if slides:
    frame.setFontSizeForSlide(1.0,1.0,16.0/9.0)
    #frame.setSize(1000,500)
    frame.setSize(1600,800)
  else:
    frame.setFontSizeForPrint(8.0,widthPoints)
    frame.setSize(1000,600)
  frame.setVisible(True)
  if pngDir:
    if slides:
      frame.paintToPng(1920,1.0,pngDir+'obj.png')
    else:
      frame.paintToPng(1080,3.5,pngDir+'obj.png')
  plotTemplatesForLegend(slides,colors,styles)

def plotTemplatesForLegend(slides,colors,styles):
  def template(color,markStyle,title):
    t = zerofloat(3)
    panel = PlotPanel()
    panel.setVLimits(-1.0,1.0)
    point = panel.addPoints(t)
    point.setLineColor(color)
    point.setMarkColor(color)
    if slides:
      #point.setLineWidth(1.5)
      pass
    else:
      point.setLineWidth(3.0)
    point.setMarkStyle(markStyle)
    frame = PlotFrame(panel)
    frame.setFontSizeForSlide(1.0,1.0,16.0/9.0)
    frame.setSize(1000,500)
    frame.setVisible(True)
    if pngDir:
      if slides:
        #frame.paintToPng(1024,2.0,pngDir+title+'.png')
        frame.paintToPng(1920,1.0,pngDir+title+'.png')
      else:
        frame.paintToPng(720.0,widthPoints/72.0,pngDir+title+'.png')
  #if slides:
  #  colors = [Color.BLACK,Color.BLUE,Color.RED,Color.RED]
  #  styles = [PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE,
  #            PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.HOLLOW_CIRCLE]
  #  titles = ['black_filled','blue_filled','red_filled','red_hollow']
  #else:
  #  pass
  titles = ['template1','template2','template3','template4']
  for i in range(len(colors)):
    template(colors[i],styles[i],titles[i])
  #template(Color.BLACK,PointsView.Mark.FILLED_CIRCLE,'black_filled')
  #template(Color.BLUE,PointsView.Mark.FILLED_CIRCLE,'blue_filled')
  #template(Color.RED,PointsView.Mark.FILLED_CIRCLE,'red_filled')
  #template(Color.RED,PointsView.Mark.HOLLOW_CIRCLE,'red_hollow')

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
  frame.setFontSizeForSlide(1.5,1.5,16.0/9.0)
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
def plot(f,cmap=gray,cmin=0,cmax=0,perc=100,sperc=None,cbar=None,cwidth=None,
  cint=None,hint=None,vint=None,grid=False,title=None):
  panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  cb = panel.addColorBar()
  if cbar:
    cb.setLabel(cbar)
  if cint:
    cb.setInterval(cint)
  if cwidth:
    cb.setWidthMinimum(cwidth)
  elif widthPoints is None:
    #cb.setWidthMinimum(140)
    cb.setWidthMinimum(190)
  elif widthPoints==175.0:
    cb.setWidthMinimum(200)
  elif widthPoints==240.0:
    #cb.setWidthMinimum(150)
    cb.setWidthMinimum(160)
  elif widthPoints==260.0:
    cb.setWidthMinimum(140)
  if hint:
    panel.setHInterval(hint)
  if vint:
    panel.setVInterval(vint)
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
  if perc<100:
    pixel.setPercentiles(100-perc,perc)
  if sperc is not None: # symmetric percentile clip (for plotting gradients)
    clips = Clips(100-sperc,sperc,f)
    clip = max(abs(clips.getClipMin()),abs(clips.getClipMax()))
    pixel.setClips(-clip,clip)
  if cmin<cmax:
    pixel.setClips(cmin,cmax)
  pixel.setInterpolation(PixelsView.Interpolation.LINEAR)
  if grid:
    gv = panel.addGrid()
    gv.setColor(Color.ORANGE)
    gv.setStyle(GridView.Style.DASH)
  frame = PlotFrame(panel)
  if widthPoints is None:
    frame.setFontSizeForSlide(1.0,1.0,16.0/9.0)
  else:
    frame.setFontSizeForPrint(8.0,widthPoints)
  if (len(f[0])==nz):
    if widthPoints is None:
      #frame.setSize(1200,770)
      #frame.setSize(1024,670)
      frame.setSize(1600,800)
    elif widthPoints==240.0 or widthPoints==175.0:
      #frame.setSize(1200,468)
      frame.setSize(1200,464)
    elif widthPoints==260.0:
      frame.setSize(1200,465)
    elif widthPoints==504.0:
      frame.setSize(1200,445)
  else:
    if widthPoints is None:
      #frame.setSize(1200,770)
      frame.setSize(1024,670)
    else:
      #frame.setSize(1200,1000)
      frame.setSize(1200,800)
  if title:
    frame.setTitle(title)
  frame.setVisible(True)
  if title and pngDir:
    if widthPoints is None:
      frame.paintToPng(1920,1.0,pngDir+title+'.png')
    else:
      frame.paintToPng(720.0,widthPoints/72.0,pngDir+title+'.png')
      #frame.paintToPng(3000.0,widthPoints/72.0,pngDir+title+'.png') # poster
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
