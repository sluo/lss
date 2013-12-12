#############################################################################
# Least-squares migration

from imports import *
from warp import *

#############################################################################

pngdatDir = None
pngdatDir = os.getenv('HOME')+'/Desktop/pngdat/'
#pngdatDir = os.getenv('HOME')+'/Desktop/pngdat2/'
#pngdatDir = os.getenv('HOME')+'/Desktop/pngdat3/'

gfile = None
pfile = None
sfile = None
#gfile = '/home/sluo/Desktop/pngdat/g_2.dat'
#pfile = '/home/sluo/Desktop/pngdat/p_2.dat'
#sfile = '/home/sluo/Desktop/pngdat/s1_2.dat'
#gfile = '/home/sluo/Desktop/save/lsm/marmousi/95p/ares2/g_4.dat'
#pfile = '/home/sluo/Desktop/save/lsm/marmousi/95p/ares2/p_4.dat'
#sfile = '/home/sluo/Desktop/save/lsm/marmousi/95p/ares2/s1_4.dat'

def main(args):
  #setupForLayered()
  setupForMarmousi()

  #initialize()
  #compareData()
  #showData()
  #WaveformInversion()
  AmplitudeInversion()
  #plotFiles()

  #plotWarpings()
  #plotDataResiduals()
  #plotLastResidual()

def setupForMarmousi():
  global sz,sx,st,nz,nx,nt,dz,dx,dt
  global kzs,kxs,kzr,kxr,ns,nr
  global fpeak,psou,niter
  global tt,t0,t1,s0,s1,u,a
  sz = Sampling(265,0.012,0.0)
  sx = Sampling(767,0.012,0.0)
  #st = Sampling(4001,0.0015,0.0)
  st = Sampling(5001,0.0012,0.0)
  nz,nx,nt = sz.count,sx.count,st.count
  dz,dx,dt = sz.delta,sx.delta,st.delta
  #kxs,kzs = [0],[0]
  #kxs,kzs = [nx/2],[0]
  #kxs,kzs = [nx-1],[0]
  #kxs,kzs = rampint(6,58,14),fillint(0,14)
  #kxs,kzs = rampint(3,33,24),fillint(0,24)
  #kxs,kzs = rampint(1,15,52),fillint(0,52)
  kxs,kzs = rampint(3,10,77),fillint(0,77)
  #kxs,kzs = rampint(3,5,153),fillint(0,153)
  kxr,kzr = rampint(0,1,nx),fillint(0,nx)
  ns,nr = len(kxs),len(kxr)
  #tt,t0,t1,s0,s1 = getMarmousi() # no error in background slowness s0
  tt,t0,t1,s0,s1 = getMarmousi(econst=-0.05) # constant error in s0
  #tt,t0,t1,s0,s1 = getMarmousi(egauss=0.25) # gaussian error in s0
  #tt,t0,t1,s0,s1 = getMarmousi(erand=0.25) # random error in s0
  #tt,t0,t1,s0,s1 = getMarmousi(erand=-0.25) # random error in s0
  #tt,t0,t1,s0,s1 = getMarmousi(erand=0.10) # random error in s0
  #tt,t0,t1,s0,s1 = getMarmousi(erand=-0.10) # random error in s0
  plots(tt,t0,t1,s0,s1)
  if sfile is not None:
    s1 = read(sfile)
  psou = min(14,ns)
  #psou = min(8,ns)
  fpeak = 10.0
  niter = 10
  sw = Stopwatch(); sw.start()
  u = zerofloat4(nz,nx,nt,psou)
  a = zerofloat4(nz,nx,nt,psou)
  report('allocate',sw)

def setupForLayered():
  global sz,sx,st,nz,nx,nt,dz,dx,dt
  global kzs,kxs,kzr,kxr,ns,nr
  global fpeak,psou,niter
  global t0,t1,s0,s1,u,a
  sz = Sampling(301,0.012,0.0)
  sx = Sampling(501,0.012,0.0)
  st = Sampling(2750,0.0015,0.0)
  nz,nx,nt = sz.count,sx.count,st.count
  dz,dx,dt = sz.delta,sx.delta,st.delta
  #kxs,kzs = [0],[0]
  #kxs,kzs = [nx/2],[0]
  kxs,kzs = rampint(1,10,51),fillint(0,51)
  #kxs,kzs = rampint(1,5,101),fillint(0,101)
  kxr,kzr = rampint(0,1,nx),fillint(0,nx)
  ns,nr = len(kxs),len(kxr)
  #tt,t0,t1,s0,s1 = getLayered1(s0mul=0.95)
  #tt,t0,t1,s0,s1 = getLayered1(s0mul=1.0)
  #tt,t0,t1,s0,s1 = getLayered2(s0mul=0.95)
  tt,t0,t1,s0,s1 = getLayered2(s0mul=1.0)
  #tt,t0,t1,s0,s1 = getLayered2(s0mul=1.05)
  plots(tt,t0,t1,s0,s1)
  if sfile is not None:
    s1 = read(sfile)
  psou = min(18,ns)
  fpeak = 10.0
  niter = 1
  sw = Stopwatch(); sw.start()
  u = zerofloat4(nz,nx,nt,psou)
  a = zerofloat4(nz,nx,nt,psou)
  report('allocate',sw)

def plots(tt,t0,t1,s0,s1):
  s0min,s0max = min(t0),max(t0)
  #s0min,s0max = min(min(s0),min(t0)),max(max(s0),max(t0))
  s1min,s1max = -max(abs(t1)),max(abs(t1))
  plot(tt,cmap=jet,cmin=s0min,cmax=s0max,cbar='Slowness (s/km)',title='s_true')
  plot(t0,cmap=jet,cmin=s0min,cmax=s0max,cbar='Slowness (s/km)',
       title='s0_true')
  plot(t1,sperc=100,cbar='Reflectivity',title='s1_true')
  plot(s0,cmap=jet,cmin=s0min,cmax=s0max,cbar='Slowness (s/km)',
       title='s0_init')
  plot(s1,sperc=100,cbar='Reflectivity',title='s1_init')

def initialize():
  global _t,_t0,_t1,_s0,_s1,_u,_a
  #_t,_t0,_t1,_s0,_s1 = getGaussian1(gmul=0.5)
  #_t,_t0,_t1,_s0,_s1 = getLayered1(s0mul=0.95)
  #_t,_t0,_t1,_s0,_s1 = getLayered1(s0mul=1.0)
  #_t,_t0,_t1,_s0,_s1 = getLayered2(s0mul=0.95)
  #_t,_t0,_t1,_s0,_s1 = getLayered2(s0mul=1.0)
  #_t,_t0,_t1,_s0,_s1 = getLayered2(s0mul=1.05)
  #_t,_t0,_t1,_s0,_s1 = getMarmousi(sigma=0.1)
  _t,_t0,_t1,_s0,_s1 = getMarmousi(sigma=0.1,s0mul=0.95)
  sw = Stopwatch(); sw.start()
  _u = zerofloat4(nz,nx,nt,psou)
  _a = zerofloat4(nz,nx,nt,psou)
  report('allocate',sw)
  s0min,s0max = min(_t0),max(_t0)
  #s0min,s0max = min(min(_s0),min(_t0)),max(max(_s0),max(_t0))
  s1min,s1max = -max(abs(_t1)),max(abs(_t1))
  plot(_t,cmap=jet,cmin=s0min,cmax=s0max,cbar='Slowness (s/km)',title='s_true')
  plot(_t0,cmap=jet,cmin=s0min,cmax=s0max,cbar='Slowness (s/km)',
       title='s0_true')
  plot(_t1,sperc=100,cbar='Reflectivity',title='s1_true')
  #plot(_s,cmap=jet,cmin=s0min,cmax=s0max,cbar='Slowness (s/km)',
  #     title='s_init')
  plot(_s0,cmap=jet,cmin=s0min,cmax=s0max,cbar='Slowness (s/km)',
       title='s0_init')
  plot(_s1,sperc=100,cbar='Reflectivity',title='s1_init')
  if sfile is not None:
    _s1 = read(sfile)
    #plot(_s1,sperc=100,cbar='Reflectivity',title='s1_file')
    plot(_s1,cmap=rwb,cmin=s1min,cmax=s1max,
      cbar='Reflectivity',title='s1_file')

#############################################################################

def plotFiles():
  if gfile is not None:
    plot(read(gfile),sperc=100,title='gfile')
  if pfile is not None:
    plot(read(pfile),sperc=100,title='pfile')
  if sfile is not None:
    plot(read(sfile),cmap=jet,cbar='Slowness (s/km)',title='sfile')

def plotLastResidual():
  iiter = '19'
  sdir = '/home/sluo/Desktop/save/lsm/marmousi/95p/dres3/'
  tres = zerofloat(20)
  rres = zerofloat(21)
  read(sdir+'res.dat',tres)
  copy(tres,rres)
  s0 = read(sdir+'s0_true.dat')
  s1 = read(sdir+'s1_true.dat')
  do = modelBornData(s0,s1)
  ds = zerofloat(nt,nr,ns)
  s0 = read(sdir+'s0_init.dat')
  s1 = read(sdir+'s1_'+iiter+'.dat')
  modelBornData(s0,s1,d=ds)
  r = zerofloat(nt,nr,ns)
  sub(ds,do,r); t = 'dres'
  #makeWarpedResidualP(ds,do,ra=r); t = 'ares'
  rres[20] = sum(mul(r,r))
  points(rres,title=t)

def plotDataResiduals():
  nm = 20 # number of iterations run
  #sdir = '/home/sluo/Desktop/save/lsm/marmousi/100p/dres2/'; ampRes = False
  #sdir = '/home/sluo/Desktop/save/lsm/marmousi/95p/dres2/'; ampRes = False
  #sdir = '/home/sluo/Desktop/save/lsm/marmousi/95p/ares5/'; ampRes = True
  sdir = '/home/sluo/Desktop/save/lsm/marmousi/random/10p/plus/dres/'
  #ampRes = True
  ampRes = False

  s0 = read(sdir+'s0_true.dat')
  s1 = read(sdir+'s1_true.dat')
  #mask(s1,18) # mask water bottom
  do = modelBornData(s0,s1)
  ds = zerofloat(nt,nr,ns)
  dw = zerofloat(nt,nr,ns)
  ra = zerofloat(nt,nr,ns)
  rd = zerofloat(nt,nr,ns)
  v = zerofloat(nt,nr,ns)
  dres = zerofloat(nm+1) # data residual
  ares = zerofloat(nm+1) # amplitude residual
  dres[0] = ares[0] = sum(mul(do,do)) 

  for im in range(1,nm+1):
    print 'im=%d'%im
    s0 = read(sdir+'s0_init.dat')
    s1 = read(sdir+'s1_'+str(im-1)+'.dat')
    #mask(s1,18) # mask water bottom
    modelBornData(s0,s1,d=ds)
    sub(ds,do,rd)
    dres[im] = sum(mul(rd,rd))
    if ampRes:
      makeWarpedResidualP(ds,do,ra=ra,dw=dw,v=v)
      ares[im] = sum(mul(ra,ra))

  points(dres,title='dres')
  if ampRes:
    points(ares,title='ares')

def plotWarpings():
  computeData = False

  nm = 20 # number of models
  if computeData:
    print 'computing data...'
    pre = '/home/sluo/Desktop/save/lsm/marmousi/95p/ares4/'
    s0,s1 = read(pre+'s0_true.dat'),read(pre+'s1_true.dat')
    do = modelBornData(s0,s1,isou=ns/2)
    plot(do,title='do')
    class Loop(Parallel.LoopInt):
      def compute(self,im):
        s0,s1 = read(pre+'s0_init.dat'),read(pre+'s1_'+str(im)+'.dat')
        ds = modelBornData(s0,s1,isou=ns/2)
        plot(ds,title='ds_'+str(im))
    Parallel.loop(nm,Loop())

  datDir = os.getenv('HOME')+'/Desktop/pngdat3/'

  do = zerofloat(nt,nr)
  read(datDir+'do.dat',do)

  print 'warping...'
  v = zerofloat(nt,nr,nm)
  class Loop(Parallel.LoopInt):
    def compute(self,im):
      ds = zerofloat(nt,nr)
      read(datDir+'ds_'+str(im)+'.dat',ds)
      warp(do,ds,v=v[im]) # wrong order gives v as a function of do time
  Parallel.loop(nm,Loop())

  print 'differentiating...'
  thresh = 10.0 # don't compare shifts if shifts are less than thresh
  for im in range(1,nm):
    vi,vm = v[im],v[im-1]
    total = 0.0
    count = 1.0
    for ir in range(nr):
      for it in range(nt):
        vii = vi[ir][it]
        vmi = vm[ir][it]
        if abs(vii)>thresh and abs(vmi)>thresh:
          e = vii-vmi
          total += e*e
          count += 1.0
    plot(vi,sperc=100,title='v_'+str(im))
    print 'dv=%f'%(total/count)

def compareData():
  dob = modelBornData(t0,t1); dob = dob[ns/2]
  #GaussianTaper.apply2(0.8,s1,s1)
  dsb = modelBornData(s0,s1); dsb = dsb[ns/2]
  ra,v = like(dsb),like(dsb)
  makeWarpedResidual(dsb,dob,ra=ra,v=v) # warped residuals
  cmin = -0.5*max(abs(dob))
  cmax =  0.5*max(abs(dob))
  plot(dob,cmin=cmin,cmax=cmax,title='dob')
  plot(dsb,cmin=cmin,cmax=cmax,title='dsb')
  plot(ra,cmin=cmin,cmax=cmax,title='ra')
  plot(v,sperc=100,title='v')
  #dof = sub(modelData(_t),modelDirectArrival(_t)); plot(dof[ns/2],title='dof')

#############################################################################
# Data

def modelBornData(s0,s1,isou=None,d=None):
  if isou is None:
    if d is None:
      d = zerofloat(nt,nr,ns)
    def compute(isou,fsou,lsou):
      ui = u[isou-fsou]
      wave = Wavefield(sz,sx,st)
      wave.modelAcousticWavefield(
        Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou]),s0,ui)
      wave.modelAcousticData(
        Wavefield.WavefieldSource(dt,s1,ui),
        Wavefield.Receiver(kzr,kxr),
        s0,d[isou])
    ChunkLoopInt(compute)
  else:
    if d is None:
      d = zerofloat(nt,nr)
    ui = u[0]
    wave = Wavefield(sz,sx,st)
    wave.modelAcousticWavefield(
      Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou]),s0,ui)
    wave.modelAcousticData(
      Wavefield.WavefieldSource(dt,s1,ui),
      Wavefield.Receiver(kzr,kxr),
      s0,d)
  return d

def modelData(s,isou=None):
  if isou is None:
    sw = Stopwatch(); sw.start()
    d = zerofloat(nt,nr,ns)
    Parallel.loop(ns,DataP(s,d))
    report('data',sw)
  else:
    d = zerofloat(nt,nr)
    wave = Wavefield(sz,sx,st)
    source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
    receiver = Wavefield.Receiver(kzr,kxr)
    wave.modelAcousticData(source,receiver,s,d)
  return d

def modelDirectArrival(s,isou=None):
  removeDirectArrival = True
  removeWaterBottom = False
  if removeDirectArrival:
    if removeWaterBottom:
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
      for iz in range(zbot+2,nz):
        r = 1.0+0.1*(iz-zbot-3.0)/(nz-zbot-10.0)
        mul(r,t[iz],t[iz])
        ref = RecursiveExponentialFilter(sigma)
        ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
        ref.apply(t[iz],t[iz])
        sigma += 1.0
      t = transpose(t)
    else:
      t = fillfloat(s[0][0],nz,nx)
    #plot(t,cmap=jet)
    return modelData(t,isou)
  else:
    return zerofloat(nt,nr,ns) if isou is None else zerofloat(nt,nr)

class DataP(Parallel.LoopInt):
  def __init__(self,s,d):
    self.s = s # slowness
    self.d = d # output array
  def compute(self,isou):
    wave = Wavefield(sz,sx,st)
    source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
    receiver = Wavefield.Receiver(kzr,kxr)
    wave.modelAcousticData(source,receiver,self.s,self.d[isou])

def showData():
  sw = Stopwatch(); sw.restart()
  do = modelBornData(t0,t1,isou=ns/2)
  ds = modelBornData(s0,s1,isou=ns/2)
  sw.stop(); print 'data total: %.2fs'%sw.time()
  rd = sub(ds,do) # data residual
  ra,rt,rc,dw,s = like(do),like(do),like(do),like(do),like(do)
  makeWarpedResidual(ds,do,ra,rt,rc,dw,s) # warped residuals
  #mul(s,dt,s)
  dmin,dmax = 1.0*min(min(do),min(ds)),1.0*max(max(do),max(ds))
  #rmin,rmax = 1.0*min(min(ra),min(rc)),1.0*max(max(ra),max(rc))
  rmin,rmax = dmin,dmax
  plot(do,cmin=dmin,cmax=dmax,cbar='Amplitude',title='observed')
  plot(ds,cmin=dmin,cmax=dmax,cbar='Amplitude',title='simulated')
  plot(dw,cmin=dmin,cmax=dmax,cbar='Amplitude',title='warped')
  plot(s,sperc=100.0,cbar='Traveltime shift (s)',title='shifts')
  plot(rd,cmin=rmin,cmax=rmax,cbar='Amplitude',title='data_residual')
  plot(ra,cmin=rmin,cmax=rmax,cbar='Amplitude',title='amplitude_residual')
  plot(rt,cmin=rmin,cmax=rmax,cbar='Amplitude',title='traveltime_residual')
  plot(rc,cmin=rmin,cmax=rmax,cbar='Amplitude',title='combined_residual')

#############################################################################
# Line Search

def updateModel(misfitFunction):
  print 'searching for step length...'
  #a,b = -1.0*max(abs(t1)),0.1*max(abs(t1)); tol = 0.10*abs(b-a)
  #a,b = -1.0*max(abs(t1)),0.1*max(abs(t1)); tol = 0.20*abs(b-a)
  #a,b = -1.0*max(abs(t1)),0.0*max(abs(t1)); tol = 0.20*abs(b-a)
  a,b = -0.5*max(abs(t1)),0.1*max(abs(t1)); tol = 0.20*abs(b-a)
  #a,b = -0.5*max(abs(t1)),0.0*max(abs(t1)); tol = 0.20*abs(b-a)
  sw = Stopwatch(); sw.restart()
  step = BrentMinFinder(misfitFunction).findMin(a,b,tol)
  print 'a =',a
  #print 'b =',b
  #print 'tol =',tol
  print 'step =',step
  report('line search',sw)
  add(mul(step,misfitFunction.p),s1,s1)

class MisfitFunction(BrentMinFinder.Function):
  def __init__(self,p,isou,do):
    self.p = div(p,max(abs(p))) # normalized conjugate ascent direction
    self.isou = isou
    self.do = do[isou]
    #plot(self.do,title='do')
  def evaluate(self,a):
    print 'evaluating'
    s1p = add(s1,mul(a,self.p))
    ds = modelBornData(s0,s1p,self.isou)
    #plot(ds,title='a='+str(a))
    r = self.residual(ds,self.do)
    return sum(mul(r,r))

class WaveformMisfitFunction(MisfitFunction):
  def residual(self,ds,do,rd=None):
    if rd is None:
      rd = like(ds)
    sub(ds,do,rd)
    return rd

class AmplitudeMisfitFunction(MisfitFunction):
  def residual(self,ds,do,ra=None):
    if ra is None:
      ra = like(ds)
    makeWarpedResidual(ds,do,ra=ra)
    return ra

#############################################################################
# Inversion

def processGradient(g):
  if nx==767:
    maskWater = True
  else:
    maskWater = False
  exponentialFilter = True # highpass using ref
  laplacianFilter = False
  highPassFilter = False
  h = copy(g)

  # Cheat
  stack = False
  mute = False
  if stack:
    for ix in range(1,nx):
      add(g[ix],g[0],g[0])
    for ix in range(1,nx):
      copy(g[0],g[ix])
    div(g,nx,g)
  if mute:
    for ix in range(nx):
      for iz in range(2*nz/3-40):
        g[ix][iz] = 0.0
      for iz in range(2*nz/3+40,nz):
        g[ix][iz] = 0.0

  if exponentialFilter:
    sigma = 0.100
    ref = RecursiveExponentialFilter(sigma/dx)
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
    ref.apply(g,h)
    #RecursiveGaussianFilter(sigma/dx).apply00(g,h)
    sub(g,h,h)
  if laplacianFilter:
    laplacian(h,h)
  if highPassFilter:
    bp = BandPassFilter(0.02,0.50,0.01,0.99)
    bp.apply(h,h)
  if maskWater:
    #mask(h,17) # mask water layer
    mask(h,15) # mask water layer
  #plot(g,sperc=100,title='before')
  #plot(h,sperc=100,title='after')
  copy(h,g)

class Inversion():
  """Abstract inversion class."""
  def __init__(self):
    self.do = modelBornData(t0,t1) # observed data
    self.ds = zerofloat(nt,nr,ns) # simulated data
    self.r = zerofloat(nt,nr,ns) # residual
    plot(self.do[ns/2],title='do')
    self.invert()
  def invert(self):
    if gfile is not None and pfile is not None:
      gm = read(gfile)
      pm = read(pfile)
      giter = int(gfile[-5:-4])
      if sfile is None or int(sfile[-5:-4])<giter:
        updateModel(self.getMisfitFunction(pm,ns/2,self.do))
      self.plots(gm,pm,s1,giter)
      fiter = giter+1
    else:
      gm = None
      pm = None
      fiter = 0
    res = zerofloat(niter)
    for iter in range(fiter,niter+fiter):
      print '\niteration',iter
      sw = Stopwatch(); sw.start()
      g = self.computeGradient(iter)
      p = conjugateDirection(g,gm,pm)
      mul(self.r,self.r,self.r); res[iter-fiter] = sum(self.r)
      if niter>1:
        updateModel(self.getMisfitFunction(p,ns/2,self.do)) # line search
      gm,pm = g,p
      self.plots(g,p,s1,iter)
      report('iteration',sw)
    points(res,title='res')
  def computeGradient(self,iter):
    def compute(isou,fsou,lsou):
      ui,ai = u[isou-fsou],a[isou-fsou]
      do,ds,r = self.do[isou],self.ds[isou],self.r[isou]
      wave = Wavefield(sz,sx,st)
      wave.modelAcousticWavefield(
        Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou]),s0,ui)
      if iter>0:
        wave.modelAcousticData(
          Wavefield.WavefieldSource(dt,s1,ui),
          Wavefield.Receiver(kzr,kxr),
          s0,ds)
      self.residual(ds,do,r)
      #GaussianTaper.apply2(r,r) # taper
      wave.modelAcousticWavefield(
        Wavefield.AdjointSource(dt,kzr,kxr,r),s0,ai)
      return makeGradient(ui,ai,d2=False)
    sw = Stopwatch(); sw.start()
    g = zerofloat(nz,nx)
    ChunkReduceInt(compute,g)
    processGradient(g)
    report('gradient',sw)
    return g
  def plots(self,g,p,m,iter=None):
    post = 'file' if iter is None else str(iter)
    plot(g,sperc=100,title='g_'+post)
    plot(p,sperc=100,title='p_'+post)
    plot(m,sperc=100,cbar='Reflectivity',title='s1_'+post)
    plot(self.ds[ns/2],cmin=min(self.do[ns/2]),cmax=max(self.do[ns/2]),
      title='ds_'+post)
  def getMisfitFunction(g,isou,do):
    pass
  def residual(ds,do,r):
    pass

def ChunkLoopInt(computeMethod,nsou=None):
  stopwatch = Stopwatch()
  if nsou is None:
    nsou = ns
  nchunk = (nsou+psou-1)/psou
  for ichunk in range(nchunk):
    stopwatch.restart()
    fsou = ichunk*psou # first source
    lsou = min(fsou+psou,nsou) # last source
    class Loop(Parallel.LoopInt):
      def compute(self,isou):
        computeMethod(isou,fsou,lsou)
    Parallel.loop(fsou,lsou,Loop())
    report('lsou=%d'%lsou,stopwatch)
def ChunkReduceInt(computeMethod,v):
  stopwatch = Stopwatch()
  nchunk = (ns+psou-1)/psou
  for ichunk in range(nchunk):
    stopwatch.restart()
    fsou = ichunk*psou # first source
    lsou = min(fsou+psou,ns) # last source
    class Reduce(Parallel.ReduceInt):
      def compute(self,isou):
        return computeMethod(isou,fsou,lsou)
      def combine(self,v1,v2):
        return add(v1,v2)
    add(Parallel.reduce(fsou,lsou,Reduce()),v,v)
    report('lsou=%d'%lsou,stopwatch)

class WaveformInversion(Inversion):
  def residual(self,ds,do,rd=None):
    if rd is None:
      rd = like(ds)
    sub(ds,do,rd)
    return rd
  def getMisfitFunction(self,g,isou,do):
    return WaveformMisfitFunction(g,isou,do)
    #return AmplitudeMisfitFunction(g,isou,do)

class AmplitudeInversion(Inversion):
  def residual(self,ds,do,ra=None):
    if ra is None:
      ra = like(ds)
    makeWarpedResidual(ds,do,ra=ra)
    return ra
  def getMisfitFunction(self,g,isou,do):
    return AmplitudeMisfitFunction(g,isou,do)

#############################################################################
# Dynamic Warping

def makeWarpedResidualP(ds,do,ra=None,rt=None,rc=None,dw=None,v=None):
  class Loop(Parallel.LoopInt):
    def compute(self,isou):
      doi = do[isou]
      dsi = ds[isou]
      rai = None if ra is None else ra[isou]
      rti = None if rt is None else rt[isou]
      rci = None if rc is None else rc[isou]
      dwi = None if dw is None else dw[isou]
      vi = None if v is None else v[isou]
      makeWarpedResidual(dsi,doi,rai,rti,rci,dwi,vi)
  Parallel.loop(len(ds),Loop())
#  nsou = len(ds)
#  reverseOrder = False
#  if v is None:
#    v = like(do)
#  if dw is None:
#    dw = like(do)
#  class Loop(Parallel.LoopInt):
#    def compute(self,isou):
#      doi = do[isou]
#      dsi = ds[isou]
#      dwi = dw[isou]
#      vi = v[isou]
#      rai = None if ra is None else ra[isou]
#      rti = None if rt is None else rt[isou]
#      rci = None if rc is None else rc[isou]
#      warp(dsi,doi,dwi,vi,rti)
#      if reverseOrder:
#        mul(-1.0,vi,vi)
#        if rai is not None:
#          sub(dwi,doi,rai)
#        if rti is not None:
#          mul(vi,timeDerivative(dsi),rti)
#      else:
#        if rai is not None:
#          sub(dsi,dwi,rai)
#  Parallel.loop(nsou,Loop())

def makeWarpedResidual(ds,do,ra=None,rt=None,rc=None,dw=None,v=None):
  reverseOrder = False
  if dw is None:
    dw = like(ds)
  if v is None:
    v = like(ds)
  if reverseOrder:
    warp(do,ds,dw,v) # wrong order
    mul(-1.0,v,v)
    if ra is not None:
      sub(dw,do,ra)
    if rt is not None:
      mul(v,timeDerivative(ds),rt)
    if rc is not None:
      sub(dw,do,rc)
      add(mul(v,timeDerivative(ds)),rc,rc)
  else:
    if (rt is None) and (rc is not None):
      rt = like(ds)
    warp(ds,do,dw,v,rt) # right order
    if ra is not None:
      sub(ds,dw,ra)
    if rc is not None:
      sub(ds,dw,rc)
      add(rt,rc,rc)
def xmakeWarpedResidual(ds,do,u=None,wd=None):
  reverseOrder = False
  if reverseOrder:
    dw,v = warp(do,ds) # wrong order
    mul(-1.0,v,v)
    rt = mul(v,timeDerivative(ds)) # traveltime residual
    ra = sub(dw,do) # warped residual
  else:
    rt = like(do)
    dw,v = warp(ds,do,rt=rt) # right order
    ra = sub(ds,dw) # warped residual
  if u is not None:
    copy(v,u)
  if wd is not None:
    copy(dw,wd)
  return rt,ra,add(rt,ra)

def xwarp(ds,do,dw=None,v=None,rt=None):
  """Smooth dynamic warping"""
  doRmsFilter = True # filter by local rms amplitudes
  doEgain = False # exponential gain from first arrivals
  maxShift = 0.5 # max shift in seconds
  sw = Stopwatch(); sw.start()
  doc = copy(do)
  if sum(ds)==0.0:
    return doc,like(doc)
  if doRmsFilter:
    ds,do = rmsFilter(ds,do,sigma=0.5*maxShift)
  if doEgain:
    ds,do = egain(ds),egain(do)
  ds,do = addRandomNoise(10.0,ds,do,sigma=1.0)

  r1max,r2max = 0.25,0.10 # strain limits (r1max in s/s, r2max in s/km)
  #d1min,d2min = 20*dt,20*dx # smoothness (d1min in seconds, d2min in km)
  d1min,d2min = 50*dt,50*dx # smoothness (d1min in seconds, d2min in km)
  warp = DynamicWarpingR(-maxShift,maxShift,st,sx)
  warp.setStrainLimits(-r1max,r1max,-r2max,r2max)
  warp.setSmoothness(d1min,d2min)
  if v is None:
    v = warp.findShifts(st,ds,st,do)
  else:
    copy(warp.findShifts(st,ds,st,do),v)
  if dw is None:
    dw = warp.applyShifts(st,doc,v)
  else:
    copy(warp.applyShifts(st,doc,v),dw)
  if rt is not None:
    rt = warp.applyShifts(st,timeDerivative(doc),v)
    mul(v,rt,rt)

  #report('warping',sw)
  #plot(ds,title='f') # used for dynamic warping
  #plot(do,title='g') # used for dynamic warping
  #plot(dw,title='h') # warped
  #plot(v,sperc=100.0,title='shifts')
  return dw,v

def warp(ds,do,dw=None,v=None,rt=None):
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
  #strainMax1,strainMax2 = 1.00,0.50
  strainMax1,strainMax2 = 0.50,0.20
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
  #plot(v,sperc=100.0,title='shifts')
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

def agc(do,ds):
  # AGC observed data
  ro = mul(do,do)
  ref = RecursiveExponentialFilter(50.0)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref.apply(ro,ro); div(ro,max(ro),ro); add(0.05,ro,ro)
  eo = div(do,ro)

  # AGC simulated data
  rs = mul(ds,ds)
  ref = RecursiveExponentialFilter(1000.0) # larger window
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref.apply(rs,rs); div(rs,max(rs),rs); add(0.05,rs,rs)
  es = div(ds,rs)

  # Divide by (scaled) rms
  """
  scale = 2.5
  print 'rms(eo) =',rms(eo)
  print 'rms(es) =',rms(es)
  div(eo,rms(eo),eo)
  div(es,rms(es)*scale,es)
  """

  #plot(ro,cmap=jet,title='ro')
  #plot(rs,cmap=jet,title='rs')
  #plot(do,title='do')
  #plot(ds,title='ds')
  #plot(eo,title='eo')
  #plot(es,title='es')
  return eo,es

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

def pickFirstArrivals(f):
  """First arrival picking.
  Reference: Wong, Han, Bancroft, and Stewart, 2009, Automatic
  time-picking of first arrivals on noisy microseismic data.
  NB: Equation 4 is wrong; use reciprocal.
  """
  smoothPicks = True
  te = 2.000 # energy window in seconds
  eps = 1.0E-8*max(abs(f)) # epsilon to stabilize division
  ne = int(te/dt)
  er = like(f) # modified energy ratio
  pp = zerofloat(nr) # picks (as floats)
  #print "picking..."; sw = Stopwatch(); sw.start()
  class Loop(Parallel.LoopInt):
    def compute(self,ir):
      # it=0
      num,den = f[ir][0]*f[ir][0],0.0
      for jt in range(0,min(ne+1,nt)):
        fi = f[ir][jt]
        den += fi*fi
      er[ir][0] = abs(f[ir][0])*den/(num+eps)
      # 0<it<nt
      for it in range(1,nt):
        fa = f[ir][it-ne-1] if it-ne-1>=0 else 0.0
        fb = f[ir][it-1]
        fc = f[ir][it]
        fd = f[ir][it+ne] if it+ne<nt else 0.0
        num += fc*fc-fa*fa
        den += fd*fd-fb*fb
        er[ir][it] = abs(f[ir][it])*den/(num+eps)
      # 0<it<nt
      i = zeroint(1)
      max(er[ir],i)
      pp[ir] = i[0]
  Parallel.loop(nr,Loop())
  #sw.stop(); print 'pick: %.2fs'%sw.time()
  #plot(er,title='energy_ratio')
  if smoothPicks:
    ref = RecursiveExponentialFilter(2.0)
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
    ref.apply(pp,pp)
  return pp

def egain(f,picks=None,sigma=None):
  """Exponential gain, starting from picked first arrivals."""
  square = False # square amplitudes
  normalize = True # normalize traces by rms amplitude
  if picks is None:
    picks = pickFirstArrivals(f)
  if sigma is None:
    #sigma = 0.25/fpeak
    sigma = 0.50/fpeak
  ref = RecursiveExponentialFilter(sigma/dt)
  ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_VALUE)
  e = zerofloat(nt,nr)
  for ir in range(nr):
    e[ir][int(picks[ir])] = 1.0
    ref.apply(e[ir],e[ir])
  mul(e,1.0/max(abs(e)),e)
  g = mul(e,f)
  if square:
    for ir in range(nr):
      mul(g[ir],g[ir],g[ir])
  if normalize:
    for ir in range(nr):
      r = 1.0/sqrt(sum(mul(g[ir],g[ir]))/nt)
      mul(g[ir],r,g[ir])
  #plot(e,title='exp_gain')
  return g

#############################################################################

def makeGradient(u,a,d2):
  g = correlate(u,a,d2) # gradient
  return g

def correlate(u,a,d2=False):
  """Zero-lag correlation."""
  class ReduceInt(Parallel.ReduceInt):
    def compute(self,it):
      if d2:
        t = add(add(u[it-1],u[it+1]),mul(-2.0,u[it])) # 2nd time derivative
        #return mul(-1.0/(dt*dt),mul(a[it],t))
        return mul(-1.0,mul(a[it],t))
      else:
        return mul(a[it],u[it])
    def combine(self,g1,g2):
      return add(g1,g2)
  return Parallel.reduce(1,nt-1,ReduceInt())

def conjugateDirection(g,gm=None,pm=None):
  """Polak-Ribiere nonlinear conjugate gradient method.
  Parameters:
    g - gradient ascent direction for current iteration.
    gm - gradient ascent direction for previous iteration.
    pm - conjugate ascent direction for previous iteration.
  Returns:
    the conjugate ascent direction for the current iteration.
  """
  if gm is None and pm is None:
    # For first iteration, use steepest descent direction.
    return g
  else:
    b = sum(mul(g,sub(g,gm)))/sum(mul(gm,gm))
    if b<0.0:
      b = 0.0
      print "  CG DIRECTION RESET"
    return add(g,mul(b,pm))


def maskWaterLayer(g,value=0.0):
  t0 = tt[0][0]
  for ix in range(nx):
    iz = 0
    while tt[ix][iz]==t0:
      g[ix][iz] = value
      iz += 1

def mask(g,l1=0):
  n1,n2 = len(g[0]),len(g)
  for i2 in range(n2):
    for i1 in range(l1):
      g[i2][i1] = 0.0

def getGaussian1(gmul=1.0):
  constantBackground = True
  tb = 0.5 # background slowness
  t = fillfloat(tb,nz,nx)
  for ix in range(nx):
    for iz in range(2*nz/3,nz):
      t[ix][iz] = 0.2
  t0,t1 = makeBornModel(t)
  GaussianTaper.apply2(t1,t1)
  if constantBackground:
    fill(tb,t0)

  g = zerofloat(nz,nx)
  g[nx/2][nz/3] = 1.0
  RecursiveGaussianFilter(nz/8.0).apply00(g,g)
  mul(g,gmul*tb/max(abs(g)),g)
  add(g,t0,t0)
  plot(g,cmap=jet)

  s0 = copy(t0)
  s1 = like(t1)
  return t,t0,t1,s0,s1

def getLayered1(s0mul=1.0):
  constantBackground = True
  tb = 0.5 # background slowness
  t = fillfloat(tb,nz,nx)
  for ix in range(nx):
    for iz in range(3*nz/5,nz):
      t[ix][iz] = 0.2
  t0,t1 = makeBornModel(t)
  GaussianTaper.apply2(t1,t1)
  if constantBackground:
    fill(tb,t0)
  s0 = copy(t0)
  mul(s0mul,s0,s0)
  s1 = like(t1)
  return t,t0,t1,s0,s1

def getLayered2(s0mul=1.0):
  constantBackground = True
  tb = 0.5 # background slowness
  t = fillfloat(tb,nz,nx)
  for ix in range(nx):
    for iz in range(nz/3,2*nz/3):
      t[ix][iz] = 0.38
    for iz in range(2*nz/3,nz):
      t[ix][iz] = 0.2
  t0,t1 = makeBornModel(t)
  GaussianTaper.apply2(t1,t1)
  if constantBackground:
    fill(tb,t0)
  s0 = copy(t0)
  mul(s0mul,s0,s0)
  s1 = like(t1)
  return t,t0,t1,s0,s1

import socket
def getMarmousi(sigma=0.1,econst=0.0,egauss=0.0,erand=0.0):
  """Marmousi model.
  Parameters:
    sigma - half-width of smoothing window for scale separation
    econst - constant error in background slowness s0
    egauss - gaussian error in background slowness s0
    erand - random error in background slowness s0
  Returns:
    tt - true model
    t0 - true background slowness
    t1 - true reflectivity
    s0 - true background slowness
    s1 - true reflectivity
  """
  p = zerofloat(751,2301)
  if socket.gethostname()=='backus.Mines.EDU':
    read("/data/sluo/marmousi/marmousi.dat",p)
  else:
    read("/data/seis/marmousi/marmousi.dat",p)
  p = copy(743,2301,8,0,p)
  div(1000.0,p,p) # slowness (s/km)
  t = fillfloat(2.0/3.0,nz,nx)
  copy(248,767,0,0,3,3,p,17,0,1,1,t)
  t0,t1 = makeBornModel(t,sigma); mask(t1,15)
  s0,s1 = mul(1.0+econst,t0),like(t1)
  if erand!=0.0:
    random = Random(0)
    r = sub(randfloat(random,nz,nx),0.5)
    RecursiveGaussianFilter(62.5).apply00(r,r)
    scale = erand*sum(s0)/(nz*nx) # 
    mul(scale/max(abs(r)),r,r)
    add(r,s0,s0)
    plot(mul(-1.0,r),sperc=100,title='s0_error')
  if egauss!=0.0:
    r = zerofloat(nz,nx);
    r[nx/2][nz/2] = 1.0
    RecursiveGaussianFilter(62.5).apply00(r,r)
    scale = egauss*sum(s0)/(nz*nx) # 
    mul(scale/max(abs(r)),r,r)
    add(r,s0,s0)
    plot(mul(-1.0,r),sperc=100,title='s0_error')
  return t,t0,t1,s0,s1

def linearRegression(x,y):
  """Fits a line to (x,y).
  Parameters:
    x - x-coordinates.
    y - y-coordinates.
  Returns:
    The sampled best fit line in an array the same size as x and y.
  """
  n = len(x)
  sumx = 0.0 # sum_n(x_i)
  sumy = 0.0 # sum_n(y_i)
  sumxx = 0.0 # sum_n(x_i*x_i)
  sumxy = 0.0 # sum_n(x_i*y_i)
  for i in range(n):
    xi = x[i]
    yi = y[i]
    sumx += xi
    sumy += yi
    sumxx += xi*xi
    sumxy += xi*yi
  beta = (sumxy-sumx*sumy/n)/(sumxx-sumx*sumx/n)
  alpha = (sumy-beta*sumx)/n
  z = zerofloat(n)
  for i in range(n):
    z[i] = alpha+beta*x[i]
  return z

def makeBornModel(s,sigma0=0.100,sigma1=None):
  """
  sigma0: smoothing for background model
  sigma1: smoothing for perturbation
  """
  s0,s1 = like(s),like(s)
  ref0 = RecursiveExponentialFilter(sigma0/dx)
  ref0.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref0.apply(s,s0)
  t = copy(s)
  if sigma1 is None:
    sigma1 = sigma0
  ref1 = RecursiveExponentialFilter(sigma1/dx)
  ref1.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
  ref1.apply(s,t)
  sub(s,t,s1)
  #r = mul(mul(2.0,s0),s1)
  #r = laplacian(mul(mul(2.0,s0),s1))
  #r = laplacian(sub(div(s1,s0),1.0))
  r = sub(mul(s,s),mul(t,t))
  return s0,r

def laplacian(x,y=None):
  n1,n2 = len(x[0]),len(x)
  y = like(x)
  class Loop(Parallel.LoopInt):
    def compute(self,i2):
      for i1 in range(1,n1-1):
        y[i2][i1] = (1.0/6.0)*(-20.0*x[i2][i1]+
          4.0*(x[i2-1][i1  ]+x[i2  ][i1+1]+x[i2  ][i1-1]+x[i2+1][i1  ])+
               x[i2+1][i1+1]+x[i2+1][i1-1]+x[i2-1][i1+1]+x[i2-1][i1-1]);
  Parallel.loop(1,n2-1,Loop())
  if y is None:
    return y
  else:
    copy(x,y)

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(f,cmap=gray,cmin=0,cmax=0,perc=100,sperc=None,cbar=None,title=None):
  panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  cb = panel.addColorBar()
  if cbar:
    cb.setLabel(cbar)
  cb.setWidthMinimum(160)
  #cb.setWidthMinimum(170)
  if len(f[0])==nz and len(f)==nx:
    pixel = panel.addPixels(sz,sx,f)
    panel.setHLabel('Distance (km)')
    panel.setVLabel('Depth (km)')
  elif len(f[0])==nt and len(f)==nr:
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
    pixel.setColorModel(rwb)
  pixel.setInterpolation(PixelsView.Interpolation.LINEAR)
  frame = PlotFrame(panel)
  frame.setFontSizeForSlide(1.5,1.5)
  #frame.setFontSizeForPrint(8.0,222.0)
  if (len(f[0])==nz):
    width = 1200
    height = int(1.25*width*nz/nx)
    frame.setSize(width,height)
  else:
    frame.setSize(1200,1000)
    #frame.setSize(1200,750)
  if title:
    frame.setTitle(title)
  frame.setVisible(True)
  if title and pngdatDir:
    frame.paintToPng(360,3.0,pngdatDir+title+'.png')
    #frame.paintToPng(720,3.08,pngdatDir+title+'.png')
    write(pngdatDir+title+'.dat',f)
  return panel

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
  if title and pngdatDir:
    frame.paintToPng(720,3.08,pngdatDir+title+'.png')
    write(pngdatDir+title+'.dat',f)

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

def checkForNaN(x):
  n1,n2 = len(x[0]),len(x)
  class Loop(Parallel.LoopInt):
    def compute(self,i2):
      for i1 in range(n1):
        if Float.isNaN(x[i2][i1]):
          raise RuntimeError('found NaN')
  Parallel.loop(n2,Loop())

def report(str,stopwatch):
  s = stopwatch.time()
  h = int(s/3600); s -= h*3600
  m = int(s/60); s -= m*60
  if h>0:
    print str+': %02d:%02d:%02d'%(h,m,s)
  else:
    print str+': %02d:%02d'%(m,s)

def smin(x1,x2=None,x3=None,x4=None):
  return minmax(ArrayMath.min,x1,x2,x3,x4)
def smax(x1,x2=None,x3=None,x4=None):
  return minmax(ArrayMath.max,x1,x2,x3,x4)
def minmax(method,x1,x2,x3,x4):
  if x2 is None:
    return method(x1)
  elif x3 is None:
    return method(method(x1),method(x2))
  elif x4 is None:
    return method(method(x1),method(x2),method(x3))
  else:
    return method(method(x1),method(x2),method(x3),method(x4))

#############################################################################
# Do everything on Swing thread.
import sys,time
class RunMain(Runnable):
  def run(self):
    start = time.time()
    if pngdatDir is not None:
      print 'cleaning '+pngdatDir.split('/')[-2]
      cleanDir(pngdatDir)
    print 'gfile = None' if gfile is None else 'gfile = '+gfile.split('/')[-1]
    print 'pfile = None' if pfile is None else 'pfile = '+pfile.split('/')[-1]
    print 'sfile = None' if sfile is None else 'sfile = '+sfile.split('/')[-1]
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02d'%(h,m,s)
if __name__=='__main__':
  SwingUtilities.invokeLater(RunMain())
