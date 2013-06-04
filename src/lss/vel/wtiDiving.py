"""
Combined waveform and traveltime inversion
"""
from imports import *

#############################################################################
sz = Sampling(265,0.012,0.0)
sx = Sampling(767,0.012,0.0)
#st = Sampling(2750,0.0015,0.0)
st = Sampling(4001,0.0015,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta

#kxs,kzs = [0],[0]
#kxs,kzs = [nx/2],[0]
#kxs,kzs = [nx-1],[0]
#kxs,kzs = rampint(8,50,16),fillint(0,16)
#kxs,kzs = rampint(8,30,26),fillint(0,26)
#kxs,kzs = rampint(3,20,39),fillint(0,39)
#kxs,kzs = rampint(1,15,52),fillint(0,52)
kxs,kzs = rampint(3,10,77),fillint(0,77)
kxr,kzr = rampint(0,1,nx),fillint(0,nx)
ns,nr = len(kxs),len(kxr)
fpeak = 5.0
niter = 5

pngdatDir = None
#pngdatDir = os.getenv('HOME')+'/Desktop/pngdat/'
#pngdatDir = os.getenv('HOME')+'/Desktop/pngdat2/'
pngdatDir = os.getenv('HOME')+'/Desktop/pngdat3/'

sfile = None
gfile = None
pfile = None
#sfile = '/home/sluo/Desktop/fwi/marmousi/2000m_100p/cres4/dat/s_9.dat'
#gfile = '/home/sluo/Desktop/fwi/marmousi/2000m_100p/cres4/dat/g_9.dat'
#pfile = '/home/sluo/Desktop/fwi/marmousi/2000m_100p/cres4/dat/p_9.dat'

psou = min(14,ns) # number of parallel sources
REMOVE_DIRECT_ARRIVAL = True
#############################################################################

def main(args):
  setModel()
  #showData()
  #goMigration()
  #WaveformInversion()
  #TraveltimeInversion()
  #CombinedInversion()
  plotDataResiduals()
  #plotLastResidual()

def setModel():
  global _t,_s
  #_t,_s = getGaussian()
  #_t,_s = getMarmousi(sigma=0.2,smul=1.0)
  #_t,_s = getMarmousi(sigma=0.5,smul=1.0)
  _t,_s = getMarmousi(sigma=2.0,smul=1.0)
  #_t,_s = getMarmousi(sigma=2.0,smul-0.95)
  #_t,_s = getMarmousi(None,1.0) # linear initial model
  #_t,_s = getMarmousi(None,0.95) # linear initial model
  #_t,_s = getMarmousi(None,0.90) # linear initial model
  g = sub(_s,_t); div(g,max(abs(g)),g)
  plot(_t,cmap=jet,cbar='Slowness (s/km)',title='s_true')
  plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
    title='s_init')
  plot(g,cmap=rwb,cmin=-0.95,cmax=0.95,title='g_true')
  if sfile is not None:
    _s = read(sfile)
    plot(_s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
      title='s_'+sfile[-5:-4])

def plotDataResiduals():
  niter = 20

  #sdir,cres = '/home/sluo/Desktop/save/fwi/marmousi/2000m_100p/dres/',False
  #sdir,cres = '/home/sluo/Desktop/save/fwi/marmousi/2000m_100p/dres2/',False
  #sdir,cres = '/home/sluo/Desktop/save/fwi/marmousi/2000m_100p/cres4a/',True
  sdir,cres = '/home/sluo/Desktop/save/fwi/marmousi/2000m_100p/cres6/',True

  do = zerofloat(nt,nr,ns)
  eo = zerofloat(nt,nr,ns)
  ds = zerofloat(nt,nr,ns)
  es = zerofloat(nt,nr,ns)
  dw = zerofloat(nt,nr,ns)
  rd = zerofloat(nt,nr,ns)
  rc = zerofloat(nt,nr,ns)
  v = zerofloat(nt,nr,ns)
  dres = zerofloat(niter+1) # data residual
  if cres:
    cres = zerofloat(niter+1) # combined residual

  t = read(sdir+'s_true.dat')
  modelData(t,d=do)
  modelDirectArrival(t,d=eo)
  sub(do,eo,do)

  dres[0] = cres[0] = sum(mul(do,do)) 
  for iiter in range(1,niter+1):
    print 'iiter=%d'%iiter
    s = read(sdir+'s_'+str(iiter-1)+'.dat')
    modelData(s,d=ds)
    modelDirectArrival(s,d=es)
    sub(ds,es,ds)
    sub(ds,do,rd)
    dres[iiter] = sum(mul(rd,rd))
    if cres:
      makeWarpedResidualP(ds,do,rc=rc,dw=dw,v=v)
      cres[iiter] = sum(mul(rc,rc))

  points(dres,title='dres')
  if cres:
    points(cres,title='cres')

def plotLastResidual():
  liter = 19 # last iteration number
  #sdir,cres = '/home/sluo/Desktop/save/fwi/marmousi/2000m_100p/dres2/',False
  sdir,cres = '/home/sluo/Desktop/save/fwi/marmousi/2000m_100p/cres6/',True

  tres = zerofloat(liter+1)
  rres = zerofloat(liter+2)
  read(sdir+'res.dat',tres)
  copy(tres,rres)
  do = zerofloat(nt,nr,ns)
  eo = zerofloat(nt,nr,ns)
  ds = zerofloat(nt,nr,ns)
  es = zerofloat(nt,nr,ns)
  dw = zerofloat(nt,nr,ns)
  rr = zerofloat(nt,nr,ns)
  v = zerofloat(nt,nr,ns)
  t = read(sdir+'s_true.dat')
  s = read(sdir+'s_'+str(liter)+'.dat')
  modelData(t,d=do)
  modelData(s,d=ds)
  modelDirectArrival(t,d=eo)
  modelDirectArrival(s,d=es)
  sub(do,eo,do)
  sub(ds,es,ds)
  if cres:
    makeWarpedResidualP(ds,do,rc=rr,dw=dw,v=v)
  else:
    sub(ds,do,rr)
  rres[liter+1] = sum(mul(rr,rr))
  points(rres,title='cres' if cres else 'dres')


#############################################################################
# Data

def modelData(s,isou=None,d=None):
  if isou is None:
    sw = Stopwatch(); sw.start()
    if d is None:
      d = zerofloat(nt,nr,ns)
    Parallel.loop(ns,DataP(s,d))
    report('data',sw)
  else:
    if d is None:
      d = zerofloat(nt,nr)
    wave = Wavefield(sz,sx,st)
    source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
    receiver = Wavefield.Receiver(kzr,kxr)
    wave.modelAcousticData(source,receiver,s,d)
  return d

def modelDirectArrival(s,isou=None,d=None):
  if REMOVE_DIRECT_ARRIVAL:
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
    #plot(t,cmap=jet)
    return modelData(t,isou,d)
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
  do = modelData(_t) # observed data
  ds = modelData(_s) # simulated data
  eo = modelDirectArrival(_t)  # direct arrival for observed data
  es = modelDirectArrival(_s)  # direct arrival for simulated data
  report('data total',sw)
  do = do[ns/2]
  ds = ds[ns/2]
  eo = eo[ns/2]
  es = es[ns/2]
  sub(do,eo,do) # observed
  sub(ds,es,ds) # simulated
  po = pickFirstArrivals(do)
  ps = pickFirstArrivals(ds)
  go = egain(do,po)
  gs = egain(ds,ps)
  ra = sub(ds,do) # amplitude residual
  s,dw,rt,rw,rc = like(do),like(do),like(do),like(do),like(do)
  makeWarpedResidual(ds,do,rw,rt,rc,dw,s) # warped residuals
  mul(s,dt,s)
  smax = 0.9*max(abs(s))
  rmin,rmax = 0.5*min(min(ra),min(rc)),0.5*max(max(ra),max(rc))
  dmin,dmax = 0.2*min(min(do),min(ds)),0.2*max(max(do),max(ds))
  gmin,gmax = 1.0*min(min(go),min(gs)),1.0*max(max(go),max(gs))
  #emin,emax = 0.5*min(min(eo),min(es)),0.5*max(max(eo),max(es))
  emin,emax = 2.0*dmin,2.0*dmax
  plot(eo,cmin=emin,cmax=emax,title='direct arrival (observed)')
  plot(es,cmin=emin,cmax=emax,title='direct arrival (simulated)')
  plot(do,sperc=99.6,cbar='Amplitude',title='observed')\
    .addPoints(mul(po,dt),rampfloat(0,dx,nr))
  plot(ds,sperc=99.6,cbar='Amplitude',title='simulated')\
    .addPoints(mul(ps,dt),rampfloat(0,dx,nr))
  plot(go,cmin=gmin,cmax=gmax,cbar='Amplitude',title='observed_gained')
  plot(gs,cmin=gmin,cmax=gmax,cbar='Amplitude',title='simulated_gained')
  plot(dw,sperc=99.6,cbar='Amplitude',title='warped')
  plot(s,cmap=rwb,cmin=-smax,cmax=smax,
    cbar='Traveltime shift (s)',title='shifts')
  plot(ra,cmin=rmin,cmax=rmax,cbar='Amplitude',title='residual')
  plot(rw,cmin=rmin,cmax=rmax,cbar='Amplitude',title='warped_residual')
  plot(rt,cmin=rmin,cmax=rmax,cbar='Amplitude',title='traveltime_residual')
  plot(rc,cmin=rmin,cmax=rmax,cbar='Amplitude',title='combined_residual')

#############################################################################
# Line Search

def updateModel(misfitFunction):
  print 'searching for step length...'
  #a,b  = -1.0*max(abs(sub(_t,_s))),0.0; tol = 0.25*abs(a)
  #a,b = -0.10,0.00; tol = 0.25*abs(b-a)
  #a,b = -0.10,0.02; tol = 0.20*abs(b-a)
  #a,b = -0.05,0.00; tol = 0.25*abs(b-a)
  a,b = -0.05,0.01; tol = 0.20*abs(b-a)
  sw = Stopwatch(); sw.restart()
  step = BrentMinFinder(misfitFunction).findMin(a,b,tol)
  print 'a =',a
  #print 'b =',b
  #print 'tol =',tol
  print 'step =',step
  report('line search',sw)
  add(mul(step,misfitFunction.p),_s,_s)

class MisfitFunction(BrentMinFinder.Function):
  """Abstract misfit function."""
  def __init__(self,p,isou,do):
    self.p = div(p,max(abs(p))) # normalized conjugate ascent direction
    self.isou = isou
    self.do = do[isou]
    self.wave = Wavefield(sz,sx,st)
  def evaluate(self,a):
    print 'evaluating'
    s = add(_s,mul(a,self.p))
    ds = modelData(s,self.isou)
    es = modelDirectArrival(s,self.isou)
    sub(ds,es,ds)
    r = self.residual(ds,self.do)
    res = sum(mul(r,r))
    #plot(r,title='residual='+str(res)+'_a='+str(a))
    return res
  def residual(self,ds,do):
    pass

class WaveformMisfitFunction(MisfitFunction):
  def residual(self,ds,do,rd=None):
    if rd is None:
      rd = like(ds)
    sub(ds,do,rd)
    return rd

class TraveltimeMisfitFunction(MisfitFunction):
  def residual(self,ds,do,rt=None):
    if rt is None:
      rt = like(ds)
    makeWarpedResidual(ds,do,rt=rt)
    return rt

class CombinedMisfitFunction(MisfitFunction):
  def residual(self,ds,do,rc=None):
    if rc is None:
      rc = like(ds)
    makeWarpedResidual(ds,do,rc=rc)
    return rc

#############################################################################
# Inversion

class Inversion():
  """Abstract inversion class."""
  def __init__(self):
    sw = Stopwatch(); sw.start()
    self.u = zeros(nz,nx,nt,psou) # source wavefields
    self.a = zeros(nz,nx,nt,psou) # receiver wavefields
    report('allocate',sw)
    self.do = sub(modelData(_t),modelDirectArrival(_t)) # observed data
    self.ds = zerofloat(nt,nr,ns) # simulated data
    self.r = zerofloat(nt,nr,ns) # residual
    self.invert()
  def invert(self):
    if gfile is not None and pfile is not None:
      gm = read(gfile)
      pm = read(pfile)
      giter = int(gfile[-5:-4])
      if sfile is None or int(sfile[-5:-4])<giter:
        updateModel(self.getMisfitFunction(pm,ns/2,self.do))
      self.plots(gm,pm,_s,giter)
      fiter = giter+1
    else:
      gm = None
      pm = None
      fiter = 0
    res = zerofloat(niter)
    for iter in range(fiter,niter+fiter):
      print '\niteration',iter
      sw = Stopwatch(); sw.start()
      g = self.computeGradient()
      p = conjugateDirection(g,gm,pm)
      mul(self.r,self.r,self.r); res[iter-fiter] = sum(self.r)
      if niter>0:
        updateModel(self.getMisfitFunction(p,ns/2,self.do))
        s = _s
      else:
        s = None
      self.plots(g,p,s,iter)
      gm,pm = g,p # rotate arrays
      report('iteration',sw)
    points(res,title='res')
  def computeGradient(self):
    do,ds,r,u,a,residual = self.do,self.ds,self.r,self.u,self.a,self.residual
    es = modelDirectArrival(_s) # direct arrival for simulated data
    g = zerofloat(nz,nx)
    #p = zerofloat(nz,nx) # preconditioner
    sw = Stopwatch(); sw.start()
    for chunkSource in range(0,(ns+psou-1)/psou):
      sw = Stopwatch(); sw.start()
      fsou = chunkSource*psou # first source
      lsou = min(fsou+psou,ns) # last source
      class Reduce(Parallel.ReduceInt):
        def compute(self,isou):
          wave = Wavefield(sz,sx,st)
          ui,ai = u[isou-fsou],a[isou-fsou]
          doi,dsi,esi,ri = do[isou],ds[isou],es[isou],r[isou]
          zero(ui) # necessary?
          zero(ai) # necessary?
          wave.modelAcousticDataAndWavefield( # forward-propagate
            Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou]),
            Wavefield.Receiver(kzr,kxr),
            _s,dsi,ui)
          sub(dsi,esi,dsi) # remove direct arrival
          checkForNaN(dsi) # throws an exception if NaN
          residual(dsi,doi,ri) # residual
          #GaussianTaper.apply(ri,ri) # taper
          wave.modelAcousticWavefield( # back-propagate
            Wavefield.AdjointSource(dt,kzr,kxr,ri),
            _s,ai)
          g,p = makeGradient(ui,ai)
          return g
        def combine(self,g1,g2):
          return add(g1,g2)
      add(Parallel.reduce(fsou,lsou,Reduce()),g,g) # stack
      report('lsou=%d'%lsou,sw)
    #div(g,p,g)
    maskWaterLayer(g)
    report('gradient',sw)
    return g
  def plots(self,g,p,s=None,iter=None):
    post = 'file' if iter is None else str(iter)
    gtitle = 'g_'+post
    ptitle = 'p_'+post
    stitle = 's_'+post
    plot(g,cmap=rwb,sperc=100.0,title=gtitle)
    plot(p,cmap=rwb,sperc=100.0,title=ptitle)
    if s is not None:
      plot(s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
        title=stitle)
  def getMisfitFunction(g,isou,do):
    pass
  def residual(ds,do):
    pass

class xInversion():
  """Abstract inversion class."""
  def __init__(self):
    self.ds = zerofloat(nt,nr)
    self.u = zerofloat(nz,nx,nt)
    self.a = zerofloat(nz,nx,nt)
    self.wave = Wavefield(sz,sx,st)
    self.do = modelData(_t) # observed data
    eo = modelDirectArrival(_t) # direct arrival
    sub(self.do,eo,self.do)
    self.invert()
  def invert(self):
    if gfile is not None and pfile is not None:
      print 'reading gfile'
      print 'reading pfile'
      gm = read(gfile)
      pm = read(pfile)
      giter = int(gfile[-5:-4])
      #gm = sub(_s,_t)
      #pm = sub(_s,_t)
      if sfile is None or int(sfile[-5:-4])<giter:
        updateModel(self.getMisfitFunction(pm,ns/2,self.do))
      self.plots(gm,pm,_s,giter)
      fiter = giter+1
    else:
      gm = None
      pm = None
      fiter = 0
    for iter in range(fiter,niter+fiter):
      print '\niteration',iter
      sw = Stopwatch(); sw.start()
      g = self.computeGradient()
      print 'sum(g) =',sum(g)
      p = conjugateDirection(g,gm,pm)
      print 'gradient: %.2fm'%(sw.time()/60.0)
      if niter>0:
        updateModel(self.getMisfitFunction(p,ns/2,self.do))
        s = _s
      else:
        s = None
      self.plots(g,p,s,iter)
      gm,pm = g,p # rotate arrays
      sw.stop(); print 'iteration: %.2fm'%(sw.time()/60.0)
  def computeGradient(self):
    g = zerofloat(nz,nx)
    p = zerofloat(nz,nx) # preconditioner
    es = modelDirectArrival(_s) # direct arrival for simulated data
    for isou in range(ns):
      num,den = self.computeGradientForOneSource(isou,es)
      add(num,g,g)
      add(den,p,p)
    #div(g,p,g)
    maskWaterLayer(g)
    #div(g,max(abs(g)),g)
    return g
  def computeGradientForOneSource(self,isou,es):
    sw0 = Stopwatch(); sw0.start()
    source = Wavefield.RickerSource(fpeak,kzs[isou],kxs[isou])
    receiver = Wavefield.Receiver(kzr,kxr)
    sw = Stopwatch(); sw.start()
    self.wave.modelAcousticDataAndWavefield(source,receiver,_s,self.ds,self.u)
    print 'forward: %.2fs'%sw.time()
    sub(self.ds,es[isou],self.ds) # subtract direct arrival
    checkForNaN(self.ds) # throws an exception if NaN
    r = self.residual(self.ds,self.do[isou]) # residual
    #GaussianTaper.apply(r,r) # taper residuals
    source = Wavefield.AdjointSource(dt,kzr,kxr,r)
    sw.restart()
    self.wave.modelAcousticWavefield(source,_s,self.a)
    print 'reverse: %.2fs'%sw.time()
    g,p = makeGradient(self.u,self.a)
    sw0.stop(); print 'source %d: %.2fs'%(isou,sw0.time())
    return g,p
  def plots(self,g,p,s=None,iter=None):
    post = 'file' if iter is None else str(iter)
    gtitle = 'g_'+post
    ptitle = 'p_'+post
    stitle = 's_'+post
    plot(g,cmap=rwb,sperc=100.0,title=gtitle)
    plot(p,cmap=rwb,sperc=100.0,title=ptitle)
    if s is not None:
      plot(s,cmap=jet,cmin=min(_t),cmax=max(_t),cbar='Slowness (s/km)',
        title=stitle)
  def getMisfitFunction(g,isou,do):
    pass
  def residual(ds,do):
    pass

class WaveformInversion(Inversion):
  def residual(self,ds,do,rd=None):
    if rd is None:
      rd = like(ds)
    sub(ds,do,rd)
    return rd
  def getMisfitFunction(self,g,isou,do):
    return WaveformMisfitFunction(g,isou,do)
    #return AmplitudeMisfitFunction(g,isou,do)

class TraveltimeInversion(Inversion):
  def residual(self,ds,do,rt=None):
    if rt is None:
      rt = like(ds)
    makeWarpedResidual(ds,do,rt=rt)
    return rt
  def getMisfitFunction(self,g,isou,do):
    return TraveltimeMisfitFunction(g,isou,do)

class CombinedInversion(Inversion):
  def residual(self,ds,do,rc=None):
    if rc is None:
      rc = like(ds)
    makeWarpedResidual(ds,do,rc=rc)
    return rc
  def getMisfitFunction(self,g,isou,do):
    #return WaveformMisfitFunction(g,isou,do)
    #return TraveltimeMisfitFunction(g,isou,do)
    return CombinedMisfitFunction(g,isou,do)

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

def makeWarpedResidual(ds,do,
  ra=None,rt=None,rc=None,dw=None,v=None,reverseOrder=True):
  #reverseOrder = True
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

def warp(ds,do,dw=None,v=None,rt=None):
  doRmsFilter = False # filter by local rms amplitudes
  doEgain = True # exponential gain from first arrivals
  doAgc = False # agc
  td = 5 # time decimation
  rd = 1 # receiver decimation
  maxShift = 0.5 # max shift in seconds
  sw = Stopwatch(); sw.start()
  nt,nr = len(do[0]),len(do)
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
  #strainMax1,strainMax2 = 0.50,0.20
  #strainMax1,strainMax2 = 0.25,0.10
  strainMax1,strainMax2 = 0.10,0.05
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
    v = like(do)
  for ir in range(nr):
    z = ir*dz
    for it in range(nt):
      t = it*dt
      v[ir][it] = li.interpolate(t,z) # interpolate shifts
  if dw is None:
    dw = like(do)
  warp.applyShifts(v,doc,dw)
  if rt is not None:
    warp.applyShifts(v,timeDerivative(doc),rt)
    mul(v,rt,rt)
  #report('warping',sw)
  #plot(ds,title='f') # used for dynamic warping
  #plot(do,title='g') # used for dynamic warping
  #plot(v,cmap=jet,title='shifts (interpolated)')
  return dw,v

#def makeWarpedResidual(ds,do,u=None,wd=None):
#  reverseOrder = False
#  if reverseOrder:
#    v,dw,_ = warp(do,ds) # wrong order
#    mul(-1.0,v,v)
#    rt = mul(v,timeDerivative(ds)) # traveltime residual
#    rw = sub(dw,do) # warped residual
#  else:
#    v,dw,rt = warp(ds,do) # right order
#    rw = sub(ds,dw) # warped residual
#  #print '  max shift =',max(abs(v))
#  #print '  max combined residual =',max(abs(add(rt,rw)))
#  if u is not None:
#    copy(v,u)
#  if wd is not None:
#    copy(dw,wd)
#  return rt,rw,add(rt,rw)
#
#def warp(ds,do):
#  doEgain = True # exponential gain from first arrivals
#  doRmsEqual = True # equalize rms value
#  doRmsFilter = False # filter by local rms amplitudes
#  doAgc = False # agc
#  td = 5 # time decimation
#  rd = 1 # receiver decimation
#  #maxShift = 0.6 # max shift in seconds
#  maxShift = 0.5 # max shift in seconds
#  sw = Stopwatch(); sw.start()
#  doc = copy(do)
#  if doRmsEqual:
#    do,ds = copy(do),copy(ds)
#    mul(rms(ds)/rms(do),do,do) # equalize rms
#  if doRmsFilter:
#    ds,do = rmsFilter(ds,do,sigma=0.5*maxShift)
#  if doEgain:
#    ds = egain(ds)
#    do = egain(do)
#  if doAgc:
#    ds,do = agc(ds,do)
#  ds,do = addRandomNoise(10.0,ds,do,sigma=1.0)
#  ds = copy(nt/td,nr/rd,0,0,td,rd,ds)
#  do = copy(nt/td,nr/rd,0,0,td,rd,do)
#  #strainMax1 = 0.50
#  #strainMax2 = 0.20
#  #strainMax1 = 0.25
#  #strainMax2 = 0.10
#  strainMax1 = 0.10
#  strainMax2 = 0.05
#  shiftMax = int(maxShift/(td*dt))
#  warp = DynamicWarping(-shiftMax,shiftMax)
#  warp.setErrorExponent(1.0)
#  warp.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.AVERAGE)
#  warp.setStrainMax(strainMax1,strainMax2)
#  warp.setShiftSmoothing(32.0/td,8.0/rd) # shift smoothing
#  warp.setErrorSmoothing(2) # number of smoothings of alignment errors
#  u = warp.findShifts(ds,do)
#  mul(td,u,u) # scale shifts to compensate for decimation
#  li = LinearInterpolator()
#  li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT)
#  li.setUniform(nt/td,td*dt,0.0,nr/rd,rd*dx,0.0,u)
#  v = zerofloat(nt,nr)
#  for ir in range(nr):
#    z = ir*dz
#    for it in range(nt):
#      t = it*dt
#      v[ir][it] = li.interpolate(t,z) # interpolate shifts
#  dw = warp.applyShifts(v,doc)
#  #rt = warp.applyShifts(v,mul(v,timeDerivative(doc)))
#  rt = mul(v,warp.applyShifts(v,timeDerivative(doc)))
#  #report('warping',sw)
#  #plot(u,cmap=jet,title='shifts')
#  #plot(v,cmap=jet,title='shifts (interpolated)')
#  #plot(ds,title='f') # used for dynamic warping
#  #plot(do,title='g') # used for dynamic warping
#  #plot(dw,title='dw')
#  #plot(rt,title='rt')
#  return v,dw,rt # shifts, warped data, traveltime residual

def rmsFilter(ds,do,sigma=0.1):
  x,y = copy(ds),copy(do)
  mul(rms(x)/rms(y),y,y) # equalize rms
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
  nt,n = len(f[0]),len(f)
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
  nt,nr = len(f[0]),len(f)
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
  nt,nr = len(f[0]),len(f)
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

def makeGradient(u,a,d2=True):
  g = correlate(u,a,d2) # gradient
  #p = correlate(u,u) # illumination preconditioner
  p = fillfloat(1.0,nz,nx)
  return g,p

def correlate(u,a,d2=False):
  """Zero-lag correlation."""
  class ReduceInt(Parallel.ReduceInt):
    def compute(self,it):
      if d2:
        t = add(add(u[it-1],u[it+1]),mul(-2.0,u[it])) # 2nd time derivative
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
  t0 = _t[0][0]
  for ix in range(nx):
    iz = 0
    while _t[ix][iz]==t0:
        g[ix][iz] = value
        iz += 1

def laplacian(x):
  n1,n2 = len(x[0]),len(x)
  y = like(x)
  class Loop(Parallel.LoopInt):
    def compute(self,i2):
      for i1 in range(1,n1-1):
        y[i2][i1] = (1.0/6.0)*(-20.0*x[i2][i1]+
          4.0*(x[i2-1][i1  ]+x[i2  ][i1+1]+x[i2  ][i1-1]+x[i2+1][i1  ])+
               x[i2+1][i1+1]+x[i2+1][i1-1]+x[i2-1][i1+1]+x[i2-1][i1-1]);
  Parallel.loop(1,n2-1,Loop())
  return y

import socket
def getMarmousi(sigma=0.5,smul=1.0):
  p = zerofloat(751,2301)
  if socket.gethostname()=='backus.Mines.EDU':
    read("/data/sluo/marmousi/marmousi.dat",p)
  else:
    read("/data/seis/marmousi/marmousi.dat",p)
  p = copy(743,2301,8,0,p)
  div(1000.0,p,p) # convert to slowness
  q = copy(p)
  if sigma is None: # linear s(z) initial model
    mx,mz = len(q),len(q[0])
    sumx = zerofloat(mz)
    for ix in range(mx):
      add(q[ix],sumx,sumx)
    div(sumx,mx,sumx)
    z = linearRegression(rampint(0,1,mz),sumx)
    #SimplePlot.asPoints(sumx).addPoints(mul(smul,z))
    for ix in range(mx):
      copy(z,q[ix])
  else:
    ref = RecursiveExponentialFilter(sigma/0.004)
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE)
    ref.apply(q,q)
  mul(smul,q,q) # scale slowness, e.g., so that simulated data arrives earlier
  v = fillfloat(2.0/3.0,nz,nx)
  c = copy(v)
  copy(248,767,0,0,3,3,p,17,0,1,1,v)
  copy(248,767,0,0,3,3,q,17,0,1,1,c)
  return v,c

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

def getGaussian(splus=0.01):
  v = rampfloat(0.5,-0.001,0.0,nz,nx)
  c = copy(v)
  t = zerofloat(nz,nx); t[nx/2][nz/2] = 1.0
  RecursiveGaussianFilter(40.0).apply00(t,t)
  div(t,max(t),t)
  mul(t,splus,t)
  add(t,v,v)
  for ix in range(nx):
    for iz in range(17):
      v[ix][iz] = 2.0/3.0
      c[ix][iz] = 2.0/3.0
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
  pixel.setInterpolation(PixelsView.Interpolation.LINEAR)
  frame = PlotFrame(panel)
  frame.setFontSizeForSlide(1.5,1.5)
  #frame.setFontSizeForPrint(8.0,222.0)
  if (len(f[0])==nz):
    frame.setSize(1200,500)
  else:
    frame.setSize(1200,1000)
    #frame.setSize(1200,600)
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

def like(x):
  return zerofloat(len(x[0]),len(x))

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

#############################################################################
# Do everything on Swing thread.
import sys,time
class RunMain(Runnable):
  def run(self):
    start = time.time()
    if pngdatDir is not None:
      print 'cleaning '+pngdatDir.split('/')[-2]
      cleanDir(pngdatDir)
    print 'sfile = None' if sfile is None else 'sfile = '+sfile.split('/')[-1]
    print 'gfile = None' if gfile is None else 'gfile = '+gfile.split('/')[-1]
    print 'pfile = None' if pfile is None else 'pfile = '+pfile.split('/')[-1]
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02d'%(h,m,s)
if __name__=='__main__':
  SwingUtilities.invokeLater(RunMain())
