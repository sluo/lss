"""
Test weighted semblance accuracy
"""
from imports import *

#############################################################################
# parameters

#dir = '/Users/sluo/Home/doc/Research/Semblance/temp/'
dir = '/Users/sluo/Desktop/'
datadir = '/Users/sluo/Home/box/idh/trunk/bench/src/sim/cdps/'

st = Sampling(1001,0.004,0.0) # sampling for sythetic data
sx = Sampling(60,0.050,0.050)
#sv = Sampling(101,0.020,1.5)
#sv = Sampling(251,0.008,1.5)
#sv = Sampling(501,0.004,1.5)
sv = Sampling(601,0.002,1.9)
nt,nx,nv = st.count,sx.count,sv.count
dt,dx,dv = st.delta,sx.delta,sv.delta
ft,fx,fv = st.first,sx.first,sv.first

vp = (2.00,3.00)
#vm = (1.98,2.70)
vm = [1.98,2.50]
#vm = (2.00,3.00)
#vm = None

fpeak = 25.0
#tsigma = 8.0 # smoothing window
tsigma = 16.0
#snr = 1.0e8
snr = 1.0e0

ngather = 1000

#############################################################################
# functions

pngDir = '/Users/sluo/Desktop/'
#pngDir = '/Users/sluo/Desktop/png/'

def main(args):
  #goSynthetic()
  #goOneEvent()
  #goAccuracy2()
  #goAccuracy()
  #goAccuracyPlot() 
  #goAccuracyNew2() 
  #testRms()

  #goAccuracyNew() 
  readAndPlot()

def readAndPlot():
  fit,lit = int(0.25*nt),int(0.75*nt)+1
  nit = lit-fit
  pc0 = zerofloat(nit)
  pw0 = zerofloat(nit)
  pc1 = zerofloat(nit)
  pw1 = zerofloat(nit)
  readImage('/Users/sluo/Desktop/sdat2/pc0.dat',pc0)
  readImage('/Users/sluo/Desktop/sdat2/pw0.dat',pw0)
  readImage('/Users/sluo/Desktop/sdat2/pc1.dat',pc1)
  readImage('/Users/sluo/Desktop/sdat2/pw1.dat',pw1)
  stt = Sampling(nit,dt,ft+dt*fit)

  def p(pc,pw,png):
    pixels = 650
    points = pixels*120.0/600.0
    inches = pixels*1.667/600.0
    width,height = pixels,1000
    sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
    sp.setSize(500,1000)
    sp.setVLabel('Time (s)')
    sp.setHLabel('Rms error (m/s)')
    sp.setHLimits(0.0,20.0)
    sp.setFontSizeForPrint(8,points)
    sp.setSize(width,height)
    ppc = sp.addPoints(stt,pc)
    ppw = sp.addPoints(stt,pw)
    ppc.setLineWidth(4.0)
    ppw.setLineWidth(4.0)
    ppc.setLineColor(Color.GRAY)
    sp.paintToPng(720,inches,'/Users/sluo/Desktop/'+png+'.png')

  p(pc0,pw0,'noisefree')
  p(pc1,pw1,'noisy')
  print 'conventional semblance, noise-free average : %f'%(sum(pc0)/nit)
  print '    weighted semblance, noise-free average : %f'%(sum(pw0)/nit)
  print 'conventional semblance,      noisy average : %f'%(sum(pc1)/nit)
  print '    weighted semblance,      noisy average : %f'%(sum(pw1)/nit)

def goAccuracyPlot():
  n = 11

  #pc0 = (7.189,7.178,7.163,7.161,7.176,7.224,7.313,7.447,7.624,7.840,8.081)
  #pw0 = (6.004,5.977,5.950,5.942,5.959,6.011,6.112,6.262,6.452,6.676,6.924)
  #mc0 = (4.707,4.796,4.896,5.002,5.114,5.241,5.382,5.558,5.790,5.983,6.170)
  #mw0 = (3.900,3.970,4.046,4.124,4.211,4.308,4.418,4.541,4.759,4.958,5.132)

  #sv = Sampling(n,0.01,2.4)
  #pc0 = (7.258,7.254,7.250,7.252,7.273,7.325,7.414,7.548,7.731,7.958,8.206)
  #pw0 = (6.035,6.019,6.006,6.007,6.031,6.092,6.195,6.343,6.531,6.749,6.984)
  #mc0 = (4.791,4.877,4.967,5.063,5.170,5.292,5.427,5.575,5.760,5.956,6.153)
  #mw0 = (3.953,4.018,4.085,4.156,4.239,4.336,4.450,4.579,4.715,4.860,5.014)
  #pc1 = (7.485,7.450,7.435,7.435,7.447,7.487,7.575,7.731,7.923,8.179,8.498)
  #pw1 = (6.804,6.718,6.644,6.591,6.581,6.618,6.711,6.863,7.079,7.328,7.618)
  #mc1 = (5.387,5.404,5.449,5.524,5.628,5.765,5.961,6.241,6.542,6.867,7.192)
  #mw1 = (4.836,4.851,4.900,4.999,5.128,5.279,5.444,5.623,5.825,6.061,6.311)

  sv = Sampling(n,1.0,20.0)
  pc0 = (8.206,7.958,7.731,7.548,7.414,7.325,7.273,7.252,7.250,7.254,7.258)
  pw0 = (6.984,6.749,6.531,6.343,6.195,6.092,6.031,6.007,6.006,6.019,6.035)
  mc0 = (6.153,5.956,5.760,5.575,5.427,5.292,5.170,5.063,4.967,4.877,4.791)
  mw0 = (5.014,4.860,4.715,4.579,4.450,4.336,4.239,4.156,4.085,4.018,3.953)
  pc1 = (8.498,8.179,7.923,7.731,7.575,7.487,7.447,7.435,7.435,7.450,7.485)
  pw1 = (7.618,7.328,7.079,6.863,6.711,6.618,6.581,6.591,6.644,6.718,6.804)
  mc1 = (7.192,6.867,6.542,6.241,5.961,5.765,5.628,5.524,5.449,5.404,5.387)
  mw1 = (6.311,6.061,5.825,5.623,5.444,5.279,5.128,4.999,4.900,4.851,4.836)
  
  ppoints(sv,pc0,pw0,mc0,mw0,'notnoisy')
  ppoints(sv,pc1,pw1,mc1,mw1,'noisy')

  #d = 0.4
  #pc2 = (4+3*d,4+3*d,4+3*d)
  #pw2 = (4+2*d,4+2*d,4+2*d)
  #mc2 = (4+1*d,4+1*d,4+1*d)
  #mw2 = (4+0*d,4+0*d,4+0*d)
  #ppoints(Sampling(3,0.01,2.4),pc2,pw2,mc2,mw2,'legend')

def ppoints(sv,pc,pw,mc,mw,png=None):
  sp = SimplePlot()
  sp.setSize(1000,585)
  sp.setVLimits(3.5,9.0)
  #sp.setHLimits(2.395,2.505)
  sp.setHLimits(19.5,30.5)
  sp.setVLabel('rms error (m/s)')
  sp.setHLabel('Primary velocity minus multiple velocity (m/s)')
  sp.setFontSizeForPrint(8.0,240.0)
  pvpc = sp.addPoints(sv,pc)
  pvpw = sp.addPoints(sv,pw)
  pvmc = sp.addPoints(sv,mc)
  pvmw = sp.addPoints(sv,mw)
  pvpc.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pvpw.setMarkStyle(PointsView.Mark.FILLED_SQUARE)
  pvmc.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
  pvmw.setMarkStyle(PointsView.Mark.HOLLOW_SQUARE)
  if png:
    sp.paintToPng(720.0,3.33,'/Users/sluo/Desktop/'+png+'.png')

def goOneEvent():
  vps = Velan.makeLinearVelocity(vp[0],vp[1],st)
  vms = Velan.makeLinearVelocity(vm[0],vm[1],st)
  p = Velan.xmakeRickerGather(fpeak,vps,st,sx)
  #p = add(p,Velan.xmakeRickerGather(fpeak,vms,st,sx)) # multiples
  #p = Velan.addRandomNoise(snr,p) # noise
  sc = Velan.velocitySpectrum(st,sx,p,sv,tsigma,False); splot(sc,'sc')
  sw = Velan.velocitySpectrum(st,sx,p,sv,tsigma,True); splot(sw,'sw')
  SimplePlot.asPixels(p)

def goAccuracy2():
  #for i in range(4):
  #  vm[1] = 2.5-i*0.01
  for i in range(10):
    vm[1] = 2.40+i*0.01
    print '\nvm[1]=%.2f'%vm[1]
    goAccuracy()

def goAccuracy():
  if not vm:
    return
  vps = Velan.makeLinearVelocity(vp[0],vp[1],st)
  vms = Velan.makeLinearVelocity(vm[0],vm[1],st)
  vmin = fillfloat(sv.first,nt)
  vmax = fillfloat(sv.last,nt)
  vhaf = Velan.makeLinearVelocity(0.5*(vp[0]+vm[0]),0.5*(vp[1]+vm[1]),st)
  def makeGather(rand):
    p = Velan.makeRickerGather(rand,fpeak,vps,st,sx)
    p = add(p,Velan.makeRickerGather(fpeak,vms,st,sx)) # multiples
    p = Velan.addRandomNoise(snr,p) # noise
    return p

  random = Random(314159) # TODO try different seed
  print ngather,'realizations'
  ps = zerofloat(nt,nx,ngather)
  for igather in range(ngather):
    ps[igather] = makeGather(random)

  print 'Semblance'
  sc = Velan.velocitySpectrum(st,sx,ps,sv,tsigma,False)
  sw = Velan.velocitySpectrum(st,sx,ps,sv,tsigma,True)

  print 'Picking'
  ppc = Velan.pickPeaks(sv,vhaf,vmax,sc) # primary,conventional
  ppw = Velan.pickPeaks(sv,vhaf,vmax,sw) # primary,weighted
  pmc = Velan.pickPeaks(sv,vmin,vhaf,sc) # multiple,conventional
  pmw = Velan.pickPeaks(sv,vmin,vhaf,sw) # multiple,weighted

  print 'Errors'
  fit,lit = int(0.25*nt),int(0.75*nt)
  nit = lit-fit
  epc = zerofloat(lit-fit,ngather)
  epw = zerofloat(lit-fit,ngather)
  emc = zerofloat(lit-fit,ngather)
  emw = zerofloat(lit-fit,ngather)
  for igather in range(ngather):
    for it in range(fit,lit):
      epc[igather][it-fit] = ppc[igather][it]-vps[it]
      epw[igather][it-fit] = ppw[igather][it]-vps[it]
      emc[igather][it-fit] = pmc[igather][it]-vms[it]
      emw[igather][it-fit] = pmw[igather][it]-vms[it]
  minpc,maxpc,rmspc = min(abs(epc)),max(abs(epc)),rms(epc),
  minpw,maxpw,rmspw = min(abs(epw)),max(abs(epw)),rms(epw),
  minmc,maxmc,rmsmc = min(abs(emc)),max(abs(emc)),rms(emc),
  minmw,maxmw,rmsmw = min(abs(emw)),max(abs(emw)),rms(emw),
  print 'primaries'
  print '  c: min=%.8f, max=%.8f, rms=%.8f'%(minpc,maxpc,rmspc)
  print '  w: min=%.8f, max=%.8f, rms=%.8f'%(minpw,maxpw,rmspw)
  print 'multiples'
  print '  c: min=%.8f, max=%.8f, rms=%.8f'%(minmc,maxmc,rmsmc)
  print '  w: min=%.8f, max=%.8f, rms=%.8f'%(minmw,maxmw,rmsmw)

  print 'Averages'
  asc = zerofloat(nt,nv) # semblance
  asw = zerofloat(nt,nv)
  appc = zerofloat(nt) # picks
  appw = zerofloat(nt)
  apmc = zerofloat(nt)
  apmw = zerofloat(nt)
  for igather in range(ngather):
    add(asc,sc[igather],asc)
    add(asw,sw[igather],asw)
    add(appc,ppc[igather],appc)
    add(appw,ppw[igather],appw)
    add(apmc,pmc[igather],apmc)
    add(apmw,pmw[igather],apmw)
  div(asc,ngather,asc)
  div(asw,ngather,asw)
  div(appc,ngather,appc)
  div(appw,ngather,appw)
  div(apmc,ngather,apmc)
  div(apmw,ngather,apmw)
  plotPicks(asc,appc,apmc,title='conventional',png='%.2fsc'%vm[1])
  plotPicks(asw,appw,apmw,title='weighted',png='%.2fsw'%vm[1])

  ## old (wrong) rms
  #sub(appc,vps,appc)
  #sub(appw,vps,appw)
  #sub(apmc,vms,apmc)
  #sub(apmw,vms,apmw)
  #def rms1(x):
  #  y = 0.0
  #  for i1 in range(fit,lit):
  #    xi = x[i1]
  #    y += xi*xi
  #  return sqrt(y/nit)
  #print '\nprimaries'
  #print '  c: rms=%.8f'%rms1(appc)
  #print '  w: rms=%.8f'%rms1(appw)
  #print 'multiples'
  #print '  c: rms=%.8f'%rms1(apmc)
  #print '  w: rms=%.8f'%rms1(apmw)

  #fti = nt/4
  #nti = nt/2
  #print 'primaries'
  #plotError(copy(nti,fti,vps),copy(nti,fti,vpc),copy(nti,fti,vpw),\
  #          'primaries','%.2fprimaries'%vm[1])
  #print 'multiples'
  #plotError(copy(nti,fti,vms),copy(nti,fti,vmc),copy(nti,fti,vmw),\
  #          'multiples','%.2fmultiples'%vm[1])

def goAccuracyNew2():
  for i in range(5):
    vm[1] = 2.40+i*0.05
    print '\nvm[1]=%.2f'%vm[1]
    goAccuracyNew()

def goAccuracyNew():
  if not vm:
    return
  vps = Velan.makeLinearVelocity(vp[0],vp[1],st)
  vms = Velan.makeLinearVelocity(vm[0],vm[1],st)
  vmin = fillfloat(sv.first,nt)
  vmax = fillfloat(sv.last,nt)
  vhaf = Velan.makeLinearVelocity(0.5*(vp[0]+vm[0]),0.5*(vp[1]+vm[1]),st)
  def makeGather(rand):
    p = Velan.makeRickerGather(rand,fpeak,vps,st,sx)
    p = add(p,Velan.makeRickerGather(fpeak,vms,st,sx)) # multiples
    p = Velan.addRandomNoise(snr,p) # noise
    return p

  random = Random(314) # TODO try different seed
  print ngather,'realizations'
  ps = zerofloat(nt,nx,ngather)
  for igather in range(ngather):
    ps[igather] = makeGather(random)

  print 'Semblance'
  sc = Velan.velocitySpectrum(st,sx,ps,sv,tsigma,False)
  sw = Velan.velocitySpectrum(st,sx,ps,sv,tsigma,True)

  print 'Picking'
  ppc = Velan.pickPeaks(sv,vhaf,vmax,sc) # primary,conventional
  ppw = Velan.pickPeaks(sv,vhaf,vmax,sw) # primary,weighted
  pmc = Velan.pickPeaks(sv,vmin,vhaf,sc) # multiple,conventional
  pmw = Velan.pickPeaks(sv,vmin,vhaf,sw) # multiple,weighted

  #print 'Errors'
  #fit,lit = int(0.25*nt),int(0.75*nt)
  #nit = lit-fit
  #epc = zerofloat(lit-fit,ngather)
  #epw = zerofloat(lit-fit,ngather)
  #emc = zerofloat(lit-fit,ngather)
  #emw = zerofloat(lit-fit,ngather)
  #for igather in range(ngather):
  #  for it in range(fit,lit):
  #    epc[igather][it-fit] = ppc[igather][it]-vps[it]
  #    epw[igather][it-fit] = ppw[igather][it]-vps[it]
  #    emc[igather][it-fit] = pmc[igather][it]-vms[it]
  #    emw[igather][it-fit] = pmw[igather][it]-vms[it]
  #rpc = zerofloat(nit)
  #rpw = zerofloat(nit)
  #rmc = zerofloat(nit)
  #rmw = zerofloat(nit)
  #for igather in range(ngather):
  #  for it in range(fit,lit):
  #    rpci = epc[igather][it-fit]
  #    rpwi = epw[igather][it-fit]
  #    rmci = emc[igather][it-fit]
  #    rmwi = emw[igather][it-fit]
  #    rpc[it] += rpci*rpci
  #    rpw[it] += rpwi*rpwi
  #    rmc[it] += rmci*rmci
  #    rmw[it] += rmwi*rmwi
  #div(rpc,nit*ngather,rpc)
  #div(rpw,nit*ngather,rpw)
  #div(rmc,nit*ngather,rmc)
  #div(rmw,nit*ngather,rmw)
  #sqrt(rpc,rpc)
  #sqrt(rpw,rpw)
  #sqrt(rmc,rmc)
  #sqrt(rmw,rmw)

  print 'Errors'
  fit,lit = int(0.25*nt),int(0.75*nt)+1
  nit = lit-fit
  rpc = zerofloat(nit)
  rpw = zerofloat(nit)
  rmc = zerofloat(nit)
  rmw = zerofloat(nit)
  rmspc,rmspw,rmsmc,rmsmw = 0.0,0.0,0.0,0.0
  for igather in range(ngather):
    for it in range(fit,lit):
      epc = 1000.0*(ppc[igather][it]-vps[it])
      epw = 1000.0*(ppw[igather][it]-vps[it])
      emc = 1000.0*(pmc[igather][it]-vms[it])
      emw = 1000.0*(pmw[igather][it]-vms[it])
      rpc[it-fit] += epc*epc
      rpw[it-fit] += epw*epw
      rmc[it-fit] += emc*emc
      rmw[it-fit] += emw*emw
      rmspc += epc*epc
      rmspw += epw*epw
      rmsmc += emc*emc
      rmsmw += emw*emw
  div(rpc,ngather,rpc)
  div(rpw,ngather,rpw)
  div(rmc,ngather,rmc)
  div(rmw,ngather,rmw)
  sqrt(rpc,rpc)
  sqrt(rpw,rpw)
  sqrt(rmc,rmc)
  sqrt(rmw,rmw)
  print 'primaries'
  print '  c: rms=%.8f'%sqrt(rmspc/nit/ngather)
  print '  w: rms=%.8f'%sqrt(rmspw/nit/ngather)
  print 'multiples'
  print '  c: rms=%.8f'%sqrt(rmsmc/nit/ngather)
  print '  w: rms=%.8f'%sqrt(rmsmw/nit/ngather)

  print ''
  print sum(rpc)
  print sum(rpw)
  print sum(rmc)
  print sum(rmw)

  pixels = 650
  points = pixels*120.0/600.0
  inches = pixels*1.667/600.0
  width,height = pixels,1000

  multiples = False
  stt = Sampling(nit,dt,ft+dt*fit)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #sp.addTitle('%.2f'%vm[1])
  sp.setSize(500,1000)
  sp.setVLabel('Time (s)')
  sp.setHLabel('Rms error (m/s)')
  sp.setHLimits(0.0,20.0)
  sp.setFontSizeForPrint(8,points)
  sp.setSize(width,height)
  p0pc = sp.addPoints(stt,rpc)
  p0pw = sp.addPoints(stt,rpw)
  p0pc.setLineWidth(4.0)
  p0pw.setLineWidth(4.0)
  p0pc.setLineColor(Color.GRAY)
  if multiples:
    p1mc = sp.addPoints(stt,rmc)
    p1mw = sp.addPoints(stt,rmw)
    p1mc.setLineWidth(4.0)
    p1mw.setLineWidth(4.0)
    p1mc.setLineColor(Color.GRAY)
    p1mc.setLineStyle(PointsView.Line.DOT)
    p1mw.setLineStyle(PointsView.Line.DOT)
  sp.paintToPng(720,inches,'/Users/sluo/Desktop/%d'%(1000.0*vm[1])+'.png')
  writeImage('/Users/sluo/Desktop/sdat/pc.dat',rpc)
  writeImage('/Users/sluo/Desktop/sdat/pw.dat',rpw)
  writeImage('/Users/sluo/Desktop/sdat/mc.dat',rmc)
  writeImage('/Users/sluo/Desktop/sdat/mw.dat',rmw)

def testRms():
  fill = 5.0
  n1,n2 = 1234,4848
  x = fillfloat(fill,n1,n2)
  #print 'x = fillfloat(%.1f,%d,%d)'%(fill,n1,n2)
  print rms(x)

def rms(x):
  n1,n2 = len(x[0]),len(x)
  y = 0.0
  for i1 in range(n1):
    for i2 in range(n2):
      xi = x[i2][i1]
      y += xi*xi
  return sqrt(y/n1/n2)

#############################################################################
# plots

def splot(x,title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(500,800)
  if title:
    sp.addTitle(title)
  pv = sp.addPixels(x)
  pv.setColorModel(ColorMap.JET)
  #pv.setClips(0.0,1.0)
  cb = sp.addColorBar()
  #cb.setWidthMinimum(80)
  if png:
    sp.paintToPng(180,5.0,pngDir+png+'.png')

def plotPicks(s,primary=None,multiple=None,title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(500,800)
  cb = sp.addColorBar()
  cb.setWidthMinimum(80)
  if title:
    sp.addTitle(title)
  pv = sp.addPixels(st,sv,s)
  pv.setColorModel(jet)
  #pv.setClips(0,0.01)
  if primary:
    sp.addPoints(st,primary)
  if multiple:
    sp.addPoints(st,multiple)
  if png:
    sp.paintToPng(180,5.0,pngDir+png+'.png')

def plotError(vt,vc,vw,title=None,png=None):
  t = sub(vt,vt)
  c = sub(vc,vt)
  w = sub(vw,vt)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(300,800)
  sp.setVLabel('Time (s)')
  sp.setHLabel('Error (km/s)')
  if title:
    sp.addTitle(title)
  def rms(x):
    n = len(x)
    y = 0.0
    for i in range(n):
      y += x[i]*x[i]
    return sqrt(y/n)
  minc = min(abs(c))
  minw = min(abs(w))
  maxc = max(abs(c))
  maxw = max(abs(w))
  rmsc = rms(c)
  rmsw = rms(w)
  print 'convention: min=%.8f, max=%.8f, rms=%.8f'%(minc,maxc,rmsc)
  print '  weighted: min=%.8f, max=%.8f, rms=%.8f'%(minw,maxw,rmsw)
  s = Sampling(len(vt),dt,1.0)
  pt = sp.addPoints(s,t)
  pc = sp.addPoints(s,c); pc.setLineColor(Color.BLUE)
  pw = sp.addPoints(s,w); pw.setLineColor(Color.RED)
  if png:
    sp.paintToPng(180,5.0,pngDir+png+'.png')


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
  return pp

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

def readImage(name,image):
  ais = ArrayInputStream(name)
  ais.readFloats(image)
  ais.close()

def writeImage(name,image):
  aos = ArrayOutputStream(name)
  aos.writeFloats(image)
  aos.close()

def readImage(name,image):
  ais = ArrayInputStream(name)
  ais.readFloats(image)
  ais.close()

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
