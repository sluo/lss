"""
Figures for slides
"""
from imports import *

#subset = None
#subset = "a"
subset = "b"

pngDir = None
#pngDir = "/Users/sluo/Desktop/"
#pngDir = "/Users/sluo/Desktop/png/"
#pngDir = "/Users/sluo/Home/doc/research/fault/tex/figures/"

if subset=='a':
  dataDir = '/data/seis/f3/suba/'
  s1 = Sampling(90,0.004,1.484)
  s2 = Sampling(221,0.025,1.25)
  s3 = Sampling(220,0.025,2.5)
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  k1,k2,k3 = 22,48,151
  #k1 = 28
  #k1 = 34
  #h1 = 26
  h1 = 32
  #h1s = [20,40,60]
  h1s = [32,42,52,62]
elif subset=='b':
  dataDir = '/data/seis/f3/subb/'
  s1 = Sampling(90,0.004,1.024)
  s2 = Sampling(221,0.025,1.25)
  s3 = Sampling(220,0.025,2.5)
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  k1,k2,k3 = 56,48,150
  #k1 = 12
  #k1 = 37
  h1 = 53
  h1s = [20,40,60]
else:
  dataDir = '/data/seis/f3/' # f3x =tranpose23(f3)
  s1 = Sampling(462,0.004,0.0)
  s2 = Sampling(951,0.025,0.0)
  s3 = Sampling(651,0.025,0.0)
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  k2,k3 = n2/2,n3/2
  #k1 = 315 # chaotic zone
  k1 = 415 # en echelon faults

azimuth,elevation = 290,50 # for 3D views

#dataDir = '/data/seis/f3/sub'+subset+'/'
#n1,n2,n3 = 90,221,220
#d1,d2,d3 = 0.004,0.025,0.025
##f1,f2,f3 = 1.204 if subset=='a' else 3.044,1.25,2.5
#f1,f2,f3 = 1.484 if subset=='a' else 1.024,1.25,2.5
#s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
#k1 =  25 if subset=='a' else 56
#k2 = 108 if subset=='a' else 48
#k3 = 161 if subset=='a' else 150
#azimuth,elevation =210,25 # for 3D views
#h1s = [32] if subset=="a" else [55] # horizons

#############################################################################

def main(args):
  #goFull()
  goSubset()

def goSubset():
  if subset==None:
    return
  def goPlot3():
    seismic = True
    blended = False
    shiftsr = True
    shiftss = False
    if seismic:
      g = read("g")
      h = read("h")
      f = read("f")
      cmin,cmax = 0.8*min(g),0.8*max(g)
      cmin,cmax = 0.8*min(g),0.8*max(g)
      plot3(g,cmin=cmin,cmax=cmax,cbar="Amplitude",name="g"+subset)
      plot3(h,cmin=cmin,cmax=cmax,cbar="Amplitude",name="h"+subset)
      plot3(f,cmin=cmin,cmax=cmax,cbar="Amplitude",name="f"+subset)
      #display(g,cmin=cmin,cmax=cmax);
    if blended:
      g,h,f = read('g'),read('h'),read('f')
      #g,h,f = read('g'),read('h'),read('tmprr') # XXX
      r = array(read('r1'),read('r2'),read('r3'))
      t1 = mul(d1*1000.0,read("t1"))
      q1 = mul(-d1*1000.0,read("q1"))
      ming,maxg = 0.8*min(g),0.8*max(g)
      
      null = -0.4938
      minq,_ = getClipsForPlot3(q1); maxq = 0.4*max(q1)
      knotch = int(256.0*abs((null-minq)/(maxq-minq)))
      amap = getNotchAlphaColorMap(jet,knotch)
      """
      minq,maxq = getClipsForPlot3(q1,0.4)
      amap = getLinearAlphaColorMap(jet)
      """

      plot3(g,cmin=ming,cmax=maxg,
            t=t1,cmap1=amap,cmin1=minq,cmax1=maxq,
            cbar="Vertical throw (ms)",name="t1"+subset)
      plot3(q1,cmap=jet,cmin=minq,cmax=maxq,nearest=False,
            cbar="Vertical throw (ms)",name="q1"+subset)

      q2 = mul(-1000*d2,read('q2'))
      q3 = mul(-1000*d3,read('q3'))
      minq2,_ = getClipsForPlot3(q2); maxq2 = 0.4*max(q2)
      minq3,_ = getClipsForPlot3(q3); maxq3 = 0.4*max(q3)
      plot3(q2,cmap=jet,cmin=minq2,cmax=maxq2,nearest=False,
            cbar="Inline throw (m)",name="q2"+subset)
      plot3(q3,cmap=jet,cmin=minq3,cmax=maxq3,nearest=False,
            cbar="Crossline throw (m)",name="q3"+subset)

    if shiftsr:
      r1 = mul(d1*1000.0,read('r1'))
      plot3(r1,cmap=jet,cbar="Vertical shift (ms)",
            name="r1"+subset)
      #r2 = mul(d2*1000.0,read('r2'))
      r2 = read('r2')
      plot3(r2,cmap=jet,cbar="Inline shift (m)",
            name="r2"+subset)
      #r3 = mul(d3*1000.0,read('r3'))
      r3 = read('r3')
      plot3(r3,cmap=jet,cbar="Crossline shift (m)",
            name="r3"+subset)

    if shiftss:
      s1 = mul(d1*1000.0,read('s1'))
      plot3(s1,cmap=jet,cbar="Vertical shift (ms)",
            name="s1"+subset)
  def plot3(f,cmap=gray,cmin=0.0,cmax=0.0,perc=100.0,
            t=None,cmap1=jet,cmin1=0.0,cmax1=0.0,perc1=0.0,
            nearest=False,cbar=None,name=None):
    o = PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT
    a = PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM
    panel = PlotPanelPixels3(o,a,s1,s2,s3,f)
    if t:
      s12 = slice12(k3,t)
      s13 = slice13(k2,t)
      s23 = slice23(k1,t)
      p12 = panel.addPixels(1,0,s1,s2,s12)
      p13 = panel.addPixels(1,1,s1,s3,s13)
      p23 = panel.addPixels(0,0,s3,s2,transpose(s23))
      p12.setInterpolation(PixelsView.Interpolation.NEAREST)
      p13.setInterpolation(PixelsView.Interpolation.NEAREST)
      p23.setInterpolation(PixelsView.Interpolation.NEAREST)
      if cmin1<cmax1:
        p12.setClips(cmin1,cmax1)
        p13.setClips(cmin1,cmax1)
        p23.setClips(cmin1,cmax1)
      p12.setColorModel(cmap1)
      p13.setColorModel(cmap1)
      p23.setColorModel(cmap1)
    if nearest:
      panel.setInterpolation(PixelsView.Interpolation.NEAREST)
    panel.setColorModel(cmap)
    panel.setSlices(k1,k2,k3)
    panel.setLabel1("Time (s)")
    panel.setLabel2("Inline (km)")
    panel.setLabel3("Crossline (km)")
    panel.setLineColor(Color.YELLOW)
    if (cbar):
      panel.addColorBar(cbar)
      panel.setColorBarWidthMinimum(130)
    panel.setBackground(Color.WHITE)
    panel.setInterval1(0.1)
    panel.setInterval2(1.0)
    panel.setInterval3(1.0)
    if cmin<cmax:
      panel.setClips(cmin,cmax)
    if perc<100.0:
      panel.setPercentiles(100-perc,perc)
    #panel.mosaic.setWidthElastic(1,100)
    #panel.mosaic.setHeightElastic(0,205)
    panel.mosaic.setHeightElastic(0,140)
    frame = PlotFrame(panel)
    frame.setFontSizeForSlide(1.0,1.0)
    #frame.setSize(1150,800)
    frame.setSize(1050,760)
    frame.setVisible(True)
    if name and pngDir:
      frame.paintToPng(720,7.0,pngDir+name+".png")
    return frame

  def go3dViews():
    seismic = True
    normals = False
    blended = False
    shiftss = False
    horizon = False
    g = read("g")
    #cmin,cmax = 0.8*min(g),0.8*max(g)
    cmin,cmax = 0.7*min(g),0.7*max(g)
    if seismic:
      h = read("h")
      f = read("f")
      displaySeismic(g,cmin,cmax,name='g')
      displaySeismic(h,cmin,cmax,name='h')
      displaySeismic(f,cmin,cmax,name='f')
    if normals:
      h = read('h')
      u1 = mul(d1*1000,read('u1'))
      u2 = mul(d2*1000,read('u2'))
      u3 = mul(d3*1000,read('u3'))
      jettrans = ColorMap.setAlpha(jet,0.40)
      displaySeismicAndShifts(h,u1,jettrans,
                              cbar='Vertical component (ms)',name='u1')
      displaySeismicAndShifts(h,u2,jettrans,
                              cbar='Inline component (m)',name='u2')
      displaySeismicAndShifts(h,u3,jettrans,
                              cbar='Crossline component (m)',name='u3')

    if blended:
      g,h,f = read('g'),read('h'),read('f')
      r = array(read('r1'),read('r2'),read('r3'))
      t1 = mul(d1*1000.0,read("t1"))
      q1 = mul(-d1*1000.0,read("q1"))

      null = -0.4938
      minq,_ = getClipsForPlot3(q1); maxq = 0.4*max(q1)
      knotch = int(256.0*abs((null-minq)/(maxq-minq)))
      jetnotch = getNotchAlphaColorMap(jet,knotch)
      displaySeismicAndShifts(g,t1,jetnotch,cmin2=minq,cmax2=maxq,
                              cbar='Vertical throw (ms)',name='t1')
      """
      jetlinear = getLinearAlphaColorMap(jet)
      displaySeismicAndShifts(g,t1,jetlinear)
      """

      jettrans = ColorMap.setAlpha(jet,0.40)
      displaySeismicAndShifts(g,q1,jettrans,cmin2=minq,cmax2=maxq,
                              cbar='Vertical throw (ms)',name='q1')

      """
      q2 = mul(-1000*d2,read('q2'))
      q3 = mul(-1000*d3,read('q3'))
      minq2,_ = getClipsForPlot3(q2); maxq2 = 0.4*max(q2)
      minq3,_ = getClipsForPlot3(q3); maxq3 = 0.4*max(q3)
      displaySeismicAndShifts(g,q2,jettrans,cmin2=minq2,cmax2=maxq2,
                              cbar='Inline throw (m)',name='q2')
      displaySeismicAndShifts(g,q3,jettrans,cmin2=minq3,cmax2=maxq3,
                              cbar='Inline throw (m)',name='q3')
      """


    if shiftss:
      #s1 = mul(d1*1000.0,read('s1'))
      s1 = mul(d1*1000,read('s1'))
      displaySeismicAndShifts(g,s1,jettrans,cmin2=min(s1),cmax2=max(s1),
                              cbar='Vertical shift (ms)',name='s1_3d')
    if horizon:

      world1,frame = displaySeismic(g,cmin,cmax)
      xyz = getHorizonVertices(h1)
      tg = QuadGroup(True,xyz)
      tg.setColor(Color.ORANGE)
      world1.addChild(tg)
      """
      world2 = displaySeismic(g,cmin,cmax)
      colors = [Color.ORANGE,Color.BLUE,Color.RED,Color.CYAN]
      r = array(read('r1'),read('r2'),read('r3'))
      x = array(read('x1'),read('x2'),read('x3'))
      for i in range(len(h1s)):
        xyz = getHorizonVertices(h1s[i],r,x)
        tg = QuadGroup(True,xyz)
        tg.setColor(colors[i])
        world2.addChild(tg)
      """
      if pngDir:
        #frame.getViewCanvas().paintToFile(pngDir+name+'.png')
        frame.paintToFile(pngDir+'horizon.png')

  def displaySeismic(f,cmin=0,cmax=0,perc=100,name=None):
    world = World()
    ipg = ImagePanelGroup(s1,s2,s3,f)
    ipg.setColorModel(ColorMap.getGray())
    ipg.setSlices(3*n1/4,n2/4,6*n3/8)
    #ipg.setSlices(3*n1/4,n2/4,5*n3/8)
    #ipg.setSlices(3*n1/4,n2/4,4*n3/8)
    #ipg.setSlices(3*n1/4,n2/4,3*n3/8)
    #ipg.setSlices(3*n1/4,n2/4,2*n3/8)
    world.addChild(ipg)
    frame = makeFrame(world,name)
    colorbar = addColorBar(frame,'Amplitude')
    colorbar.setWidthMinimum(250)
    ipg.addColorMapListener(colorbar)
    if cmax>cmin:
      ipg.setClips(cmin,cmax)
    if perc<100:
      ipg.setPercentiles(100-perc,perc)
    if pngDir and name:
      #frame.getViewCanvas().paintToFile(pngDir+name+'.png')
      frame.paintToFile(pngDir+name+'.png')
      #colorbar.paintToPng(colorbar.getWidth(),1.0,pngDir+'c.png')
      pass
    return world,frame

  def displaySeismicAndShifts(f,s,cmap2=jet,cmin2=0.0,cmax2=0.0,\
                              cbar=None,name=None):
    world = World()

    #addImageToWorld(world,f,gray)
    #cmap = getLinearAlphaColorMap(jet)
    #cmap.setValueRange(0.0,max(s))
    #addImageToWorld(world,s,cmap)

    ipg = ImagePanelGroup2(s1,s2,s3,f,s)
    ipg.setColorModel1(ColorMap.getGray())
    ipg.setColorModel2(cmap2)
    #ipg.setSlices(k1,k2,k3)
    ipg.setSlices(3*n1/4,n2/4,3*n3/4)
    cmin1,cmax1 = 0.8*min(f),0.8*max(f)
    #cmin2,cmax2 = 0.0,0.4*max(s)
    ipg.setClips1(cmin1,cmax1)
    world.addChild(ipg)
    frame = makeFrame(world,name)
    if cbar:
      colorbar = addColorBar(frame,cbar)
      colorbar.setWidthMinimum(250)
      ipg.addColorMap2Listener(colorbar)
    if cmax2>cmin2:
      ipg.setClips2(cmin2,cmax2)
    else:
      ipg.setClips2(min(s),max(s))
    if pngDir and name:
      #frame.getViewCanvas().paintToFile(pngDir+name+'.png')
      frame.paintToFile(pngDir+name+'.png')
      pass

  #goPlot3()
  go3dViews()


def getHorizonVertices(k1,r=None,x=None):
  """Computes horizon vertex coordinates for use in QuadGroup."""
  # Coordinates: x = x(x)
  if r:
    r1,r2,r3 = r[0],r[1],r[2]
  else:
    #r1,r2,r3 = read('r1'),read('r2'),read('r3')
    r = array(read('r1'),read('r2'),read('r3'))
  if x:
    x1,x2,x3 = x[0],x[1],x[2]
  else:
    x1,x2,x3 = read('x1'),read('x2'),read('x3')
    """
    x2 = rampfloat(0.0,0.0,1.0,0.0,n1,n2,n3)
    x3 = rampfloat(0.0,0.0,0.0,1.0,n1,n2,n3)
    x1 = Faulter.applyShiftsRLinear(x1,s) # x1(u) = x1(u-s(u))
    x2 = Faulter.applyShiftsRLinear(x2,s) # x2(u) = x2(u-s(u))
    x3 = Faulter.applyShiftsRLinear(x3,s) # x3(u) = x3(u-s(u))
    """
  # Fault location: w = w(x)
  t1,t2,t3,z1,z2,z3 = getGriddedThrows()
  e = Faulter.getGriddedThrows(0.0,t1,t2,t3,z1,z2,z3)
  z1,z2,z3 = e[3],e[4],e[5] # live value coords
  w = zerofloat(n1,n2,n3)
  for i in range(len(z1)):
    i1,i2,i3 = int(z1[i]),int(z2[i]),int(z3[i])
    w[i3][i2][i1] = 1.0
  w = FlattenerUtil.applyShiftsR(w,r) # w(u) = w(u-r(u))

  # Scale by sampling intervals
  x1 = add(f1,mul(d1,x1))
  x2 = add(f2,mul(d2,x2))
  x3 = add(f3,mul(d3,x3))

  # Slice
  x1 = slice23(k1,x1) # x1(u=k1)
  x2 = slice23(k1,x2) # x2(u=k1)
  x3 = slice23(k1,x3) # x3(u=k1)
  w = slice23(k1,w)   # w(u=k1)
  return FlattenerUtil.makeQuadVertices(x1,x2,x3,w)
  

def goFull():
  if subset:
    return
  def plot3(f,cmin=0.0,cmax=0.0,perc=100.0,cbar=None,name=None):
    o = PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT
    a = PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM
    panel = PlotPanelPixels3(o,a,s1,s2,s3,f)
    panel.setColorModel(gray)
    panel.setSlices(k1,k2,k3)
    panel.setLabel1("Time (s)")
    panel.setLabel2("Inline (km)")
    panel.setLabel3("Crossline (km)")
    panel.setLineColor(Color.YELLOW)
    if (cbar):
      panel.addColorBar(cbar)
      panel.setColorBarWidthMinimum(90)
    panel.setBackground(Color.WHITE)
    panel.setInterval1(1.0)
    panel.setInterval2(5.0)
    panel.setInterval3(5.0)
    if cmin<cmax:
      panel.setClips(cmin,cmax)
    if perc<100.0:
      panel.setPercentiles(100-perc,perc)
    #panel.mosaic.setWidthElastic(0,130)
    #panel.mosaic.setHeightElastic(0,130)
    #panel.mosaic.setWidthElastic(1,120)
    #panel.mosaic.setHeightElastic(0,130)
    #panel.mosaic.setWidthMinimum(1,380)
    #panel.mosaic.setHeightMinimum(0,380)
    frame = PlotFrame(panel)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.setSize(1024,768)
    frame.setVisible(True)
    if name and pngDir:
      frame.paintToPng(720,7.0,pngDir+name+".png")
    return frame

  g = read('f3x')
  #g = zeros(n1,n2,n3)
  display(g,perc=99.5)
  plot3(g,perc=99.0,name='f3')



#############################################################################
# Plots for CWP report.

def goCwp():

  teaser = True
  plot3 = True

  def goPlot3():
    seismic = False
    blended = True
    shiftss = False
    if seismic:
      g = read("g")
      h = read("h")
      f = read("f")
      cmin,cmax = 0.9*min(g),0.9*max(g)
      plot3(g,cmin=cmin,cmax=cmax,cbar="Amplitude",name="g"+subset)
      plot3(h,cmin=cmin,cmax=cmax,cbar="Amplitude",name="h"+subset)
      plot3(f,cmin=cmin,cmax=cmax,cbar="Amplitude",name="f"+subset)
      #display(g,cmin=cmin,cmax=cmax);
    if blended:
      g,h,f = read('g'),read('h'),read('f')
      r = array(read('r1'),read('r2'),read('r3'))
      t1 = mul(d1*1000.0,read("t1"))
      q1 = mul(-d1*1000.0,read("q1"))
      ming,maxg = 0.9*min(g),0.9*max(g)
      
      null = -0.4938
      minq,_ = getClipsForPlot3(q1); maxq = 0.4*max(q1)
      knotch = int(256.0*abs((null-minq)/(maxq-minq)))
      amap = getNotchAlphaColorMap(jet,knotch)
      """
      minq,maxq = getClipsForPlot3(q1,0.4)
      amap = getLinearAlphaColorMap(jet)
      """
      plot3(g,cmin=ming,cmax=maxg,
            t=t1,cmap1=amap,cmin1=minq,cmax1=maxq,
            cbar="Vertical component of throw (ms)",name="t"+subset)
      plot3(q1,cmap=jet,cmin=minq,cmax=maxq,nearest=False,
            cbar="Vertical component of throw (ms)",name="q"+subset)

      #plot3(h,t1,cmin=cmin,cmax=cmax,
      #      cbar="Vertical component of throw (ms)",name="t"+subset)
      #plot3(f,FlattenerUtil.applyShiftsR(t1,r),cmin=cmin,cmax=cmax,
      #      cbar="Vertical component of throw (ms)",name="t"+subset)
      plot3(h,cmin=ming,cmax=maxg,cbar="Amplitude",name="h"+subset)
      plot3(f,cmin=ming,cmax=maxg,cbar="Amplitude",name="f"+subset)

    if shiftss:
      s1 = mul(d1*1000.0,read('s1'))
      plot3(s1,cmap=jet,cbar="Vertical component of composite shift (ms)",
            name="s1"+subset)

  def goTeaser():
    def p(x,cmap=gray,cmin=0.0,cmax=0.0,name=None):
      pan = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
      pan.setHLabel("Crossline (km)")
      pan.setVLabel("Time (s)")
      #cb = pan.addColorBar("Amplitude")
      #cb.setWidthMinimum(150)
      pix = pan.addPixels(s1,s2,x)
      pix.setColorModel(cmap)
      if cmin<cmax:
        pix.setClips(cmin,cmax)
      pix.setInterpolation(PixelsView.Interpolation.LINEAR)
      frame = PlotFrame(pan)
      frame.setFontSizeForSlide(1.0,1.0)
      frame.setSize(1200,800)
      if name:
        frame.setTitle(name)
      frame.setVisible(True)
      if pngDir:
        frame.paintToPng(720,2.17,pngDir+name+'.png')
    g = slice12(k3,read("g"))
    h = slice12(k3,read("h"))
    f = slice12(k3,read("f"))
    cmin,cmax = 0.9*min(g),0.9*max(g)
    p(g,cmin=cmin,cmax=cmax,name="g")
    p(h,cmin=cmin,cmax=cmax,name="h")
    p(f,cmin=cmin,cmax=cmax,name="f")
  def plot3(f,cmap=gray,cmin=0.0,cmax=0.0,perc=100.0,
            t=None,cmap1=jet,cmin1=0.0,cmax1=0.0,perc1=0.0,
            nearest=False,cbar=None,name=None):
    o = PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT
    a = PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM
    panel = PlotPanelPixels3(o,a,s1,s2,s3,f)
    if t:
      s12 = slice12(k3,t)
      s13 = slice13(k2,t)
      s23 = slice23(k1,t)
      p12 = panel.addPixels(1,0,s1,s2,s12)
      p13 = panel.addPixels(1,1,s1,s3,s13)
      p23 = panel.addPixels(0,0,s3,s2,transpose(s23))
      p12.setInterpolation(PixelsView.Interpolation.NEAREST)
      p13.setInterpolation(PixelsView.Interpolation.NEAREST)
      p23.setInterpolation(PixelsView.Interpolation.NEAREST)
      if cmin1<cmax1:
        p12.setClips(cmin1,cmax1)
        p13.setClips(cmin1,cmax1)
        p23.setClips(cmin1,cmax1)
      p12.setColorModel(cmap1)
      p13.setColorModel(cmap1)
      p23.setColorModel(cmap1)
    if nearest:
      panel.setInterpolation(PixelsView.Interpolation.NEAREST)
    panel.setColorModel(cmap)
    panel.setSlices(k1,k2,k3)
    panel.setLabel1("Time (s)")
    panel.setLabel2("Inline (km)")
    panel.setLabel3("Crossline (km)")
    panel.setLineColor(Color.YELLOW)
    if (cbar):
      panel.addColorBar(cbar)
      panel.setColorBarWidthMinimum(90)
    panel.setBackground(Color.WHITE)
    panel.setInterval1(0.1)
    panel.setInterval2(0.5)
    panel.setInterval3(0.5)
    if cmin<cmax:
      panel.setClips(cmin,cmax)
    if perc<100.0:
      panel.setPercentiles(100-perc,perc)
    #panel.mosaic.setWidthElastic(1,100)
    #panel.mosaic.setHeightElastic(0,205)
    panel.mosaic.setHeightElastic(0,230)
    frame = PlotFrame(panel)
    frame.setFontSizeForSlide(1.0,1.0)
    #frame.setSize(1150,800)
    frame.setSize(1200,795)
    frame.setVisible(True)
    if name and pngDir:
      frame.paintToPng(720,7.0,pngDir+name+".png")
    return frame
  if plot3:
    goPlot3()
  if teaser:
    goTeaser()

#############################################################################
# Plots for SEG abstract
      
def goSeg():

  teaser = False
  plot3 = True

  def goPlot3():
    seismic = False
    blended = True
    shiftss = True
    if seismic:
      g = read("g")
      h = read("h")
      f = read("f")
      cmin,cmax = 0.9*min(g),0.9*max(g)
      plot3(g,cmin=cmin,cmax=cmax,cbar="Amplitude",name="g"+subset)
      plot3(h,cmin=cmin,cmax=cmax,cbar="Amplitude",name="h"+subset)
      plot3(f,cmin=cmin,cmax=cmax,cbar="Amplitude",name="f"+subset)
      #display(g,cmin=cmin,cmax=cmax);
    if blended:
      g,h,f = read('g'),read('h'),read('f')
      r = array(read('r1'),read('r2'),read('r3'))
      t1 = mul(d1*1000.0,read("t1"))
      q1 = mul(-d1*1000.0,read("q1"))
      ming,maxg = 0.9*min(g),0.9*max(g)
      
      null = -0.4938
      minq,_ = getClipsForPlot3(q1); maxq = 0.4*max(q1)
      knotch = int(256.0*abs((null-minq)/(maxq-minq)))
      amap = getNotchAlphaColorMap(jet,knotch)
      """
      """
      minq,maxq = getClipsForPlot3(q1,0.4)
      amap = getLinearAlphaColorMap(jet)

      plot3(g,cmin=ming,cmax=maxg,
            t=t1,cmap1=amap,cmin1=minq,cmax1=maxq,
            cbar="Vertical component of throw (ms)",name="t"+subset)

      plot3(h,cmin=ming,cmax=maxg,
            t=t1,cmap1=amap,cmin1=minq,cmax1=maxq,
            cbar="Vertical component of throw (ms)",name="h"+subset)
      #plot3(h,cmin=ming,cmax=maxg,cbar="Amplitude",name="h"+subset)

      #plot3(f,FlattenerUtil.applyShiftsR(t1,r),cmin=cmin,cmax=cmax,
      #      cbar="Vertical component of throw (ms)",name="t"+subset)
      plot3(f,cmin=ming,cmax=maxg,cbar="Amplitude",name="f"+subset)

      #plot3(q1,cmap=jet,cmin=minq,cmax=maxq,nearest=False,
      #      cbar="Vertical component of throw (ms)",name="q"+subset)

    if shiftss:
      s1 = mul(d1*1000.0,read('s1'))
      plot3(s1,cmap=jet,cbar="Vertical component of composite shift (ms)",
            name="s1"+subset)

  def goTeaser():
    def p(x,cmap=gray,cmin=0.0,cmax=0.0,name=None):
      pan = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
      pan.setHLabel("Crossline (km)")
      pan.setVLabel("Time (s)")
      #cb = pan.addColorBar("Amplitude")
      #cb.setWidthMinimum(150)
      pix = pan.addPixels(s1,s2,x)
      pix.setColorModel(cmap)
      if cmin<cmax:
        pix.setClips(cmin,cmax)
      pix.setInterpolation(PixelsView.Interpolation.LINEAR)
      frame = PlotFrame(pan)
      frame.setFontSizeForPrint(7.0,150)
      frame.setSize(1200,800)
      if name:
        frame.setTitle(name)
      frame.setVisible(True)
      if pngDir:
        frame.paintToPng(720,2.17,pngDir+name+'.png')
    g = slice12(k3,read("g"))
    h = slice12(k3,read("h"))
    f = slice12(k3,read("f"))
    cmin,cmax = 0.9*min(g),0.9*max(g)
    p(g,cmin=cmin,cmax=cmax,name="g")
    p(h,cmin=cmin,cmax=cmax,name="h")
    p(f,cmin=cmin,cmax=cmax,name="f")

  def plot3(f,cmap=gray,cmin=0.0,cmax=0.0,perc=100.0,
            t=None,cmap1=jet,cmin1=0.0,cmax1=0.0,perc1=0.0,
            nearest=False,cbar=None,name=None):
    o = PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT
    a = PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM
    panel = PlotPanelPixels3(o,a,s1,s2,s3,f)
    if t:
      s12 = slice12(k3,t)
      s13 = slice13(k2,t)
      s23 = slice23(k1,t)
      p12 = panel.addPixels(1,0,s1,s2,s12)
      p13 = panel.addPixels(1,1,s1,s3,s13)
      p23 = panel.addPixels(0,0,s3,s2,transpose(s23))
      p12.setInterpolation(PixelsView.Interpolation.NEAREST)
      p13.setInterpolation(PixelsView.Interpolation.NEAREST)
      p23.setInterpolation(PixelsView.Interpolation.NEAREST)
      if cmin1<cmax1:
        p12.setClips(cmin1,cmax1)
        p13.setClips(cmin1,cmax1)
        p23.setClips(cmin1,cmax1)
      p12.setColorModel(cmap1)
      p13.setColorModel(cmap1)
      p23.setColorModel(cmap1)
    if nearest:
      panel.setInterpolation(PixelsView.Interpolation.NEAREST)
    panel.setColorModel(cmap)
    panel.setSlices(k1,k2,k3)
    panel.setLabel1("Time (s)")
    panel.setLabel2("Inline (km)")
    panel.setLabel3("Crossline (km)")
    panel.setLineColor(Color.YELLOW)
    if (cbar):
      panel.addColorBar(cbar)
      panel.setColorBarWidthMinimum(125)
    panel.setBackground(Color.WHITE)
    panel.setInterval1(0.2)
    panel.setInterval2(1.0)
    panel.setInterval3(1.0)
    if cmin<cmax:
      panel.setClips(cmin,cmax)
    if perc<100.0:
      panel.setPercentiles(100-perc,perc)
    #panel.mosaic.setWidthElastic(1,100)
    #panel.mosaic.setHeightElastic(0,205)
    panel.mosaic.setHeightElastic(0,215)
    frame = PlotFrame(panel)
    frame.setFontSizeForSlide(1.0,1.0)
    #frame.setSize(1150,800)
    frame.setSize(1200,795)
    frame.setVisible(True)
    if name and pngDir:
      frame.paintToPng(720,4.0,pngDir+name+".png")
    return frame
  if plot3:
    goPlot3()
  if teaser:
    goTeaser()

#############################################################################

def goHorizon():
  g,r1,r2,r3 = read("g"),read("r1"),read("r2"),read("r3")
  #world = World(); addImageToWorld(world,g).setSlices(k1,k2,k3)
  for h1 in h1s:
    world = World(); addImageToWorld(world,g).setSlices(k1,k2,k3)
    xyz = makeVertices(h1,[r1,r2,r3])
    tg = QuadGroup(True,xyz)
    tg.setColor(Color.YELLOW)
    world.addChild(tg)
    makeFrame(world)

def makeVertices(k1,r):
  """Make vertices for QuadGroup."""
  # Coordinates: x = x(x)
  x1 = rampfloat(0.0,1.0,0.0,0.0,n1,n2,n3)
  x2 = rampfloat(0.0,0.0,1.0,0.0,n1,n2,n3)
  x3 = rampfloat(0.0,0.0,0.0,1.0,n1,n2,n3)
  x1 = Faulter.applyShiftsRLinear(x1,r) # x1(u) = x1(u-r(u))
  x2 = Faulter.applyShiftsRLinear(x2,r) # x2(u) = x2(u-r(u))
  x3 = Faulter.applyShiftsRLinear(x3,r) # x3(u) = x3(u-r(u))
  # Fault location: w = w(x)
  t1,t2,t3,z1,z2,z3 = getGriddedThrows()
  e = Faulter.getGriddedThrows(0.0,t1,t2,t3,z1,z2,z3)
  z1,z2,z3 = e[3],e[4],e[5] # live value coords
  w = zerofloat(n1,n2,n3)
  for i in range(len(z1)):
    i1,i2,i3 = int(z1[i]),int(z2[i]),int(z3[i])
    w[i3][i2][i1] = 1.0
  w = Faulter.applyShiftsR(w,r) # w(u) = w(u-r(u))
  # Slice
  x1 = slice23(k1,x1) # x1(u=k1)
  x2 = slice23(k1,x2) # x2(u=k1)
  x3 = slice23(k1,x3) # x3(u=k1)
  w = slice23(k1,w)   # w(u=k1)
  return FlattenerUtil.makeQuadVertices(x1,x2,x3,w)

def getGriddedThrows():
  null = -0.12345
  t1,t2,t3 = read("t1"),read("t2"),read("t3")
  s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  f1 = SimpleGridder3.getGriddedSamples(null,s1,s2,s3,t1)
  f2 = SimpleGridder3.getGriddedSamples(null,s1,s2,s3,t2)
  f3 = SimpleGridder3.getGriddedSamples(null,s1,s2,s3,t3)
  return mul(-1.0,f1[0]),mul(-1.0,f2[0]),mul(-1.0,f3[0]),f1[1],f1[2],f1[3]

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,cmap=gray,cmin=0.0,cmax=0.0,cbar=None,name=None):
  pan = panel(cbar)
  pix = pan.addPixels(x)
  pix.setColorModel(cmap)
  if cmin<cmax:
    pix.setClips(cmin,cmax)
  pix.setInterpolation(PixelsView.Interpolation.LINEAR)
  frame(pan,name)

def panel(cbar=None):
  p = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p.setHLabel("index i2")
  p.setVLabel("index i1")
  cb = p.addColorBar()
  if cbar:
    cb.setLabel(cbar)
  return p

def frame(panel,name=None):
  frame = PlotFrame(panel)
  #frame.setBackground(Color(204,204,204,255))
  frame.setFontSizeForSlide(1.0,1.0)
  frame.setSize(1200,800)
  if name:
    frame.setTitle(name)
  frame.setVisible(True)
  if name and pngDir:
    frame.paintToPng(720,widthinches,pngDir+name+'.png')

def getClipsForPlot3(f,r=1.0):
  s12 = slice12(k3,f)
  s13 = slice13(k2,f)
  s23 = slice23(k1,f)
  #cmin = 0.0
  cmin = r*min(min(s12),min(min(s13),min(s23)))
  cmax = r*max(max(s12),max(max(s13),max(s23)))
  return cmin,cmax

def getLinearAlphaColorMap(cmap):
  alpha = rampfloat(0.0,1.0/255.0,256)
  return ColorMap.setAlpha(cmap,alpha)

def getNotchAlphaColorMap(cmap,knotch):
  alpha = fillfloat(1.0,256)
  alpha[knotch] = 0.0
  if knotch>0:
    alpha[knotch-1] = 0.0
  if knotch<255:
    alpha[knotch+1] = 0.0
  return ColorMap.setAlpha(cmap,alpha)

def plot3(f,cmap=gray,cmin=0.0,cmax=0.0,perc=100.0,
          t=None,cmap1=jet,cmin1=0.0,cmax1=0.0,perc1=0.0,
          nearest=False,cbar=None,name=None):
  o = PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT
  a = PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM
  panel = PlotPanelPixels3(o,a,s1,s2,s3,f)
  if t:
    s12 = slice12(k3,t)
    s13 = slice13(k2,t)
    s23 = slice23(k1,t)
    p12 = panel.addPixels(1,0,s1,s2,s12)
    p13 = panel.addPixels(1,1,s1,s3,s13)
    p23 = panel.addPixels(0,0,s3,s2,transpose(s23))
    p12.setInterpolation(PixelsView.Interpolation.NEAREST)
    p13.setInterpolation(PixelsView.Interpolation.NEAREST)
    p23.setInterpolation(PixelsView.Interpolation.NEAREST)
    if cmin1<cmax1:
      p12.setClips(cmin1,cmax1)
      p13.setClips(cmin1,cmax1)
      p23.setClips(cmin1,cmax1)
    p12.setColorModel(cmap1)
    p13.setColorModel(cmap1)
    p23.setColorModel(cmap1)
  if nearest:
    panel.setInterpolation(PixelsView.Interpolation.NEAREST)
  panel.setColorModel(cmap)
  panel.setSlices(k1,k2,k3)
  panel.setLabel1("Time (s)")
  panel.setLabel2("Inline (km)")
  panel.setLabel3("Crossline (km)")
  panel.setLineColor(Color.YELLOW)
  if (cbar):
    panel.addColorBar(cbar)
    panel.setColorBarWidthMinimum(90)
  panel.setBackground(Color(255,254,255,255))
  panel.setInterval1(0.1)
  panel.setInterval2(0.5)
  panel.setInterval3(0.5)
  if cmin<cmax:
    panel.setClips(cmin,cmax)
  if perc<100.0:
    panel.setPercentiles(100-perc,perc)
  #panel.mosaic.setWidthElastic(1,100)
  #panel.mosaic.setHeightElastic(0,205)
  panel.mosaic.setHeightElastic(0,230)
  frame = PlotFrame(panel)
  frame.setFontSizeForSlide(1.0,1.0)
  #frame.setSize(1150,800)
  frame.setSize(1200,795)
  frame.setVisible(True)
  if name and pngDir:
    frame.paintToPng(720,7.0,pngDir+name+".png")
  return frame

def read(name,image=None,dir=None):
  if not image:
    image = zerofloat(n1,n2,n3)
  if not dir:
    dir = dataDir
  fileName = dir+name+".dat"
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def write(name,image):
  fileName = dataDir+name+'.dat'
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
    
def normal(x,y):
  div(x,max(abs(x)),y)

def transpose(x):
  n1,n2 = len(x[0]),len(x)
  y = zerofloat(n2,n1)
  for i1 in range(n1):
    for i2 in range(n2):
      y[i1][i2] = x[i2][i1]
  return y

def transpose23(x):
  m1,m2,m3 = len(x[0][0]),len(x[0]),len(x)
  y = zeros(m1,m3,m2)
  for i3 in range(m3):
    for i2 in range(m2):
      copy(x[i3][i2],y[i2][i3])
  return y

def array(x1,x2,x3=None,x4=None):
  if x3 and x4:
    return jarray.array([x1,x2,x3,x4],Class.forName('[[[F'))
  elif x3:
    return jarray.array([x1,x2,x3],Class.forName('[[[F'))
  else:
    return jarray.array([x1,x2],Class.forName('[[[F'))

def zeros(n1,n2,n3):
  return zerofloat(n1,n2,n3)

#############################################################################
# graphics

def display(image,cmap=gray,cmin=0,cmax=0,perc=100,name=None):
  world = World()
  addImageToWorld(world,image,cmap,cmin,cmax,perc)
  makeFrame(world,name)

def display2(image1,image2,cmap1=gray,cmap2=gray,name=None):
  world = World()
  addImageToWorld(world,image1,cmap1)
  addImageToWorld(world,image2,cmap2)
  makeFrame(world,name)

def addImageToWorld(world,image,cmap=gray,cmin=0,cmax=0,perc=100):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  ipg.setColorModel(cmap)
  ipg.setSlices(k1,k2,k3)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
  if perc<100:
    ipg.setPercentiles(100-perc,perc)
  world.addChild(ipg)
  return ipg

def addImage2ToWorld(world,image1,image2):
  ipg = ImagePanelGroup2(s1,s2,s3,image1,image2)
  ipg.setColorModel1(ColorMap.getGray())
  #ipg.setColorModel2(ColorMap.getJet(0.3))
  ipg.setColorModel2(ColorMap.getHue(0.0,20.0,0.3))
  world.addChild(ipg)
  return ipg

def addColorBar(frame,label):
  cbar = ColorBar(label)
  cbar.setFont(cbar.getFont().deriveFont(64.0))
  frame.add(cbar,BorderLayout.EAST)
  #frame.viewCanvas.setBackground(frame.getBackground())
  return cbar

def makeFrame(world,name=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  frame = SimpleFrame(world)
  #frame.setBackground(Color(204,204,204,255))
  frame.setBackground(Color.WHITE)
  if name:
    frame.setTitle(name)
  view = frame.getOrbitView()
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.55)
  view.setAzimuth(azimuth)
  view.setElevation(elevation)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3-1.0,f2,f1,l3,l2,l1)))
  frame.viewCanvas.setBackground(frame.getBackground())
  #frame.setSize(1020,750)
  frame.setSize(1800,1200)
  frame.setVisible(True)
  return frame

def slice12(k3,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n2)
  SimpleFloat3(f).get12(n1,n2,0,0,k3,s)
  return s

def slice13(k2,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n3)
  SimpleFloat3(f).get13(n1,n3,0,k2,0,s)
  return s

def slice23(k1,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n2,n3)
  SimpleFloat3(f).get23(n2,n3,k1,0,0,s)
  return s

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
    print '%02d:%02d:%02ds'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
