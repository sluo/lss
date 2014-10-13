##############################################################################
# Figures for Eni data

from imports import *

#############################################################################

#sx = Sampling(1921,0.00625,0.0)
#sz = Sampling(181,0.00625,0.0)
sx = Sampling(1751,0.00625,0.7)
sz = Sampling(145,0.00625,0.2)
st = Sampling(3751,0.0004,0.0)
sr = Sampling(197,0.00625,-1.225)
nz,nx,nt,nr = sz.count,sx.count,st.count,sr.count
dz,dx,dt,dr = sz.delta,sx.delta,st.delta,sr.delta
fz,fx,ft,fr = sz.first,sx.first,st.first,sr.first

home = os.getenv('HOME')
#dresDir = home+'/Desktop/subg/nonlinear/vfile/dres00/'
#aresDir = home+'/Desktop/subg/nonlinear/vfile/ares00/'
dresDir = home+'/Desktop/subg/nonlinear/vz/dres00/'
#aresDir = home+'/Desktop/subg/nonlinear/vz/ares00/'
#aresDir = home+'/Desktop/subg/nonlinear/vz/ares01/'
#aresDir = home+'/Desktop/subg/nonlinear/vz/ares02/'
#aresDir = home+'/Desktop/subg/nonlinear/vz/ares03/'
#aresDir = home+'/Desktop/subg/nonlinear/vz/ares04/'
aresDir = home+'/Desktop/subg/nonlinear/vz/ares05/' # SD, shift obs (best)
#aresDir = home+'/Desktop/subg/nonlinear/vz/ares05a/' # CG, shift pre
#aresDir = home+'/Desktop/subg/nonlinear/vz/ares05b/' # CG, shift obs
#aresDir = home+'/Desktop/subg/nonlinear/vz/ares05c/' # SD, shift pre
#aresDir = home+'/Desktop/subg/nonlinear/vz/ares06/'
#aresDir = home+'/Desktop/pngdat3/'
#aresDir = home+'/Desktop/pngdat2/'

savDir = None
#savDir = '/Users/sluo/Desktop/png/'
#savDir = '/Users/sluo/Dropbox/png/'

plotVelocity = True # convert slowness to velocity

widthPoints = None # slides
#widthPoints = 190.0 # 1/3 page
#widthPoints = 225.0 # 1 column
#widthPoints = 470.0 # 2 column (full page)

#rclip = 0.015 # clip for migration images
rclip = 0.020 # clip for migration images ###
#rclip = 0.025 # clip for migration images

#cmapVel = ColorMap.GRAY # colormap for slowness model
cmapVel = ColorMap.JET # colormap for slowness model
cmapMig = ColorMap.GRAY # colormap for migration images
#cmapMig = ColorMap.GRAY_YELLOW_RED # colormap for migration images

#############################################################################

def main(args):
  plotImages()
  #plotImageSubsets()
  plotSlownesses()
  #plotData() # ares05 only

def plotData():
  #x = '1.225000'
  #x = '1.850000' # (shown in CWP report)
  #x = '2.475000'
  #x = '3.100000'
  x = '5.550000' # (smallest shifts)
  """
  do = zerofloat(nt,nr)
  dp = zerofloat(nt,nr)
  dw = zerofloat(nt,nr)
  w = zerofloat(nt,nr)
  readImage(aresDir+'data10/do_x='+x+'.dat',do)
  readImage(aresDir+'data10/dp_x='+x+'.dat',dp)
  readImage(aresDir+'data10/dw_x='+x+'.dat',dw)
  readImage(aresDir+'data10/w_x='+x+'.dat',w)
  mul(1000.0*dt,w,w) # convert shift in samples to ms
  """
  rdir = home+'/Desktop/subg/nonlinear/vz/ares05/warp3d/'
  do = zerofloat(nt,nr,432); readImage(rdir+'do.dat',do)
  dp = zerofloat(nt,nr,432); readImage(rdir+'dp.dat',dp)
  dw = zerofloat(nt,nr,432); readImage(rdir+'dw.dat',dw)
  w = zerofloat(nt,nr,432); readImage(rdir+'w.dat',w); mul(1000.0*dt,w,w)
  ks = int((float(x)-1.225)/0.025)
  do = do[ks]; dp = dp[ks]; dw = dw[ks]; w = w[ks];
  clip = 8.0
  tmin,tmax = 0.3,1.4
  plot(do,cmin=-clip,cmax=clip,vmin=tmin,vmax=tmax,
    cbar='Amplitude',title='do')
  plot(dw,cmin=-clip,cmax=clip,vmin=tmin,vmax=tmax,
    cbar='Amplitude',title='dw')
  plot(dp,cmin=-clip,cmax=clip,vmin=tmin,vmax=tmax,
    cbar='Amplitude',title='dp')
  plot(w,cmap=jet,vmin=tmin,vmax=tmax,cint=5.0,
    cbar='Time shift (ms)',title='w',
    cmin=-17.5,cmax=17.5)
    #sperc=100.0)

def plotSlownesses():
  vtFile = '/Users/sluo/Desktop/subg/nonlinear/vfile/dres00/s.dat'
  vzFile = '/Users/sluo/Desktop/subg/nonlinear/vz/ares00/s.dat'
  vt = readImage(vtFile)
  vz = readImage(vzFile)
  cbar,cmap,clip,cint = 'Slowness (s/km)',rwb,0.04,0.04
  if plotVelocity:
    cbar,cmap,clip,cint = 'Velocity (km/s)',jet,0.3,None
    div(1.0,vt,vt)
    div(1.0,vz,vz)
  vd = sub(vz,vt)
  cmin = min(min(vt),min(vz))
  cmax = max(max(vt),max(vz))
  plot(vt,cmap=cmapVel,cmin=cmin,cmax=cmax,
    cbar=cbar,title='s_true')
  plot(vz,cmap=cmapVel,cmin=cmin,cmax=cmax,
    cbar=cbar,title='s_z')
  plot(vd,cmap=cmap,cmin=-clip,cmax=clip,cint=cint,title='s_z-s_true',
    cbar='Velocity difference (km/s)' if plotVelocity\
      else 'Slowness difference (s/km)')

def getImages():
  s = readImage(dresDir+'s.dat')
  #ra = readImage(aresDir+'r1.dat')
  #rb = readImage(dresDir+'r1.dat')
  ra = readImage(aresDir+'r9.dat')
  rb = readImage(dresDir+'r9.dat')
  return s,ra,rb

def plotImages():
  s,ra,rb = getImages()
  #plot(ra,cmap=cmapMig,sperc=99.8,title='ra')
  #plot(rb,cmap=cmapMig,sperc=99.8,title='rb')
  #plot(ra,cmap=cmapMig,cmin=-rclip,cmax=rclip,title='ra')
  #plot(rb,cmap=cmapMig,cmin=-rclip,cmax=rclip,title='rb')
  plot(ra,cmap=cmapMig,cmin=-rclip,cmax=rclip,cbar='Reflectivity',title='ra')
  plot(rb,cmap=cmapMig,cmin=-rclip,cmax=rclip,cbar='Reflectivity',title='rb')
  plot(s,cmap=cmapVel,cbar='Slowness (s/km)',title='s')
  plot(div(1.0,s),cmap=cmapVel,cbar='Velocity (km/s)',title='s')

def plotImageSubsets():
  xmin,xmax = 0.5,3.0
  zmin,zmax = 0.5,0.8
  s,ra,rb = getImages()
  nxm = 1+int((xmax-xmin)/dx)
  nzm = 1+int((zmax-zmin)/dz)
  sxm = Sampling(nxm,dx,xmin)
  szm = Sampling(nzm,dz,zmin)
  ca = copy(nzm,nxm,int((zmin-fz)/dz),int((xmin)/dx),ra)
  cb = copy(nzm,nxm,int((zmin-fz)/dz),int((xmin)/dx),rb)
  def subPlot(f,title):
    widthPoints = None # slides
    #widthPoints = 190.0 # 1/3 page
    #widthPoints = 225.0 # single column
    #widthPoints = 280.0 # half page (for 2 column figure)
    panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
    panel.setHLabel('Distance (km)')
    panel.setVLabel('Depth (km)')
    panel.setVInterval(0.1)
    pixel = panel.addPixels(szm,sxm,f)
    pixel.setColorModel(cmapMig)
    pixel.setClips(-rclip,rclip)
    frame = PlotFrame(panel)
    if widthPoints is None:
      #frame.setFontSizeForSlide(0.4,1.0,16.0/9.0)
      #frame.setFontSizeForSlide(1.0,0.8,16.0/9.0)
      frame.setFontSize(54.0)
    else:
      frame.setFontSizeForPrint(8.0,widthPoints)
    if widthPoints is None:
      #frame.setSize(1115,895)
      frame.setSize(985,800)
    elif widthPoints==225.0:
      frame.setSize(1115,735)
    elif widthPoints==190.0:
      frame.setSize(1115,740)
    else:
      frame.setSize(1115,720)
    frame.setTitle(title)
    frame.setVisible(True)
    if title and savDir:
      if widthPoints is None:
        frame.paintToPng(1020.0,1.0,savDir+title+'.png')
      else:
        frame.paintToPng(720.0,widthPoints/72.0,savDir+title+'.png')
  subPlot(ca,'ca')
  subPlot(cb,'cb')

#############################################################################

gray = ColorMap.GRAY
gyr = ColorMap.GRAY_YELLOW_RED
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(f,cmap=gray,cmin=0,cmax=0,perc=100,sperc=None,cbar=None,cwidth=None,
  cint=None,hint=None,vint=None,vmin=None,vmax=None,grid=False,title=None):
  panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  if cbar:
    cb = panel.addColorBar(cbar)
    cb.setWidthMinimum(150)
  if cbar and cwidth:
    cb.setWidthMinimum(cwidth)
  elif cbar and widthPoints is None:
    if len(f[0])==nz:
      cb.setWidthMinimum(190)
    else:
      cb.setWidthMinimum(160)
  elif cbar and widthPoints==225.0:
    cb.setWidthMinimum(130)
  elif cbar and widthPoints==190.0:
    cb.setWidthMinimum(180)
  if cint:
    cb.setInterval(cint)
  if hint:
    panel.setHInterval(hint)
  if vint:
    panel.setVInterval(vint)
  if len(f[0])==nz and len(f)==nx:
    #pixel = panel.addPixels(sz,sx,f)
    pixel = panel.addPixels(sz,Sampling(nx,0.00625,0.0),f)
    panel.setHLabel('Distance (km)')
    panel.setVLabel('Depth (km)')
  elif len(f[0])==nt and len(f)==nx:
    pixel = panel.addPixels(st,sx,f)
    panel.setHLabel('Distance (km)')
    panel.setVLabel('Time (s)')
  elif len(f[0])==nt and len(f)==nr:
    pixel = panel.addPixels(st,sr,f)
    panel.setHLabel('Offset (km)')
    panel.setVLabel('Time (s)')
  else:
    pixel = panel.addPixels(f)
  if vmin is not None and vmax is not None:
    panel.setVLimits(vmin,vmax)
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
  if len(f[0])==nz:
    if widthPoints is None:
      frame.setFontSizeForSlide(1.0,1.0,16.0/9.0)
      frame.setSize(1600,800)
    else:
      frame.setFontSizeForPrint(8.0,widthPoints)
      if widthPoints==225.0:
        frame.setSize(1115,476)
        frame.setFontSizeForPrint(7.0,widthPoints)
      elif widthPoints==240.0 or widthPoints==175.0:
        frame.setSize(1200,464)
      elif widthPoints==260.0:
        frame.setSize(1200,465)
      elif widthPoints==470.0:
        frame.setSize(1100,550)
  else:
    if widthPoints is None:
      frame.setSize(980,800)
      #frame.setFontSizeForSlide(1.0,1.0,16.0/9.0)
      frame.setFontSize(45.0)
    else:
      frame.setSize(1200,900)
      frame.setFontSizeForPrint(8.0,widthPoints)
  if title:
    frame.setTitle(title)
  frame.setVisible(True)
  if title and savDir:
    if widthPoints is None:
      #frame.paintToPng(1024,2.0,savDir+title+'.png')
      frame.paintToPng(1920,1.0,savDir+title+'.png')
    else:
      frame.paintToPng(720.0,widthPoints/72.0,savDir+title+'.png')
      #frame.paintToPng(3000.0,widthPoints/72.0,savDir+title+'.png') # poster
  return panel

def readImage(fname,image=None):
  imageIsNone = image is None
  if imageIsNone:
    image = zerofloat(181,1921)
  ais = ArrayInputStream(fname)
  ais.readFloats(image)
  ais.close()
  if (imageIsNone) and (nz!=181 or nx!=1921):
    return copy(nz,nx,int(fz/dz),int(fx/dx),image)
  else:
    return image

def cleanDir(dir):
  os.chdir(dir)
  for f in os.listdir(dir):
    os.remove(f)

#############################################################################
# Do everything on Swing thread.
class RunMain(Runnable):
  def run(self):
    if savDir is not None:
      print 'cleaning savDir'
      cleanDir(savDir)
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
