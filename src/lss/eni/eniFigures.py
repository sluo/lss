##############################################################################
# Figures for Eni data

from imports import *

#############################################################################

#sx = Sampling(3201,0.00625,0.0)
#sz = Sampling(181,0.00625,0.0)
sx = Sampling(2881,0.00625,1.0)
sz = Sampling(165,0.00625,0.1)
st = Sampling(3751,0.0004,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
fz,fx,ft = sz.first,sx.first,st.first

savDir = None
savDir = '/Users/sluo/Desktop/png/'

#widthPoints = None # slides
widthPoints = 225.0 # 1 column
#widthPoints = 470.0 # 2 column (full page)

rclip = 0.015 # clip for migration images
#rclip = 0.020 # clip for migration images

#cmapVel = ColorMap.GRAY # colormap for slowness model
cmapVel = ColorMap.JET # colormap for slowness model
#cmapMig = ColorMap.GRAY # colormap for migration images
cmapMig = ColorMap.GRAY_YELLOW_RED # colormap for migration images

#############################################################################

def main(args):
  vz = True
  plotFiles(vz)
  plotSubset(vz)

def readFiles(vz):
  if vz:
    dira,dirb = 'iter005vz/','iter525vz/'
  else:
    dira,dirb = 'iter005/','iter525/'
  s = read('/Users/sluo/Dropbox/save/eni/subc/'+dira+'s.dat');
  ra = read('/Users/sluo/Dropbox/save/eni/subc/'+dira+'r0.dat');
  rb = read('/Users/sluo/Dropbox/save/eni/subc/'+dirb+'r5.dat');
  return s,ra,rb

def plotFiles(vz):
  s,ra,rb = readFiles(vz)
  plot(ra,cmap=cmapMig,cmin=-rclip,cmax=rclip,title='ra')
  plot(rb,cmap=cmapMig,cmin=-rclip,cmax=rclip,title='rb')
  plot(s,cmap=cmapVel,cbar='Slowness (s/km)',title='s')
  #plot(div(1.0,s),cmap=cmapVel,cbar='Velocity (km/s)',title='s')

def plotSubset(vz):
  s,ra,rb = readFiles(vz)
  xmin,xmax = 6.5,10.5
  zmin,zmax = 0.7,1.0
  nxm = 1+int((xmax-xmin)/dx)
  nzm = 1+int((zmax-zmin)/dz)
  sxm = Sampling(nxm,dx,xmin)
  szm = Sampling(nzm,dz,zmin)
  ca = copy(nzm,nxm,int((zmin-fz)/dz),int((xmin)/dx),ra)
  cb = copy(nzm,nxm,int((zmin-fz)/dz),int((xmin)/dx),rb)
  def subPlot(f,title):
    widthPoints = 225.0 # single column
    #widthPoints = 280.0 # half page (for 2 column figure)
    panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
    panel.setHLabel('Distance (km)')
    panel.setVLabel('Depth (km)')
    pixel = panel.addPixels(szm,sxm,f)
    pixel.setColorModel(cmapMig)
    pixel.setClips(-rclip,rclip)
    frame = PlotFrame(panel)
    frame.setFontSizeForPrint(8.0,widthPoints)
    #frame.setSize(1115,475)
    frame.setSize(1115,580)
    frame.setTitle(title)
    frame.setVisible(True)
    if title and savDir:
      frame.paintToPng(720.0,widthPoints/72.0,savDir+title+'.png')
  subPlot(ca,'ca')
  subPlot(cb,'cb')

#def plotSubset(f,zmin,zmax,xmin,xmax,width,height,cmin,cmax,title):
#  g = copy(n1,n2,j1,j2,f)
#  panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
#  panel.setHLabel('Distance (km)')
#  panel.setVLabel('Depth (km)')
#
#
#  widthPoints = 260.0
#  j1,n1 = int(zmin/dz),int((zmax-zmin)/dz)
#  j2,n2 = int(xmin/dx),int((xmax-xmin)/dx)
#  s1 = Sampling(n1,dz,zmin)
#  s2 = Sampling(n2,dx,xmin)
#  g = copy(n1,n2,j1,j2,f)
#  panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
#  panel.setHLabel('Distance (km)')
#  panel.setVLabel('Depth (km)')
#  panel.setHInterval(0.5)
#  panel.setVInterval(0.5)
#  cb = panel.addColorBar('Reflectivity')
#  cb.setWidthMinimum(150)
#  cb.setInterval(0.05)
#  pixel = panel.addPixels(s1,s2,g)
#  pixel.setColorModel(gray)
#  pixel.setClips(cmin,cmax)
#  pixel.setInterpolation(PixelsView.Interpolation.LINEAR)
#  frame = PlotFrame(panel)
#  frame.setFontSizeForPrint(8.0,widthPoints)
#  #frame.setSize(1024,670)
#  frame.setSize(width,height)
#  frame.setTitle(title)
#  frame.setVisible(True)
#  if title and savDir:
#    frame.paintToPng(720.0,widthPoints/72.0,savDir+title+'.png')

#############################################################################

gray = ColorMap.GRAY
gyr = ColorMap.GRAY_YELLOW_RED
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(f,cmap=gray,cmin=0,cmax=0,perc=100,sperc=None,cbar=None,cwidth=None,
  cint=None,hint=None,vint=None,grid=False,title=None):
  panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  if cbar:
    cb = panel.addColorBar(cbar)
  if cint:
    cb.setInterval(cint)
  if cwidth:
    cb.setWidthMinimum(cwidth)
  elif widthPoints is None:
    cb.setWidthMinimum(140)
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
    #pixel = panel.addPixels(sz,sx,f)
    pixel = panel.addPixels(sz,Sampling(nx,0.00625,0.0),f)
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
    frame.setFontSizeForSlide(1.0,1.0)
  else:
    frame.setFontSizeForPrint(8.0,widthPoints)
  if (len(f[0])==nz):
    if widthPoints is None:
      #frame.setSize(1200,770)
      frame.setSize(1024,670)
    elif widthPoints==225.0:
      frame.setSize(1115,440)
    elif widthPoints==240.0 or widthPoints==175.0:
      #frame.setSize(1200,468)
      frame.setSize(1200,464)
    elif widthPoints==260.0:
      frame.setSize(1200,465)
    elif widthPoints==470.0:
      #frame.setSize(1200,445)
      frame.setSize(1200,475)
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
  if title and savDir:
    if widthPoints is None:
      frame.paintToPng(1024,2.0,savDir+title+'.png')
    else:
      frame.paintToPng(720.0,widthPoints/72.0,savDir+title+'.png')
      #frame.paintToPng(3000.0,widthPoints/72.0,savDir+title+'.png') # poster
  return panel

def read(name,image=None):
  if not image:
    image = zerofloat(181,3201)
  fileName = name
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  if nz!=181 or nx!=3201:
    return copy(nz,nx,int(fz/dz),int(fx/dx),image)
  else:
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
class RunMain(Runnable):
  def run(self):
    if savDir is not None:
      print 'cleaning savDir'
      cleanDir(savDir)
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
