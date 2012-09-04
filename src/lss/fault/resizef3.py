"""
Resize Dave's F3 subset
"""
from imports import *
from java.awt.image import IndexColorModel

readDir = "/Users/sluo/Dropbox/Share/"
writeDir = "/data/seis/f3/subb/"
n1,n2,n3 = 120,221,220

#############################################################################

def main(args):
  g,t1,t2,t3 = read('g'),read('t1'),read('t2'),read('t3')
  g = copy(90,n2,n3,15,0,0,g)
  t1 = copy(90,n2,n3,15,0,0,t1)
  t2 = copy(90,n2,n3,15,0,0,t2)
  t3 = copy(90,n2,n3,15,0,0,t3)
  write('g',g)
  write('t1',t1)
  write('t2',t2)
  write('t3',t3)

def read(name,image=None):
  if not image:
    image = zerofloat(n1,n2,n3)
  fileName = readDir+'b'+name+".dat"
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def write(name,image):
  fileName = writeDir+name+'.dat'
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
    print '%02d:%02d:%02ds'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
