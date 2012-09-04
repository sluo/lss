package lss.dev;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

// TESTING
import java.awt.*;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.mosaic.*;

/**
 * Interpolation of missing data by projection onto convex sets.
 * @author Simon Luo
 * @version 2012.02.15
 */
public class PocsInterpolator {

  public PocsInterpolator(int n1, int n2) {
    _n1 = n1;
    _n2 = n2;
    int nfft1 = FftReal.nfftFast(n1*PAD1);
    int nfft2 = FftComplex.nfftFast(n2*PAD2);
    _nf1 = nfft1/2+1;
    _nf2 = nfft2;
    _fr = new FftReal(nfft1);
    _fc = new FftComplex(nfft2);
  }

  public float[][] interpolate(
    float fnull, float[][] nulls, float[][] f)
  {
    int n1c = 2*_nf1;
    int n2c = _nf2;
    float[][] t = new float[n2c][n1c];
    copy(f,t);

    float[][] z = copy(t);
    forwardFft(t,z);
    float cmin = 0.5f*min(z);
    float cmax = 0.5f*max(z);

    for (float thresh=S; thresh>0.01f; thresh*=R) {
      forwardFft(t,t);
      plot(t,cmin,cmax,"before_"+thresh);
      cthreshold(thresh,t);
      plot(t,cmin,cmax,"after_"+thresh);
      inverseFft(t,t);
      for (int i2=0; i2<n2c; ++i2) {
        for (int i1=0; i1<n1c; ++i1) {
          if (i2>=_n2 || i1>=_n1) {
            t[i2][i1] = 0.0f;
          } else if (!(nulls[i2][i1]==fnull)) {
            t[i2][i1] = f[i2][i1]; // replace known values
          }
        }
      }
      System.out.println("threshold="+thresh);
    }
    return copy(_n1,_n2,0,0,t);
  }

  // threshold fk spectrum
  private static void cthreshold(float r, float[][] c) {
    float[][] cab = cabs(c);
    int n1 = cab[0].length;
    int n2 = cab.length;
    float cmin = r*max(cab);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (cab[i2][i1]<cmin) {
          c[i2][2*i1  ] = 0.0f;
          c[i2][2*i1+1] = 0.0f;
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////

  private static final float S = 0.5f; // starting threshold
  private static final float R = 0.5f; // rate of threshold decrease;
  private static final int PAD1 = 2;   // fft padding in 1st dimension
  private static final int PAD2 = 1;   // fft padding in 2nd dimension

  private int _n1;
  private int _n2;
  private int _nf1;
  private int _nf2;
  private float[][] _f;
  private boolean[][] _null;
  private FftReal _fr;
  private FftComplex _fc;

  private float[][] forwardFft(float[][] x, float[][] y) {
    _fr.realToComplex1(1,_nf2,x,y);
    _fc.complexToComplex2(1,_nf1,y,y);
    return y;
  }

  private float[][] inverseFft(float[][] x, float[][] y) {
    _fc.complexToComplex2(-1,_nf1,x,y);
    _fr.complexToReal1(-1,_nf2,y,y);
    _fc.scale(_nf1,_nf2,y);
    _fr.scale(_nf1,_nf2,y);
    return y;
  }

  //////////////////////////////////////////////////////////////////////////
  // TESTING

  private static final boolean PNG = false;

  private static void plot(
    float[][] x, float cmin, float cmax, String name)
  {
    PlotPanel pan = panel();
    PixelsView pix = pan.addPixels(x);
    //pix.setColorModel(ColorMap.JET);
    pix.setColorModel(ColorMap.RED_WHITE_BLUE);
    pix.setInterpolation(PixelsView.Interpolation.NEAREST);
    if (cmin<cmax)
      pix.setClips(cmin,cmax);
    frame(pan,name);
  }

  private static PlotPanel panel() {
    PlotPanel p = new PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT);
    ColorBar cb = p.addColorBar();
    cb.setWidthMinimum(80);
    return p;
  }

  private static void frame(PlotPanel panel, String name) {
    PlotFrame frame = new PlotFrame(panel);
    frame.setBackground(new Color(204,204,204,255));
    frame.setFontSizeForSlide(1.5,1.5);
    frame.setSize(1024,700);
    frame.setTitle(name);
    frame.setVisible(true);
    if (PNG)
      frame.paintToPng(360,3.0,"/Users/sluo/Desktop/"+name+".png");
  }

}
