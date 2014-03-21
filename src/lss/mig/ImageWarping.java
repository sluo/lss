package lss.mig;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import lss.mod.Receiver;

public class ImageWarping {

  public ImageWarping(
    int bstrain1, int bstrain2, int bstrain3,
    double smooth1, double smooth2, double smooth3,
    double maxShift, double d1)
  {
    this(
      min(1.0,1.0/(bstrain1-0.5)),
      min(1.0,1.0/(bstrain2-0.5)),
      min(1.0,1.0/(bstrain3-0.5)),
      smooth1,smooth2,smooth3,maxShift,d1
    );
  }

  public ImageWarping(
    double strain1, double strain2, double strain3,
    double smooth1, double smooth2, double smooth3,
    double maxShift, double d1)
  {
    this(strain1,strain2,strain3,smooth1,smooth2,smooth3,maxShift,d1,1);
  }

  public ImageWarping(
    double strain1, double strain2, double strain3,
    double smooth1, double smooth2, double smooth3,
    double maxShift, double d1, int t1)
  {
    Check.argument(t1>=1,"t1>=1");
    int shiftMax = (int)(maxShift/(t1*d1));
    _warp = new DynamicWarping(-shiftMax,shiftMax);
    _warp.setShiftSmoothing(smooth1/t1,smooth2,smooth3); // shift smoothing
    _warp.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.AVERAGE);
    _warp.setErrorSmoothing(2); // number of smoothings of alignment errors
    if (strain3<=0.0) {
      _warp.setStrainMax(strain1,strain2);
      _3d = false;
    } else {
      _warp.setStrainMax(strain1,strain2,strain3);
      _3d = true;
    }
    _t1 = t1;
    _sigmaRms = 1.0*maxShift/(t1*d1);
  }

  public Receiver[] warp(Receiver[] rp, Receiver[] ro) {
    float[][] d = rp[0].getData();
    int nt = d[0].length;
    int nr = d.length;
    int ns = rp.length;
    float[][][] u = new float[ns][nr][nt];
    return warp(rp,ro,u);
  }

  public Receiver[] warp(
    final Receiver[] rp, final Receiver[] ro, final float[][][] u)
  {
    final int nt = u[0][0].length;
    final int nr = u[0].length;
    final int ns = u.length;
    final int ntm = nt/_t1;
    float[][][] ep = new float[ns][][];
    float[][][] eo = new float[ns][][];
    float[][][] uu = new float[ns][nr][ntm];
    for (int is=0; is<ns; ++is) {
      ep[is] = copy(ntm,nr,0,0,_t1,1,rp[is].getData());
      eo[is] = copy(ntm,nr,0,0,_t1,1,ro[is].getData());
    }
    filterAndFindShifts(ep,eo,uu);
    if (_t1>1) { // if decimated, interpolate shifts
      mul((float)_t1,uu,uu); // scale shifts to compensate for decimation
      final LinearInterpolator li = new LinearInterpolator();
      li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
      li.setUniform(ntm,_t1,0.0,nr,1.0,0.0,ns,1.0,0.0,uu);
      Parallel.loop(ns,new Parallel.LoopInt() {
      public void compute(int is) {
        for (int ir=0; ir<nr; ++ir) {
          for (int it=0; it<nt; ++it) {
            u[is][ir][it] = li.interpolate(it,ir,is);
          }
        }
      }});
    }
    final Receiver[] rw = new Receiver[ns];
    Parallel.loop(ns,new Parallel.LoopInt() {
    public void compute(int is) {
      rw[is] = new Receiver(ro[is].getXIndices(),ro[is].getZIndices(),nt);
      _warp.applyShifts(u[is],ro[is].getData(),rw[is].getData());
    }});
    return rw;
  }

  public void findShifts(
    final float[][][] f, final float[][][] g, final float[][][] u)
  {
    zero(u);
    if (_3d) {
      _warp.findShifts(f,g,u);
    } else {
      int ns = f.length;
      Parallel.loop(ns,new Parallel.LoopInt() {
      public void compute(int is) {
        findShifts(f[is],g[is],u[is]);
      }});
    }
  }

  public void findShifts(float[][] f, float[][] g, float[][] u) {
    zero(u);
    _warp.findShifts(f,g,u);
  }

  public float[][] findShifts(float[][] f, float[][] g) {
    return _warp.findShifts(f,g);
  }

  public float[][][] applyShifts(float[][][] u, float[][][] g) {
    return _warp.applyShifts(u,g);
  }

  public void applyShifts(float[][][] u, float[][][] g, float[][][] h) {
    _warp.applyShifts(u,g,h);
  }

  public float[][] applyShifts(float[][] u, float[][] g) {
    return _warp.applyShifts(u,g);
  }

  public void applyShifts(float[][] u, float[][] g, float[][] h) {
    _warp.applyShifts(u,g,h);
  }

  public void setWindowSizeAndOverlap(int l2, int l3, double f2, double f3) {
    _warp.setWindowSizeAndOverlap(l2,l3,f2,f3);
  }

  public void setNoise(boolean noise) {
    _noise = noise;
  }

  //////////////////////////////////////////////////////////////////////////

  private int _t1; // decimate along in direction 1, e.g., time
  private boolean _3d; // true for 3D warping
  private double _sigmaRms; // sigma for RMS filtering
  private DynamicWarping _warp;
  private boolean _noise = true; // add random noise before finding shifts

  // public for now
  public void filterAndFindShifts(
    final float[][][] rp, final float[][][] ro, final float[][][] u)
  {
    rmsFilter(rp,ro);
    if (_noise)
      addRandomNoise(10.0f,rp,ro);
    findShifts(rp,ro,u);
  }

  private void addRandomNoise(float snr, float[][][] x, float[][][] y) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    //float xrms = sqrt(sum(mul(x,x))/n1/n2/n3);
    //float yrms = sqrt(sum(mul(y,y))/n1/n2/n3);
    float xrms = rms(x);
    float yrms = rms(y);
    float arms = 0.5f*(xrms+yrms); // average rms of signal
    java.util.Random random = new java.util.Random(12345);
    float[][][] s = sub(randfloat(random,n1,n2,n3),0.5f);
    new RecursiveGaussianFilter(1.0).apply000(s,s); // bandlimited noise
    float srms = rms(s); // rms of noise
    mul(arms/(srms*snr),s,s);
    add(s,x,x);
    add(s,y,y);
  }

  private void rmsFilter(float[][][] x, float[][][] y) {
    int ns = x.length;
    plot(x[ns/2],"x before");
    plot(y[ns/2],"y before");

    equalize(x,y);
    float[][][] xx = abs(x);
    float[][][] yy = abs(y);

    /*
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(_sigmaRms);
    if (_3d) {
      rgf.apply000(xx,xx);
      rgf.apply000(yy,yy);
    } else {
      rgf.apply0XX(xx,xx);
      rgf.applyX0X(xx,xx);
      rgf.apply0XX(yy,yy);
      rgf.applyX0X(yy,yy);
    }
    */
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(_sigmaRms);
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE);
    if (_3d) {
      ref.apply(xx,xx);
      ref.apply(yy,yy);
    } else {
      ref.apply2(xx,xx);
      ref.apply1(xx,xx);
      ref.apply2(yy,yy);
      ref.apply1(yy,yy);
    }

    plot(xx[ns/2],"xx smoothed");
    plot(yy[ns/2],"yy smoothed");

    // numerator = 2*xx*yy
    float[][][] num = mul(2.0f,xx);
    mul(yy,num,num);
    plot(num[ns/2],"num");

    // denominator = xx^2+yy^2
    mul(xx,xx,xx);
    mul(yy,yy,yy);
    add(xx,yy,xx); // den = xx^2+yy^2
    add(1.0E-3f*max(xx),xx,xx);
    plot(xx[ns/2],"den");

    // scale = numerator/denominator
    div(num,xx,num);
    //pow(num,0.50f,num); // increase small values
    mul(1.0f/max(num),num,num); // normalized
    mul(num,x,x);
    mul(num,y,y);
    equalize(x,y);

    plot(num[ns/2],"scale");
    plot(x[ns/2],"x after");
    plot(y[ns/2],"y after");
  }

  private static void plot(float[][] x, String title) {
    /*
    edu.mines.jtk.mosaic.SimplePlot sp =
      edu.mines.jtk.mosaic.SimplePlot.asPixels(x);
    sp.setTitle(title);
    sp.addColorBar();
    */
  }

  private static void equalize(float[][][] x, float[][][] y) {
    float rmsx = rms(x);
    float rmsy = rms(y);
    mul(rms(y)/rms(x),x,x);
  }

  private static float rms(final float[][][] x) {
    final int n1 = x[0][0].length;
    final int n2 = x[0].length;
    final int n3 = x.length;
    float r = 0.0f;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float xi = x[i3][i2][i1];
          r += xi*xi;
        }
      }
    }
    return sqrt(r/n1/n2/n3);
  }

}
