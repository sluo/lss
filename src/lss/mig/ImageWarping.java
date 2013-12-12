package lss.mig;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import lss.mod.Receiver;

public class ImageWarping {

  public ImageWarping(
    double strainT, double strainR, double strainS,
    double smoothT, double smoothR, double smoothS,
    double maxShift, double dt)
  {
    this(strainT,strainR,strainS,smoothT,smoothR,smoothS,maxShift,dt,1);
  }

  public ImageWarping(
    double strainT, double strainR, double strainS,
    double smoothT, double smoothR, double smoothS,
    double maxShift, double dt, int td)
  {
    Check.argument(td>=1,"td>=1");
    int shiftMax = (int)(maxShift/(td*dt));
    _warp = new DynamicWarping(-shiftMax,shiftMax);
    _warp.setShiftSmoothing(smoothT/td,smoothR,smoothS); // shift smoothing
    _warp.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.AVERAGE);
    _warp.setErrorSmoothing(2); // number of smoothings of alignment errors
    if (strainS<=0.0) {
      _warp.setStrainMax(strainT,strainR);
      _3d = false;
    } else {
      _warp.setStrainMax(strainT,strainR,strainS);
      _3d = true;
    }
    _td = td;
    _sigmaRms = 1.0*maxShift/(td*dt);
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
    final int ntm = nt/_td;
    float[][][] ep = new float[ns][][];
    float[][][] eo = new float[ns][][];
    float[][][] uu = new float[ns][nr][ntm];
    for (int is=0; is<ns; ++is) {
      ep[is] = copy(ntm,nr,0,0,_td,1,rp[is].getData());
      eo[is] = copy(ntm,nr,0,0,_td,1,ro[is].getData());
    }
    findShifts(ep,eo,uu);
    if (_td>1) { // if decimated, interpolate shifts
      mul((float)_td,uu,uu); // scale shifts to compensate for decimation
      final LinearInterpolator li = new LinearInterpolator();
      li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
      li.setUniform(ntm,_td,0.0,nr,1.0,0.0,ns,1.0,0.0,uu);
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

  //////////////////////////////////////////////////////////////////////////

  private int _td; // time decimation
  private boolean _3d; // true for 3D warping
  private double _sigmaRms; // sigma for RMS filtering
  private DynamicWarping _warp;

  private void findShifts(
    final float[][][] rp, final float[][][] ro, final float[][][] u)
  {
    rmsFilter(rp,ro);
    addRandomNoise(10.0f,rp,ro);
    if (_3d) {
      _warp.findShifts(rp,ro,u);
    } else {
      int ns = rp.length;
      Parallel.loop(ns,new Parallel.LoopInt() {
      public void compute(int is) {
        _warp.findShifts(rp[is],ro[is],u[is]);
      }});
    }
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
    plot(x[0],"x before");
    plot(y[0],"y before");

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
      ref.apply(xx,xx);
    } else {
      ref.apply2(xx,xx);
      ref.apply1(xx,xx);
      ref.apply2(yy,yy);
      ref.apply1(yy,yy);
    }

    plot(xx[0],"xx smoothed");
    plot(yy[0],"yy smoothed");

    // numerator = 2*xx*yy
    float[][][] num = mul(2.0f,xx);
    mul(yy,num,num);
    plot(num[0],"num");

    // denominator = xx^2+yy^2
    mul(xx,xx,xx);
    mul(yy,yy,yy);
    add(xx,yy,xx); // den = xx^2+yy^2
    add(1.0E-3f*max(xx),xx,xx);
    plot(xx[0],"den");

    // scale = numerator/denominator
    div(num,xx,num);
    //pow(num,0.50f,num); // increase small values
    mul(1.0f/max(num),num,num); // normalized
    mul(num,x,x);
    mul(num,y,y);
    equalize(x,y);

    plot(num[0],"scale");
    plot(x[0],"x after");
    plot(y[0],"y after");
  }

  private static void plot(float[][] x, String title) {
    //edu.mines.jtk.mosaic.SimplePlot sp =
    //  edu.mines.jtk.mosaic.SimplePlot.asPixels(x);
    //sp.setTitle(title);
    //sp.addColorBar();
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
