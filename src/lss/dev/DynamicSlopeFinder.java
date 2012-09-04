package lss.dev;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;

/**
 * Estimates local slopes of features in 2D and 3D images.
 * Slopes are estimated by first comparing references traces to other
 * traces in an image and estimating the shifts between them using
 * dynamic warping. The shifts can then be easily converted to
 * slopes or horizons.
 * @author Simon Luo
 * @version 2012.08.16
 */
public class DynamicSlopeFinder {

  /**
   * Constructs a slope finder with default parameters.
   */
  public DynamicSlopeFinder() {
    this(SIGMA,SHIFT_MAX);
  }

  /**
   * Constructs a slope finder with specified half-width.
   * The half-width determines the size of the Gaussian window used
   * when averaging slopes estimated for different reference traces.
   * Increasing the half-width may improve slope estimates at faults.
   * @param sigma half-width of Gaussian window used when averaging slopes.
   */
  public DynamicSlopeFinder(double sigma) {
    this(sigma,SHIFT_MAX);
  }

  /**
   * Constructs a slope finder with specified half-width and maximum shift.
   * The half-width determines the size of the Gaussian window used
   * when averaging slopes estimated for different reference traces.
   * The maximum shift should be larger than the half the largest possible 
   * vertical shift expected within the aforementioned Gaussian window.
   * @param sigma half-width of Gaussian window used when averaging slopes.
   * @param shiftMax maximum shift allowed.
   */
  public DynamicSlopeFinder(double sigma, int shiftMax) {
    _dw = new DynamicWarping(-shiftMax,shiftMax);
    _dw.setStrainMax(STRAIN_MAX1,STRAIN_MAX2);
    _dw.setShiftSmoothing(1.0,1.0); // shift smoothing
    _dw.setErrorSmoothing(2); // number of smoothings of alignment errors
    _sigma = sigma;
  }

  /**
   * Sets strain limits for dynamic warping.
   * @param strainMax1 strain limit in the 1st dimension.
   * @param strainMax2 strain limit in the 2nd dimension.
   */
  public void setStrainMax(double strainMax1, double strainMax2) {
    _dw.setStrainMax(strainMax1,strainMax2);
  }

  // Default parameters
  private static final int STEP = 2; // step between reference traces
  private static final int SHIFT_MAX = 32; // maximum shift
  private static final double STRAIN_MAX1 = 0.50; // 50% strain
  private static final double STRAIN_MAX2 = 1.00; // 100% strain
  private static final double SIGMA = 32.0; // Gaussian half-width

  private DynamicWarping _dw;
  private double _sigma;

  //////////////////////////////////////////////////////////////////////////
  // shifts

  /**
   * Finds shifts between a 2D image and reference image.
   * Shifts are found by using dynamic warping to warp the input
   * image to a reference image consisting of the reference trace
   * repeated across the entire image.
   * @param k2 index of reference trace.
   * @param f input array of image samples.
   * @param s output array of shifts.
   */
  public void findShiftsFromReference(int k2, float[][] f, float[][] s) {
    findShiftsFromReference(f[k2],f,s);
  }

  /**
   * Finds shifts between a 2D image and reference image.
   * Shifts are found by using dynamic warping to warp the input
   * image to a reference image consisting of the reference trace
   * repeated across the entire image.
   * @param r input reference trace.
   * @param f input array of image samples.
   * @param s output array of shifts.
   */
  public void findShiftsFromReference(float[] r, float[][] f, float[][] s) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] t = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2)
      copy(r,t[i2]);
    _dw.findShifts(t,f,s);
  }

  /**
   * Applies shifts.
   * @param s input array of shifts.
   * @param f input array of image samples.
   * @return array of shifted image samples.
   */
  public float[][] applyShifts(float[][] s, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] g = new float[n2][n1];
    applyShifts(s,f,g);
    return g;
  }

  /**
   * Applies shifts.
   * @param s input array of shifts.
   * @param f input array of image samples.
   * @param g output array of shifted image samples.
   */
  public void applyShifts(float[][] s, float[][] f, float[][] g) {
    final int n1 = s[0].length;
    final int n2 = s.length;
    final float[][] sf = s;
    final float[][] ff = f;
    final float[][] gf = g;
    final Parallel.Unsafe<SincInterpolator> siu =
      new Parallel.Unsafe<SincInterpolator>();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      SincInterpolator si = siu.get();
      if (si==null) {
        si = new SincInterpolator();
        si.setUniformSampling(n1,1.0,0.0);
        siu.set(si);
      }
      si.setUniformSamples(ff[i2]);
      for (int i1=0; i1<n1; ++i1) {
        gf[i2][i1] = si.interpolate(i1+sf[i2][i1]);
      }
    }});
  }

  //////////////////////////////////////////////////////////////////////////
  // horizons

  /**
   * Finds horizons for a 2D image for a single reference trace.
   * @param k2 index of reference trace.
   * @param f input array[n2][n1] of image samples.
   * @return output array[n1][n2] of horizons.
   */
  public float[][] findHorizonsFromReference(int k2, float[][] f) {
    return findHorizonsFromReference(f[k2],f);
  }

  /**
   * Finds horizons for a 2D image for a single reference trace.
   * @param r input reference trace.
   * @param f input array[n2][n1] of image samples.
   * @return output array[n1][n2] of horizons.
   */
  public float[][] findHorizonsFromReference(float[] r, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] s = new float[n2][n1];
    findShiftsFromReference(r,f,s);
    float[][] h = transpose(s);
    for (int i1=0; i1<n1; ++i1)
      add(i1,h[i1],h[i1]);
    return h;
  }

  //////////////////////////////////////////////////////////////////////////
  // slopes

  /**
   * Finds slopes for a 2D image for a single reference trace.
   * @param k2 index of reference trace.
   * @param f input array of image samples.
   * @param p2 output array of slopes.
   */
  public void findSlopesFromReference(
    int k2, float[][] f, float[][] p2)
  {
    findSlopesFromReference(f[k2],f,p2);
  }

  /**
   * Finds slopes for a 2D image for a single reference trace.
   * @param r input reference trace.
   * @param f input array of image samples.
   * @param p2 output array of slopes.
   */
  public void findSlopesFromReference(
    float[] r, float[][] f, float[][] p2)
  {
    findShiftsFromReference(r,f,p2);
    new RecursiveGaussianFilter(1.0).apply01(p2,p2);
    //new RecursiveGaussianFilter(1.0).applyX1(p,p);
  }

  /**
   * Finds slopes for a 2D image.
   * @param f input array of image samples.
   * @param p2 output array of slopes.
   */
  public void findSlopes(float[][] f, float[][] p2) {
    final int n1 = f[0].length;
    final int n2 = f.length;
    final float[][] ff = f;

    final RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(_sigma);
    float[][][] nd = reduce(0,n2,STEP,new ReduceInt<float[][][]>() {
      public float[][][] compute(int k2) {
        float[][][] nd = new float[2][n2][n1]; // slopes and weights
        int w2 = (int)(3.0*_sigma);
        int f2 = max(k2-w2,0);
        int l2 = min(k2+w2,n2-1);
        int m2 = 1+l2-f2;
        float[][] g = new float[m2][];
        float[][] p = new float[m2][];
        for (int i2=f2,j2=0 ; i2<=l2; ++i2,++j2) {
          g[j2] = ff[i2];
          p[j2] = nd[0][i2];
        }
        findSlopesFromReference(ff[k2],g,p);
        fill(1.0f,nd[1][k2]);
        rgf.applyX0(nd[1],nd[1]);
        mul(nd[1],nd[0],nd[0]);
        return nd;
      }
      public float[][][] combine(float[][][] a, float[][][] b) {
        return add(a,b);
      }
    });
    div(nd[0],nd[1],p2);
  }

  /**
   * Finds slopes for a 3D image.
   * @param f input array of image samples.
   * @param p2 output array of slopes along 2nd dimension.
   * @param p3 output array of slopes along 3rd dimension.
   */
  public void findSlopes(float[][][] f, float[][][] p2, float[][][] p3) {
    findSlopes2(f,p2);
    findSlopes3(f,p3);
  }

  /**
   * Finds slopes for a 3D image along the 2nd dimension.
   * @param f input array of image samples.
   * @param p2 output array of slopes along 2nd dimension.
   */
  public void findSlopes2(float[][][] f, float[][][] p2) {
    findSlopes(f,p2);
  }

  /**
   * Finds slopes for a 3D image along the 3rd dimension.
   * @param f input array of image samples.
   * @param p3 output array of slopes along 3rd dimension.
   */
  public void findSlopes3(float[][][] f, float[][][] p3) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] ft = new float[n2][n3][];
    float[][][] p3t = new float[n2][n3][];
    // Transpose 2nd and 3rd dimensions, then find slopes.
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        ft[i2][i3] = f[i3][i2];
        p3t[i2][i3] = p3[i3][i2];
      }
    }
    findSlopes(ft,p3t);
  }

  private void findSlopes(float[][][] f, float[][][] p) {
    final int n3 = f.length;
    final float[][][] ff = f;
    final float[][][] pf = p;
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      //System.out.println("i3="+i3);
      findSlopes(ff[i3],pf[i3]);
    }});
  }
}
