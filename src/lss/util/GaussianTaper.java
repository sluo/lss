package lss.util;

import edu.mines.jtk.dsp.RecursiveGaussianFilter;
import static edu.mines.jtk.util.ArrayMath.*;

public class GaussianTaper {

  /*
  public static void apply1(float[][] x, float[][] y) {
    apply1(R,x,y);
  }

  public static void apply2(float[][] x, float[][] y) {
    apply2(R,x,y);
  }
  */

  public static void apply1(float r, float[][] x, float[][] y) {
    int n1 = y[0].length;
    int n2 = y.length;
    float[] t = new float[n1];
    t[n1/2] = 1.0f;
    new RecursiveGaussianFilter(r*n1).apply0(t,t);
    div(t,max(t),t);
    for (int i2=0; i2<n2; ++i2)
      mul(t,x[i2],y[i2]);
  }

  public static void apply2(float r, float[][] x, float[][] y) {
    int n1 = y[0].length;
    int n2 = y.length;
    float[] t = new float[n2];
    t[n2/2] = 1.0f;
    new RecursiveGaussianFilter(r*n2).apply0(t,t);
    div(t,max(t),t);
    for (int i2=0; i2<n2; ++i2)
      mul(t[i2],x[i2],y[i2]);
  }

  public static float[] getValuesForLength(int n) {
    float[] t = new float[n];
    t[n/2] = 1.0f;
    new RecursiveGaussianFilter(R*n).apply0(t,t);
    div(t,max(t),t);
    return t;
  }

  private static final float R = 0.25f; // sigma=R*n

  ////////////////////////////////////////////////////////////////////////////

  public GaussianTaper(double sigma) {
    _rgf = new RecursiveGaussianFilter(sigma);
  }

  public void apply(float[] x, float[] y) {
    int n = y.length;
    float[] w = getWeights(n,_rgf);
    mul(w,x,y);
  }

  public void apply1(float[][] x, float[][] y) {
    int n1 = y[0].length;
    int n2 = y.length;
    float[] w = getWeights(n1,_rgf);
    for (int i2=0; i2<n2; ++i2) {
      mul(w,x[i2],y[i2]);
    }
  }

  public void apply2(float[][] x, float[][] y) {
    int n1 = y[0].length;
    int n2 = y.length;
    float[] w = getWeights(n2,_rgf);
    for (int i2=0; i2<n2; ++i2) {
      mul(w[i2],x[i2],y[i2]);
    }
  }

  public float[] getWeightsForLength(int n) {
    return getWeights(n,_rgf);
  }

  private static float[] getWeights(int n, RecursiveGaussianFilter rgf) {
    float[] w = fillfloat(1.0f,n);
    rgf.apply0(w,w);
    mul(1.0f/max(w),w,w);
    return w;
  }

  private RecursiveGaussianFilter _rgf;
}
