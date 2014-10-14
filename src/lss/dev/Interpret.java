package lss.dev;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

// TESTING
import edu.mines.jtk.mosaic.*;

public class Interpret {

  public float[][] computeErrors(float[] f, float[] g) {
    int n = f.length;
    float[][] e = new float[n][n];
    //Error error = new AbsoluteError();
    Error error = new GaussianAbsoluteError(12.0);
    error.compute(f,g,e);
    return e;
  }
  
  //////////////////////////////////////////////////////////////////////////

  private interface Error {
    public void compute(float[] f, float[] g, float[][] e);
  }
  private class AbsoluteError implements Error {
    public void compute(float[] f, float[] g, float[][] e) {
      int n = f.length;
      for (int i=0; i<n; ++i)
        for (int j=0; j<n; ++j)
          e[i][j] = abs(f[i]-g[j]);
    }
  }

  private abstract class WindowedError implements Error {
    public float sigma;
    public WindowedError(double sigma) {
      this.sigma = (float)sigma;
    }
    public abstract void compute(float[] f, float[] g, float[][] e);
  }
  private class GaussianAbsoluteError extends WindowedError {
    public GaussianAbsoluteError(double sigma) {
      super(sigma);
    }
    public void compute(float[] f, float[] g, float[][] e) {
      int n = f.length;
      float[][] w = new float[n][n];
      for (int i=0; i<n; ++i)
        w[i][i] = 1.0f;
      new RecursiveGaussianFilter(super.sigma).apply0X(w,w);
      for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
          float[] fw = mul(f,w[i]);
          float[] gw = mul(g,w[j]);
          int d = j-i;
          int ad = abs(d);
          float s = 0.0f;
          for (int k=ad; k<n-ad; ++k)
            s += abs(f[k]-g[k+d]);
          e[i][j] = s;
        }
      }
    }
  }

  private static void computeErrors(
    Error error, float[] f, float[] g, float[][] e)
  { 
    error.compute(f,g,e);
  }
}
