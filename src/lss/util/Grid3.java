package lss.util;

import edu.mines.jtk.dsp.Sampling;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;

// FOR TESTING
import java.util.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;

/**
 * Gridding by nearest-neighbor interpolation or discrete Sibson 
 * interpolation of scattered samples f(x1,x2).
 * @author Simon Luo
 * @version 2011.12.08
 */
public class Grid3 {

  /**
   * Constructs a gridder.
   * @param f sample values f(x1,x2,x3).
   * @param x1 sample coordinates x1
   * @param x2 sample coordinates x2
   * @param x3 sample coordinates x3
   */
  public Grid3(
    float[] f, float[] x1, float[] x2, float[] x3)
  {
    _f = f;
    _x = new float[][]{x1,x2,x3};
    _tree = new KdTree(_x);
  }

  /**
   * Sets the known samples.
   * @param f sample values f(x1,x2,x3).
   * @param x1 sample coordinates x1
   * @param x2 sample coordinates x2
   * @param x3 sample coordinates x3
   */
  public void setScattered(
    float[] f, float[] x1, float[] x2, float[] x3)
  {
    _f = f;
    _x = new float[][]{x1,x2,x3};
    _tree = new KdTree(_x);
  }

  /**
   * Computes gridded sample values using 
   * nearest-neighbor interpolation.
   * @param s1 sampling of x1.
   * @param s2 sampling of x2.
   * @param s3 sampling of x3.
   * @return array of gridded sample values.
   */
  public float[][][] gridNearest(
    Sampling s1, Sampling s2, Sampling s3)
  {
    return gridNearest(s1,s2,s3,null);
  }

  /**
   * Computes gridded sample values using 
   * nearest-neighbor interpolation.
   * @param s1 sampling of x1.
   * @param s2 sampling of x2.
   * @param s3 sampling of x3.
   * @param d array of distances to nearest sample
   * @return array of gridded sample values.
   */
  public float[][][] gridNearest(
    Sampling s1, Sampling s2, Sampling s3, float[][][] d)
  {
    return _parallel?gridNearestP(s1,s2,s3,d):gridNearestS(s1,s2,s3,d);
  }

  /**
   * Computes gridded sample values using 
   * discrete Sibson interpolation.
   * @param s1 sampling of x1.
   * @param s2 sampling of x2.
   * @param s3 sampling of x3.
   * @return array of gridded sample values.
   */
  public float[][][] gridSibson(
    Sampling s1, Sampling s2, Sampling s3) {
    return _parallel?gridSibsonP(s1,s2,s3):gridSibsonS(s1,s2,s3);
  }

  /**
   * Enables or disables or parallel processing.
   * @param parallel true for parallel processing
   */
  public static void setParallel(boolean parallel) {
    _parallel = parallel;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final int NREDUCE = 8;
  private static boolean _parallel = true;

  private float[] _f;   // values
  private float[][] _x; // coordinates
  private KdTree _tree;

  private float[][][] gridNearestS(
    Sampling s1, Sampling s2, Sampling s3, float[][][] d)
  {
    int n1 = s1.getCount(); 
    int n2 = s2.getCount();
    int n3 = s3.getCount();
    final float[][][] y = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      gridNearestSlice3(i3,y,s1,s2,s3,d);
    }
    return y;
  }

  private float[][][] gridNearestP(
    final Sampling s1, final Sampling s2, final Sampling s3,
    final float[][][] d)
  {
    final int n1 = s1.getCount(); 
    final int n2 = s2.getCount();
    final int n3 = s3.getCount();
    final float[][][] y = new float[n3][n2][n1];
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      gridNearestSlice3(i3,y,s1,s2,s3,d);
    }});
    return y;
  }

  private void gridNearestSlice3(
    int i3, float[][][] y, Sampling s1, Sampling s2, Sampling s3,
    float[][][] d)
  {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float x1 = (float)s1.getValue(i1);
        float x2 = (float)s2.getValue(i2);
        float x3 = (float)s3.getValue(i3);
        float[] xq = new float[]{x1,x2,x3};
        int nearest = _tree.findNearest(xq);
        y[i3][i2][i1] = _f[nearest];
        if (d!=null) d[i3][i2][i1] = sqrt(distance(nearest,_x,xq));
      }
    }
  }

  private float[][][] gridSibsonS(
    Sampling s1, Sampling s2, Sampling s3)
  {
    int n3 = s3.getCount();
    float[][][][] yc = accumulateSibson(0,n3,s1,s2,s3);
    return div(yc[0],yc[1]);
  }

  public float[][][] gridSibsonP(
    final Sampling s1, final Sampling s2, final Sampling s3)
  {
    final int n3 = s3.getCount();
    final int k3 = n3/NREDUCE;
    float[][][][] yc = reduce(NREDUCE,new ReduceInt<float[][][][]>() {
      public float[][][][] compute(int i) {
        int min3 = i*k3;
        int max3 = (i<NREDUCE-1)?min3+k3:n3;
        return accumulateSibson(min3,max3,s1,s2,s3);
      }
      public float[][][][] combine(float[][][][] yc1, float[][][][] yc2) {
        add(yc1[0],yc2[0],yc2[0]);
        add(yc1[1],yc2[1],yc2[1]);
        return yc2;
      }
    });
    return div(yc[0],yc[1]);
  }

  public float[][][] sdiv(float[][][] num, float[][][] den) {
    int n1 = num[0][0].length;
    int n2 = num[0].length;
    int n3 = num.length;
    float[][][] y = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float deni = den[i3][i2][i1];
          y[i3][i2][i1] = (deni==0.0f)?0.0f:num[i3][i2][i1]/deni;
        }
      }
    }
    return y;
  }

  private float[][][][] accumulateSibson(
    int p3, int q3, Sampling s1, Sampling s2, Sampling s3)
  {
    int n1 = s1.getCount(); 
    int n2 = s2.getCount();
    int n3 = s3.getCount();
    double d1 = s1.getDelta();
    double d2 = s2.getDelta();
    double d3 = s3.getDelta();
    double f1 = s1.getFirst();
    double f2 = s2.getFirst();
    double f3 = s3.getFirst();
    float[] xq = new float[3];
    float[][][] c = new float[n3][n2][n1];
    float[][][] y = new float[n3][n2][n1];
    for (int i3=p3; i3<q3; ++i3) {
      float v3 = (float)s3.getValue(i3);
      xq[2] = v3;
      for (int i2=0; i2<n2; ++i2) {
        float v2 = (float)s2.getValue(i2);
        xq[1] = v2;
        for (int i1=0; i1<n1; ++i1) {
          float v1 = (float)s1.getValue(i1);
          xq[0] = v1;

          // Find distance to nearest neighbor
          int nearest = _tree.findNearest(xq);
          float rs = distance(nearest,_x,xq);
          float f = _f[nearest];

          // Fill circle and increment counter
          float r = sqrt(rs);
          int min1 = s1.indexOfNearest(v1-r);
          int max1 = s1.indexOfNearest(v1+r);
          int min2 = s2.indexOfNearest(v2-r);
          int max2 = s2.indexOfNearest(v2+r);
          int min3 = s3.indexOfNearest(v3-r);
          int max3 = s3.indexOfNearest(v3+r);
          for (int j3=min3; j3<=max3; ++j3) {
            float w3 = (float)s3.getValue(j3);
            float r3 = w3-v3;
            float r3s = r3*r3;
            for (int j2=min2; j2<=max2; ++j2) {
              float w2 = (float)s2.getValue(j2);
              float r2 = w2-v2;
              float r2s = r2*r2;
              for (int j1=min1; j1<=max1; ++j1) {
                float w1 = (float)s1.getValue(j1);
                float r1 = w1-v1;
                float r1s = r1*r1;
                if (r1s+r2s+r3s<rs) {
                  y[j3][j2][j1] += f;
                  c[j3][j2][j1] += 1.0f;
                }
              }
            }
          }
        }
      }
    }
    return new float[][][][]{y,c};
  }

  private static float distance(int i, float[][] x, float[] xq) {
    int k = x.length;
    double ds = 0.0f;
    for (int ik=0; ik<k; ++ik) {
      double d = x[ik][i]-xq[ik];
      ds += d*d;
    }
    return (float)ds;
  }

}
