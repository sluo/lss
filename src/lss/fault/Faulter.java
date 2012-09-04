package lss.fault;

import java.util.ArrayList;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import dnp.*;
import lss.flat.*;

/**
 * Adds faults to 2-D images.
 * @author Simon Luo
 * @version 2012.01.10
 */
public class Faulter {

  // fault growth
  private static final float GROWTH = 0.0002f;

  public static float[][][] addOneFault(float r, float[][] f) {
    int n2 = f.length;
    float c2 = 0.5f*n2+r;
    return addOneFault(r,c2,f);
  }

  public static float[][][] addTwoFaults(
    float r1, float r2, float[][] f)
  {
    int n2 = f.length;
    float c21 = 0.33f*n2+r1;
    float c22 = 0.66f*n2+r2;
    float[][][] y = addOneFault(r1,c21,f);
    // TODO how to combine r
    return addOneFault(r2,c22,y[0]);
  }

  public static float[][][] addOneFault(
    float rad, float c2, float[][] f)
  {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][][] r = new float[2][n2][n1];
    float[][][] s = new float[2][n2][n1];
    float[][] r1 = r[0], r2 = r[1];
    float[][] s1 = s[0], s2 = s[1];
    float rads = rad*rad;
    float theta = 0.0f; // starting angle
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float th = theta+i1*GROWTH;
        float i1s = i1*i1;
        float i2mc2 = i2-c2;
        float k1s = rads-i2mc2*i2mc2;
        if (i1s<k1s) {
          float v1 = (float)i1;
          float v2 = c2-i2;
          float w1 = v2*sin(th)+v1*cos(th);
          float w2 = v2*cos(th)-v1*sin(th);
          r1[i2][i1] = i1-w1;
          r2[i2][i1] = i2-(c2-w2);
        }
      }
    }
    ArrayList<Float> r1list = new ArrayList<Float>();
    ArrayList<Float> r2list = new ArrayList<Float>();
    ArrayList<Float> x1list = new ArrayList<Float>();
    ArrayList<Float> x2list = new ArrayList<Float>();
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        if (r1[i2][i1]!=0.0f) {
          r1list.add(r1[i2][i1]);
          r2list.add(r2[i2][i1]);
          x1list.add((float)i1);
          x2list.add((float)i2);
          r1list.add(0.0f);
          r2list.add(0.0f);
          x1list.add((float)i1);
          x2list.add((float)i2-1);
          break;
        }
      }
    }
    int size = r1list.size();
    System.out.println("size="+size);
    float[][] q = new float[4][size];
    for (int i=0; i<size; ++i) {
      q[0][i] = r1list.get(i).floatValue();
      q[1][i] = r2list.get(i).floatValue();
      q[2][i] = x1list.get(i).floatValue();
      q[3][i] = x2list.get(i).floatValue();
    }
    convertR2S(r,s);
    float[][] g = applyShiftsS(f,s);
    return new float[][][]{g,r1,r2,s1,s2,q};
  }

  public static float[][] applyShiftsR(float[][] f, float[][][] r) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] r1 = r[0], r2 = r[1];
    SincInterpolator si = new SincInterpolator();
    si.setUniform(n1,1.0,0.0,n2,1.0,0.0,f);
    float[][] g = zerofloat(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        g[i2][i1] = si.interpolate(i1-r1[i2][i1],i2-r2[i2][i1]);
      }
    }
    return g;
  }

  public static float[][][] applyShiftsR(float[][][] f, float[][][][] r) {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] r1 = r[0], r2 = r[1], r3 = r[2];
    final SincInterpolator si = new SincInterpolator();
    si.setUniform(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,f);
    final float[][][] g = zerofloat(n1,n2,n3);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          g[i3][i2][i1] = si.interpolate(i1-r1[i3][i2][i1],
                                         i2-r2[i3][i2][i1],
                                         i3-r3[i3][i2][i1]);
    }});
    return g;
  }

  public static float[][][] applyShiftsRLinear(
    float[][][] f, float[][][][] r)
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] r1 = r[0], r2 = r[1], r3 = r[2];
    final LinearInterpolator li = new LinearInterpolator();
    li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li.setUniform(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,f);
    final float[][][] g = zerofloat(n1,n2,n3);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          g[i3][i2][i1] = li.interpolate(i1-r1[i3][i2][i1],
                                         i2-r2[i3][i2][i1],
                                         i3-r3[i3][i2][i1]);
    }});
    return g;
  }

  public static float[][] applyShiftsS(float[][] f, float[][][] s) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] s1 = s[0], s2 = s[1];
    float[][] g = new float[n2][n1];
    SincInterpolator si = new SincInterpolator();
    //si.setExtrapolation(SincInterpolator.Extrapolation.ZERO);
    si.setUniform(n1,1.0,0.0,n2,1.0,0.0,f);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        g[i2][i1] = si.interpolate(i1+s1[i2][i1],i2+s2[i2][i1]);
      }
    }
    return g;
  }

  public static float[][][] applyShiftsS(
    final float[][][] f, final float[][][][] s)
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] g = zerofloat(n1,n2,n3);
    final float[][][] s1 = s[0], s2 = s[1], s3 = s[2];
    final SincInterpolator si = new SincInterpolator();
    si.setUniform(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,f);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          g[i3][i2][i1] = si.interpolate(i1+s1[i3][i2][i1],
                                         i2+s2[i3][i2][i1],
                                         i3+s3[i3][i2][i1]);
    }});
    return g;
  }

  public static float[][][] resample(final int m, final float[][][] x) {
    final int n1 = x[0][0].length;
    final int n2 = x[0].length;
    final int n3 = x.length;
    final int m1 = m*n1;
    final int m2 = m*n2;
    final int m3 = m*n3;
    final float s = 1.0f/m;
    final float[][][] y = new float[m3][m2][m1];
    final LinearInterpolator li = new LinearInterpolator();
    li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li.setUniform(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,x);
    Parallel.loop(m3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<m2; ++i2)
        for (int i1=0; i1<m1; ++i1)
          y[i3][i2][i1] = li.interpolate(s*i1,s*i2,s*i3);
    }});
    return y;
  }

  public static float[][] getGriddedThrows(
    float fnull, 
    float[] t1, float[] t2, float[] t3,
    float[] x1, float[] x2, float[] x3)
  {
    int n = t1.length;
    int nnull = 0;
    for (int i=0; i<n; ++i)
      if (t1[i]==fnull && t2[i]==fnull && t3[i]==fnull)
      //if (t1[i]==fnull)
        ++nnull;
    int nlive = n-nnull;
    float[] s1 = new float[nlive];
    float[] s2 = new float[nlive];
    float[] s3 = new float[nlive];
    float[] y1 = new float[nlive];
    float[] y2 = new float[nlive];
    float[] y3 = new float[nlive];
    float[] z1 = new float[nnull];
    float[] z2 = new float[nnull];
    float[] z3 = new float[nnull];
    for (int i=0, iz=0, iy=0; i<n; ++i) {
      if (t1[i]==fnull && t2[i]==fnull && t3[i]==fnull) {
      //if (t1[i]==fnull) {
        z1[iz] = x1[i];
        z2[iz] = x2[i];
        z3[iz] = x3[i];
        ++iz;
      } else {
        s1[iy] = t1[i];
        s2[iy] = t2[i];
        s3[iy] = t3[i];
        y1[iy] = x1[i];
        y2[iy] = x2[i];
        y3[iy] = x3[i];
        ++iy;
      }
    }
    return new float[][]{s1,s2,s3,y1,y2,y3,z1,z2,z3};
  }

  /**
   * Adjusts distance for use with BlendedGridder3
   * @param adj amount to adjust distances by.
   * @param d array {r1,r2,r3} of distances.
   */
  public static void adjustDistances(
    final float adj, final float[][][] d)
  {
    final int n1 = d[0][0].length;
    final int n2 = d[0].length;
    final int n3 = d.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          d[i3][i2][i1] = max(0.0f,d[i3][i2][i1]-adj);
    }});
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Invert shifts using fixed-point iterations.
  private static void convertS2R(float[][][] s, float[][][] r) {
    int n1 = r[0][0].length;
    int n2 = r[0].length;
    float[][] r1 = r[0], r2 = r[1];
    float[][] s1 = s[0], s2 = s[1];
    LinearInterpolator li1 = new LinearInterpolator();
    li1.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li1.setUniform(n1,1.0,0.0,n2,1.0,0.0,s1);
    LinearInterpolator li2 = new LinearInterpolator();
    li2.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li2.setUniform(n1,1.0,0.0,n2,1.0,0.0,s2);
    for (int i2=0; i2<n2; ++i2) {
      double u2 = i2;
      for (int i1=0; i1<n1; ++i1) {
        double u1 = i1;
        double r1i = s1[i2][i1];
        double r2i = s2[i2][i1];
        double r1p,r2p,dr1,dr2;
        for (int iter=0; iter<100; ++iter) {
          r1p = r1i;
          r2p = r2i;
          double x1 = u1-r1i;
          double x2 = u2-r2i;
          r1i = li1.interpolate(x1,x2);
          r2i = li2.interpolate(x1,x2);
          dr1 = r1i-r1p;
          dr2 = r2i-r2p;
          if (dr1*dr1+dr2*dr2<0.0001) 
            break;
        }
        r1[i2][i1] = (float)r1i;
        r2[i2][i1] = (float)r2i;
      }
    }
  }
  private static void convertR2S(float[][][] r, float[][][] s) {
    int n1 = r[0][0].length;
    int n2 = r[0].length;
    float[][] r1 = r[0], r2 = r[1];
    float[][] s1 = s[0], s2 = s[1];
    LinearInterpolator li1 = new LinearInterpolator();
    li1.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li1.setUniform(n1,1.0,0.0,n2,1.0,0.0,r1);
    LinearInterpolator li2 = new LinearInterpolator();
    li2.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li2.setUniform(n1,1.0,0.0,n2,1.0,0.0,r2);
    for (int i2=0; i2<n2; ++i2) {
      double x2 = i2;
      for (int i1=0; i1<n1; ++i1) {
        double x1 = i1;
        double s1i = r1[i2][i1];
        double s2i = r2[i2][i1];
        double s1p,s2p,ds1,ds2;
        for (int iter=0; iter<100; ++iter) {
          s1p = s1i;
          s2p = s2i;
          double u1 = x1+s1i;
          double u2 = x2+s2i;
          s1i = li1.interpolate(u1,u2);
          s2i = li2.interpolate(u1,u2);
          ds1 = s1i-s1p;
          ds2 = s2i-s2p;
          if (ds1*ds1+ds2*ds2<0.0001) 
            break;
        }
        s1[i2][i1] = (float)s1i;
        s2[i2][i1] = (float)s2i;
      }
    }
  }

}
