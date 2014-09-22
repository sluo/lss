package lss.flat;

import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import dnp.*;
import lss.util.*;

/**
 * Utilies for flattening.
 * @author Simon Luo
 * @version 2012.02.01
 */
public class FlattenerUtil {

  /**
   * Computes vertices for use in QuadGroup.
   * @param x1 sample coordinates x1
   * @param x2 sample coordinates x2
   * @param x3 sample coordinates x3
   * @param s fault locations
   * @return array of packed vertex coordinates
   */
  public static float[] makeQuadVertices(
    float[][] x1, float[][] x2, float[][] x3, float[][] s)
  {
    int n2 = x1[0].length-1;
    int n3 = x1.length-1;
    float threshold = 0.1f;

    /*
    float[] xyz = new float[3*4*n2*n3];
    for (int i3=0, i=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
    */
    int trim = 4; // trim edges
    float[] xyz = new float[3*4*(n2-2*trim)*(n3-2*trim)];
    for (int i3=trim, i=0; i3<n3-trim; i3++) {
      for (int i2=trim; i2<n2-trim; i2++) {


        //          d---c
        // a quad:  |   |
        //          a---b
        float xa = x3[i3  ][i2  ];
        float xb = x3[i3+1][i2  ];
        float xc = x3[i3+1][i2+1];
        float xd = x3[i3  ][i2+1];
        float ya = x2[i3  ][i2  ];
        float yb = x2[i3+1][i2  ];
        float yc = x2[i3+1][i2+1];
        float yd = x2[i3  ][i2+1];
        float za = x1[i3  ][i2  ];
        float zb = x1[i3+1][i2  ];
        float zc = x1[i3+1][i2+1];
        float zd = x1[i3  ][i2+1];

        // move vertex if at a fault
        boolean a = s[i3  ][i2  ]>threshold;
        boolean b = s[i3+1][i2  ]>threshold;
        boolean c = s[i3+1][i2+1]>threshold;
        boolean d = s[i3  ][i2+1]>threshold;
        float xm = 0.25f*(xa+xb+xc+xd); // midpoints
        float ym = 0.25f*(ya+yb+yc+yd);
        float zm = 0.25f*(za+zb+zc+zd);
        if (a) {
          if (b || c || d) {
            continue;
          } else {
            xa = xm;
            ya = ym;
            za = zm;
          }
        }
        if (b) {
          if (a || c || d) {
            continue;
          } else {
            xb = xm;
            yb = ym;
            zb = zm;
          }
        }
        if (c) {
          if (a || b || d) {
            continue;
          } else {
            xc = xm;
            yc = ym;
            zc = zm;
          }
        }
        if (d) {
          if (a || b || c) {
            continue;
          } else {
            xd = xm;
            yd = ym;
            zd = zm;
          }
        }

        xyz[i++] = xa; xyz[i++] = ya;  xyz[i++] = za;
        xyz[i++] = xb; xyz[i++] = yb;  xyz[i++] = zb;
        xyz[i++] = xc; xyz[i++] = yc;  xyz[i++] = zc;
        xyz[i++] = xd; xyz[i++] = yd;  xyz[i++] = zd;
      }
    }
    return xyz;
  }

  public static float[][] makeQuadVerticesAndRgbFloats(
    float[][] x1, float[][] x2, float[][] x3, float[][] s)
  {
    float[] xyz = makeQuadVertices(x1,x2,x3,s);
    int n = xyz.length;
    int nv = n/3;
    float[] v = new float[nv];
    ColorMap cmap = new ColorMap(-max(x1),-min(x1),ColorMap.JET);
    for (int iv=0, ixyz=2; iv<nv; ++iv, ixyz+=3)
      v[iv] = -xyz[ixyz];
    float[] rgb = cmap.getRgbFloats(v);
    return new float[][]{xyz,rgb};
  }
  
  // make vertices for TriangleGroup
  public static float[] makeTriangleVertices(
    float[][] x1, float[][] x2, float[][] x3, float[][] ep)
  {
    int n2 = x1[0].length;
    int n3 = x1.length;
    float epmin = 0.66f; //minimum planarity
    int nv = 0;
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        if (ep[i3][i2]>epmin) nv++;
      }
    }
    float[] xyz = new float[3*6*nv-36];
    for (int i3=0,i=0; i3<n3-1; i3++) {
      for (int i2=0; i2<n2-1; i2++) {
        if (ep[i3][i2]<=epmin)
          continue;
        float x00 = x3[i3  ][i2  ]; float x10 = x3[i3+1][i2  ];
        float y00 = x2[i3  ][i2  ]; float y01 = x2[i3  ][i2+1];
        float z00 = x1[i3  ][i2  ]; float z01 = x1[i3  ][i2+1];
        float z10 = x1[i3+1][i2  ]; float z11 = x1[i3+1][i2+1];
        xyz[i++] = x00;  xyz[i++] = y00;  xyz[i++] = z00;
        xyz[i++] = x00;  xyz[i++] = y01;  xyz[i++] = z01;
        xyz[i++] = x10;  xyz[i++] = y00;  xyz[i++] = z10;
        xyz[i++] = x10;  xyz[i++] = y00;  xyz[i++] = z10;
        xyz[i++] = x00;  xyz[i++] = y01;  xyz[i++] = z01;
        xyz[i++] = x10;  xyz[i++] = y01;  xyz[i++] = z11;
      }
    }
    return xyz;
  }

  /**
   * Applies the specified shifts using sinc interpolation.
   * The returned array is a sampling of g(u) = f(u-r(u)).
   * @param f input array to which shifts are to be applied.
   * @param r array {r1,r2} of shifts.
   * @return array with shifts applied.
   */
  public static float[][] applyShiftsR(
    final float[][] f, final float[][][] r)
  {
    final int n1 = f[0].length;
    final int n2 = f.length;
    final float[][] r1 = r[0], r2 = r[1];
    final SincInterpolator si = new SincInterpolator();
    final float[][] g = zerofloat(n1,n2);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        g[i2][i1] = si.interpolate(
          n1,1.0,0.0,n2,1.0,0.0,f,
          i1-r1[i2][i1],i2-r2[i2][i1]);
      }
    }});
    return g;
  }

  /**
   * Applies the specified shifts using linear interpolation.
   * The returned array is a sampling of g(u) = f(u-r(u)).
   * @param f input array to which shifts are to be applied.
   * @param r array {r1,r2} of shifts.
   * @return array with shifts applied.
   */
  public static float[][] applyShiftsRLinear(float[][] f, float[][][] r) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] r1 = r[0], r2 = r[1];
    LinearInterpolator li = new LinearInterpolator();
    li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li.setUniform(n1,1.0,0.0,n2,1.0,0.0,f);
    float[][] g = zerofloat(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        g[i2][i1] = li.interpolate(i1-r1[i2][i1],i2-r2[i2][i1]);
      }
    }
    return g;
  }

  /**
   * Applies the specified shifts using sinc interpolation.
   * The returned array is a sampling of g(u) = f(u-r(u)).
   * @param f input array to which shifts are to be applied.
   * @param r array {r1,r2,r3} of shifts.
   * @return array with shifts applied.
   */
  public static float[][][] applyShiftsR(
    final float[][][] f, final float[][][][] r)
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] r1 = r[0], r2 = r[1], r3 = r[2];
    final SincInterpolator si = new SincInterpolator();
    final float[][][] g = zerofloat(n1,n2,n3);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          g[i3][i2][i1] = si.interpolate(
            n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,f, 
            i1-r1[i3][i2][i1],i2-r2[i3][i2][i1],i3-r3[i3][i2][i1]);
    }});
    return g;
  }

  /**
   * Applies the specified shifts using linear interpolation.
   * The returned array is a sampling of g(u) = f(u-r(u)).
   * @param f input array to which shifts are to be applied.
   * @param r array {r1,r2,r3} of shifts.
   * @return array with shifts applied.
   */
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

  /**
   * Applies the specified shifts using sinc interpolation.
   * The returned array is a sampling of g(x) = f(x+s(x)).
   * @param f input array to which shifts are to be applied.
   * @param s array {s1,s2} of shifts.
   * @return array with shifts applied.
   */
  public static float[][] applyShiftsS(
    final float[][] f, final float[][][] s)
  {
    final int n1 = f[0].length;
    final int n2 = f.length;
    final float[][] s1 = s[0], s2 = s[1];
    final float[][] g = new float[n2][n1];
    final SincInterpolator si = new SincInterpolator();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        g[i2][i1] = si.interpolate(
          n1,1.0,0.0,n2,1.0,0.0,f,
          i1+s1[i2][i1],i2+s2[i2][i1]);
      }
    }});
    return g;
  }

  /**
   * Applies the specified shifts using linear interpolation.
   * The returned array is a sampling of g(x) = f(x+s(x)).
   * @param f input array to which shifts are to be applied.
   * @param s array {s1,s2} of shifts.
   * @return array with shifts applied.
   */
  public static float[][] applyShiftsSLinear(
    float[][] f, float[][][] s)
  {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] s1 = s[0], s2 = s[1];
    float[][] g = new float[n2][n1];
    LinearInterpolator li = new LinearInterpolator();
    li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li.setUniform(n1,1.0,0.0,n2,1.0,0.0,f);
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        g[i2][i1] = li.interpolate(i1+s1[i2][i1],i2+s2[i2][i1]);
    return g;
  }

  /**
   * Applies the specified shifts using sinc interpolation.
   * The returned array is a sampling of g(x) = f(x+s(x)).
   * @param f input array to which shifts are to be applied.
   * @param s array {s1,s2,s3} of shifts.
   * @return array with shifts applied.
   */
  public static float[][][] applyShiftsS(
    final float[][][] f, final float[][][][] s)
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] g = zerofloat(n1,n2,n3);
    final float[][][] s1 = s[0], s2 = s[1], s3 = s[2];
    final SincInterpolator si = new SincInterpolator();
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          g[i3][i2][i1] = si.interpolate(
            n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,f,
            i1+s1[i3][i2][i1],i2+s2[i3][i2][i1],i3+s3[i3][i2][i1]);
    }});

    return g;
  }

  /**
   * Applies the specified shifts using linear interpolation.
   * The returned array is a sampling of g(u) = f(u-r(u)).
   * @param f input array to which shifts are to be applied.
   * @param s array {s1,s2,s3} of shifts.
   * @return array with shifts applied.
   */
  public static float[][][] applyShiftsSLinear(
    float[][][] f, float[][][][] s)
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] s1 = s[0], s2 = s[1], s3 = s[2];
    final LinearInterpolator li = new LinearInterpolator();
    li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li.setUniform(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,f);
    final float[][][] g = zerofloat(n1,n2,n3);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          g[i3][i2][i1] = li.interpolate(i1+s1[i3][i2][i1],
                                         i2+s2[i3][i2][i1],
                                         i3+s3[i3][i2][i1]);
    }});
    return g;
  }

  /**
   * Gets inverse of the specified shifts.
   * @param s array {s1,s2} of shifts.
   * @return array {r1,r2} of inverse shifts.
   */
  public static float[][][] getShiftsR(float[][][] s) {
    int n1 = s[0][0].length;
    int n2 = s[0].length;
    float[][] s1 = s[0], s2 = s[1];
    float[][] r1 = new float[n2][n1];
    float[][] r2 = new float[n2][n1];
    convertS2R(s1,s2,r1,r2);
    return new float[][][]{r1,r2};
  }

  /**
   * Gets inverse of the specified shifts.
   * @param s array {s1,s2,s3} of shifts.
   * @return array {r1,r2,r3} of inverse shifts.
   */
  public static float[][][][] getShiftsR(float[][][][] s) {
    int n1 = s[0][0][0].length;
    int n2 = s[0][0].length;
    int n3 = s[0].length;
    float[][][] s1 = s[0], s2 = s[1], s3 = s[2];
    float[][][] r1 = new float[n3][n2][n1];
    float[][][] r2 = new float[n3][n2][n1];
    float[][][] r3 = new float[n3][n2][n1];
    convertS2R(s1,s2,s3,r1,r2,r3);
    return new float[][][][]{r1,r2,r3};
  }

  /**
   * Gets inverse of the specified shifts.
   * @param r array {r1,r2} of shifts.
   * @return array {s1,s2} of inverse shifts.
   */
  public static float[][][] getShiftsS(float[][][] r) {
    int n1 = r[0][0].length;
    int n2 = r[0].length;
    float[][] r1 = r[0], r2 = r[1];
    float[][] s1 = new float[n2][n1];
    float[][] s2 = new float[n2][n1];
    convertR2S(r1,r2,s1,s2);
    return new float[][][]{s1,s2};
  }

  /**
   * Gets inverse of the specified shifts.
   * @param r array {r1,r2,r3} of shifts.
   * @return array {s1,s2,s3} of inverse shifts.
   */
  public static float[][][][] getShiftsS(float[][][][] r) {
    int n1 = r[0][0][0].length;
    int n2 = r[0][0].length;
    int n3 = r[0].length;
    float[][][] r1 = r[0], r2 = r[1], r3 = r[2];
    float[][][] s1 = new float[n3][n2][n1];
    float[][][] s2 = new float[n3][n2][n1];
    float[][][] s3 = new float[n3][n2][n1];
    convertR2S(r1,r2,r3,s1,s2,s3);
    return new float[][][][]{s1,s2,s3};
  }

  /**
   * Gets partial derivatives of the specified shifts.
   * @param s array {s1,s2} of shifts.
   * @return array {s11,s12,s21,s22} of partial derivatives.
   */
  public static float[][][] getDerivatives(float[][][] s) {
    float[][][] s1 = getDerivatives(s[0]);
    float[][][] s2 = getDerivatives(s[1]);
    return new float[][][]{s1[0],s1[1],s2[0],s2[1]};
  }

  /**
   * Gets partial derivatives of the specified shifts.
   * @param s shifts.
   * @return array {s1,s2} of partial derivatives.
   */
  public static float[][][] getDerivatives(float[][] s) {
    int n1 = s[0].length;
    int n2 = s.length;
    float[][][] t = new float[2][n2][n1];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float s00 = s[i2  ][i1  ];
        float s01 = s[i2  ][i1-1];
        float s10 = s[i2-1][i1  ];
        float s11 = s[i2-1][i1-1];
        float sa = s00-s11;
        float sb = s01-s10;
        t[0][i2][i1] = 0.5f*(sa-sb);
        t[1][i2][i1] = 0.5f*(sa+sb);
      }
    }
    LinearInterpolator li = new LinearInterpolator();
    li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    for (int i=0; i<2; ++i) {
      float[][] ti = copy(n1-1,n2-1,1,1,t[i]);
      li.setUniform(n1-1,1.0,0.5,n2-1,1.0,0.5,ti);
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          t[i][i2][i1] = li.interpolate(i1,i2);
        }
      }
    }
    return t;
  }

  /**
   * Gets partial derivatives of the specified shifts.
   * @param s array {s1,s2,s3} of shifts.
   * @return array {s11,s12,s13,s21,s22,s23,s31,s32,s33} of derivatives.
   */
  public static float[][][][] getDerivatives(final float[][][][] s) {
    final int n1 = s[0][0][0].length;
    final int n2 = s[0][0].length;
    final int n3 = s[0].length;
    final float[][][] s1 = s[0], s2 = s[1], s3 = s[2];
    final float[][][][] t = new float[9][n3][n2][n1];
    final float[][][] s11 = t[0], s12 = t[1], s13 = t[2];
    final float[][][] s21 = t[3], s22 = t[4], s23 = t[5];
    final float[][][] s31 = t[6], s32 = t[7], s33 = t[8];
    Parallel.loop(1,n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=1; i2<n2; ++i2) {
        float[] s100 = s1[i3  ][i2  ];
        float[] s101 = s1[i3  ][i2-1];
        float[] s110 = s1[i3-1][i2  ];
        float[] s111 = s1[i3-1][i2-1];
        float[] s200 = s2[i3  ][i2  ];
        float[] s201 = s2[i3  ][i2-1];
        float[] s210 = s2[i3-1][i2  ];
        float[] s211 = s2[i3-1][i2-1];
        float[] s300 = s3[i3  ][i2  ];
        float[] s301 = s3[i3  ][i2-1];
        float[] s310 = s3[i3-1][i2  ];
        float[] s311 = s3[i3-1][i2-1];
        for (int i1=1, i1m=0; i1<n1; ++i1, ++i1m) {
          float s1000 = s100[i1 ];
          float s1001 = s100[i1m];
          float s1010 = s101[i1 ];
          float s1011 = s101[i1m];
          float s1100 = s110[i1 ];
          float s1101 = s110[i1m];
          float s1110 = s111[i1 ];
          float s1111 = s111[i1m];
          float s2000 = s200[i1 ];
          float s2001 = s200[i1m];
          float s2010 = s201[i1 ];
          float s2011 = s201[i1m];
          float s2100 = s210[i1 ];
          float s2101 = s210[i1m];
          float s2110 = s211[i1 ];
          float s2111 = s211[i1m];
          float s3000 = s300[i1 ];
          float s3001 = s300[i1m];
          float s3010 = s301[i1 ];
          float s3011 = s301[i1m];
          float s3100 = s310[i1 ];
          float s3101 = s310[i1m];
          float s3110 = s311[i1 ];
          float s3111 = s311[i1m];
          float s1a = s1000-s1111;
          float s1b = s1001-s1110;
          float s1c = s1010-s1101;
          float s1d = s1100-s1011;
          float s2a = s2000-s2111;
          float s2b = s2001-s2110;
          float s2c = s2010-s2101;
          float s2d = s2100-s2011;
          float s3a = s3000-s3111;
          float s3b = s3001-s3110;
          float s3c = s3010-s3101;
          float s3d = s3100-s3011;
          s11[i3][i2][i1] = 0.25f*(s1a-s1b+s1c+s1d);
          s12[i3][i2][i1] = 0.25f*(s1a+s1b-s1c+s1d);
          s13[i3][i2][i1] = 0.25f*(s1a+s1b+s1c-s1d);
          s21[i3][i2][i1] = 0.25f*(s2a-s2b+s2c+s2d);
          s22[i3][i2][i1] = 0.25f*(s2a+s2b-s2c+s2d);
          s23[i3][i2][i1] = 0.25f*(s2a+s2b+s2c-s2d);
          s31[i3][i2][i1] = 0.25f*(s3a-s3b+s3c+s3d);
          s32[i3][i2][i1] = 0.25f*(s3a+s3b-s3c+s3d);
          s33[i3][i2][i1] = 0.25f*(s3a+s3b+s3c-s3d);
        }
      }
    }});
    final float[][][] tc = new float[n3-1][n2-1][n1-1];
    final LinearInterpolator li = new LinearInterpolator();
    li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li.setUniform(n1-1,1.0,0.5,n2-1,1.0,0.5,n3-1,1.0,0.5,tc);
    for (int i=0; i<9; ++i) {
      final float[][][] ti = t[i];
      copy(n1-1,n2-1,n3-1,1,1,1,ti,0,0,0,tc);
      Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2)
          for (int i1=0; i1<n1; ++i1)
            ti[i3][i2][i1] = li.interpolate(i1,i2,i3);
      }});
      zero(tc);
    }
    return t;
  }

  /**
   * Gets the frame corresponding to the specified shifts.
   * The G-frame, as defined by Mallet, is equivalent to the 
   * columns of the Jacobian matrix, and is specified at
   * each sample by three vectors.
   * @param r array {r1,r2} of shifts.
   * @param x1 array {x11,x12} of frame vectors.
   * @param x2 array {x21,x22} of frame vectors.
   */
  public static void getFrame(
    final float[][][] r,
    final float[][][] x1, final float[][][] x2)
  {
    final int n1 = r[0][0].length;
    final int n2 = r[0].length;
    final float[][][] dr = getDerivatives(r);
    final float[][] r11 = dr[0], r12 = dr[1];
    final float[][] r21 = dr[2], r22 = dr[3];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (x1!=null) {
          x1[0][i2][i1] = 1.0f-r11[i2][i1];
          x1[1][i2][i1] = -r21[i2][i1];
        }
        if (x2!=null) {
          x2[0][i2][i1] = -r12[i2][i1];
          x2[1][i2][i1] = 1.0f-r22[i2][i1];
        }
      }
    }
  }

  /**
   * Gets the frame corresponding to the specified shifts.
   * The G-frame, as defined by Mallet, is equivalent to the 
   * columns of the Jacobian matrix, and is specified at
   * each sample by three vectors.
   * @param r array {r1,r2,r3} of shifts.
   * @param x1 array {x11,x12,x13} of frame vectors.
   * @param x2 array {x21,x22,x23} of frame vectors.
   * @param x3 array {x31,x32,x33} of frame vectors.
   */
  public static void getFrame(
    final float[][][][] r,
    final float[][][][] x1, final float[][][][] x2, final float[][][][] x3)
  {
    final int n1 = r[0][0][0].length;
    final int n2 = r[0][0].length;
    final int n3 = r[0].length;
    final float[][][][] dr = getDerivatives(r);
    final float[][][] r11 = dr[0], r12 = dr[1], r13 = dr[2];
    final float[][][] r21 = dr[3], r22 = dr[4], r23 = dr[5];
    final float[][][] r31 = dr[6], r32 = dr[7], r33 = dr[8];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          /*
          if (x1!=null) {
            x1[0][i3][i2][i1] = 1.0f-r11[i3][i2][i1];
            x1[1][i3][i2][i1] = -r12[i3][i2][i1];
            x1[2][i3][i2][i1] = -r13[i3][i2][i1];
          }
          if (x2!=null) {
            x2[0][i3][i2][i1] = -r21[i3][i2][i1];
            x2[1][i3][i2][i1] = 1.0f-r22[i3][i2][i1];
            x2[2][i3][i2][i1] = -r23[i3][i2][i1];
          }
          if (x3!=null) {
            x3[0][i3][i2][i1] = -r31[i3][i2][i1];
            x3[1][i3][i2][i1] = -r32[i3][i2][i1];
            x3[2][i3][i2][i1] = 1.0f-r33[i3][i2][i1];
          }
          */
          if (x1!=null) {
            x1[0][i3][i2][i1] = 1.0f-r11[i3][i2][i1];
            x1[1][i3][i2][i1] = -r21[i3][i2][i1];
            x1[2][i3][i2][i1] = -r31[i3][i2][i1];
          }
          if (x2!=null) {
            x2[0][i3][i2][i1] = -r12[i3][i2][i1];
            x2[1][i3][i2][i1] = 1.0f-r22[i3][i2][i1];
            x2[2][i3][i2][i1] = -r32[i3][i2][i1];
          }
          if (x3!=null) {
            x3[0][i3][i2][i1] = -r13[i3][i2][i1];
            x3[1][i3][i2][i1] = -r23[i3][i2][i1];
            x3[2][i3][i2][i1] = 1.0f-r33[i3][i2][i1];
          }
        }
      }
    }});
  }

  /**
   * Gets metric tensors for the specified shifts.
   * @param r array {r1,r2} of shifts.
   * @return array {a11,12,a22} of tensors.
   */
  public static SimpleTensors2 getMetricTensors(float[][][] r) {
    final int n1 = r[0][0].length;
    final int n2 = r[0].length;
    final float[][][] a = new float[3][n2][n1];
    final float[][][] x1 = new float[2][n2][n1];
    final float[][][] x2 = new float[2][n2][n1];
    getFrame(r,x1,x2);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float x11i = x1[0][i2][i1];
        float x12i = x1[1][i2][i1];
        float x21i = x2[0][i2][i1];
        float x22i = x2[1][i2][i1];
        a[0][i2][i1] = x11i*x11i+x12i*x12i;
        a[1][i2][i1] = x11i*x21i+x12i*x22i;
        a[2][i2][i1] = x21i*x21i+x22i*x22i;
      }
    }
    return new SimpleTensors2(a);
  }

  /**
   * Gets metric tensors for the specified shifts.
   * @param r array {r1,r2,r3} of shifts.
   * @return array {a11,12,a13,a22,a23,a33} of tensors.
   */
  public static SimpleTensors3 getMetricTensors(float[][][][] r) {
    final int n1 = r[0][0][0].length;
    final int n2 = r[0][0].length;
    final int n3 = r[0].length;
    final float[][][][] a = new float[6][n3][n2][n1];
    final float[][][][] x1 = new float[3][n3][n2][n1];
    final float[][][][] x2 = new float[3][n3][n2][n1];
    final float[][][][] x3 = new float[3][n3][n2][n1];
    getFrame(r,x1,x2,x3);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x11i = x1[0][i3][i2][i1];
          float x12i = x1[1][i3][i2][i1];
          float x13i = x1[2][i3][i2][i1];
          float x21i = x2[0][i3][i2][i1];
          float x22i = x2[1][i3][i2][i1];
          float x23i = x2[2][i3][i2][i1];
          float x31i = x3[0][i3][i2][i1];
          float x32i = x3[1][i3][i2][i1];
          float x33i = x3[2][i3][i2][i1];
          a[0][i3][i2][i1] = x11i*x11i+x12i*x12i+x13i*x13i;
          a[1][i3][i2][i1] = x11i*x21i+x12i*x22i+x13i*x23i;
          a[2][i3][i2][i1] = x11i*x31i+x12i*x32i+x13i*x33i;
          a[3][i3][i2][i1] = x21i*x21i+x22i*x22i+x23i*x23i;
          a[4][i3][i2][i1] = x21i*x31i+x22i*x32i+x23i*x33i;
          a[5][i3][i2][i1] = x31i*x31i+x32i*x32i+x33i*x33i;
        }
      }
    }});
    return new SimpleTensors3(a);
  }

  /**
   * Gets curvature tensors for the specified shifts.
   * @param r array {r1,r2,r3} of shifts.
   * @param u array {u1,u2,u3} of normal vectors.
   * @return array {a11,12,a13,a22,a23,a33} of tensors.
   */
  public static SimpleTensors3 getCurvatureTensors(
    float[][][][] r, float[][][][] u)
  {
    final int n1 = r[0][0][0].length;
    final int n2 = r[0][0].length;
    final int n3 = r[0].length;
    final float[][][][] a = new float[6][n3][n2][n1];
    final float[][][][] x1 = new float[3][n3][n2][n1];
    final float[][][][] x2 = new float[3][n3][n2][n1];
    final float[][][][] x3 = new float[3][n3][n2][n1];
    getFrame(r,x1,x2,x3);
    final float[][][][] du = getDerivatives(u);
    final float[][][] u11 = du[0], u12 = du[1], u13 = du[2];
    final float[][][] u21 = du[3], u22 = du[4], u23 = du[5];
    final float[][][] u31 = du[6], u32 = du[7], u33 = du[8];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x11i = x1[0][i3][i2][i1];
          float x12i = x1[1][i3][i2][i1];
          float x13i = x1[2][i3][i2][i1];
          float x21i = x2[0][i3][i2][i1];
          float x22i = x2[1][i3][i2][i1];
          float x23i = x2[2][i3][i2][i1];
          float x31i = x3[0][i3][i2][i1];
          float x32i = x3[1][i3][i2][i1];
          float x33i = x3[2][i3][i2][i1];
          float u11i = u11[i3][i2][i1];
          float u12i = u12[i3][i2][i1];
          float u13i = u13[i3][i2][i1];
          float u21i = u21[i3][i2][i1];
          float u22i = u22[i3][i2][i1];
          float u23i = u23[i3][i2][i1];
          float u31i = u31[i3][i2][i1];
          float u32i = u32[i3][i2][i1];
          float u33i = u33[i3][i2][i1];
          a[0][i3][i2][i1] = -(x11i*u11i+x12i*u12i+x13i*u13i);
          a[1][i3][i2][i1] = -(x11i*u21i+x12i*u22i+x13i*u23i);
          a[2][i3][i2][i1] = -(x11i*u31i+x12i*u32i+x13i*u33i);
          a[3][i3][i2][i1] = -(x21i*u21i+x22i*u22i+x23i*u23i);
          a[4][i3][i2][i1] = -(x21i*u31i+x22i*u32i+x23i*u33i);
          a[5][i3][i2][i1] = -(x31i*u31i+x32i*u32i+x33i*u33i);
        }
      }
    }});
    return new SimpleTensors3(a);
  }


  /**
   * Computes cross products of two vector fields.
   * @param x array {x1,x2,x3} of vectors.
   * @param y array {y1,y2,y3} of vectors.
   * @param z array {z1,z2,z3} of cross products.
   */
  public static void cross(
    float[][][][] x, float[][][][] y, float[][][][] z)
  {
    final int n1 = x[0][0][0].length;
    final int n2 = x[0][0].length;
    final int n3 = x[0].length;
    final float[][][] x1 = x[0], x2 = x[1], x3 = x[2];
    final float[][][] y1 = y[0], y2 = y[1], y3 = y[2];
    final float[][][] z1 = z[0], z2 = z[1], z3 = z[2];
    Parallel.loop(n3, new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x1i = x1[i3][i2][i1];
          float x2i = x2[i3][i2][i1];
          float x3i = x3[i3][i2][i1];
          float y1i = y1[i3][i2][i1];
          float y2i = y2[i3][i2][i1];
          float y3i = y3[i3][i2][i1];
          z1[i3][i2][i1] = x3i*y2i-x2i*y3i;
          z2[i3][i2][i1] = x1i*y3i-x3i*y1i;
          z3[i3][i2][i1] = x2i*y1i-x1i*y2i;
        }
      }
    }});
  }

  /**
   * Normalizes a vector field to produce a unit vector field.
   * @param x array {x1,x2} of vectors.
   * @param y array {y1,y2} of unit vectors.
   */
  public static void normalize(float[][][] x, float[][][] y) {
    final int n1 = x[0][0].length;
    final int n2 = x[0].length;
    final float[][] x1 = x[0], x2 = x[1];
    final float[][] y1 = y[0], y2 = y[1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float x1i = x1[i2][i1];
        float x2i = x2[i2][i1];
        float s = 1.0f/sqrt(x1i*x1i+x2i*x2i);
        y1[i2][i1] = s*x1i;
        y2[i2][i1] = s*x2i;
      }
    }
  }

  /**
   * Normalizes a vector field to produce a unit vector field.
   * @param x array {x1,x2,x3} of vectors.
   * @param y array {y1,y2,y3} of unit vectors.
   */
  public static void normalize(float[][][][] x, float[][][][] y) {
    final int n1 = x[0][0][0].length;
    final int n2 = x[0][0].length;
    final int n3 = x[0].length;
    final float[][][] x1 = x[0], x2 = x[1], x3 = x[2];
    final float[][][] y1 = y[0], y2 = y[1], y3 = y[2];
    Parallel.loop(n3, new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x1i = x1[i3][i2][i1];
          float x2i = x2[i3][i2][i1];
          float x3i = x3[i3][i2][i1];
          float s = 1.0f/sqrt(x1i*x1i+x2i*x2i+x3i*x3i);
          y1[i3][i2][i1] = s*x1i;
          y2[i3][i2][i1] = s*x2i;
          y3[i3][i2][i1] = s*x3i;
        }
      }
    }});
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final int MAX_ITER = 200; // max fixed-point iterations

  // Invert shifts using fixed-point iterations.
  private static void convertS2R(
    float[][] s1, float[][] s2, 
    float[][] r1, float[][] r2)
  {
    int n1 = s1[0].length;
    int n2 = s1.length;
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
        for (int iter=0; iter<MAX_ITER; ++iter) {
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
  private static void convertS2R(
    final float[][][] s1, final float[][][] s2, final float[][][] s3, 
    final float[][][] r1, final float[][][] r2, final float[][][] r3)
  {
    final int n1 = s1[0][0].length;
    final int n2 = s1[0].length;
    final int n3 = s1.length;
    final LinearInterpolator li1 = new LinearInterpolator();
    li1.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li1.setUniform(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,s1);
    final LinearInterpolator li2 = new LinearInterpolator();
    li2.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li2.setUniform(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,s2);
    final LinearInterpolator li3 = new LinearInterpolator();
    li3.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li3.setUniform(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,s3);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      double u3 = i3;
      for (int i2=0; i2<n2; ++i2) {
        double u2 = i2;
        for (int i1=0; i1<n1; ++i1) {
          double u1 = i1;
          double r1i = s1[i3][i2][i1];
          double r2i = s2[i3][i2][i1];
          double r3i = s3[i3][i2][i1];
          double r1p,r2p,r3p,dr1,dr2,dr3;
          for (int iter=0; iter<MAX_ITER; ++iter) {
            r1p = r1i;
            r2p = r2i;
            r3p = r3i;
            double x1 = u1-r1i;
            double x2 = u2-r2i;
            double x3 = u3-r3i;
            r1i = li1.interpolate(x1,x2,x3);
            r2i = li2.interpolate(x1,x2,x3);
            r3i = li3.interpolate(x1,x2,x3);
            dr1 = r1i-r1p;
            dr2 = r2i-r2p;
            dr3 = r3i-r3p;
            if (dr1*dr1+dr2*dr2+dr3*dr3<0.001) 
              break;
          }
          r1[i3][i2][i1] = (float)r1i;
          r2[i3][i2][i1] = (float)r2i;
          r3[i3][i2][i1] = (float)r3i;
        }
      }
    }});
  }
  private static void convertR2S(
    float[][] r1, float[][] r2, 
    float[][] s1, float[][] s2)
  {
    int n1 = r1[0].length;
    int n2 = r1.length;
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
        for (int iter=0; iter<MAX_ITER; ++iter) {
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

  private static void convertR2S(
    final float[][][] r1, final float[][][] r2, final float[][][] r3,
    final float[][][] s1, final float[][][] s2, final float[][][] s3)
  {
    final int n1 = r1[0][0].length;
    final int n2 = r1[0].length;
    final int n3 = r1.length;
    final LinearInterpolator li1 = new LinearInterpolator();
    li1.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li1.setUniform(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,r1);
    final LinearInterpolator li2 = new LinearInterpolator();
    li2.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li2.setUniform(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,r2);
    final LinearInterpolator li3 = new LinearInterpolator();
    li3.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li3.setUniform(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,r3);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      double x3 = i3;
      for (int i2=0; i2<n2; ++i2) {
        double x2 = i2;
        for (int i1=0; i1<n1; ++i1) {
          double x1 = i1;
          double s1i = r1[i3][i2][i1];
          double s2i = r2[i3][i2][i1];
          double s3i = r3[i3][i2][i1];
          double s1p,s2p,s3p,ds1,ds2,ds3;
          for (int iter=0; iter<MAX_ITER; ++iter) {
            s1p = s1i;
            s2p = s2i;
            s3p = s3i;
            double u1 = x1+s1i;
            double u2 = x2+s2i;
            double u3 = x3+s3i;
            s1i = li1.interpolate(u1,u2,u3);
            s2i = li2.interpolate(u1,u2,u3);
            s3i = li3.interpolate(u1,u2,u3);
            ds1 = s1i-s1p;
            ds2 = s2i-s2p;
            ds3 = s3i-s3p;
            if (ds1*ds1+ds2*ds2+ds3*ds3<0.001) 
              break;
          }
          s1[i3][i2][i1] = (float)s1i;
          s2[i3][i2][i1] = (float)s2i;
          s3[i3][i2][i1] = (float)s3i;
        }
      }
    }});
  }

}
