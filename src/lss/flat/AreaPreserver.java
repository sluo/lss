package lss.flat;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import dnp.*;
import lss.*;

/**
 * Computes additional horizontal shifts for 
 * area-preserving flattening.
 * @author Simon Luo 
 */
public class AreaPreserver {

  public AreaPreserver(double sigma1, double sigma2) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
  }

  public float[][] findShifts(float[][] s1, float[][] el) {
    int n1 = s1[0].length;
    int n2 = s1.length;
    float[][] r = new float[n2][n1];
    float[][] s = new float[n2][n1];
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    VecArrayFloat2 vs = new VecArrayFloat2(s);
    Smoother2 smoother = new Smoother2(_sigma1,_sigma2,el);
    A2 a2 = new A2(_epsilon,smoother,s1,el);
    CgSolver cg = new CgSolver(_small,_niter);
    makeRhs(s1,el,r);
    smoother.applyTranspose(r);
    cg.solve(a2,vr,vs);
    smoother.apply(s);
    return s;
  }

  public float[][] estimateDeformation(
    final float[][] s1, final float[][] el)
  {
    final int n1 = s1[0].length;
    final int n2 = s1.length;
    final float[][] d = new float[n2][n1];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute (int i2) {
      float[] t = copy(s1[i2]);
      float[][] r1 = copy(s1);
      for (int j2=0; j2<n2; ++j2)
        sub(r1[j2],t,r1[j2]);
      float[][] r2 = findShifts(r1,el);
      float[][] c = getCurl(r1,r2);
      copy(c[i2],d[i2]);
    }});
    return d;
  }

  public static float[][] getCurl(float[][] s1, float[][] s2) {
    int n1 = s1[0].length;
    int n2 = s1.length;
    float[][][] ds = FlattenerVS.getDerivatives(new float[][][]{s1,s2});
    float[][] s12 = ds[1], s21 = ds[2];
    float[][] c = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        c[i2][i1] = s21[i2][i1]-s12[i2][i1];
    return c;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma1 = 6.0f; // half-width of smoother in 1st dimension
  private float _sigma2 = 6.0f; // half-width of smoother in 2nd dimension
  private float _epsilon = 0.000f; // damping for stability?
  private float _small = 0.01f; // stop CG iterations if residuals are small
  private int _niter = 500; // maximum number of CG iterations

  private static void transpose(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        y[i1][i2] = x[i2][i1];
      }
    }
  }

  private static class A2 implements CgSolver.A {
    A2(float epsilon, Smoother2 smoother, float[][] s1, float[][] el) {
      _epsilon = epsilon;
      _smoother = smoother;
      _s1 = s1;
      _el = el;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      VecArrayFloat2 v2z = v2x.clone();
      v2y.zero();
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = v2z.getArray();
      _smoother.apply(z);
      applyLhs(_s1,_el,z,y);
      _smoother.applyTranspose(y);
      if (_epsilon>0.0f)
        v2y.add(1.0,v2x,_epsilon*_epsilon);
    }
    private float _epsilon;
    private Smoother2 _smoother;
    private float[][] _s1;
    private float[][] _el;
  }

  private static void makeRhs(
    float[][] s, float[][] el, float[][] y)
  {
    int n1 = y[0].length;
    int n2 = y.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float eli = el[i2][i1];
        float s00 = s[i2  ][i1  ];
        float s01 = s[i2  ][i1-1];
        float s10 = s[i2-1][i1  ];
        float s11 = s[i2-1][i1-1];
        float sa = s00-s11;
        float sb = s01-s10;
        float s1 = 0.5f*(sa-sb);
        float s2 = 0.5f*(sa+sb);
        float b1 = -s1*eli;
        float y1 = b1*s2;
        float y2 = b1*(1.0f-s1);
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  private static void applyLhs(
    float[][] s, float[][] el, float[][] x, float[][] y)
  {
    int n1 = y[0].length;
    int n2 = y.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float eli = el[i2][i1];

        float x00 = x[i2  ][i1  ];
        float x01 = x[i2  ][i1-1];
        float x10 = x[i2-1][i1  ];
        float x11 = x[i2-1][i1-1];
        float xa = x00-x11;
        float xb = x01-x10;
        float x1 = 0.5f*(xa-xb);
        float x2 = 0.5f*(xa+xb);

        float s00 = s[i2  ][i1  ];
        float s01 = s[i2  ][i1-1];
        float s10 = s[i2-1][i1  ];
        float s11 = s[i2-1][i1-1];
        float sa = s00-s11;
        float sb = s01-s10;
        float s1 = 0.5f*(sa-sb);
        float s2 = 0.5f*(sa+sb);
        float s1m1 = 1.0f-s1;

        float b1 = (s2*x1+s1m1*x2)*eli;

        float y1 = b1*s2;
        float y2 = b1*s1m1;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }

  }

}
