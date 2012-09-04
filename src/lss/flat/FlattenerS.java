package lss.flat;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import dnp.*;

/**
 * FlattenerCg with some tweaks.
 * @author Simon Luo and Dave Hale
 * @version 2012.03.25
 */
public class FlattenerS {

  public FlattenerS() {
  }

  public FlattenerS(double sigma1, double sigma2) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
  }

  public float[][] findShifts(float[][] p2) {
    return findShifts(p2,null);
  }

  public float[][] findShifts(float[][] p2, float[][] el) {
    int n1 = p2[0].length;
    int n2 = p2.length;
    float[][] r = new float[n2][n1]; // right-hand side
    float[][] s = new float[n2][n1]; // the shifts
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    VecArrayFloat2 vs = new VecArrayFloat2(s);
    Smoother2 s2 = new Smoother2(_sigma1,_sigma2,el);
    A2 a2 = new A2(_epsilon,s2,p2,el);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(p2,el,r);
    s2.applyTranspose(r);
    cs.solve(a2,vr,vs);
    s2.apply(s);
    invertShifts(s);
    return s;
  }

  public float[][][] findShifts(
    float[][][] p2, float[][][] p3, float[][][] ep) 
  {
    int n1 = p2[0][0].length;
    int n2 = p2[0].length;
    int n3 = p2.length;
    float[][][] r = new float[n3][n2][n1]; // right-hand side
    float[][][] s = new float[n3][n2][n1]; // the shifts
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    VecArrayFloat3 vs = new VecArrayFloat3(s);
    Smoother3 s3 = new Smoother3(_sigma1,_sigma2,_sigma2,ep);
    A3 a3 = new A3(_epsilon,s3,p2,p3,ep);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(p2,p3,ep,r);
    s3.applyTranspose(r);
    cs.solve(a3,vr,vs);
    s3.apply(s);
    invertShifts(s);
    return s;
  }

  public static float[][] applyShifts(float[][] f, float[][] s) {
    int n1 = f[0].length;
    int n2 = f.length;
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(n1,1.0,0.0);
    float[] r = rampfloat(0.0f,1.0f,n1);
    float[] t = zerofloat(n1);
    float[][] g = zerofloat(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      sub(r,s[i2],t);
      si.setUniformSamples(f[i2]);
      si.interpolate(n1,t,g[i2]);
    }
    return g;
  }

  public static float[][] applyShiftsLinear(float[][] f, float[][] s) {
    int n1 = f[0].length;
    int n2 = f.length;
    LinearInterpolator li = new LinearInterpolator();
    li.setUniformSampling(n1,1.0,0.0);
    float[] r = rampfloat(0.0f,1.0f,n1);
    float[] t = zerofloat(n1);
    float[][] g = zerofloat(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      sub(r,s[i2],t);
      li.setUniformSamples(f[i2]);
      li.interpolate(n1,t,g[i2]);
    }
    return g;
  }

  public static float[][][] applyShifts(float[][][] f, float[][][] s) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(n1,1.0,0.0);
    float[] r = rampfloat(0.0f,1.0f,n1);
    float[] t = zerofloat(n1);
    float[][][] g = zerofloat(n1,n2,n3);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        sub(r,s[i3][i2],t);
        si.setUniformSamples(f[i3][i2]);
        si.interpolate(n1,t,g[i3][i2]);
      }
    }
    return g;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final boolean PARALLEL = true; // true if multi-threaded
  private static final float PRESERVE = 0.000f; // thickness preservation

  private float _sigma1 = 6.0f; // half-width of smoother in 1st dimension
  private float _sigma2 = 6.0f; // half-width of smoother in 2nd dimension
  private float _epsilon = 0.000f; // damping for stability?
  private float _small = 0.001f; // stop CG iterations if residuals are small
  private int _niter = 200; // maximum number of CG iterations

  // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(float epsilon, Smoother2 s2, float[][] p2, float[][] el) {
      _epsilon = epsilon;
      _s2 = s2;
      _p2 = p2;
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
      _s2.apply(z);
      applyLhs(_p2,_el,z,y);
      _s2.applyTranspose(y);
      if (_epsilon>0.0f)
        v2y.add(1.0,v2x,_epsilon*_epsilon);
    }
    private float _epsilon;
    private Smoother2 _s2;
    private float[][] _p2;
    private float[][] _el;
  }
  private static class A3 implements CgSolver.A {
    A3(float epsilon, Smoother3 s3,
       float[][][] p2, float[][][] p3, float[][][] ep) 
    {
      _epsilon = epsilon;
      _s3 = s3;
      _p2 = p2;
      _p3 = p3;
      _ep = ep;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      VecArrayFloat3 v3z = v3x.clone();
      v3y.zero();
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      float[][][] z = v3z.getArray();
      _s3.apply(z);
      applyLhs(_p2,_p3,_ep,z,y);
      _s3.applyTranspose(y);
      if (_epsilon>0.0f)
        v3y.add(1.0,v3x,_epsilon*_epsilon);
    }
    private float _epsilon;
    private Smoother3 _s3;
    private float[][][] _p2;
    private float[][][] _p3;
    private float[][][] _ep;
  }

  private static void makeRhs(float[][] p2, float[][] el, float[][] y) {
    int n1 = y[0].length;
    int n2 = y.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float eli = (el!=null)?el[i2][i1]:1.0f;
        float p2i = p2[i2][i1];
        float b12 = p2i*eli;
        float b22 = eli;
        float x2 = -0.5f*p2i*eli;
        float y1 = b12*x2;
        float y2 = b22*x2;
        float ya = y1+y2;
        float yb = y1-y2;
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }
  private static void applyLhs(
    float[][] p2, float[][] el, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float eli = (el!=null)?el[i2][i1]:1.0f;
        float p2i = p2[i2][i1];
        float els = eli*eli;
        float p2s = p2i*p2i;
        float d11 = (p2s+PRESERVE)*els;
        float d12 = p2i*els;
        float d22 = els;
        float x00 = x[i2  ][i1  ];
        float x01 = x[i2  ][i1-1];
        float x10 = x[i2-1][i1  ];
        float x11 = x[i2-1][i1-1];
        float xa = x00-x11;
        float xb = x01-x10;
        float x1 = 0.25f*(xa-xb);
        float x2 = 0.25f*(xa+xb);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = y1+y2;
        float yb = y1-y2;
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  private static void makeRhs(
    float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] y) 
  {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    int n3 = y.length;
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float epi = (ep!=null)?ep[i3][i2][i1]:1.0f;
          float p2i = p2[i3][i2][i1];
          float p3i = p3[i3][i2][i1];
          float b12 = p2i*epi;
          float b13 = p3i*epi;
          float b22 = epi;
          float b33 = epi;
          float x2 = -0.25f*p2i*epi;
          float x3 = -0.25f*p3i*epi;
          float y1 = b12*x2+b13*x3;
          float y2 = b22*x2;
          float y3 = b33*x3;
          float ya = y1+y2+y3;
          float yb = y1-y2+y3;
          float yc = y1+y2-y3;
          float yd = y1-y2-y3;
          y[i3  ][i2  ][i1  ] += ya;
          y[i3  ][i2  ][i1-1] -= yd;
          y[i3  ][i2-1][i1  ] += yb;
          y[i3  ][i2-1][i1-1] -= yc;
          y[i3-1][i2  ][i1  ] += yc;
          y[i3-1][i2  ][i1-1] -= yb;
          y[i3-1][i2-1][i1  ] += yd;
          y[i3-1][i2-1][i1-1] -= ya;
        }
      }
    }
  }
  private static void applyLhs(
    float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] x, float[][][] y) 
  {
    if (PARALLEL) {
      applyLhsP(p2,p3,ep,x,y);
    } else {
      applyLhsS(p2,p3,ep,x,y);
    }
  }

  private static void applyLhsS(
    float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] x, float[][][] y) 
  {
    int n3 = x.length;
    for (int i3=1; i3<n3; ++i3)
      applyLhsSlice3(i3,p2,p3,ep,x,y);
  }

  private static void applyLhsP(
    final float[][][] p2, final float[][][] p3, final float[][][] ep, 
    final float[][][] x, final float[][][] y) 
  {
    final int n3 = x.length;

    // i3 = 1,3,5,...
    Parallel.loop(1,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,p2,p3,ep,x,y);
    }});

    // i3 = 2,4,6,...
    Parallel.loop(2,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,p2,p3,ep,x,y);
    }});
  }

  private static void applyLhsSlice3(
    int i3,
    float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] x, float[][][] y) 
  {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    for (int i2=1; i2<n2; ++i2) {
      float[] x00 = x[i3  ][i2  ];
      float[] x01 = x[i3  ][i2-1];
      float[] x10 = x[i3-1][i2  ];
      float[] x11 = x[i3-1][i2-1];
      float[] y00 = y[i3  ][i2  ];
      float[] y01 = y[i3  ][i2-1];
      float[] y10 = y[i3-1][i2  ];
      float[] y11 = y[i3-1][i2-1];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        float epi = (ep!=null)?ep[i3][i2][i1]:1.0f;
        float p2i = p2[i3][i2][i1];
        float p3i = p3[i3][i2][i1];
        float eps = epi*epi;
        float p2s = p2i*p2i;
        float p3s = p3i*p3i;
        float d11 = (p2s+p3s+PRESERVE)*eps;
        float d12 = p2i*eps;
        float d13 = p3i*eps;
        float d22 = eps;
        float d33 = eps;
        float x000 = x00[i1 ];
        float x001 = x00[i1m];
        float x010 = x01[i1 ];
        float x011 = x01[i1m];
        float x100 = x10[i1 ];
        float x101 = x10[i1m];
        float x110 = x11[i1 ];
        float x111 = x11[i1m];
        float xa = x000-x111;
        float xb = x001-x110;
        float xc = x010-x101;
        float xd = x100-x011;
        float x1 = 0.0625f*(xa-xb+xc+xd);
        float x2 = 0.0625f*(xa+xb-xc+xd);
        float x3 = 0.0625f*(xa+xb+xc-xd);
        float y1 = d11*x1+d12*x2+d13*x3;
        float y2 = d12*x1+d22*x2       ;
        float y3 = d13*x1       +d33*x3;
        float ya = y1+y2+y3;
        float yb = y1-y2+y3;
        float yc = y1+y2-y3;
        float yd = y1-y2-y3;
        y00[i1 ] += ya;
        y00[i1m] -= yd;
        y01[i1 ] += yb;
        y01[i1m] -= yc;
        y10[i1 ] += yc;
        y10[i1m] -= yb;
        y11[i1 ] += yd;
        y11[i1m] -= ya;
      }
    }
  }

  // Post-processing and inversion of computed shifts.
  private static void cleanShifts(float[] s) {
    int n1 = s.length;
    for (int i1=1; i1<n1; ++i1) {
      if (s[i1]<=s[i1-1]-0.99f)
        s[i1] = s[i1-1]-0.99f;
    }
  }
  private static void invertShifts(
    InverseInterpolator ii, float[] u, float[] t, float[] s) 
  {
    cleanShifts(s);
    int n1 = s.length;
    for (int i1=0; i1<n1; ++i1)
      s[i1] += u[i1];
    ii.invert(s,t);
    float tmin = -n1;
    float tmax = n1+n1;
    for (int i1=0; i1<n1; ++i1) {
      if (t[i1]<tmin) t[i1] = tmin;
      if (t[i1]>tmax) t[i1] = tmax;
      s[i1] = u[i1]-t[i1];
    }
  }
  private static void invertShifts(float[][] s) {
    int n1 = s[0].length;
    int n2 = s.length;
    float[] u = rampfloat(0.0f,1.0f,n1);
    float[] t = zerofloat(n1);
    InverseInterpolator ii = new InverseInterpolator(n1,n1);
    for (int i2=0; i2<n2; ++i2)
      invertShifts(ii,u,t,s[i2]);
  }
  private static void invertShifts(float[][][] s) {
    int n1 = s[0][0].length;
    int n2 = s[0].length;
    int n3 = s.length;
    float[] u = rampfloat(0.0f,1.0f,n1);
    float[] t = zerofloat(n1);
    InverseInterpolator ii = new InverseInterpolator(n1,n1);
    for (int i3=0; i3<n3; ++i3)
      for (int i2=0; i2<n2; ++i2)
        invertShifts(ii,u,t,s[i3][i2]);
  }
}
