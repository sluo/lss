package lss.flat;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import lss.opt.*;

/**
 * Estimates shift vectors to flatten features in 2D and 3D images.
 * In 2D, the shift vector field has two components r1(x1,x3) and r3(x1,x3).
 * In 3D, the vector field has three components r1(x1,x2,x3), r2(x1,x2,x3),
 * and r3(x1,x2,x3).
 * @author Simon Luo and Dave Hale
 * @version 2012.02.01
 */
public class FlattenerR {

  /**
   * Constructs a flattener.
   * @param sigma1 smoother half-width for 1st dimension.
   * @param sigma2 smoother half-width for 2nd dimension.
   */
  public FlattenerR(double sigma1, double sigma2) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
  }

  /**
   * Constructs a flattener.
   * @param sigma1 smoother half-width for 1st dimension.
   * @param sigma2 smoother half-width for 2nd dimension.
   * @param d1 sampling interval for 1st dimension.
   * @param d2 sampling interval for 2nd dimension.
   */
  public FlattenerR(
    double sigma1, double sigma2, double d1, double d2)
  {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
    _d1 = d1;
    _d2 = d2;
  }

  /**
   * Constructs a flattener.
   * @param sigma1 smoother half-width for 1st dimension.
   * @param sigma2 smoother half-width for 2nd and 3rd dimensions.
   * @param d1 sampling interval for 1st dimension.
   * @param d2 sampling interval for 2nd dimension.
   * @param d3 sampling interval for 2nd dimension.
   */
  public FlattenerR(
    double sigma1, double sigma2, double d1, double d2, double d3)
  {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
    _d1 = d1;
    _d2 = d2;
    _d3 = d3;
  }

  /**
   * Estimates shift vectors for a 2D image.
   * @param p array of parameters {u1,u2,el,a}.
   * @return array of shifts {r1,r2}.
   */
  public float[][][] findShifts(float[][][] p) {
    return findShifts(p,null);
  }

  /**
   * Estimates shift vectors for a 2D image.
   * @param p array of parameters {u1,u2,el,a}.
   * @param res squared norm of the residuals.
   * @return array of shifts {r1,r2}.
   */
  public float[][][] findShifts(float[][][] p, float[][] res) {
    int n1 = p[0][0].length;
    int n2 = p[0].length;
    adjustNormals(_d1,_d2,p);
    float[][][] p0 = copy(p);
    float[][][] r = new float[2][n2][n1];
    float[][][] b = new float[2][n2][n1];
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    Smoother2 s2 = new Smoother2(_sigma1,_sigma2,p[2]);
    CgSolver cg = new CgSolver(_small,_inner);
    A2 ma = new A2(_epsilon,s2,p);
    for (int outer=0, inner=1; outer<_outer && inner>0; ++outer) {
      if (outer>0) {
        float[][][] t = copy(r);
        s2.apply(t);
        copy(p0,p);
        updateParameters(t,p);
        vb.zero();
      }
      makeRhs(_d1,_d2,p,b);
      s2.applyTranspose(b);
      inner = cg.solve(ma,vb,vr).niter;
    }
    if (res!=null)
      getResiduals(ma,vb,vr,res);
    s2.apply(r);
    return r;
  }

  private static void getResiduals(
    CgSolver.A a, Vec b, Vec x, float[][] res)
  {
    zero(res);
    int n1 = res[0].length;
    int n2 = res.length;
    VecArrayFloat3 q = (VecArrayFloat3)b.clone();
    VecArrayFloat3 r = (VecArrayFloat3)b.clone();
    a.apply(x,q); // q = Ax
    r.add(1.0,q,-1.0); // r = b-Ax
    float[][][] t = r.getArray();
    int n = t.length;
    for (int i=0; i<n; ++i) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float ti = t[i][i2][i1];
          res[i2][i1] += ti*ti;
        }
      }
    }
  }

  /**
   * Estimates shift vectors for a 3D image.
   * @param p array of parameters {u1,u2,u3,ep,a}.
   * @return array of shifts {r1,r2,r3}.
   */
  public float[][][][] findShifts(float[][][][] p) {
    int n1 = p[0][0][0].length;
    int n2 = p[0][0].length;
    int n3 = p[0].length;
    adjustNormals(_d1,_d2,_d3,p);
    float[][][][] p0 = scopy(p);
    float[][][][] r = new float[3][n3][n2][n1];
    float[][][][] b = new float[3][n3][n2][n1];
    VecArrayFloat4 vr = new VecArrayFloat4(r);
    VecArrayFloat4 vb = new VecArrayFloat4(b);
    Smoother3 s3 = new Smoother3(_sigma1,_sigma2,_sigma2,p[3]);
    CgSolver cg = new CgSolver(_small,_inner);
    A3 ma = new A3(_epsilon,s3,p);
    for (int outer=0, inner=1; outer<_outer && inner>0; ++outer) {
      if (outer>0) {
        float[][][][] t = scopy(r);
        s3.apply(t);
        scopy(p0,p);
        updateParameters(t,p);
        vb.zero();
      }
      makeRhs(_d1,_d2,_d3,p,b);
      s3.applyTranspose(b);
      inner = cg.solve(ma,vb,vr).niter;
    }
    s3.apply(r);
    return r;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float TP = 0.000f; // thickness preservation weight

  private double _d1 = 1.0; // sampling interval in 1st dimension
  private double _d2 = 1.0; // sampling interval in 2nd dimension
  private double _d3 = 1.0; // sampling interval in 3rd dimension
  private float _sigma1 = 8.0f; // half-width of smoother in 1st dimension
  private float _sigma2 = 8.0f; // half-width of smoother in 2nd dimension
  private float _epsilon = 0.000f; // damping for stability?
  private float _small = 0.001f; // stop CG iterations if residuals are small
  private int _inner = 10; // maximum number of inner CG iterations
  private int _outer = 20; // maximum number of outer iterations

  private static class A2 implements CgSolver.A {
    A2(float epsilon, Smoother2 smoother, float[][][] p) {
      _epsilon = epsilon;
      _smoother = smoother;
      _p = p;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      VecArrayFloat3 v3z = v3x.clone();
      v3y.zero();
      float[][][] x = v3z.getArray();
      float[][][] y = v3y.getArray();
      _smoother.apply(x);
      applyLhs(_p,x,y);
      _smoother.applyTranspose(y);
      if (_epsilon>0.0f)
        v3y.add(1.0,v3x,_epsilon*_epsilon);
    }
    private float _epsilon = 0.0f;
    private Smoother2 _smoother;
    private float[][][] _p;
  }
  private static class A3 implements CgSolver.A {
    A3(float epsilon, Smoother3 smoother, float[][][][] p) {
      _epsilon = epsilon;
      _smoother = smoother;
      _p = p;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat4 v4x = (VecArrayFloat4)vx;
      VecArrayFloat4 v4y = (VecArrayFloat4)vy;
      VecArrayFloat4 v4z = v4x.clone();
      v4y.zero();
      float[][][][] x = v4z.getArray();
      float[][][][] y = v4y.getArray();
      _smoother.apply(x);
      applyLhs(_p,x,y);
      _smoother.applyTranspose(y);
      if (_epsilon>0.0f)
        v4y.add(1.0,v4x,_epsilon*_epsilon);
    }
    private float _epsilon = 0.0f;
    private Smoother3 _smoother;
    private float[][][][] _p;
  }

  // 2D RHS
  private static void makeRhs(
    double d1, double d2, float[][][] p, float[][][] y)
  {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    float[][] y1 = y[0]; float[][] y2 = y[1];
    float[][] u1 = p[0]; float[][] u2 = p[1];
    float[][] el = p[2]; float[][] a = p[3];
    float d1d2 = (float)(d1/d2);
    float d2d1 = (float)(d2/d1);
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        float u2i = u2[i2][i1];
        float eli = el[i2][i1];
        float ai = a[i2][i1];
        float a1m = 1.0f-ai;
        float u1m = 1.0f-u1i;

        //float b1 = ai*u1m;
        //float b2 = u2i-ai*u1m*u2i;
        //float b3 = -ai*u2i;
        //float b4 = ai*u1m;
        //b1 *= W1*eli;
        //b2 *= W2*eli;
        //b3 *= W3*eli;
        //b4 *= W4*eli;
        //float y11 = b1;
        //float y12 = b2*u1i;
        //float y21 = b3;
        //float y22 = b4;
        //float y1a = 0.5f*(y11+y12);
        //float y1b = 0.5f*(y11-y12);
        //float y2a = 0.5f*(y21+y22);
        //float y2b = 0.5f*(y21-y22);

        float b1 = eli*TP*(ai*u1m);
        float b2 = eli*(ai*u2i+a1m*u2i/u1i);
        float b3 = eli*TP*(-ai*u2i);
        float b4 = eli*(ai*u1m);
        float y1a = 0.5f*(b1+b2);
        float y1b = 0.5f*(b1-b2);
        float y2a = 0.5f*(b3+b4);
        float y2b = 0.5f*(b3-b4);

        y1[i2  ][i1  ] += y1a;
        y1[i2  ][i1-1] -= y1b;
        y1[i2-1][i1  ] += y1b;
        y1[i2-1][i1-1] -= y1a;
        y2[i2  ][i1  ] += y2a;
        y2[i2  ][i1-1] -= y2b;
        y2[i2-1][i1  ] += y2b;
        y2[i2-1][i1-1] -= y2a;
      }
    }
  }

  // 2D LHS
  private static void applyLhs(
    float[][][] p, float[][][] x, float[][][] y)
  {
    int n1 = y[0][0].length;
    int n2 = y[0].length;
    float[][] el = p[2];
    float[][] x1 = x[0]; float[][] x2 = x[1];
    float[][] y1 = y[0]; float[][] y2 = y[1];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float eli = el[i2][i1];
        float x100 = x1[i2  ][i1  ];
        float x101 = x1[i2  ][i1-1];
        float x110 = x1[i2-1][i1  ];
        float x111 = x1[i2-1][i1-1];
        float x200 = x2[i2  ][i1  ];
        float x201 = x2[i2  ][i1-1];
        float x210 = x2[i2-1][i1  ];
        float x211 = x2[i2-1][i1-1];
        float x1a = x100-x111;
        float x1b = x101-x110;
        float x2a = x200-x211;
        float x2b = x201-x210;
        float x11 = 0.5f*(x1a-x1b);
        float x12 = 0.5f*(x1a+x1b);
        float x21 = 0.5f*(x2a-x2b);
        float x22 = 0.5f*(x2a+x2b);
        float b1 = x11*eli*TP; // y11;
        float b2 = x12*eli;    // y12;
        float b3 = x21*eli*TP; // y21;
        float b4 = x22*eli;    // y22;
        float y1a = 0.5f*(b1+b2);
        float y1b = 0.5f*(b1-b2);
        float y2a = 0.5f*(b3+b4);
        float y2b = 0.5f*(b3-b4);
        y1[i2  ][i1  ] += y1a;
        y1[i2  ][i1-1] -= y1b;
        y1[i2-1][i1  ] += y1b;
        y1[i2-1][i1-1] -= y1a;
        y2[i2  ][i1  ] += y2a;
        y2[i2  ][i1-1] -= y2b;
        y2[i2-1][i1  ] += y2b;
        y2[i2-1][i1-1] -= y2a;
      }
    }
  }

  private static void makeRhs(
    final double d1, final double d2, final double d3, 
    final float[][][][] p, final float[][][][] y)
  { 
    final int n3 = y[0].length;
    // i3 = 1,3,5,...
    Parallel.loop(1,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      makeRhsSlice3(i3,d1,d2,d3,p,y);
    }});
    // i3 = 2,4,6,...
    Parallel.loop(2,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      makeRhsSlice3(i3,d1,d2,d3,p,y);
    }});
  }

  private static void makeRhsSlice3(
    int i3, double d1, double d2, double d3,
    float[][][][] p, float[][][][] y)
  { 
    int n1 = y[0][0][0].length;
    int n2 = y[0][0].length;
    float d1d2 = (float)(d1/d2);
    float d1d3 = (float)(d1/d3);
    float d2d1 = (float)(d2/d1);
    float d2d3 = (float)(d2/d3);
    float d3d1 = (float)(d3/d1);
    float d3d2 = (float)(d3/d2);
    float[][][] y1 = y[0]; float[][][] y2 = y[1]; float[][][] y3 = y[2];
    float[][][] u1 = p[0]; float[][][] u2 = p[1]; float[][][] u3 = p[2];
    float[][][] ep = p[3]; float[][][] a = p[4];
    for (int i2=1; i2<n2; ++i2) {
      float[] y100 = y1[i3  ][i2  ];
      float[] y101 = y1[i3  ][i2-1];
      float[] y110 = y1[i3-1][i2  ];
      float[] y111 = y1[i3-1][i2-1];
      float[] y200 = y2[i3  ][i2  ];
      float[] y201 = y2[i3  ][i2-1];
      float[] y210 = y2[i3-1][i2  ];
      float[] y211 = y2[i3-1][i2-1];
      float[] y300 = y3[i3  ][i2  ];
      float[] y301 = y3[i3  ][i2-1];
      float[] y310 = y3[i3-1][i2  ];
      float[] y311 = y3[i3-1][i2-1];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];
        float epi = ep[i3][i2][i1];
        float ai = a[i3][i2][i1];
        float ou1p = 1.0f/(1.0f+u1i);
        float u11 = u1i*u1i;
        float u22 = u2i*u2i;
        float u23 = u2i*u3i;
        float u33 = u3i*u3i;
        float am = 1.0f-ai;
        float b1 = epi*TP*(ai-ai*u1i);      // y11
        float b2 = epi*(ai*u2i+am*u2i/u1i); // y12
        float b3 = epi*(ai*u3i+am*u3i/u1i); // y13
        float b4 = epi*TP*(-ai*u2i);        // y21
        float b5 = epi*(ai*u22*ou1p);       // y22
        float b6 = epi*(ai*u23*ou1p);       // y23
        float b7 = epi*TP*(-ai*u3i);        // y31
        float b8 = epi*(b6);                // y32
        float b9 = epi*(ai*u33*ou1p);       // y33
        float y1a = 0.25f*(b1+b2+b3);
        float y1b = 0.25f*(b1-b2+b3);
        float y1c = 0.25f*(b1+b2-b3);
        float y1d = 0.25f*(b1-b2-b3);
        float y2a = 0.25f*(b4+b5+b6);
        float y2b = 0.25f*(b4-b5+b6);
        float y2c = 0.25f*(b4+b5-b6);
        float y2d = 0.25f*(b4-b5-b6);
        float y3a = 0.25f*(b7+b8+b9);
        float y3b = 0.25f*(b7-b8+b9);
        float y3c = 0.25f*(b7+b8-b9);
        float y3d = 0.25f*(b7-b8-b9);
        y100[i1 ] += y1a;
        y100[i1m] -= y1d;
        y101[i1 ] += y1b;
        y101[i1m] -= y1c;
        y110[i1 ] += y1c;
        y110[i1m] -= y1b;
        y111[i1 ] += y1d;
        y111[i1m] -= y1a;
        y200[i1 ] += y2a;
        y200[i1m] -= y2d;
        y201[i1 ] += y2b;
        y201[i1m] -= y2c;
        y210[i1 ] += y2c;
        y210[i1m] -= y2b;
        y211[i1 ] += y2d;
        y211[i1m] -= y2a;
        y300[i1 ] += y3a;
        y300[i1m] -= y3d;
        y301[i1 ] += y3b;
        y301[i1m] -= y3c;
        y310[i1 ] += y3c;
        y310[i1m] -= y3b;
        y311[i1 ] += y3d;
        y311[i1m] -= y3a;
      }
    }
  }

  private static void applyLhs(
    final float[][][][] p, final float[][][][] x, final float[][][][] y)
  { 
    final int n3 = y[0].length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,p,x,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyLhsSlice3(i3,p,x,y);
    }});
  }

  // 3D LHS
  private static void applyLhsSlice3(
    int i3, float[][][][] p, float[][][][] x, float[][][][] y)
  {
    int n1 = y[0][0][0].length;
    int n2 = y[0][0].length;
    float[][][] ep = p[3];
    float[][][] x1 = x[0]; float[][][] x2 = x[1]; float[][][] x3 = x[2];
    float[][][] y1 = y[0]; float[][][] y2 = y[1]; float[][][] y3 = y[2];
    for (int i2=1; i2<n2; ++i2) {
      float[] x100 = x1[i3  ][i2  ];
      float[] x101 = x1[i3  ][i2-1];
      float[] x110 = x1[i3-1][i2  ];
      float[] x111 = x1[i3-1][i2-1];
      float[] x200 = x2[i3  ][i2  ];
      float[] x201 = x2[i3  ][i2-1];
      float[] x210 = x2[i3-1][i2  ];
      float[] x211 = x2[i3-1][i2-1];
      float[] x300 = x3[i3  ][i2  ];
      float[] x301 = x3[i3  ][i2-1];
      float[] x310 = x3[i3-1][i2  ];
      float[] x311 = x3[i3-1][i2-1];
      float[] y100 = y1[i3  ][i2  ];
      float[] y101 = y1[i3  ][i2-1];
      float[] y110 = y1[i3-1][i2  ];
      float[] y111 = y1[i3-1][i2-1];
      float[] y200 = y2[i3  ][i2  ];
      float[] y201 = y2[i3  ][i2-1];
      float[] y210 = y2[i3-1][i2  ];
      float[] y211 = y2[i3-1][i2-1];
      float[] y300 = y3[i3  ][i2  ];
      float[] y301 = y3[i3  ][i2-1];
      float[] y310 = y3[i3-1][i2  ];
      float[] y311 = y3[i3-1][i2-1];
      for (int i1=1,i1m=0; i1<n1; ++i1,++i1m) {
        float epi = ep[i3][i2][i1];
        float x1000 = x100[i1 ];
        float x1001 = x100[i1m];
        float x1010 = x101[i1 ];
        float x1011 = x101[i1m];
        float x1100 = x110[i1 ];
        float x1101 = x110[i1m];
        float x1110 = x111[i1 ];
        float x1111 = x111[i1m];
        float x2000 = x200[i1 ];
        float x2001 = x200[i1m];
        float x2010 = x201[i1 ];
        float x2011 = x201[i1m];
        float x2100 = x210[i1 ];
        float x2101 = x210[i1m];
        float x2110 = x211[i1 ];
        float x2111 = x211[i1m];
        float x3000 = x300[i1 ];
        float x3001 = x300[i1m];
        float x3010 = x301[i1 ];
        float x3011 = x301[i1m];
        float x3100 = x310[i1 ];
        float x3101 = x310[i1m];
        float x3110 = x311[i1 ];
        float x3111 = x311[i1m];
        float x1a = x1000-x1111;
        float x1b = x1001-x1110;
        float x1c = x1010-x1101;
        float x1d = x1100-x1011;
        float x2a = x2000-x2111;
        float x2b = x2001-x2110;
        float x2c = x2010-x2101;
        float x2d = x2100-x2011;
        float x3a = x3000-x3111;
        float x3b = x3001-x3110;
        float x3c = x3010-x3101;
        float x3d = x3100-x3011;
        float x11 = 0.0625f*(x1a-x1b+x1c+x1d);
        float x12 = 0.0625f*(x1a+x1b-x1c+x1d);
        float x13 = 0.0625f*(x1a+x1b+x1c-x1d);
        float x21 = 0.0625f*(x2a-x2b+x2c+x2d);
        float x22 = 0.0625f*(x2a+x2b-x2c+x2d);
        float x23 = 0.0625f*(x2a+x2b+x2c-x2d);
        float x31 = 0.0625f*(x3a-x3b+x3c+x3d);
        float x32 = 0.0625f*(x3a+x3b-x3c+x3d);
        float x33 = 0.0625f*(x3a+x3b+x3c-x3d);
        float b1 = x11*epi*TP; // y11
        float b2 = x12*epi;    // y12
        float b3 = x13*epi;    // y13
        float b4 = x21*epi*TP; // y21
        float b5 = x22*epi;    // y22
        float b6 = x23*epi;    // y23
        float b7 = x31*epi*TP; // y31
        float b8 = x32*epi;    // y32
        float b9 = x33*epi;    // y33
        float y1a = b1+b2+b3;
        float y1b = b1-b2+b3;
        float y1c = b1+b2-b3;
        float y1d = b1-b2-b3;
        float y2a = b4+b5+b6;
        float y2b = b4-b5+b6;
        float y2c = b4+b5-b6;
        float y2d = b4-b5-b6;
        float y3a = b7+b8+b9;
        float y3b = b7-b8+b9;
        float y3c = b7+b8-b9;
        float y3d = b7-b8-b9;
        y100[i1 ] += y1a;
        y100[i1m] -= y1d;
        y101[i1 ] += y1b;
        y101[i1m] -= y1c;
        y110[i1 ] += y1c;
        y110[i1m] -= y1b;
        y111[i1 ] += y1d;
        y111[i1m] -= y1a;
        y200[i1 ] += y2a;
        y200[i1m] -= y2d;
        y201[i1 ] += y2b;
        y201[i1m] -= y2c;
        y210[i1 ] += y2c;
        y210[i1m] -= y2b;
        y211[i1 ] += y2d;
        y211[i1m] -= y2a;
        y300[i1 ] += y3a;
        y300[i1m] -= y3d;
        y301[i1 ] += y3b;
        y301[i1m] -= y3c;
        y310[i1 ] += y3c;
        y310[i1m] -= y3b;
        y311[i1 ] += y3d;
        y311[i1m] -= y3a;
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // utilities

  private static final float MIN_U1 = 0.1f; // minimum u1 value

  private static float[][][][] scopy(float[][][][] x) {
    int n1 = x[0][0][0].length;
    int n2 = x[0][0].length;
    int n3 = x[0].length;
    int n4 = x.length;
    float[][][][] y = new float[n4][n3][n2][n1];
    scopy(x,y);
    return y;
  }

  private static void scopy(float[][][][] x, float[][][][] y) {
    int n4 = x.length;
    for (int i4=0; i4<n4; i4++)
      copy(x[i4],y[i4]);
  }

  private static void updateParameters(float[][][] r, float[][][] p) {
    int n1 = p[0][0].length;
    int n2 = p[0].length;
    int np = p.length;
    //for (int i=0; ip<n-2; ++i) { // shift normal vectors only
    for (int i=0; i<np-1; ++i) { // shift normal vectors and el
      float[][] t = FlattenerUtil.applyShiftsRLinear(p[i],r);
      copy(t,p[i]);
    }
    normalize(p[0],p[1]);
  }

  private static void updateParameters(float[][][][] r, float[][][][] p) {
    int n1 = p[0][0][0].length;
    int n2 = p[0][0].length;
    int n3 = p[0].length;
    int np = p.length;
    //for (int i=0; i<np-2; ++i) { // shift normal vectors only
    for (int i=0; i<np-1; ++i) { // shift normal vectors and ep
      float[][][] t = FlattenerUtil.applyShiftsRLinear(p[i],r);
      copy(t,p[i]);
    }
    normalize(p[0],p[1],p[2]);
  }

  // scale and limit normals
  private static void adjustNormals(
    double d1, double d2, float[][][] p)
  {
    if (d1!=1.0 || d2!=1.0)
      scaleNormals(d1,d2,p); 
    limitNormals(p); 
  }
  private static void adjustNormals(
    double d1, double d2, double d3, float[][][][] p)
  {
    if (d1!=1.0 || d2!=1.0 || d3!=1.0)
      scaleNormals(d1,d2,d3,p); 
    limitNormals(p); 
  }

  // scale normals by sampling interval
  private static void scaleNormals(
    double d1, double d2, float[][][] p)
  {
    int n1 = p[0][0].length;
    int n2 = p[0].length;
    float[][] u1 = p[0];
    float[][] u2 = p[1];
    float od1 = (float)(1.0/d1);
    float od2 = (float)(1.0/d2);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float z1i = u1[i2][i1]*od1;
        float z2i = u2[i2][i1]*od2;
        float sca = 1.0f/sqrt(z1i*z1i+z2i*z2i);
        u1[i2][i1] = z1i*sca;
        u2[i2][i1] = z2i*sca;
      }
    }
  }
  private static void scaleNormals(
    double d1, double d2, double d3, float[][][][] p)
  {
    final int n1 = p[0][0][0].length;
    final int n2 = p[0][0].length;
    final int n3 = p[0].length;
    final float[][][] u1 = p[0];
    final float[][][] u2 = p[1];
    final float[][][] u3 = p[2];
    final float od1 = (float)(1.0/d1);
    final float od2 = (float)(1.0/d2);
    final float od3 = (float)(1.0/d3);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float z1i = u1[i3][i2][i1]*od1;
          float z2i = u2[i3][i2][i1]*od2;
          float z3i = u3[i3][i2][i1]*od3;
          float sca = 1.0f/sqrt(z1i*z1i+z2i*z2i+z3i*z3i);
          u1[i3][i2][i1] = z1i*sca;
          u2[i3][i2][i1] = z2i*sca;
          u3[i3][i2][i1] = z3i*sca;
        }
      }
    }});
  }

  // limit minimum u1
  private static void limitNormals(float[][][] p) {
    int n1 = p[0][0].length;
    int n2 = p[0].length;
    float[][] u1 = p[0];
    float[][] u2 = p[1];
    float min = MIN_U1;
    float mis = min*min;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i2][i1];
        if (u1i<min) { 
          float u2i = u2[i2][i1];
          float u2s = u2i*u2i;
          float sca = sqrt((1.0f-mis)/u2s);
          u1[i2][i1] = min;
          u2[i2][i1] *= sca;
        }
      }
    }
  }
  private static void limitNormals(float[][][][] p) {
    final int n1 = p[0][0][0].length;
    final int n2 = p[0][0].length;
    final int n3 = p[0].length;
    final float[][][] u1 = p[0];
    final float[][][] u2 = p[1];
    final float[][][] u3 = p[2];
    final float min = MIN_U1;
    final float mis = min*min;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u1i = u1[i3][i2][i1];
          if (u1i<min) { 
            float u2i = u2[i3][i2][i1];
            float u3i = u3[i3][i2][i1];
            float u2s = u2i*u2i;
            float u3s = u3i*u3i;
            float sca = sqrt((1.0f-mis)/(u2s+u3s));
            u1[i3][i2][i1] = min;
            u2[i3][i2][i1] *= sca;
            u3[i3][i2][i1] *= sca;
          }
        }
      }
    }});
  }

  // Normalize vectors
  private static void normalize(float[][] a1, float[][] a2) {
    int n1 = a1[0].length;
    int n2 = a1.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float a1i = a1[i2][i1];
        float a2i = a2[i2][i1];
        float sca = 1.0f/sqrt(a1i*a1i+a2i*a2i);
        a1[i2][i1] *= sca;
        a2[i2][i1] *= sca;
      }
    }
  }
  private static void normalize(
    final float[][][] a1, final float[][][] a2, final float[][][] a3)
  {
    final int n1 = a1[0][0].length;
    final int n2 = a1[0].length;
    final int n3 = a1.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float a1i = a1[i3][i2][i1];
          float a2i = a2[i3][i2][i1];
          float a3i = a3[i3][i2][i1];
          float sca = 1.0f/sqrt(a1i*a1i+a2i*a2i+a3i*a3i);
          a1[i3][i2][i1] *= sca;
          a2[i3][i2][i1] *= sca;
          a3[i3][i2][i1] *= sca;
        }
      }
    }});
  }

}
