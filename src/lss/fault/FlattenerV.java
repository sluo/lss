package lss.fault;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import lss.opt.*;

import ipf.*;

/**
 * Image flattening using vector shifts.
 * @author Simon Luo
 * @version 2014.10.23
 */
public class FlattenerV {

  /**
   * Constructs a flattener with default parameters.
   */
  public FlattenerV() {
  }

  /**
   * Constructs a flattener.
   * @param sigma1 smoothing preconditioner half-width in 1st dimension.
   * @param sigma2 smoothing preconditioner half-width in 2nd dimension.
   * @param sigma3 smoothing preconditioner half-width in 3nd dimension.
   */
  public FlattenerV(double sigma1, double sigma2, double sigma3) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
    _sigma3 = (float)sigma3;
  }

  /**
   * Estimates shift vectors for a 3D image.
   * @param u1 1st component of normal vectors.
   * @param u2 2nd component of normal vectors.
   * @param u3 3rd component of normal vectors.
   * @param we weights for equations.
   * @return array of shifts {r1,r2,r3}.
   */
  public float[][][][] findShifts(
    float[][][] u1, float[][][] u2, float[][][] u3,
    float[][][] we)
  {
    return findShifts(u1,u2,u3,we,null);
  }

  /**
   * Estimates shift vectors for a 3D image.
   * @param u1 1st component of normal vectors.
   * @param u2 2nd component of normal vectors.
   * @param u3 3rd component of normal vectors.
   * @param we weights for equations.
   * @param fs fault skins.
   * @return array of shifts {r1,r2,r3}.
   */
  public float[][][][] findShifts(
    float[][][] u1, float[][][] u2, float[][][] u3,
    float[][][] we, FaultSkin[] fs)
  {
    int n1 = u1[0][0].length;
    int n2 = u1[0].length;
    int n3 = u1.length;
    float[][][] ws = (fs==null)?null:weightsFromFaultSkins(fs,n1,n2,n3);
    float[][][][] pc = (ws==null)?new float[][][][]{u1,u2,u3,ws}:
      new float[][][][]{u1,u2,u3,mul(we,ws),ws};
      //new float[][][][]{u1,u2,u3,we,ws};
    float[][][][] p = scopy(pc);
    float[][][][] b = new float[3][n3][n2][n1];
    float[][][][] r = new float[3][n3][n2][n1];
    VecArrayFloat4 vr = new VecArrayFloat4(r);
    VecArrayFloat4 vb = new VecArrayFloat4(b);
    CgSolver cg = new CgSolver(_small,_inner);
    A3 a3 = new A3(p);
    M3 m3 = new M3(_sigma1,_sigma2,_sigma3,(ws==null)?null:p[4]);
    for (int outer=0; outer<_outer; ++outer) {
      if (outer>0) {
        updateParameters(r,pc,p);
        vb.zero();
      }
      makeRhs(p,b);
      int inner = cg.solve(a3,m3,vb,vr).niter;
      if (inner==0) break;
    }
    return r;
  }

  private static float[][][] weightsFromFaultSkins(
    FaultSkin[] fs, int n1, int n2, int n3)
  {
    float[][][] wf = new float[n3][n2][n1];
    for (FaultSkin skin : fs) {
      for (FaultCell cell : skin) {
        //int i1 = cell.i1;
        //int i2 = cell.i2;
        //int i3 = cell.i3;
        //int i2m = cell.i2m;
        //int i3m = cell.i3m;
        //int i2p = cell.i2p;
        //int i3p = cell.i3p;
        //wf[i3m][i2m][i1] = 0.0f;
        //wf[i3 ][i2m][i1] = 0.0f;
        //wf[i3p][i2m][i1] = 0.0f;
        //wf[i3m][i2 ][i1] = 0.0f;
        //wf[i3 ][i2 ][i1] = 0.0f;
        //wf[i3p][i2 ][i1] = 0.0f;
        //wf[i3m][i2p][i1] = 0.0f;
        //wf[i3 ][i2p][i1] = 0.0f;
        //wf[i3p][i2p][i1] = 0.0f;
        int i3 = round(cell.getX3());
        int i2 = round(cell.getX2());
        int i1 = round(cell.getX1());
        wf[i3][i2][i1] = 1.0f;
      }
    }
    new RecursiveGaussianFilter(1.0).apply000(wf,wf);
    mul(1.0f/max(wf),wf,wf);
    sub(1.0f,wf,wf);
    return wf;
  }

  /**
   * Applies shifts using sinc interpolation.
   * @param r input array {r1,r2,r3} of shifts.
   * @param f input image.
   * @param g output shifted image.
   */
  public static void applyShifts(
    float[][][][] r, float[][][] f, float[][][] g)
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] ff = f;
    final float[][][] gf = g;
    final float[][][] r1 = r[0], r2 = r[1], r3 = r[2];
    final SincInterpolator si = new SincInterpolator();
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          gf[i3][i2][i1] = si.interpolate(
            n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,
            ff,i1-r1[i3][i2][i1],i2-r2[i3][i2][i1],i3-r3[i3][i2][i1]);
    }});
  }

  ///////////////////////////////////////////////////////////////////////////

  private static final float W0 = 1.000f; // thickness preservation
  private static final float W1 = 0.010f; // thickness preservation
  private static final float W2 = 0.001f; // thickness preservation

  private int _inner = 5; // max inner iterations
  private int _outer = 5; // max outer iterations
  private double _small = 0.01; // convergence criteria for CG solver
  private float _sigma1 = 6.0f; // smoother half-width in 1st dimension
  private float _sigma2 = 6.0f; // smoother half-width in 2nd dimension
  private float _sigma3 = 6.0f; // smoother half-width in 3rd dimension

  ////////////////////////////////////////////////////////////////////////////
  // Linear operator A

  private static class A3 implements CgSolver.A {
    A3(float[][][][] p) {
      _p = p;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat4 v4x = (VecArrayFloat4)vx;
      VecArrayFloat4 v4y = (VecArrayFloat4)vy;
      v4y.zero();
      float[][][][] x = v4x.getArray();
      float[][][][] y = v4y.getArray();
      applyLhs(_p,x,y);
    }
    private float[][][][] _p;
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

  private static void applyLhsSlice3(
    int i3, float[][][][] p, float[][][][] x, float[][][][] y)
  {
    int n1 = y[0][0][0].length;
    int n2 = y[0][0].length;
    int i3m = i3-1;
    float[][][] u1 = p[0];
    float[][][] u2 = p[1];
    float[][][] u3 = p[2];
    float[][][] ew = p[3];
    float[][] x1i = x[0][i3], x1m = x[0][i3m];
    float[][] x2i = x[1][i3], x2m = x[1][i3m];
    float[][] x3i = x[2][i3], x3m = x[2][i3m];
    float[][] y1i = y[0][i3], y1m = y[0][i3m];
    float[][] y2i = y[1][i3], y2m = y[1][i3m];
    float[][] y3i = y[2][i3], y3m = y[2][i3m];
    for (int i2=1, i2m=0; i2<n2; ++i2, ++i2m) {
      float[] x1ii = x1i[i2], x1im = x1i[i2m];
      float[] x1mi = x1m[i2], x1mm = x1m[i2m];
      float[] x2ii = x2i[i2], x2im = x2i[i2m];
      float[] x2mi = x2m[i2], x2mm = x2m[i2m];
      float[] x3ii = x3i[i2], x3im = x3i[i2m];
      float[] x3mi = x3m[i2], x3mm = x3m[i2m];
      float[] y1ii = y1i[i2], y1im = y1i[i2m];
      float[] y1mi = y1m[i2], y1mm = y1m[i2m];
      float[] y2ii = y2i[i2], y2im = y2i[i2m];
      float[] y2mi = y2m[i2], y2mm = y2m[i2m];
      float[] y3ii = y3i[i2], y3im = y3i[i2m];
      float[] y3mi = y3m[i2], y3mm = y3m[i2m];
      for (int i1=1, i1m=0; i1<n1; ++i1, ++i1m) {
        float x1iii = x1ii[i1], x1iim = x1ii[i1m];
        float x1imi = x1im[i1], x1imm = x1im[i1m];
        float x1mii = x1mi[i1], x1mim = x1mi[i1m];
        float x1mmi = x1mm[i1], x1mmm = x1mm[i1m];
        float x2iii = x2ii[i1], x2iim = x2ii[i1m];
        float x2imi = x2im[i1], x2imm = x2im[i1m];
        float x2mii = x2mi[i1], x2mim = x2mi[i1m];
        float x2mmi = x2mm[i1], x2mmm = x2mm[i1m];
        float x3iii = x3ii[i1], x3iim = x3ii[i1m];
        float x3imi = x3im[i1], x3imm = x3im[i1m];
        float x3mii = x3mi[i1], x3mim = x3mi[i1m];
        float x3mmi = x3mm[i1], x3mmm = x3mm[i1m];

        // Coefficients
        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];
        float ewi = ew[i3][i2][i1]*0.0625f; // 0.0625 = 0.5^4
        float u1s = u1i*u1i;
        float u1p = u1i+1.0f;
        float c22 = u1i+u2i*u2i/u1p; // alpha
        float c33 = u1i+u3i*u3i/u1p; // beta
        float c23 = -u2i*u3i/u1p; // gamma

        // Gather (1)
        float x1a = x1iii-x1mmm;
        float x1b = x1imi-x1mim;
        float x1c = x1mii-x1imm;
        float x1d = x1mmi-x1iim;
        float x2a = x2iii-x2mmm;
        float x2b = x2imi-x2mim;
        float x2c = x2mii-x2imm;
        float x2d = x2mmi-x2iim;
        float x3a = x3iii-x3mmm;
        float x3b = x3imi-x3mim;
        float x3c = x3mii-x3imm;
        float x3d = x3mmi-x3iim;

        // Gather (2) derivatives
        float x11 = x1a+x1b+x1c+x1d;
        float x12 = x1a-x1b+x1c-x1d;
        float x13 = x1a+x1b-x1c-x1d;
        float x21 = x2a+x2b+x2c+x2d;
        float x22 = x2a-x2b+x2c-x2d;
        float x23 = x2a+x2b-x2c-x2d;
        float x31 = x3a+x3b+x3c+x3d;
        float x32 = x3a-x3b+x3c-x3d;
        float x33 = x3a+x3b-x3c-x3d;

        /*
        // Gather (3) scaling
        float b1 = x11*u1i+x21*u2i+x31*u3i;
        float b2 = x12*u1i+x22*u2i+x32*u3i;
        float b3 = x13*u1i+x23*u2i+x33*u3i;

        // Weights are symmetric
        b1 *= ewi*W2;
        b2 *= ewi*W0;
        b3 *= ewi*W0;

        // Scatter (3)
        float y11 = b1*u1i;
        float y12 = b2*u1i;
        float y13 = b3*u1i;
        float y21 = b1*u2i;
        float y22 = b2*u2i;
        float y23 = b3*u2i;
        float y31 = b1*u3i;
        float y32 = b2*u3i;
        float y33 = b3*u3i;
        */

        // Gather (3) scaling
        float b1 =  x11*u1i+x21*u2i+x31*u3i;
        float b2 = -x11*u2i+x21*c33+x31*c23;
        float b3 = -x11*u3i+x21*c23+x31*c22;
        float b4 = -x12*u2i+x22*c33+x32*c23;
        float b5 = -x12*u3i+x22*c23+x32*c22;
        float b6 = -x13*u2i+x23*c33+x33*c23;
        float b7 = -x13*u3i+x23*c23+x33*c22;
        float b8 =  x12*u1i+x22*u2i+x32*u3i;
        float b9 =  x13*u1i+x23*u2i+x33*u3i;

        // Weights are symmetric
        b1 *= ewi*W2;
        b2 *= ewi*W2;
        b3 *= ewi*W2;
        b4 *= ewi*W1;
        b5 *= ewi*W1;
        b6 *= ewi*W1;
        b7 *= ewi*W1;
        b8 *= ewi*W0;
        b9 *= ewi*W0;

        // Scatter (3)
        float y11 =  b1*u1i-b2*u2i-b3*u3i;
        float y12 = -b4*u2i-b5*u3i+b8*u1i;
        float y13 = -b6*u2i-b7*u3i+b9*u1i;
        float y21 =  b1*u2i+b2*c33+b3*c23;
        float y22 =  b4*c33+b5*c23+b8*u2i;
        float y23 =  b6*c33+b7*c23+b9*u2i;
        float y31 =  b1*u3i+b2*c23+b3*c22;
        float y32 =  b4*c23+b5*c22+b8*u3i;
        float y33 =  b6*c23+b7*c22+b9*u3i;

        // Scatter (2)
        float y1a = y11+y12+y13;
        float y1b = y11-y12+y13;
        float y1c = y11+y12-y13;
        float y1d = y11-y12-y13;
        float y2a = y21+y22+y23;
        float y2b = y21-y22+y23;
        float y2c = y21+y22-y23;
        float y2d = y21-y22-y23;
        float y3a = y31+y32+y33;
        float y3b = y31-y32+y33;
        float y3c = y31+y32-y33;
        float y3d = y31-y32-y33;

        // Scatter (1)
        y1ii[i1  ] += y1a;
        y1mm[i1-1] -= y1a;
        y1im[i1  ] += y1b;
        y1mi[i1-1] -= y1b;
        y1mi[i1  ] += y1c;
        y1im[i1-1] -= y1c;
        y1mm[i1  ] += y1d;
        y1ii[i1-1] -= y1d;
        y2ii[i1  ] += y2a;
        y2mm[i1-1] -= y2a;
        y2im[i1  ] += y2b;
        y2mi[i1-1] -= y2b;
        y2mi[i1  ] += y2c;
        y2im[i1-1] -= y2c;
        y2mm[i1  ] += y2d;
        y2ii[i1-1] -= y2d;
        y3ii[i1  ] += y3a;
        y3mm[i1-1] -= y3a;
        y3im[i1  ] += y3b;
        y3mi[i1-1] -= y3b;
        y3mi[i1  ] += y3c;
        y3im[i1-1] -= y3c;
        y3mm[i1  ] += y3d;
        y3ii[i1-1] -= y3d;

      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // RHS

  private static void makeRhs(
    final float[][][][] p, final float[][][][] y)
  { 
    final int n3 = y[0].length;
    Parallel.loop(1,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      makeRhsSlice3(i3,p,y);
    }});
    Parallel.loop(2,n3,2,new Parallel.LoopInt() {
    public void compute(int i3) {
      makeRhsSlice3(i3,p,y);
    }});
  }

  private static void makeRhsSlice3(
    int i3, float[][][][] p, float[][][][] y)
  {
    int n1 = y[0][0][0].length;
    int n2 = y[0][0].length;
    int i3m = i3-1;
    float[][][] u1 = p[0];
    float[][][] u2 = p[1];
    float[][][] u3 = p[2];
    float[][][] ew = p[3];
    float[][] y1i = y[0][i3], y1m = y[0][i3m];
    float[][] y2i = y[1][i3], y2m = y[1][i3m];
    float[][] y3i = y[2][i3], y3m = y[2][i3m];
    for (int i2=1, i2m=0; i2<n2; ++i2, ++i2m) {
      float[] y1ii = y1i[i2], y1im = y1i[i2m];
      float[] y1mi = y1m[i2], y1mm = y1m[i2m];
      float[] y2ii = y2i[i2], y2im = y2i[i2m];
      float[] y2mi = y2m[i2], y2mm = y2m[i2m];
      float[] y3ii = y3i[i2], y3im = y3i[i2m];
      float[] y3mi = y3m[i2], y3mm = y3m[i2m];
      for (int i1=1, i1m=0; i1<n1; ++i1, ++i1m) {

        // Coefficients
        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];
        float ewi = ew[i3][i2][i1]*0.25f; // 0.25 = 0.5^2
        float u1s = u1i*u1i;
        float u1p = u1i+1.0f;
        float c22 = u1i+u2i*u2i/u1p; // alpha
        float c33 = u1i+u3i*u3i/u1p; // beta
        float c23 = -u2i*u3i/u1p; // gamma

        /*
        // RHS vector b
        float b1 = u1i-1.0f;
        float b2 = u2i;
        float b3 = u3i;

        // Weights are symmetric
        b1 *= ewi*W2;
        b2 *= ewi*W0;
        b3 *= ewi*W0;

        // Scatter (3)
        float y11 = b1*u1i;
        float y12 = b2*u1i;
        float y13 = b3*u1i;
        float y21 = b1*u2i;
        float y22 = b2*u2i;
        float y23 = b3*u2i;
        float y31 = b1*u3i;
        float y32 = b2*u3i;
        float y33 = b3*u3i;
        */

        // RHS vector b
        float b1 =  u1i-1.0f;
        float b2 = -u2i;
        float b3 = -u3i;
        float b4 =  c33-1.0f;
        float b5 =  c23;
        float b6 =  c23;
        float b7 =  c22-1.0f;
        float b8 =  u2i;
        float b9 =  u3i;

        // Weights are symmetric
        b1 *= ewi*W2;
        b2 *= ewi*W2;
        b3 *= ewi*W2;
        b4 *= ewi*W1;
        b5 *= ewi*W1;
        b6 *= ewi*W1;
        b7 *= ewi*W1;
        b8 *= ewi*W0;
        b9 *= ewi*W0;

        // Scatter (3)
        float y11 =  b1*u1i-b2*u2i-b3*u3i;
        float y12 = -b4*u2i-b5*u3i+b8*u1i;
        float y13 = -b6*u2i-b7*u3i+b9*u1i;
        float y21 =  b1*u2i+b2*c33+b3*c23;
        float y22 =  b4*c33+b5*c23+b8*u2i;
        float y23 =  b6*c33+b7*c23+b9*u2i;
        float y31 =  b1*u3i+b2*c23+b3*c22;
        float y32 =  b4*c23+b5*c22+b8*u3i;
        float y33 =  b6*c23+b7*c22+b9*u3i;

        // Scatter (2)
        float y1a = y11+y12+y13;
        float y1b = y11-y12+y13;
        float y1c = y11+y12-y13;
        float y1d = y11-y12-y13;
        float y2a = y21+y22+y23;
        float y2b = y21-y22+y23;
        float y2c = y21+y22-y23;
        float y2d = y21-y22-y23;
        float y3a = y31+y32+y33;
        float y3b = y31-y32+y33;
        float y3c = y31+y32-y33;
        float y3d = y31-y32-y33;

        // Scatter (1)
        y1ii[i1  ] += y1a;
        y1mm[i1-1] -= y1a;
        y1im[i1  ] += y1b;
        y1mi[i1-1] -= y1b;
        y1mi[i1  ] += y1c;
        y1im[i1-1] -= y1c;
        y1mm[i1  ] += y1d;
        y1ii[i1-1] -= y1d;
        y2ii[i1  ] += y2a;
        y2mm[i1-1] -= y2a;
        y2im[i1  ] += y2b;
        y2mi[i1-1] -= y2b;
        y2mi[i1  ] += y2c;
        y2im[i1-1] -= y2c;
        y2mm[i1  ] += y2d;
        y2ii[i1-1] -= y2d;
        y3ii[i1  ] += y3a;
        y3mm[i1-1] -= y3a;
        y3im[i1  ] += y3b;
        y3mi[i1-1] -= y3b;
        y3mi[i1  ] += y3c;
        y3im[i1-1] -= y3c;
        y3mm[i1  ] += y3d;
        y3ii[i1-1] -= y3d;

      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // Preconditioner M

  private static class M3 implements CgSolver.A {
    M3(float sigma1, float sigma2, float sigma3, float[][][] ew) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _sigma3 = sigma3;
      _ew = ew;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat4 v4x = (VecArrayFloat4)vx;
      VecArrayFloat4 v4y = (VecArrayFloat4)vy; v4y.zero();
      float[][][][] x = v4x.getArray();
      float[][][][] y = v4y.getArray();
      for (int i4=0; i4<3; ++i4) {
        float[][][] xi = x[i4];
        float[][][] yi = y[i4];
        copy(xi,yi);
        removeAverage(yi);
        smooth3(_sigma3,_ew,yi);
        smooth2(_sigma2,_ew,yi);
        smooth1(sqrt(2.0f)*_sigma1,_ew,yi);
        smooth2(_sigma2,_ew,yi);
        smooth3(_sigma3,_ew,yi);
        removeAverage(yi);
      }
    }
    float[][][] _ew; // equation weights
    float _sigma1,_sigma2,_sigma3;
  }

  public static void constrain(
    int[][] k1, int[][] k2, int[][] k3, float[][][] x)
  {
    if (k1!=null && k2!=null &&k3!=null) {
      int nc = k1.length;
      for (int ic=0; ic<nc; ++ic) {
        int nk = k1[ic].length;
        float sum = 0.0f;
        for (int ik=0; ik<nk; ++ik) {
          int i1 = k1[ic][ik];
          int i2 = k2[ic][ik];
          int i3 = k3[ic][ik];
          sum += x[i3][i2][i1];
        }
        float avg = sum/(float)nk;
        for (int ik=0; ik<nk; ++ik) {
          int i1 = k1[ic][ik];
          int i2 = k2[ic][ik];
          int i3 = k3[ic][ik];
          x[i3][i2][i1] = avg;
        }
      }
    }
  }

  private static void smooth1(
    final float sigma, final float[][][] s, final float[][][] x)
  {
    final int n3 = x.length;
    final int n2 = x[0].length;
    Parallel.loop(n3, new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] x3 = x[i3];
      float[][] s3 = (s!=null)?s[i3]:null;
      smooth1(sigma,s3,x3);
    }});
  }
 
  private static void smooth1(float sigma, float[][] s, float[][] x) {
    if (sigma<1.0f)
      return;
    int n2 = x.length;
    int n1 = x[0].length;
    float c = 0.5f*sigma*sigma;
    float[] st = fillfloat(1.0f,n1);
    float[] xt = zerofloat(n1);
    float[] yt = zerofloat(n1);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i2=0; i2<n2; ++i2) {
      if (s!=null) {
        for (int i1=0; i1<n1; ++i1)
          st[i1] = s[i2][i1];
      }
      for (int i1=0; i1<n1; ++i1)
        xt[i1] = x[i2][i1];
      lsf.apply(c,st,xt,yt);
      for (int i1=0; i1<n1; ++i1)
        x[i2][i1] = yt[i1];
    }
  }

  private static void smooth2(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] s3 = (s!=null)?s[i3]:null;
      float[][] x3 = x[i3];
      smooth2(sigma,s3,x3);
    }});
  }

  private static void smooth2(float sigma, float[][] s, float[][] x) {
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] st = fillfloat(1.0f,n2);
    float[] xt = zerofloat(n2);
    float[] yt = zerofloat(n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i1=0; i1<n1; ++i1) {
      if (s!=null) {
        for (int i2=0; i2<n2; ++i2)
          st[i2] = s[i2][i1];
      }
      for (int i2=0; i2<n2; ++i2)
        xt[i2] = x[i2][i1];
      lsf.apply(c,st,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }
  }

  private static void smooth3(
    final float sigma, final float[][][] s, final float[][][] x) 
  {
    final int n2 = x[0].length;
    final int n3 = x.length;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] s2 = (s!=null)?new float[n3][]:null;
      float[][] x2 = new float[n3][];
      for (int i3=0; i3<n3; ++i3) {
        if (s!=null)
          s2[i3] = s[i3][i2];
        x2[i3] = x[i3][i2];
      }
      smooth2(sigma,s2,x2);
    }});
  }

  private static void removeAverage(float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float nh = (float)(n2*n3);
    for (int i1=0; i1<n1; ++i1) {
      float sumx = 0.0f;
      for (int i3=0; i3<n3; ++i3)  
        for (int i2=0; i2<n2; ++i2)  
          sumx += x[i3][i2][i1];
      float avgx = sumx/nh;
      for (int i3=0; i3<n3; ++i3) 
        for (int i2=0; i2<n2; ++i2) 
          x[i3][i2][i1] -= avgx; 
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  
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

  private static void updateParameters(
    float[][][][] r, float[][][][] p, float[][][][] q)
  {
    int n1 = p[0][0][0].length;
    int n2 = p[0][0].length;
    int n3 = p[0].length;
    int n4 = p.length;
    for (int i=0; i<n4; ++i)
      applyShiftsLinear(r,p[i],q[i]);
    normalize(q[0],q[1],q[2]);
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

  private static void applyShiftsLinear(
    float[][][][] r, float[][][] f, float[][][] g)
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] gf = g;
    final float[][][] r1 = r[0], r2 = r[1], r3 = r[2];
    final LinearInterpolator li = new LinearInterpolator();
    li.setExtrapolation(LinearInterpolator.Extrapolation.CONSTANT);
    li.setUniform(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,f);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          gf[i3][i2][i1] = li.interpolate(
            i1-r1[i3][i2][i1],i2-r2[i3][i2][i1],i3-r3[i3][i2][i1]);
    }});
  }

  ///////////////////////////////////////////////////////////////////////////
  // Adjoint tests

  public static void main(String[] args) {
    int n1 = 10;
    int n2 = 11;
    int n3 = 12;
    java.util.Random random = new java.util.Random(12345);
    float[][][][] p = rfloat(random,n1,n2,n3,4);
    float[][][][] xa = rfloat(random,n1,n2,n3,3);
    float[][][][] xb = rfloat(random,n1,n2,n3,3);
    float[][][][] ya = zfloat(n1,n2,n3,3);
    float[][][][] yb = zfloat(n1,n2,n3,3);
    applyLhs(p,xa,ya);
    applyLhs(p,xb,yb);
    System.out.println(dot(xa,yb));
    System.out.println(dot(xb,ya));
  }

  private static float[][][][] rfloat(
    java.util.Random random, int n1, int n2, int n3, int n4)
  {
    float[][][][] y = new float[n4][][][];
    for (int i4=0; i4<n4; ++i4) {
      y[i4] = randfloat(random,n1,n2,n3);
    }
    return y;
  }

  private static float[][][][] zfloat(int n1, int n2, int n3, int n4) {
    float[][][][] y = new float[n4][][][];
    for (int i4=0; i4<n4; ++i4) {
      y[i4] = zerofloat(n1,n2,n3);
    }
    return y;
  }

  private static float dot(float[][][][] x, float[][][][] y) {
    int n4 = x.length;
    float s = 0.0f;
    for (int i4=0; i4<n4; ++i4)
      s += sum(mul(x[i4],y[i4]));
    return s;
  }

}
