package lss.dev;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Forward, backward, and adjoint 1D wave operator.
 * <p>
 * The forward operator models data from a model, while the backward
 * and adjoint operators migrate data to produce an image of the model.
 * In this case, the model is not a velocity model, but rather the
 * initial condition of the wavefield at time t = 0.
 * <p>
 * See Ji, J., 2009, An exact adjoint operation pair in time extrapolation
 * and its application in least-squares reverse-time migration.
 *
 * @author Simon Luo
 * @version 2013.07.11
 */
public class WaveOperator {

//      // Forward.
//      for (int iz=1; iz<nz-1; ++iz)
//        up[iz] = a[iz]*(ui[iz+1]-2.0f*ui[iz]+ui[iz-1])-um[iz]+2.0f*ui[iz];
//      up[0   ] = a[0   ]*(ui[1   ]-2.0f*ui[0   ])-um[0   ]+2.0f*ui[0   ];
//      up[nz-1] = a[nz-1]*(ui[nz-2]-2.0f*ui[nz-1])-um[nz-1]+2.0f*ui[nz-1];
//
//      // Adjoint.
//      ui[0   ] += a[1   ]*um[1   ]-2.0f*a[0   ]*um[0   ]+2.0f*um[0   ];
//      ui[nz-1] += a[nz-2]*um[nz-2]-2.0f*a[nz-1]*um[nz-1]+2.0f*um[nz-1];
//      for (int iz=0; iz<nz; ++iz)
//        up[iz] -= um[iz];
//      for (int iz=1; iz<nz-1; ++iz)
//        ui[iz] += a[iz-1]*um[iz-1]-2.0f*a[iz]*um[iz]+a[iz+1]*um[iz+1]
//          +2.0f*um[iz];

  public WaveOperator(Sampling sz, Sampling st, float[] s) {
    _nz = sz.getCount();
    _nt = st.getCount();
    _dz = sz.getDelta();
    _dt = st.getDelta();
    _a = div((float)(_dt*_dt/(_dz*_dz)),mul(s,s));
  }

  /**
   * Applies the forward operator.
   */
  public float[] applyForward(float[] m) {
    return applyForward(m,null);
  }
  public float[] applyForward(float[] m, float[][] u) {
    int nz = _nz;
    int nt = _nt;
    if (u==null)
      u = new float[nt][nz];
    // Set the initial values at t=0.
    copy(m,u[0]);
    return applyForward(u);
  }
  public float[] applyForward(float[][] u) {
    int nz = _nz;
    int nt = _nt;
    float[] a = _a;
    float[] um,ui,up;
    for (int it=0; it<nt; ++it) {
      um = (it>0)?u[it-1]:new float[nz];
      ui = u[it];
      up = (it<nt-1)?u[it+1]:new float[nz];
      for (int iz=1; iz<nz-1; ++iz)
      //  up[iz] = a[iz]*(ui[iz+1]+ui[iz-1])+2.0f*(1.0f-a[iz])*ui[iz]-um[iz];
      //up[0   ] = a[0   ]*ui[1   ]+2.0f*(1.0f-a[0   ])*ui[0   ]-um[0   ];
      //up[nz-1] = a[nz-1]*ui[nz-2]+2.0f*(1.0f-a[nz-1])*ui[nz-1]-um[nz-1];
        up[iz] += a[iz]*(ui[iz+1]+ui[iz-1])+2.0f*(1.0f-a[iz])*ui[iz]-um[iz];
      up[0   ] += a[0   ]*ui[1   ]+2.0f*(1.0f-a[0   ])*ui[0   ]-um[0   ];
      up[nz-1] += a[nz-1]*ui[nz-2]+2.0f*(1.0f-a[nz-1])*ui[nz-1]-um[nz-1];
    }
    //System.out.println(sum(u));
    //SimplePlot.asPixels(u);
    float[] d = new float[nt];
    for (int it=0; it<nt; ++it)
      d[it] = u[it][0];
    return d;
  }

  /**
   * Applies the forward operator (old version).
   */
  public float[] applyForwardX(float[] m) {
    int nz = _nz;
    int nt = _nt;
    float[] a = _a;
    float[] um,ui,up;
    float[] ut = new float[nz];
    float[][] u = new float[nt][nz];
    // Set the initial values.
    copy(m,u[0]);
    for (int it=0; it<nt; ++it) {
      um = (it>0)?u[it-1]:new float[nz];
      ui = u[it];
      up = (it<nt-1)?u[it+1]:new float[nz];
      applyForwardT(a,ui,ut);
      sub(ut,um,up);
    }
    //System.out.println(sum(u));
    //SimplePlot.asPixels(u);
    float[] d = new float[nt];
    for (int it=0; it<nt; ++it)
      d[it] = u[it][0];
    return d;
  }

  /**
   * Applies the forward operator (translated from Ji's (2009) F90 code).
   */
  public float[] applyForwardJ(float[] m) {
    int nz = _nz;
    int nt = _nt;
    float[] a = _a;
    float[] d = new float[nt]; // output data at iz=0
    float[] um = new float[nz];
    float[] ui = new float[nz];
    float[] up = new float[nz];
    float[] ut = new float[nz];
    copy(m,ui); // it=0
    d[0] = m[0];
    for (int it=1; it<nt; ++it) {
      for (int iz=1; iz<nz-1; ++iz)
        up[iz] = a[iz]*(ui[iz+1]-2.0f*ui[iz]+ui[iz-1])-um[iz]+2.0f*ui[iz];
      up[0   ] = a[0   ]*(ui[1   ]-2.0f*ui[0   ])-um[0   ]+2.0f*ui[0   ];
      up[nz-1] = a[nz-1]*(ui[nz-2]-2.0f*ui[nz-1])-um[nz-1]+2.0f*ui[nz-1];
      d[it] = up[0];
      copy(ui,um);
      copy(up,ui);
    }
    return d;
  }

  /**
   * Applies the time-reversed forward operator, i.e.,
   * backward reverse-time propagation.
   */ 
  public float[] applyBackward(float[] d) {
    return applyBackward(d,null);
  }
  public float[] applyBackward(float[] d, float[][] u) {
    int nz = _nz;
    int nt = _nt;
    float[] a = _a;
    float[] um,ui,up;
    float[] ut = new float[nz];
    if (u==null)
      u = new float[nt][nz];
    // Set boundary values.
    for (int it=0; it<nt; ++it)
      u[it][0] = d[it];
    for (int it=nt-1; it>0; --it) {
      um = u[it-1];
      ui = u[it  ];
      up = (it<nt-1)?u[it+1]:new float[nz];
      applyForwardT(a,ui,ut);

      // The code below was originally sub(ut,up,um), but this fails
      // for the same reason Ji had to ignore iz=0 (see comment below
      // in applyBackwardJ). The two lines below implement Ji's
      // equation 10 if you replace element (2,2) in the second 4x4
      // matrix and element (3,3) in the third 4x4 matrix with
      // identity matrices.
      sub(ut,up,ut);
      add(ut,um,um);
    }
    //System.out.println(sum(u));
    //SimplePlot.asPixels(u);
    return u[0];
  }

  /**
   * Applies the time-reversed forward operator
   * (translated from Ji's (2009) F90 code).
   */ 
  public float[] applyBackwardJ(float[] d) {
    int nz = _nz;
    int nt = _nt;
    float[] a = _a;
    float[] um = new float[nz];
    float[] ui = new float[nz];
    float[] up = new float[nz];
    up[0] = d[nt-1];
    ui[0] = d[nt-2];
    for (int it=nt-3; it>=0; --it) {
      um[0] = d[it];
      for (int iz=1; iz<nz-1; ++iz)
        um[iz] = a[iz]*(ui[iz+1]-2.0f*ui[iz]+ui[iz-1])-up[iz]+2.0f*ui[iz];

      // Ji chooses not to compute up[0] here, because the code breaks
      // if he does. This does not agree with his equation 10.
      um[nz-1] = a[nz-1]*(ui[nz-2]-2.0f*ui[nz-1])-up[nz-1]+2.0f*ui[nz-1];

      copy(ui,up);
      copy(um,ui);
      zero(um);
    }
    return ui;
  }

  /**
   * Applies the adjoint of the forward operator.
   */ 
  public float[] applyAdjoint(float[] d) {
    return applyAdjoint(d,null);
  }
  public float[] applyAdjoint(float[] d, float[][] u) {
    int nz = _nz;
    int nt = _nt;
    float[] a = _a;
    float[] um,ui,up;
    float[] ut = new float[nz];
    if (u==null)
      u = new float[nt][nz];
    // Set boundary values.
    for (int it=0; it<nt; ++it)
      u[it][0] = d[it];
    for (int it=nt-1; it>=0; --it) {
      up = (it<nt-1)?u[it+1]:new float[nz];
      ui = u[it];
      um = (it>0)?u[it-1]:new float[nz];

      um[0] -= up[0];
      ui[0] += 2.0f*(1.0f-a[0])*up[0];
      ui[1] += a[0]*up[0];
      up[0] = 0.0f;
      um[nz-1] -= up[nz-1];
      ui[nz-2] += a[nz-1]*up[nz-1];
      ui[nz-1] += 2.0f*(1.0f-a[nz-1])*up[nz-1];
      up[nz-1] = 0.0f;
      for (int iz=1; iz<nz-1; ++iz) {
        um[iz  ] -= up[iz];
        ui[iz-1] += a[iz]*up[iz];
        ui[iz  ] += 2.0f*(1.0f-a[iz])*up[iz];
        ui[iz+1] += a[iz]*up[iz];
      }
      /*
      sub(um,up,um);
      ui[0   ] += 2.0f*(1.0f-a[0   ])*up[0   ]+a[1   ]*up[1   ];
      ui[nz-1] += 2.0f*(1.0f-a[nz-1])*up[nz-1]+a[nz-2]*up[nz-2];
      for (int iz=1; iz<nz-1; ++iz)
        ui[iz] += 2.0f*(1.0f-a[iz])*up[iz]+a[iz+1]*up[iz+1]+a[iz-1]*up[iz-1];
      */

    }
    //System.out.println(sum(u));
    //SimplePlot.asPixels(u);
    return u[0];
  }

  /**
   * Applies the adjoint operator (old version).
   */
  public float[] applyAdjointX(float[] d) {
    int nz = _nz;
    int nt = _nt;
    float[] a = _a;
    float[] um,ui,up;
    float[] ut = new float[nz];
    float[][] u = new float[nt][nz];
    // Set boundary values.
    for (int it=0; it<nt; ++it)
      u[it][0] = d[it];
    for (int it=nt-1; it>=0; --it) {
      up = (it<nt-1)?u[it+1]:new float[nz];
      ui = u[it];
      um = (it>0)?u[it-1]:new float[nz];
      sub(um,up,um);
      applyAndAddAdjointT(a,up,ui);
    }
    //System.out.println(sum(u));
    //SimplePlot.asPixels(u);
    return u[0];
  }

  /**
   * Applies the adjoint operator (translated from Ji's (2009) F90 code).
   */
  public float[] applyAdjointJ(float[] d) { // Ji
    int nz = _nz;
    int nt = _nt;
    float[] a = _a;
    float[] m = new float[nz]; // output model at it=0
    float[] um = new float[nz];
    float[] ui = new float[nz];
    float[] up = new float[nz];
    up[0] = d[nt-1];
    ui[0] = d[nt-2];
    um[0] = d[nt-3];
    for (int it=nt-2; it>=0; --it) {
      for (int iz=0; iz<nz; ++iz)
        um[iz] -= up[iz];
      for (int iz=1; iz<nz-1; ++iz)
        ui[iz] += a[iz-1]*up[iz-1]-2.0f*a[iz]*up[iz]+a[iz+1]*up[iz+1]
          +2.0f*up[iz];
      ui[0   ] += a[1   ]*up[1   ]-2.0f*a[0   ]*up[0   ]+2.0f*up[0   ];
      ui[nz-1] += a[nz-2]*up[nz-2]-2.0f*a[nz-1]*up[nz-1]+2.0f*up[nz-1];
      copy(ui,up);
      copy(um,ui);
      if (it>1) {
        um[0] = d[it-2];
        for (int iz=1; iz<nz; ++iz)
          um[iz] = 0.0f;
      }
    }
    return up;
  }

  //////////////////////////////////////////////////////////////////////////
  // private

  private int _nz;
  private int _nt;
  private double _dz;
  private double _dt;
  private float[] _a;

  // Applies the T operator described by Ji (2009).
  private static float[] applyForwardT(float[] a, float[] x) {
    float[] y = new float[x.length];
    applyForwardT(a,x,y);
    return y;
  }
  private static void applyForwardT(float[] a, float[] x, float[] y) {
    int n = x.length;
    y[0] = 2.0f*(1.0f-a[0])*x[0]+a[0]*x[1];
    for (int i=1; i<n-1; ++i) {
      y[i] = a[i]*(x[i-1]+x[i+1])+2.0f*(1.0f-a[i])*x[i];
    }
    y[n-1] = 2.0f*(1.0f-a[n-1])*x[n-1]+a[n-1]*x[n-2];
  }

  // Applies the adjoint T described by Ji (2009).
  private static float[] applyAdjointT(float[] a, float[] x) {
    float[] y = new float[x.length];
    applyAdjointT(a,x,y);
    return y;
  }
  private static void applyAdjointT(float[] a, float[] x, float[] y) {
    int n = x.length;
    y[0] = 2.0f*(1.0f-a[0])*x[0]+a[1]*x[1];
    for (int i=1; i<n-1; ++i)
      y[i] = a[i-1]*x[i-1]+2.0f*(1.0f-a[i])*x[i]+a[i+1]*x[i+1];
    y[n-1] = 2.0f*(1.0f-a[n-1])*x[n-1]+a[n-2]*x[n-2];
  }
  private static void applyAndAddAdjointT(float[] a, float[] x, float[] y) {
    float[] z = applyAdjointT(a,x);
    add(z,y,y);
  }

  //////////////////////////////////////////////////////////////////////////
  // testing
  
  public static void main(String[] args) {
    //propagationDemo();
    adjointTest();
    //adjointTestT();
  }

  private static void propagationDemo() {
    // Construct wave operator.
    Sampling sz = new Sampling(501,0.0120,0.0);
    Sampling st = new Sampling(1001,0.0015,0.0);
    int nz = sz.getCount();
    int nt = st.getCount();
    double dz = sz.getDelta();
    double dt = st.getDelta();
    float[] s = fillfloat(0.25f,nz); // slowness
    WaveOperator wave = new WaveOperator(sz,st,s);

    // Make model.
    float[] m = new float[nz]; m[nz/2] = 1.0f;
    new RecursiveGaussianFilter(0.200/0.012).apply2(m,m);
    mul(1.0f/max(abs(m)),m,m);

    // Allocate Wavefields.
    float[][] uf = new float[nt][nz]; // forward
    float[][] ub = new float[nt][nz]; // backward
    float[][] ua = new float[nt][nz]; // adjoint

    // Apply forward operator.
    float[] d = wave.applyForward(m,uf);

    // Apply backward operator.
    float[] mb = wave.applyBackward(d,ub);

    // Apply adjoint operator.
    float[] ma = wave.applyAdjoint(d,ua);

    // Plot stuff.
    SimplePlot.asPoints(m);
    SimplePlot.asPoints(mb);
    SimplePlot.asPoints(ma);
    //SimplePlot.asPoints(d);
    SimplePlot.asPixels(uf);
    SimplePlot.asPixels(ub);
    SimplePlot.asPixels(ua);
  }

  private static void adjointTest() {
    Sampling sz = new Sampling(1001,0.0120,0.0);
    Sampling st = new Sampling(1001,0.0015,0.0);
    int nz = sz.getCount();
    int nt = st.getCount();
    double dz = sz.getDelta();
    double dt = st.getDelta();
    float[] s = fillfloat(0.25f,nz); // slowness
    WaveOperator wave = new WaveOperator(sz,st,s);
    //Random random = new Random(0123);
    Random random = new Random();
    float[] m = sub(randfloat(random,nz),0.5f);
    float[] d = sub(randfloat(random,nt),0.5f);
    float sum1 = dot(wave.applyForward(m),d);
    float sum2 = dot(wave.applyAdjoint(d),m);
    //float sum1 = dot(wave.applyForwardX(m),d); // old
    //float sum2 = dot(wave.applyAdjointX(d),m); // version
    //float sum1 = dot(wave.applyForwardJ(m),d); // Ji's
    //float sum2 = dot(wave.applyAdjointJ(d),m); // code
    System.out.println("adjoint test:");
    System.out.println(sum1);
    System.out.println(sum2);
  }

  private static void adjointTestT() {
    int n = 10;
    Random random = new Random(01234);
    float[] x = randfloat(random,n);
    float[] y = randfloat(random,n);
    float[] s = randfloat(random,n);
    System.out.println("adjoint test:");
    System.out.println(dot(applyForwardT(s,x),y));
    System.out.println(dot(applyAdjointT(s,y),x));
  }

  private static float dot(float[] x, float[] y) {
    int n = x.length;
    float sum = 0.0f;
    for (int i=0; i<n; ++i)
      sum += x[i]*y[i];
    return sum;
  }

}
