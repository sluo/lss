package lss.dev;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class WaveOperator {

  public WaveOperator(Sampling sz, Sampling st, float[] s) {
    _nz = sz.getCount();
    _nt = st.getCount();
    _dz = sz.getDelta();
    _dt = st.getDelta();
    _a = div((float)(_dt*_dt/(_dz*_dz)),mul(s,s));
  }

  public float[] applyForward(float[] m) {
    int nz = _nz;
    int nt = _nt;
    float[] a = _a;
    float[][] u = new float[nt][nz];
    float[] ut = new float[nz];
    float[] um,ui,up;
    // Set boundary values at it=0
    copy(m,u[0]);
    // End case it=1
    applyForwardT(a,u[0],u[1]);
    // All other cases 2<=it<nt
    for (int it=1; it<nt; ++it) {
      // Set references
      um = u[it-1];
      ui = u[it  ];
      up = (it<nt-1)?u[it+1]:new float[nz];
      // Step
      applyForwardT(a,ui,ut);
      sub(ut,um,up);
    }
    //SimplePlot.asPixels(u);
    // Return data at iz=0
    float[] d = new float[nt];
    for (int it=0; it<nt; ++it)
      d[it] = u[it][0];
    return d;
  }

  public float[] applyBackward(float[] d) {
    int nz = _nz;
    int nt = _nt;
    float[] a = _a;
    float[][] u = new float[nt][nz];
    float[] ut = new float[nz];
    float[] um,ui,up;
    for (int it=nt-1; it>0; --it) {
      // Set references
      um = u[it-1];
      ui = u[it  ];
      up = (it<nt-1)?u[it+1]:new float[nz];
      // Set boundary values
      um[0] = d[it-1];
      ui[0] = d[it  ];
      up[0] = (it<nt-1)?d[it+1]:0.0f;
      // Step
      applyForwardT(a,ui,ut);
      sub(ut,up,um);
    }
    //SimplePlot.asPixels(u);
    return u[0];
  }

  public float[] applyAdjoint(float[] d) {
    int nz = _nz;
    int nt = _nt;
    float[] a = _a;
    float[][] u = new float[nt][nz];
    float[] ut = new float[nz];
    float[] um,ui,up;
    for (int it=nt-1; it>0; --it) {
      // Set references
      um = u[it-1];
      ui = u[it  ];
      up = (it<nt-1)?u[it+1]:new float[nz];
      // Set boundary values
      um[0] = d[it-1];
      ui[0] = d[it  ];
      up[0] = (it<nt-1)?d[it+1]:0.0f;
      // Step
      sub(um,up,um);
      applyAndAddAjointT(a,up,ui);
    }
    //SimplePlot.asPixels(u);
    return u[0];
  }

  //////////////////////////////////////////////////////////////////////////
  // private

  private int _nz;
  private int _nt;
  private double _dz;
  private double _dt;
  private float[] _a;

  private static float[] applyForwardT(float[] a, float[] x) {
    float[] y = new float[x.length];
    applyForwardT(a,x,y);
    return y;
  }
  private static void applyForwardT(float[] a, float[] x, float[] y) {
    int n = x.length;
    float ai;
    ai = a[0];
    y[0] = 2.0f*(1.0f-ai)*x[0]+ai*x[1];
    for (int i=1; i<n-1; ++i) {
      ai = a[i];
      y[i] = ai*(x[i-1]+x[i+1])+2.0f*(1.0f-ai)*x[i];
    }
    ai = a[n-1];
    y[n-1] = 2.0f*(1.0f-ai)*x[n-1]+ai*x[n-2];
  }

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
  private static void applyAndAddAjointT(float[] a, float[] x, float[] y) {
    float[] z = applyAdjointT(a,x);
    add(z,y,y);
  }

  //////////////////////////////////////////////////////////////////////////

  public static void propagateForward(float[] a, float[][] u) {
    int nx = u[0].length;
    int nt = u.length;
    float[] um,ui,up;
    float[] ut = new float[nx];
    // Nothing required for it=0.
    applyForwardT(a,u[0],u[1]); // it=1
    for (int it=1; it<nt-1; ++it) { // 2<=it<nt
      um = u[it-1];
      ui = u[it  ];
      up = u[it+1];
      applyForwardT(a,ui,ut);
      sub(ut,um,up);
    }
  }

  // Translated from Jun Ji's F90 code.
  public static void propagateForwardJ(float[] a, float[][] u) {
    int nz = u[0].length;
    int nt = u.length;
    float[] um = new float[nz];
    float[] ui = new float[nz];
    float[] up = new float[nz];
    for (int iz=0; iz<nz; ++iz)
      ui[iz] = u[0][iz];
    for (int it=1; it<nt; ++it) {
      for (int iz=1; iz<nz-1; ++iz) {
        up[iz] = a[iz]*(ui[iz+1]-2.0f*ui[iz]+ui[iz-1])-um[iz]+2.0f*ui[iz];
      }
      up[0   ] = a[0   ]*(ui[1   ]-2.0f*ui[0   ])-um[0   ]+2.0f*ui[0   ];
      up[nz-1] = a[nz-1]*(ui[nz-2]-2.0f*ui[nz-1])-um[nz-1]+2.0f*ui[nz-1];
      copy(ui,um); // rotate
      copy(up,ui); // arrays
      copy(up,u[it]); // copy to output array
    }
  }

  public static void propagateBackwardJ(float[] a, float[][] u) {
    int nz = u[0].length;
    int nt = u.length;

//    float[] um = new float[nz];
//    float[] ui = new float[nz];
//    float[] up = new float[nz];
//    float[][] p = copy(u);
//    //um[0] = p[nt-1][0];
//    //ui[0] = p[nt-2][0];
//    //for (int it=nt-3; it>=0; --it) {
//    um[0] = p[nt-1][0];
//    for (int it=nt-2; it>=0; --it) {
//      up[0] = p[it][0];
//      for (int iz=1; iz<nz-1; ++iz) {
//        up[iz] = a[iz]*(ui[iz+1]-2.0f*ui[iz]+ui[iz-1])-um[iz]+2.0f*ui[iz];
//      }
//      up[nz-1] = a[nz-1]*(ui[nz-2]-2.0f*ui[nz-1])-um[nz-1]+2.0f*ui[nz-1];
//      copy(ui,um); // rotate
//      copy(up,ui); // arrays
//      copy(up,u[it]); // copy to output array
//    }
//    System.out.println(sum(u));

    float[][] p = copy(u);
    float[] ut = new float[nz];
    float[] up = new float[nz];
    float[] ui = u[nt-1];
    float[] um = u[nt-2];
    for (int it=nt-1; it>0; --it) {
      um = u[it-1];
      applyForwardT(a,ui,ut);
      sub(ut,up,um);
      um[0] = p[it-1][0]; // set boundaries
      up = ui;
      ui = um;
    }
    System.out.println(sum(u));

  }

  public static void propagateAdjointJ(float[] a, float[][] u) {
    int nz = u[0].length;
    int nt = u.length;

//    float[] up = new float[nz];
//    for (int it=nt-1; it>0; --it) {
//      float[] ui = u[it  ];
//      float[] um = u[it-1];
//      mul(-1.0f,up,um);
//      for (int iz=1; iz<nz-1; ++iz)
//        ui[iz] += a[iz-1]*up[iz-1]-2.0f*(1.0f-a[iz])*up[iz]+a[iz+1]*up[iz+1];
//      ui[nz-1] += a[nz-2]*up[nz-2]-2.0f*(1.0f-a[nz-1])*up[nz-1];
//      up = ui;
//    }
//    System.out.println(sum(u));

    float[][] p = copy(u);
    float[] up = copy(p[nt-1]);
    float[] ui = copy(p[nt-2]);
    float[] um = copy(p[nt-3]);
    for (int it=nt-2; it>=0; --it) {
      sub(um,up,um);
      //for (int iz=1; iz<nz-1; ++iz)
      //  ui[iz] += a[iz-1]*up[iz-1]+2.0f*(1.0f-a[iz])*up[iz]+a[iz+1]*up[iz+1];
      //ui[0   ] += a[1   ]*up[1   ]+2.0f*(1.0f-a[0   ])*up[0   ];
      //ui[nz-1] += a[nz-2]*up[nz-2]+2.0f*(1.0f-a[nz-1])*up[nz-1];
      applyAndAddAjointT(a,up,ui);
      copy(ui,up);
      copy(um,ui);
      copy(um,u[it]);
      if (it>=2)
        copy(p[it-2],um);
    }
    System.out.println(sum(u));

//    float[] up = copy(p[nt-1]);
//    float[] ui = copy(p[nt-2]);
//    float[] um = copy(p[nt-3]);
//    float[] up = new float[nz];
//    for (int it=nt-2; it>=0; --it) {
//      um = p[it-1];
//      sub(um,up,um);
//      //for (int iz=1; iz<nz-1; ++iz)
//      //  ui[iz] += a[iz-1]*up[iz-1]+2.0f*(1.0f-a[iz])*up[iz]+a[iz+1]*up[iz+1];
//      //ui[0   ] += a[1   ]*up[1   ]+2.0f*(1.0f-a[0   ])*up[0   ];
//      //ui[nz-1] += a[nz-2]*up[nz-2]+2.0f*(1.0f-a[nz-1])*up[nz-1];
//      applyAndAddAjointT(a,up,ui);
//      copy(ui,up);
//      copy(um,ui);
//      copy(um,u[it]);
//      //if (it>=2)
//      //  copy(p[it-2],um);
//    }



  }

  public static void propagateAdjoint(float[] a, float[][] u) {
    int nx = u[0].length;
    int nt = u.length;
    float[] ut = new float[nx];
    float[] um = new float[nx];
    float[] ui = new float[nx];
    float[] up = new float[nx];
    float[][] p = copy(u);
    for (int it=nt-2; it>0; --it) {
      um = u[it-1];
      //ui = u[it  ];
      //up = u[it+1];
      mul(-1.0f,up,um);
      applyAndAddAjointT(a,up,ut);
      add(ut,ui,ui);

      copy(ui,up);
      add(um,p[it-1],ui);
    }
  }
  public static void xpropagateAdjoint(float[] a, float[][] u) {
    int nx = u[0].length;
    int nt = u.length;
    float[] um,ui,up;
    float[] ut = new float[nx];
    for (int it=nt-2; it>0; --it) {
      um = u[it-1];
      ui = u[it  ];
      up = u[it+1];
      um = mul(-1.0f,up);
      applyAndAddAjointT(a,up,ut);
      add(ut,ui,ui);
    }
  }

  public static void xpropagateBackward(float[] a, float[][] u) {
    int nx = u[0].length;
    int nt = u.length;
    float[] ut = new float[nx];
    float[] um = new float[nx];
    float[] ui = new float[nx];
    float[] up = new float[nx];
    //float[] um = u[nt-3];
    //float[] ui = u[nt-2];
    //float[] up = u[nt-1];
    float[][] p = copy(u);
//    for (int it=nt-2; it>0; --it) {
//      zero(um);
//      zero(ui);
//      zero(up);
//      zero(ut);
//      copy(u[it-1],um);
//      copy(u[it  ],ui);
//      copy(u[it+1],up);
//      applyForwardT(a,ui,ut);
//      sub(ut,up,um);
//      copy(ui,up);
//      add(um,p[it-1],ui);
//      copy(um,u[it-1]);
//      copy(ui,u[it  ]);
//      copy(up,u[it+1]);
//    }
    for (int it=nt-2; it>0; --it) {
      um = u[it-1];
      //ui = u[it  ];
      //up = u[it+1];
      applyForwardT(a,ui,ut);
      sub(ut,up,um);
      copy(ui,up);
      add(um,p[it-1],ui);
    }
  }
  public static void propagateBackward(float[] a, float[][] u) {
    int nx = u[0].length;
    int nt = u.length;
    //float[] um,ui,up;
    float[] um = new float[nx];
    float[] ui = new float[nx];
    float[] up = new float[nx];
    float[] ut = new float[nx];
    float[][] p = copy(u);
    // Nothing required for it=nt-1.
    SimplePlot.asPixels(u);
    applyForwardT(a,u[nt-1],u[nt-2]); // it=nt-2
    for (int it=nt-2; it>0; --it) {
      //um = u[it-1];
      //ui = u[it  ];
      //up = u[it+1];
      copy(u[it-1],um);
      copy(u[it  ],ui);
      copy(u[it+1],up);
      applyForwardT(a,ui,ut);
      sub(ut,up,um);
      copy(um,u[it-1]);
      copy(ui,u[it  ]);
      copy(up,u[it+1]);
    }
    add(p[0],u[0],u[0]);
  }

  //////////////////////////////////////////////////////////////////////////
  // testing
  
  public static void main(String[] args) {
    propagationTest();
    //adjointTest();
    //adjointTestT();
  }

  private static void propagationTest() {
    Sampling sx = new Sampling(501,0.0120,0.0);
    Sampling st = new Sampling(1001,0.0015,0.0);
    int nx = sx.getCount();
    int nt = st.getCount();
    double dx = sx.getDelta();
    double dt = st.getDelta();
    float[] s = fillfloat(0.25f,nx); // slowness
    float[] a = div((float)(dt*dt/(dx*dx)),mul(s,s));
    float[][] u = new float[nt][nx];
    u[0][nx/2] = 1.0f;
    //new RecursiveGaussianFilter(0.050/0.012).apply0(u[0],u[0]);
    new RecursiveGaussianFilter(0.200/0.012).apply20(u,u);
    mul(1.0f/max(abs(u)),u,u);

    //propagateForward(a,u);
    propagateForwardJ(a,u); // Ji's code
    //SimplePlot.asPixels(u);
    //SimplePlot.asPoints(u[0]);

    //WaveOperator wave = new WaveOperator(sx,st,s);
    //float[] d = wave.applyForward(u[0]);

    //SimplePlot.asPixels(u);
    //SimplePlot.asPoints(u[0]);

    float[][] v = copy(u);
    for (int it=0; it<nt; ++it)
      for (int ix=1; ix<nx; ++ix)
        v[it][ix] = 0.0f;
    //float[] d = new float[nt];
    //for (int it=0; it<nt; ++it)
    //  d[it] = u[it][0];

    //propagateBackward(a,v);
    propagateBackwardJ(a,v);
    //propagateAdjoint(a,v);
    //propagateAdjointJ(a,v);
    SimplePlot.asPixels(v);
    SimplePlot.asPoints(v[0]);

    //float[] m = wave.applyBackward(d);
    //SimplePlot.asPoints(m);
   
  }

  private static void adjointTest() {
    Sampling sz = new Sampling(501,0.0120,0.0);
    Sampling st = new Sampling(1001,0.0015,0.0);
    int nz = sz.getCount();
    int nt = st.getCount();
    double dz = sz.getDelta();
    double dt = st.getDelta();
    float[] s = fillfloat(0.25f,nz); // slowness
    WaveOperator wave = new WaveOperator(sz,st,s);
    float[] m = sub(randfloat(nz),0.5f);
    float[] d = sub(randfloat(nt),0.5f);
    System.out.println(dot(wave.applyForward(m),d));
    System.out.println(dot(wave.applyAdjoint(d),m));
  }

  private static void adjointTestT() {
    int n = 10;
    float[] x = randfloat(n);
    float[] y = randfloat(n);
    float[] s = randfloat(n);
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
