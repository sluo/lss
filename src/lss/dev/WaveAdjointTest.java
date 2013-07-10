package lss.dev;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class WaveAdjointTest {

  public WaveAdjointTest(
    Sampling sz, Sampling sx, Sampling st)
  {
    _dx = sx.getDelta();
    _dt = st.getDelta();
    _nx = sx.getCount();
    _nt = st.getCount();
  }

  public static void propagateForward(float[] a, float[][] u) {
    int nx = u[0].length;
    int nt = u.length;
    float[] um,ui,up;
    float[] ut = new float[nx];
    for (int it=1; it<nt-1; ++it) {
      um = u[it-1];
      ui = u[it  ];
      up = u[it+1];
      applyForward(a,ui,ut);
      sub(ut,um,up);
    }
  }

  public static void propagateAdjoint(float[] a, float[][] u) {
    int nx = u[0].length;
    int nt = u.length;
    float[] um,ui,up;
    float[] ut = new float[nx];
    for (int it=nt-2; it>0; --it) {
      um = u[it-1];
      ui = u[it  ];
      up = u[it+1];
      um = mul(-1.0f,up);
      applyAdjoint(a,up,ut);
      add(ut,ui,ui);
    }
  }

  public static void propagateBackward(float[] a, float[][] u) {
    int nx = u[0].length;
    int nt = u.length;
    float[] um,ui,up;
    float[] ut = new float[nx];
    for (int it=nt-2; it>0; --it) {
      um = u[it-1];
      ui = u[it  ];
      up = u[it+1];
      applyForward(a,ui,ut);
      sub(ut,up,um);
    }
  }

  private static float[] applyForward(float[] a, float[] x) {
    float[] y = new float[x.length];
    applyForward(a,x,y);
    return y;
  }
  private static void applyForward(float[] a, float[] x, float[] y) {
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

  private static float[] applyAdjoint(float[] a, float[] x) {
    float[] y = new float[x.length];
    applyAdjoint(a,x,y);
    return y;
  }
  private static void applyAdjoint(float[] a, float[] x, float[] y) {
    int n = x.length;
    float am,ai,ap;
    ai = a[0];
    ap = a[1];
    y[0] = 2.0f*(1.0f-ai)*x[0]+ap*x[1];
    for (int i=1; i<n-1; ++i) {
      am = a[i-1];
      ai = a[i  ];
      ap = a[i+1];
      y[i] = am*x[i-1]+2.0f*(1.0f-ai)*x[i]+ap*x[i+1];
    }
    am = a[n-2];
    ai = a[n-1];
    y[n-1] = am*x[n-2]+2.0f*(1.0f-ai)*x[n-1];
  }

  //////////////////////////////////////////////////////////////////////////
  // private

  private int _it; // current time index
  private int _nt;
  private int _nx;
  private double _dt;
  private double _dx;

  //////////////////////////////////////////////////////////////////////////
  // testing
  
  public static void main(String[] args) {
    //adjointTest();
    propagationTest();
  }

  private static void propagationTest() {
    Sampling sx = new Sampling(501,0.0120,0.0);
    Sampling st = new Sampling(2001,0.0015,0.0);
    int nx = sx.getCount();
    int nt = st.getCount();
    double dx = sx.getDelta();
    double dt = st.getDelta();
    float[] s = fillfloat(0.25f,nx); // slowness
    float[] a = div((float)(dt*dt/(dx*dx)),mul(s,s));
    float[][] u = new float[nt][nx];
    u[0][nx/2] = 1.0f;
    //new RecursiveGaussianFilter(0.050/0.012).apply0(u[0],u[0]);
    new RecursiveGaussianFilter(0.100/0.012).apply20(u,u);
    propagateForward(a,u);
    SimplePlot.asPixels(u);
    SimplePlot.asPoints(u[0]);

    float[][] v = copy(u);
    for (int it=0; it<nt; ++it)
      for (int ix=1; ix<nx; ++ix)
        v[it][ix] = 0.0f;
    //propagateAdjoint(a,v);
    propagateBackward(a,v);
    SimplePlot.asPixels(v);
    SimplePlot.asPoints(v[0]);
   
  }

  private static void adjointTest() {
    int n = 10;
    float[] x = randfloat(n);
    float[] y = randfloat(n);
    float[] s = randfloat(n);
    System.out.println(dot(applyForward(s,x),y));
    System.out.println(dot(applyAdjoint(s,y),x));
  }

  private static float dot(float[] x, float[] y) {
    int n = x.length;
    float sum = 0.0f;
    for (int i=0; i<n; ++i)
      sum += x[i]*y[i];
    return sum;
  }
}
