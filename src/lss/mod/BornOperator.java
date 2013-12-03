package lss.mod;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Born modeling with an optional time-shift operator.
 * Note that, to obtain Born modeled data that matches data computed
 * by finite-differencing of the acoustic wave equation, the source
 * function for Born modeling must be phase-rotated by 0.5*PI,
 * twice-differentiated, and negated.
 * @author Simon Luo, Colorado School of Mines
 * @version 2013.11.20
 */
public class BornOperator {

  /**
   * Constructs a Born modeling operator.
   * @param s slowness model.
   * @param dx spatial sampling interval.
   * @param dt time sampling interval.
   * @param nabsorb width in samples of absorbing boundary.
   */
  public BornOperator(float[][] s, double dx, double dt, int nabsorb) {
    _wave = new WaveOperator(s,dx,dt,nabsorb);
    _nx = s[0].length;
    _nz = s.length;
    _nabsorb = nabsorb;
    _nzp = _nz+2*nabsorb;
    _nxp = _nx+2*nabsorb;
    _dx = (float)dx;
    _dt = (float)dt;
    _s = s;
  }

  public void setAdjoint(boolean adjoint) {
    _wave.setAdjoint(adjoint);
  }

  public void setSlowness(float[][] s) {
    copy(s,_s);
    _wave.setSlowness(s);
  }

  public WaveOperator getWaveOperator() {
    return _wave;
  }

  ////////////////////////////////////////////////////////////////////////////
  // forward

  /**
   * Applies the forward operator after computing background wavefield.
   * @param source input source for computing background wavefield.
   * @param b input array for storing background wavefield.
   * @param rx input reflectivity image.
   * @param ts input time shifts.
   * @param receiver output receiver.
   * @param u output wavefield.
   */
  public void applyForward(
    Source source, float[][][] b, float[][] rx, float[][] ts,
    Receiver receiver, float[][][] u)
  {
    computeBackgroundWavefield(source,b);
    applyForward(b,rx,ts,receiver,u);
  }

  public void applyForward(
    Source source, float[][][] b, float[][] rx, float[][] ts,
    Receiver receiver)
  {
    applyForward(source,b,rx,ts,receiver,null);
  }

  public void applyForward(
    Source source, float[][][] b, float[][] rx, float[][] ts,
    float[][][] u)
  {
    applyForward(source,b,rx,ts,null,u);
  }

  public void applyForward(
    Source source, float[][][] b, float[][] rx,
    Receiver receiver, float[][][] u)
  {
    applyForward(source,b,rx,null,receiver,u);
  }

  public void applyForward(
    Source source, float[][][] b, float[][] rx,
    Receiver receiver)
  {
    applyForward(source,b,rx,null,receiver,null);
  }

  public void applyForward(
    Source source, float[][][] b, float[][] rx,
    float[][][] u)
  {
    applyForward(source,b,rx,null,null,u);
  }

  /**
   * Applies the forward operator using precomputed background wavefield.
   * @param b input array containing precomputed background wavefield.
   * @param rx input reflectivity image.
   * @param ts input time shifts.
   * @param receiver output receiver.
   * @param u output wavefield.
   */
  public void applyForward(
    float[][][] b, float[][] rx, float[][] ts,
    Receiver receiver, float[][][] u)
  {
    Check.argument(b[0][0].length-rx[0].length==2*_nabsorb,"consistent nx");
    Check.argument(b[0].length-rx.length==2*_nabsorb,"consistent nz");
    if (u!=null) {
      if (receiver!=null)
        Check.argument(u.length==receiver.getNt(),"consistent nt");
      Check.argument(b[0][0].length==u[0][0].length,"consistent nx");
      Check.argument(b[0].length==u[0].length,"consistent nz");
    }
    _wave.applyForward(new Source.WavefieldSource(b,rx),receiver,u);
    applyForwardShifts(ts,receiver,receiver);
  }

  public void applyForward(
    float[][][] b, float[][] rx, float[][] ts,
    Receiver receiver)
  {
    applyForward(b,rx,ts,receiver,null);
  }

  public void applyForward(
    float[][][] b, float[][] rx, float[][] ts,
    float[][][] u)
  {
    applyForward(b,rx,ts,null,u);
  }

  public void applyForward(
    float[][][] b, float[][] rx,
    Receiver receiver, float[][][] u)
  {
    applyForward(b,rx,null,receiver,u);
  }

  public void applyForward(
    float[][][] b, float[][] rx,
    Receiver receiver)
  {
    applyForward(b,rx,null,receiver,null);
  }

  public void applyForward(
    float[][][] b, float[][] rx,
    float[][][] u)
  {
    applyForward(b,rx,null,null,u);
  }

  ////////////////////////////////////////////////////////////////////////////
  // adjoint

  /**
   * Applies the adjoint operator after computing the background wavefield.
   * @param source input source for computing background wavefield.
   * @param b input array for storing background wavefield.
   * @param a input array for storing adjoint wavefield.
   * @param receiver input receiver containing data to be migrated.
   * @param ts input time shifts.
   * @param ry output reflectivity image.
   */
  public void applyAdjoint(
    Source source, float[][][] b, float[][][] a,
    Receiver receiver, float[][] ts,
    float[][] ry)
  {
    computeBackgroundWavefield(source,b);
    applyAdjoint(b,a,receiver,ts,ry);
  }

  public void applyAdjoint(
    Source source, float[][][] b,
    float[][][] a, Receiver receiver, float[][] ry)
  {
    applyAdjoint(source,b,a,receiver,null,ry);
  }

  /**
   * Applies the adjoint operator after computing the background wavefield.
   * @param b input array containing precomputed background wavefield.
   * @param a input array for storing adjoint wavefield.
   * @param receiver input receiver containing data to be migrated.
   * @param ts input time shifts.
   * @param ry output reflectivity image.
   */
  public void applyAdjoint(
    float[][][] b, float[][][] a, Receiver receiver, float[][] ts,
    float[][] ry)
  {
    Check.argument(b[0][0].length-ry[0].length==2*_nabsorb,"consistent nx");
    Check.argument(b[0].length-ry.length==2*_nabsorb,"consistent nz");
    Receiver rc = receiver.clone();
    applyAdjointShifts(ts,receiver,rc);
    _wave.applyAdjoint(new Source.ReceiverSource(rc),a);
    WaveOperator.collapse(b,a,_nabsorb,ry);
  }

  public void applyAdjoint(
    float[][][] b, float[][][] a, Receiver receiver,
    float[][] ry)
  {
    applyAdjoint(b,a,receiver,null,ry);
  }

  ////////////////////////////////////////////////////////////////////////////
  // hessian

  /**
   * Applies the Hessian operator after computing the background wavefield.
   * @param source input source for computing background wavefield.
   * @param receiver input receiver containing data to be migrated.
   * @param b input array for storing background wavefield.
   * @param a input array for storing adjoint wavefield.
   * @param rx input reflectivity image.
   * @param ts input time shifts.
   * @param ry output reflectivity image.
   */
  public void applyHessian(
    Source source, Receiver receiver, float[][][] b,
    float[][][] a, float[][] rx, float[][] ts,
    float[][] ry)
  {
    computeBackgroundWavefield(source,b);
    applyHessian(receiver,b,a,rx,ts,ry);
  }

  public void applyHessian(
    Source source, Receiver receiver, float[][][] b,
    float[][][] a, float[][] rx,
    float[][] ry)
  {
    applyHessian(source,receiver,b,a,rx,null,ry);
  }

  /**
   * Applies the Hessian operator using precomputed background wavefield.
   * @param receiver input receiver containing data to be migrated.
   * @param b input array for storing background wavefield.
   * @param a input array for storing adjoint wavefield.
   * @param rx input reflectivity image.
   * @param ts time shifts.
   * @param ry output reflectivity image.
   */
  public void applyHessian(
    Receiver receiver, float[][][] b, float[][][] a,
    float[][] rx, float[][] ts,
    float[][] ry)
  {
    applyForward(b,rx,ts,receiver);
    applyAdjoint(b,a,receiver,ts,ry);
  }

  public void applyHessian(
    Receiver receiver, float[][][] b, float[][][] a, float[][] rx,
    float[][] ry)
  {
    applyHessian(receiver,b,a,rx,null,ry);
  }

  ////////////////////////////////////////////////////////////////////////////
  // illumination map

  public void applyForIllumination(
    Source source, float[][][] b, float[][] m)
  {
    int nx = b[0][0].length;
    int nz = b[0].length;
    int nt = b.length;
    Check.argument(nx-m[0].length==2*_nabsorb,"consistent nx");
    Check.argument(nz-m.length==2*_nabsorb,"consistent nz");
    _wave.applyForward(source,b);
    WaveOperator.collapse(b,b,_nabsorb,m);
    mul(1.0f/nx/nz/nt,m,m);
  }

  ////////////////////////////////////////////////////////////////////////////
  // shifts

  private static void applyForwardShifts(
    float[][] ts, Receiver rcx, Receiver rcy)
  {
    if (ts==null) return;
    //System.out.println("applying forward shifts...");
    float[][] d = rcx.getData();
    float[][] e = applyShifts(ts,d,false);
    rcy.setData(e);
  }

  private static void applyAdjointShifts(
    float[][] ts, Receiver rcx, Receiver rcy)
  {
    if (ts==null) return;
    //System.out.println("applying adjoint shifts...");
    float[][] d = rcx.getData();
    float[][] e = applyShifts(ts,d,true);
    rcy.setData(e);
  }

  private static float[][] applyShifts(
    float[][] w, float[][] d, final boolean adjoint)
  {
    final int nt = w[0].length;
    final int nr = w.length;
    final float[] rf = rampfloat(0.0f,1.0f,nt);
    final float[][] wf = w;
    final float[][] df = d;
    final float[][] ef = new float[nr][nt]; // shifted
    final SincInterp si = new SincInterp();
    Parallel.loop(nr,new Parallel.LoopInt() {
    public void compute(int ir) {
      float[] p = add(rf,wf[ir]);
      if (!adjoint) {
        si.interpolate(nt,1.0,0.0,df[ir],nt,p,ef[ir]);
      } else {
        si.accumulate(nt,p,df[ir],nt,1.0,0.0,ef[ir]);
      }
    }});
    return ef;
  }

  //////////////////////////////////////////////////////////////////////////
  // private

  private WaveOperator _wave;
  private int _nabsorb;
  private float _dx,_dt;
  private int _nx,_nz,_nxp,_nzp;
  private float[][] _s; // background slowness

  private void computeBackgroundWavefield(Source source, float[][][] b) {
    _wave.applyForward(source,b);
    //scaleLaplacian(b); // equivalent to negative 2nd time derivative
  }

  // 20th order stencil coefficients from Farhad.
  private static final int FD_ORDER = 20;
  private static final float C00 = -0.32148051f*10.0f*2.0f;
  private static final float C01 =  0.19265816f*10.0f;
  private static final float C02 = -0.43052632f*1.0f; 
  private static final float C03 =  0.15871000f*1.0f;
  private static final float C04 = -0.68711400f*0.1f;
  private static final float C05 =  0.31406935f*0.1f;
  private static final float C06 = -0.14454222f*0.1f;
  private static final float C07 =  0.65305182f*0.01f;
  private static final float C08 = -0.28531535f*0.01f;
  private static final float C09 =  0.11937032f*0.01f;
  private static final float C10 = -0.47508613f*0.001f;
  private static final float FF = 400.0f; // fudge factor
  private void scaleLaplacian(final float[][][] u) {
    final int nx = u[0][0].length;
    final int nz = u[0].length;
    final int nt = u.length;
    final int b = _nabsorb-FD_ORDER/2;
    final int ixa = FD_ORDER/2;
    final int ixb = ixa+b;
    final int ixc = ixb+nx-2*_nabsorb;
    final int ixd = ixc+b;
    final int iza = FD_ORDER/2;
    final int izb = iza+b;
    final int izc = izb+nz-2*_nabsorb;
    final int izd = izc+b;
    final float[][] s = extendModel(_s,_nabsorb);
    float[][] ua = new float[nz][nx];
    float[][] ub;
    for (int it=0; it<nt; ++it) {
      final float[][] uif = u[it];
      final float[][] uaf = ua;
      Parallel.loop(iza,izd,new Parallel.LoopInt() {
      public void compute(int iz) {
        float[] uii = uif[iz];
        float[] uim01 = uif[iz-1 ], uip01 = uif[iz+1 ];
        float[] uim02 = uif[iz-2 ], uip02 = uif[iz+2 ];
        float[] uim03 = uif[iz-3 ], uip03 = uif[iz+3 ];
        float[] uim04 = uif[iz-4 ], uip04 = uif[iz+4 ];
        float[] uim05 = uif[iz-5 ], uip05 = uif[iz+5 ];
        float[] uim06 = uif[iz-6 ], uip06 = uif[iz+6 ];
        float[] uim07 = uif[iz-7 ], uip07 = uif[iz+7 ];
        float[] uim08 = uif[iz-8 ], uip08 = uif[iz+8 ];
        float[] uim09 = uif[iz-9 ], uip09 = uif[iz+9 ];
        float[] uim10 = uif[iz-10], uip10 = uif[iz+10];
        for (int ix=ixa; ix<ixd; ++ix) {
          float f = FF*_dt*_dt*_dx;
          float d = s[iz][ix]*_dx;
          float r = -f/(d*d);
          //float q = 0.5f*r;
          //float a = 0.5461f; // (Jo et. al., 1996)
          //float b = 1.0f-a;
          float a = 1.0f;
          uaf[iz][ix] = (
            a*r*(
            C00*(uii[ix])+
            C01*(uim01[ix]+uii[ix-1 ]+uii[ix+1 ]+uip01[ix])+
            C02*(uim02[ix]+uii[ix-2 ]+uii[ix+2 ]+uip02[ix])+
            C03*(uim03[ix]+uii[ix-3 ]+uii[ix+3 ]+uip03[ix])+
            C04*(uim04[ix]+uii[ix-4 ]+uii[ix+4 ]+uip04[ix])+
            C05*(uim05[ix]+uii[ix-5 ]+uii[ix+5 ]+uip05[ix])+
            C06*(uim06[ix]+uii[ix-6 ]+uii[ix+6 ]+uip06[ix])+
            C07*(uim07[ix]+uii[ix-7 ]+uii[ix+7 ]+uip07[ix])+
            C08*(uim08[ix]+uii[ix-8 ]+uii[ix+8 ]+uip08[ix])+
            C09*(uim09[ix]+uii[ix-9 ]+uii[ix+9 ]+uip09[ix])+
            C10*(uim10[ix]+uii[ix-10]+uii[ix+10]+uip10[ix]))
            //+
            //b*q*(
            //C00*(uii[ix])+
            //C01*(uim01[ix-1 ]+uim01[ix+1 ]+uip01[ix-1 ]+uip01[ix+1 ])+
            //C02*(uim02[ix-2 ]+uim02[ix+2 ]+uip02[ix-2 ]+uip02[ix+2 ])+
            //C03*(uim03[ix-3 ]+uim03[ix+3 ]+uip03[ix-3 ]+uip03[ix+3 ])+
            //C04*(uim04[ix-4 ]+uim04[ix+4 ]+uip04[ix-4 ]+uip04[ix+4 ])+
            //C05*(uim05[ix-5 ]+uim05[ix+5 ]+uip05[ix-5 ]+uip05[ix+5 ])+
            //C06*(uim06[ix-6 ]+uim06[ix+6 ]+uip06[ix-6 ]+uip06[ix+6 ])+
            //C07*(uim07[ix-7 ]+uim07[ix+7 ]+uip07[ix-7 ]+uip07[ix+7 ])+
            //C08*(uim08[ix-8 ]+uim08[ix+8 ]+uip08[ix-8 ]+uip08[ix+8 ])+
            //C09*(uim09[ix-9 ]+uim09[ix+9 ]+uip09[ix-9 ]+uip09[ix+9 ])+
            //C10*(uim10[ix-10]+uim10[ix+10]+uip10[ix-10]+uip10[ix+10]))
          );
        }
      }});
      ub = u[it]; // rotate
      u[it] = ua; // arrays
      ua = ub;
    }
  }
  private static float[][] extendModel(float[][] c, int b) {
    int nx = c[0].length;
    int nz = c.length;
    float[][] v = new float[nz+2*b][nx+2*b];
    copy(nx,nz,0,0,c,b,b,v);
    for (int iz=b; iz<nz+b; ++iz) {
      for (int ix=0, jx=nx+b; ix<b; ++ix, ++jx) {
        v[iz][ix] = v[iz][b];
        v[iz][jx] = v[iz][nx+b-1];
      }
    }
    for (int iz=0, jz=nz+b; iz<b; ++iz, ++jz) {
      copy(v[b],v[iz]);
      copy(v[nz+b-1],v[jz]);
    }
    return v;
  }

}
