package lss.vel;

import java.util.Random;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

// testing
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;

// TODO
// Avoid copy in adjointAbsorb?
// Speed up adjointStep.

public class WaveOperator {

  // count number of same digits
  public static int compareDigits(float xa, float xb) {
    int significantDigits = 20;
    boolean matches = false;
    while (!matches && significantDigits>0) {
      Almost almost = new Almost(significantDigits);
      matches = almost.equal(xa, xb);
      if (!matches) {
        --significantDigits;
      }
    }
    return significantDigits;
  }

  public static float[][][][] randfloat(
  java.util.Random random, int n1, int n2, int n3, int n4) {
    float[][][][] r = new float[n4][][][];
    for (int i4=0; i4<n4; ++i4)
      r[i4] = edu.mines.jtk.util.ArrayMath.randfloat(random,n1,n2,n3);
    return r;
  }

  public static float dot(float[][] u, float[][] a) {
    int nx = u[0].length;
    int nz = u.length;
    float sum = 0.0f;
      for (int iz=0; iz<nz; ++iz)
        for (int ix=0; ix<nx; ++ix)
          sum += u[iz][ix]*a[iz][ix];
    return sum;
  }
  public static float dot(float[][][] u, float[][][] a) {
    return dot(u,a,0);
  }
  public static float dot(float[][][] u, float[][][] a, int nabsorb) {
    int nx = u[0][0].length;
    int nz = u[0].length;
    int nt = u.length;
    float sum = 0.0f;
    for (int it=0; it<nt; ++it)
      for (int iz=nabsorb; iz<nz-nabsorb; ++iz)
        for (int ix=nabsorb; ix<nx-nabsorb; ++ix)
          sum += u[it][iz][ix]*a[it][iz][ix];
    return sum;
  }

  //////////////////////////////////////////////////////////////////////////
  // static

  public static float[][] crop(float[][] x, int nabsorb) {
    int nx = x[0].length;
    int nz = x.length;
    return copy(nx-2*nabsorb,nz-2*nabsorb,nabsorb,nabsorb,x);
  }

  public static void zeroBoundary(float[][] x, int nabsorb) {
    int nx = x[0].length;
    int nz = x.length;
    for (int iz=0, jz=nz-nabsorb; iz<nabsorb; ++iz, ++jz) {
      for (int ix=0; ix<nx; ++ix) {
        x[iz][ix] = 0.0f;
      }
    }
    for (int iz=nabsorb; iz<nz-nabsorb; ++iz) {
      for (int ix=0, jx=nx-nabsorb; ix<nabsorb; ++ix, ++jx) {
        x[iz][ix] = 0.0f;
      }
    }
  }
  public static void zeroBoundary(float[][][] x, int nabsorb) {
    int nt = x.length;
    for (int it=0; it<nt; ++it) {
      zeroBoundary(x[it],nabsorb);
    }
  }

  public static float[][] collapse(
  float[][][] u, float[][][] a, int nabsorb) {
    int nx = u[0][0].length;
    int nz = u[0].length;
    float[][] r = new float[nz-2*nabsorb][nx-2*nabsorb];
    collapse(u,a,nabsorb,r);
    return r;
  }
  public static void collapse(
  float[][][] u, float[][][] a, int nabsorb, float[][] r) {
    zero(r);
    int nx = u[0][0].length;
    int nz = u[0].length;
    int nt = u.length;
    for (int it=0; it<nt; ++it)
      for (int iz=nabsorb; iz<nz-nabsorb; ++iz)
        for (int ix=nabsorb; ix<nx-nabsorb; ++ix)
          r[iz-nabsorb][ix-nabsorb] += u[it][iz][ix]*a[it][iz][ix];
  }

  //////////////////////////////////////////////////////////////////////////

  public WaveOperator(
  float[][] s, double dx, double dt, int nabsorb) {
    int nx = s[0].length;
    int nz = s.length;
    Check.argument(nabsorb>=FD_ORDER,"nabsorb>=FD_ORDER");
    _nabsorb = nabsorb;
    _b = nabsorb-FD_ORDER/2;
    _nz = nz+2*nabsorb;
    _nx = nx+2*nabsorb;
    _dx = (float)dx;
    _dt = (float)dt;
    _ixa = FD_ORDER/2;
    _ixb = _ixa+_b;
    _ixc = _ixb+nx;
    _ixd = _ixc+_b;
    _iza = FD_ORDER/2;
    _izb = _iza+_b;
    _izc = _izb+nz;
    _izd = _izc+_b;
    _w = makeWeights();
    setSlowness(s);
  }

  public void setAdjoint(boolean adjoint) {
    _adjoint = adjoint;
  }

  public void setSlowness(float[][] s) {
    float scale = (_dt*_dt/(_dx*_dx));
    _s = extendModel(s,_nabsorb);
    _r = copy(_s);
    mul(_r,_r,_r);
    div(scale,_r,_r);
  }

  public void applyForward(Source source, float[][][] u) {
    apply(true,source,null,u);
  }

  public void applyForward(Source source, Receiver receiver) {
    apply(true,source,receiver,null);
  }

  public void applyForward(Source source, Receiver receiver, float[][][] u) {
    apply(true,source,receiver,u);
  }

  public void applyAdjoint(Source source, float[][][] u) {
    apply(false,source,null,u);
  }

  public void applyAdjoint(Source source, Receiver receiver) {
    apply(false,source,receiver,null);
  }

  public void applyAdjoint(Source source, Receiver receiver, float[][][] u) {
    apply(false,source,receiver,u);
  }

  //////////////////////////////////////////////////////////////////////////
  // private


  private int _b,_nabsorb; // absorbing boundary
  private int _nx,_nz;
  private float _dx,_dt;
  private float[][] _s; // slowness
  private float[][] _r; // dt^2/(dx^2*s^2)
  private float[][] _w; // weights for boundary
  private int _ixa,_ixb,_ixc,_ixd;
  private int _iza,_izb,_izc,_izd;

  // Adjoint flag; if true, use for adjoint code during back propagation.
  private boolean _adjoint = true;

  private void apply(
  boolean forward, Source source, Receiver receiver, float[][][] u) {
    int nt = (u==null)?receiver.getNt():u.length;
    float[][] um,ui,up;
    float[][] ut = new float[_nz][_nx];
    if (u==null) {
      um = new float[_nz][_nx];
      ui = new float[_nz][_nx];
    } else {
      zero(u); // off for adjoint test?
      um = (forward)?u[0]:u[nt-1];
      ui = (forward)?u[1]:u[nt-2];
    }
    int fit = (forward)?2:nt-3;
    int pit = (forward)?1:-1;
    source.add(um,fit-2*pit,_nabsorb);
    source.add(ui,fit-pit  ,_nabsorb);
    if (receiver!=null) {
      receiver.setData(um,fit-2*pit,_nabsorb);
      receiver.setData(ui,fit-pit  ,_nabsorb);
    }
    for (int it=fit, count=0; count<nt-2; it+=pit, ++count) {

      // Set next time step.
      if (u==null) {
        up = ut;
        zero(ut); // necessary if using += for injection
      } else {
        up = u[it];
      }

      // Step.
      if (forward) { // forward propagate
        forwardStep(um,ui,up);
        if (source!=null)
          source.add(up,it,_nabsorb);
        forwardAbsorb(um,ui,up);
      } else if (_adjoint) { // back propagate using adjoint code
        adjointAbsorb(up,ui,um);
        if (source!=null)
          source.add(up,it,_nabsorb);
        adjointStep(up,ui,um);
      } else { // back propagate using forward code
        forwardStep(um,ui,up);
        if (source!=null)
          source.add(up,it,_nabsorb);
        forwardAbsorb(um,ui,up);
      }

      // Set data.
      if (receiver!=null)
        receiver.setData(up,it,_nabsorb);

      // Rotate arrays.
      ut = um; um = ui; ui = up;
    }
  }

  /* Indexing example for FD_ABSORB/2=2, nabsorb=3, nx=3:

          |ixa|ixb        |ixc|ixd
  |---|---|---|---|---|---|---|---|---|
  | S | S | B | M | M | M | B | S | S |
  |---|---|---|---|---|---|---|---|---|
          |   |           |   |

    S - extra stencil cells for finite-difference stencil
    B - boundary cells for absorbing boundary
    M - model cells
  */

//  // 21-point Laplacian stencil from
//  // Patra, M. and M. Karttunen, 2005, Stencils
//  // with Isotropic Error for Differential Operators.
//  private static final int FD_ORDER = 4;
//  private static final float C0 = -9.0f/2.0f;
//  private static final float C1 =  16.0f/15.0f;
//  private static final float C2 =  2.0f/15.0f;
//  private static final float C3 = -1.0f/15.0f;
//  private static final float C4 = -1.0f/120.0f;
//
//  private void forwardStep(
//  final float[][] um, final float[][] ui, final float[][] up) {
//    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
//    public void compute(int iz) {
//      int izm1 = iz-1;
//      int izp1 = iz+1;
//      int izm2 = iz-2;
//      int izp2 = iz+2;
//      for (int ix=_ixa; ix<_ixd; ++ix) {
//        int ixm1 = ix-1;
//        int ixp1 = ix+1;
//        int ixm2 = ix-2;
//        int ixp2 = ix+2;
//        float r = _r[iz][ix];
//        up[iz][ix] += (2.0f+C0*r)*ui[iz][ix]-um[iz][ix]+r*(
//          C1*(ui[izm1][ix  ]+ui[izp1][ix  ]+ui[iz  ][ixm1]+ui[iz  ][ixp1])+
//          C2*(ui[izm1][ixm1]+ui[izp1][ixm1]+ui[izm1][ixp1]+ui[izp1][ixp1])+
//          C3*(ui[izm2][ix  ]+ui[izp2][ix  ]+ui[iz  ][ixm2]+ui[iz  ][ixp2])+
//          C4*(ui[izm2][ixm2]+ui[izp2][ixm2]+ui[izm2][ixp2]+ui[izp2][ixp2]));
//        //up[iz][ix] -= um[iz  ][ix  ];
//        //up[iz][ix] += ui[iz  ][ix  ]*(2.0f+C0*r);
//        //up[iz][ix] += ui[izm1][ix  ]*C1*r;
//        //up[iz][ix] += ui[izp1][ix  ]*C1*r;
//        //up[iz][ix] += ui[iz  ][ixm1]*C1*r;
//        //up[iz][ix] += ui[iz  ][ixp1]*C1*r;
//        //up[iz][ix] += ui[izm1][ixm1]*C2*r;
//        //up[iz][ix] += ui[izp1][ixm1]*C2*r;
//        //up[iz][ix] += ui[izm1][ixp1]*C2*r;
//        //up[iz][ix] += ui[izp1][ixp1]*C2*r;
//        //up[iz][ix] += ui[izm2][ix  ]*C3*r;
//        //up[iz][ix] += ui[izp2][ix  ]*C3*r;
//        //up[iz][ix] += ui[iz  ][ixm2]*C3*r;
//        //up[iz][ix] += ui[iz  ][ixp2]*C3*r;
//        //up[iz][ix] += ui[izm2][ixm2]*C4*r;
//        //up[iz][ix] += ui[izp2][ixm2]*C4*r;
//        //up[iz][ix] += ui[izm2][ixp2]*C4*r;
//        //up[iz][ix] += ui[izp2][ixp2]*C4*r;
//      }
//    }});
//  }
//
//  private void adjointStep(
//  final float[][] um, final float[][] ui, final float[][] up) {
//
//    //for (int iz=_iza; iz<_izd; ++iz) {
//    //  int izm1 = iz-1;
//    //  int izp1 = iz+1;
//    //  int izm2 = iz-2;
//    //  int izp2 = iz+2;
//    //  for (int ix=_ixa; ix<_ixd; ++ix) {
//    //    int ixm1 = ix-1;
//    //    int ixp1 = ix+1;
//    //    int ixm2 = ix-2;
//    //    int ixp2 = ix+2;
//    //    float r = _r[iz][ix];
//    //    um[iz  ][ix  ] -= up[iz][ix];
//    //    ui[iz  ][ix  ] += up[iz][ix]*(2.0f+C0*r);
//    //    ui[izm1][ix  ] += up[iz][ix]*C1*r;
//    //    ui[izp1][ix  ] += up[iz][ix]*C1*r;
//    //    ui[iz  ][ixm1] += up[iz][ix]*C1*r;
//    //    ui[iz  ][ixp1] += up[iz][ix]*C1*r;
//    //    ui[izm1][ixm1] += up[iz][ix]*C2*r;
//    //    ui[izp1][ixm1] += up[iz][ix]*C2*r;
//    //    ui[izm1][ixp1] += up[iz][ix]*C2*r;
//    //    ui[izp1][ixp1] += up[iz][ix]*C2*r;
//    //    ui[izm2][ix  ] += up[iz][ix]*C3*r;
//    //    ui[izp2][ix  ] += up[iz][ix]*C3*r;
//    //    ui[iz  ][ixm2] += up[iz][ix]*C3*r;
//    //    ui[iz  ][ixp2] += up[iz][ix]*C3*r;
//    //    ui[izm2][ixm2] += up[iz][ix]*C4*r;
//    //    ui[izp2][ixm2] += up[iz][ix]*C4*r;
//    //    ui[izm2][ixp2] += up[iz][ix]*C4*r;
//    //    ui[izp2][ixp2] += up[iz][ix]*C4*r;
//    //  }
//    //}
//
//    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
//    public void compute(int iz) {
//      for (int ix=_ixa; ix<_ixd; ++ix) {
//        float r = _r[iz][ix];
//        ui[iz-2][ix  ] += up[iz][ix]*C3*r;
//        ui[iz-2][ix-2] += up[iz][ix]*C4*r;
//        ui[iz-2][ix+2] += up[iz][ix]*C4*r;
//      }
//    }});
//    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
//    public void compute(int iz) {
//      for (int ix=_ixa; ix<_ixd; ++ix) {
//        float r = _r[iz][ix];
//        ui[iz-1][ix  ] += up[iz][ix]*C1*r;
//        ui[iz-1][ix-1] += up[iz][ix]*C2*r;
//        ui[iz-1][ix+1] += up[iz][ix]*C2*r;
//      }
//    }});
//    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
//    public void compute(int iz) {
//      for (int ix=_ixa; ix<_ixd; ++ix) {
//        float r = _r[iz][ix];
//        um[iz  ][ix  ] -= up[iz][ix];
//        ui[iz  ][ix  ] += up[iz][ix]*(2.0f+C0*r);
//        ui[iz  ][ix-1] += up[iz][ix]*C1*r;
//        ui[iz  ][ix+1] += up[iz][ix]*C1*r;
//        ui[iz  ][ix-2] += up[iz][ix]*C3*r;
//        ui[iz  ][ix+2] += up[iz][ix]*C3*r;
//      }
//    }});
//    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
//    public void compute(int iz) {
//      for (int ix=_ixa; ix<_ixd; ++ix) {
//        float r = _r[iz][ix];
//        ui[iz+1][ix  ] += up[iz][ix]*C1*r;
//        ui[iz+1][ix-1] += up[iz][ix]*C2*r;
//        ui[iz+1][ix+1] += up[iz][ix]*C2*r;
//      }
//    }});
//    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
//    public void compute(int iz) {
//      for (int ix=_ixa; ix<_ixd; ++ix) {
//        float r = _r[iz][ix];
//        ui[iz+2][ix  ] += up[iz][ix]*C3*r;
//        ui[iz+2][ix-2] += up[iz][ix]*C4*r;
//        ui[iz+2][ix+2] += up[iz][ix]*C4*r;
//      }
//    }});
//
//    //Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
//    //public void compute(int iz) {
//    //  int izm1 = iz-1;
//    //  int izp1 = iz+1;
//    //  int izm2 = iz-2;
//    //  int izp2 = iz+2;
//    //  for (int ix=_ixa; ix<_ixd; ++ix) {
//    //    int ixm1 = ix-1;
//    //    int ixp1 = ix+1;
//    //    int ixm2 = ix-2;
//    //    int ixp2 = ix+2;
//    //    float r = _r[iz][ix];
//    //    um[iz][ix] -= up[iz  ][ix  ];
//    //    ui[iz][ix] += up[iz  ][ix  ]*(2.0f+C0*r);
//    //    ui[iz][ix] += up[izp1][ix  ]*C1*r;
//    //    ui[iz][ix] += up[izm1][ix  ]*C1*r;
//    //    ui[iz][ix] += up[iz  ][ixp1]*C1*r;
//    //    ui[iz][ix] += up[iz  ][ixm1]*C1*r;
//    //    ui[iz][ix] += up[izp1][ixp1]*C2*r;
//    //    ui[iz][ix] += up[izm1][ixp1]*C2*r;
//    //    ui[iz][ix] += up[izp1][ixm1]*C2*r;
//    //    ui[iz][ix] += up[izm1][ixm1]*C2*r;
//    //    ui[iz][ix] += up[izp2][ix  ]*C3*r;
//    //    ui[iz][ix] += up[izm2][ix  ]*C3*r;
//    //    ui[iz][ix] += up[iz  ][ixp2]*C3*r;
//    //    ui[iz][ix] += up[iz  ][ixm2]*C3*r;
//    //    ui[iz][ix] += up[izp2][ixp2]*C4*r;
//    //    ui[iz][ix] += up[izm2][ixp2]*C4*r;
//    //    ui[iz][ix] += up[izp2][ixm2]*C4*r;
//    //    ui[iz][ix] += up[izm2][ixm2]*C4*r;
//    //}});
//    //for (int iz=_iza; iz<_izd; ++iz) {
//    //  // TODO
//    //  //ui[iz  ][ix+1  ] += up[iz][ix    ]*C1*r;
//    //  //ui[iz  ][_ixd  ] += up[iz][_ixd-1]*C1*_r[iz][_ixd-1];
//    //  //ui[iz  ][_ixa  ] -= up[iz][_ixa-1]*C1*_r[iz][_ixa-1];
//    // 
//    //  //ui[iz  ][ix-1  ] += up[iz][ix  ]*C1*r;
//    //  //ui[iz  ][_ixa-1] += up[iz][_ixa]*C1*_r[iz][_ixa];
//    //  //ui[iz  ][_ixd-1] -= up[iz][_ixd]*C1*_r[iz][_ixd];
//    //}
//
//  }

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

  private void forwardStep(
  final float[][] um, final float[][] ui, final float[][] up) {
    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
    public void compute(int iz) {
      float[] upi = up[iz];
      float[] umi = um[iz];
      float[] uii = ui[iz];
      float[] uim01 = ui[iz-1 ], uip01 = ui[iz+1 ];
      float[] uim02 = ui[iz-2 ], uip02 = ui[iz+2 ];
      float[] uim03 = ui[iz-3 ], uip03 = ui[iz+3 ];
      float[] uim04 = ui[iz-4 ], uip04 = ui[iz+4 ];
      float[] uim05 = ui[iz-5 ], uip05 = ui[iz+5 ];
      float[] uim06 = ui[iz-6 ], uip06 = ui[iz+6 ];
      float[] uim07 = ui[iz-7 ], uip07 = ui[iz+7 ];
      float[] uim08 = ui[iz-8 ], uip08 = ui[iz+8 ];
      float[] uim09 = ui[iz-9 ], uip09 = ui[iz+9 ];
      float[] uim10 = ui[iz-10], uip10 = ui[iz+10];
      for (int ix=_ixa; ix<_ixd; ++ix) {
        float r = _r[iz][ix];
        //float q = 0.5f*r;
        //float a = 0.5461f; // (Jo et. al., 1996)
        float a = 1.0f;
        upi[ix] += (
          a*((2.0f+C00*r)*uii[ix]-umi[ix]+r*(
          C01*(uim01[ix]+uii[ix-1 ]+uii[ix+1 ]+uip01[ix])+
          C02*(uim02[ix]+uii[ix-2 ]+uii[ix+2 ]+uip02[ix])+
          C03*(uim03[ix]+uii[ix-3 ]+uii[ix+3 ]+uip03[ix])+
          C04*(uim04[ix]+uii[ix-4 ]+uii[ix+4 ]+uip04[ix])+
          C05*(uim05[ix]+uii[ix-5 ]+uii[ix+5 ]+uip05[ix])+
          C06*(uim06[ix]+uii[ix-6 ]+uii[ix+6 ]+uip06[ix])+
          C07*(uim07[ix]+uii[ix-7 ]+uii[ix+7 ]+uip07[ix])+
          C08*(uim08[ix]+uii[ix-8 ]+uii[ix+8 ]+uip08[ix])+
          C09*(uim09[ix]+uii[ix-9 ]+uii[ix+9 ]+uip09[ix])+
          C10*(uim10[ix]+uii[ix-10]+uii[ix+10]+uip10[ix])))
          //+
          //(1.0f-a)*((2.0f+C00*q)*uii[ix]-umi[ix]+q*(
          //C01*(uim01[ix-1 ]+uim01[ix+1 ]+uip01[ix-1 ]+uip01[ix+1 ])+
          //C02*(uim02[ix-2 ]+uim02[ix+2 ]+uip02[ix-2 ]+uip02[ix+2 ])+
          //C03*(uim03[ix-3 ]+uim03[ix+3 ]+uip03[ix-3 ]+uip03[ix+3 ])+
          //C04*(uim04[ix-4 ]+uim04[ix+4 ]+uip04[ix-4 ]+uip04[ix+4 ])+
          //C05*(uim05[ix-5 ]+uim05[ix+5 ]+uip05[ix-5 ]+uip05[ix+5 ])+
          //C06*(uim06[ix-6 ]+uim06[ix+6 ]+uip06[ix-6 ]+uip06[ix+6 ])+
          //C07*(uim07[ix-7 ]+uim07[ix+7 ]+uip07[ix-7 ]+uip07[ix+7 ])+
          //C08*(uim08[ix-8 ]+uim08[ix+8 ]+uip08[ix-8 ]+uip08[ix+8 ])+
          //C09*(uim09[ix-9 ]+uim09[ix+9 ]+uip09[ix-9 ]+uip09[ix+9 ])+
          //C10*(uim10[ix-10]+uim10[ix+10]+uip10[ix-10]+uip10[ix+10])))
        );
      }
    }});
  }
  private void adjointStep(
  final float[][] um, final float[][] ui, final float[][] up) {
    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
    public void compute(int iz) {
      float[] umi = um[iz];
      float[] uii = ui[iz];
      float[] upi = up[iz];
      for (int ix=_ixa; ix<_ixd; ++ix) {
        float r = _r[iz][ix];
        float upr = upi[ix]*r;
        umi[ix   ] -= upi[ix];
        uii[ix-10] += upr*C10;
        uii[ix-9 ] += upr*C09;
        uii[ix-8 ] += upr*C08;
        uii[ix-7 ] += upr*C07;
        uii[ix-6 ] += upr*C06;
        uii[ix-5 ] += upr*C05;
        uii[ix-4 ] += upr*C04;
        uii[ix-3 ] += upr*C03;
        uii[ix-2 ] += upr*C02;
        uii[ix-1 ] += upr*C01;
        uii[ix   ] += upi[ix]*(2.0f+C00*r);
        uii[ix+1 ] += upr*C01;
        uii[ix+2 ] += upr*C02;
        uii[ix+3 ] += upr*C03;
        uii[ix+4 ] += upr*C04;
        uii[ix+5 ] += upr*C05;
        uii[ix+6 ] += upr*C06;
        uii[ix+7 ] += upr*C07;
        uii[ix+8 ] += upr*C08;
        uii[ix+9 ] += upr*C09;
        uii[ix+10] += upr*C10;
      }
    }});
    Parallel.loop(_ixa,_ixd,new Parallel.LoopInt() {
    public void compute(int ix) {
      for (int iz=_iza; iz<_izd; ++iz) {
        float r = _r[iz][ix];
        float upr = up[iz][ix]*r;
        ui[iz-10][ix   ] += upr*C10;
        ui[iz-9 ][ix   ] += upr*C09;
        ui[iz-8 ][ix   ] += upr*C08;
        ui[iz-7 ][ix   ] += upr*C07;
        ui[iz-6 ][ix   ] += upr*C06;
        ui[iz-5 ][ix   ] += upr*C05;
        ui[iz-4 ][ix   ] += upr*C04;
        ui[iz-3 ][ix   ] += upr*C03;
        ui[iz-2 ][ix   ] += upr*C02;
        ui[iz-1 ][ix   ] += upr*C01;
        ui[iz+1 ][ix   ] += upr*C01;
        ui[iz+2 ][ix   ] += upr*C02;
        ui[iz+3 ][ix   ] += upr*C03;
        ui[iz+4 ][ix   ] += upr*C04;
        ui[iz+5 ][ix   ] += upr*C05;
        ui[iz+6 ][ix   ] += upr*C06;
        ui[iz+7 ][ix   ] += upr*C07;
        ui[iz+8 ][ix   ] += upr*C08;
        ui[iz+9 ][ix   ] += upr*C09;
        ui[iz+10][ix   ] += upr*C10;
      }
    }});
  }
  private void xadjointStep( // FIXME: not exact adjoint!
  final float[][] um, final float[][] ui, final float[][] up) {
    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
    public void compute(int iz) {
      for (int ix=_ixa; ix<_ixd; ++ix) {
        float r = _r[iz][ix];
        um[iz][ix] += (2.0f+C00*r)*ui[iz][ix]-up[iz][ix]+
          C01*(_r[iz-1 ][ix   ]*ui[iz-1 ][ix   ]+
               _r[iz   ][ix-1 ]*ui[iz   ][ix-1 ]+
               _r[iz   ][ix+1 ]*ui[iz   ][ix+1 ]+
               _r[iz+1 ][ix   ]*ui[iz+1 ][ix   ])+
          C02*(_r[iz-2 ][ix   ]*ui[iz-2 ][ix   ]+
               _r[iz   ][ix-2 ]*ui[iz   ][ix-2 ]+
               _r[iz   ][ix+2 ]*ui[iz   ][ix+2 ]+
               _r[iz+2 ][ix   ]*ui[iz+2 ][ix   ])+
          C03*(_r[iz-3 ][ix   ]*ui[iz-3 ][ix   ]+
               _r[iz   ][ix-3 ]*ui[iz   ][ix-3 ]+
               _r[iz   ][ix+3 ]*ui[iz   ][ix+3 ]+
               _r[iz+3 ][ix   ]*ui[iz+3 ][ix   ])+
          C04*(_r[iz-4 ][ix   ]*ui[iz-4 ][ix   ]+
               _r[iz   ][ix-4 ]*ui[iz   ][ix-4 ]+
               _r[iz   ][ix+4 ]*ui[iz   ][ix+4 ]+
               _r[iz+4 ][ix   ]*ui[iz+4 ][ix   ])+
          C05*(_r[iz-5 ][ix   ]*ui[iz-5 ][ix   ]+
               _r[iz   ][ix-5 ]*ui[iz   ][ix-5 ]+
               _r[iz   ][ix+5 ]*ui[iz   ][ix+5 ]+
               _r[iz+5 ][ix   ]*ui[iz+5 ][ix   ])+
          C06*(_r[iz-6 ][ix   ]*ui[iz-6 ][ix   ]+
               _r[iz   ][ix-6 ]*ui[iz   ][ix-6 ]+
               _r[iz   ][ix+6 ]*ui[iz   ][ix+6 ]+
               _r[iz+6 ][ix   ]*ui[iz+6 ][ix   ])+
          C07*(_r[iz-7 ][ix   ]*ui[iz-7 ][ix   ]+
               _r[iz   ][ix-7 ]*ui[iz   ][ix-7 ]+
               _r[iz   ][ix+7 ]*ui[iz   ][ix+7 ]+
               _r[iz+7 ][ix   ]*ui[iz+7 ][ix   ])+
          C08*(_r[iz-8 ][ix   ]*ui[iz-8 ][ix   ]+
               _r[iz   ][ix-8 ]*ui[iz   ][ix-8 ]+
               _r[iz   ][ix+8 ]*ui[iz   ][ix+8 ]+
               _r[iz+8 ][ix   ]*ui[iz+8 ][ix   ])+
          C09*(_r[iz-9 ][ix   ]*ui[iz-9 ][ix   ]+
               _r[iz   ][ix-9 ]*ui[iz   ][ix-9 ]+
               _r[iz   ][ix+9 ]*ui[iz   ][ix+9 ]+
               _r[iz+9 ][ix   ]*ui[iz+9 ][ix   ])+
          C10*(_r[iz-10][ix   ]*ui[iz-10][ix   ]+
               _r[iz   ][ix-10]*ui[iz   ][ix-10]+
               _r[iz   ][ix+10]*ui[iz   ][ix+10]+
               _r[iz+10][ix   ]*ui[iz+10][ix   ]);
      }
    }});
  }

  private static final float SQRT2 = sqrt(2.0f);

  // Liu, Y. and M. K. Sen, 2010, A hybrid scheme for absorbing
  // edge reflections in numerical modeling of wave propagation.
  public void forwardAbsorb(
  final float[][] um, final float[][] ui, final float[][] up) {
    final float oxt = (float)(1.0/(_dx*_dt));
    final float oxx = (float)(1.0/(_dx*_dx));
    final float ott = (float)(1.0/(_dt*_dt));
    final float[][] cp = up;
    //final float[][] cp = copy(up);

    for (int ix=_ixa; ix<_ixb-1; ++ix) {
      for (int iz=_izb; iz<_izc; ++iz) {
        float si = _s[iz][ix];
        float m =  1.0f-_w[iz][ix];
        float a =  0.50f*oxt;
        float b = -0.50f*ott*si;
        float c =  0.25f*oxx/si;
        float d =  1.0f/(a-b)*m;
        //up[iz][ix] -= m*up[iz][ix];
        //up[iz][ix] += um[iz-1][ix  ]*c*d;
        //up[iz][ix] += um[iz+1][ix  ]*c*d;
        //up[iz][ix] += cp[iz+1][ix+1]*c*d;
        //up[iz][ix] += cp[iz-1][ix+1]*c*d;
        //up[iz][ix] += um[iz  ][ix  ]*(a+b-2.0f*c)*d;
        //up[iz][ix] += cp[iz  ][ix+1]*(a+b-2.0f*c)*d;
        //up[iz][ix] += ui[iz  ][ix  ]*(-2.0f*b)*d;
        //up[iz][ix] += ui[iz  ][ix+1]*(-2.0f*b)*d;
        //up[iz][ix] += um[iz  ][ix+1]*(b-a)*d;
        up[iz][ix] += -m*up[iz][ix]+d*(
          (um[iz-1][ix  ]+um[iz+1][ix  ]+cp[iz+1][ix+1]+cp[iz-1][ix+1])*c+
          (um[iz  ][ix  ]+cp[iz  ][ix+1])*(a+b-2.0f*c)+
          (ui[iz  ][ix  ]+ui[iz  ][ix+1])*(-2.0f*b)+
          (um[iz  ][ix+1])*(b-a)
        );
      }
    }
    for (int ix=_ixd-1; ix>=_ixc+1; --ix) {
      for (int iz=_izb; iz<_izc; ++iz) {
        float si = _s[iz][ix];
        float m =  1.0f-_w[iz][ix];
        float a = -0.50f*oxt;
        float b =  0.50f*ott*si;
        float c = -0.25f*oxx/si;
        float d =  1.0f/(a-b)*m;
        up[iz][ix] += -m*up[iz][ix]+d*(
          (um[iz-1][ix  ]+um[iz+1][ix  ]+cp[iz+1][ix-1]+cp[iz-1][ix-1])*c+
          (um[iz  ][ix  ]+cp[iz  ][ix-1])*(a+b-2.0f*c)+
          (ui[iz  ][ix  ]+ui[iz  ][ix-1])*(-2.0f*b)+
          (um[iz  ][ix-1])*(b-a)
        );
      }
    }
    for (int iz=_iza; iz<_izb-1; ++iz) {
      for (int ix=_ixb; ix<_ixc; ++ix) {
        float si = _s[iz][ix];
        float m =  1.0f-_w[iz][ix];
        float a = -0.50f*oxt;
        float b =  0.50f*ott*si;
        float c = -0.25f*oxx/si;
        float d =  1.0f/(a-b)*m;
        up[iz][ix] += -m*up[iz][ix]+d*(
          (um[iz  ][ix-1]+um[iz  ][ix+1]+cp[iz+1][ix+1]+cp[iz+1][ix-1])*c+
          (um[iz  ][ix  ]+cp[iz+1][ix  ])*(a+b-2.0f*c)+
          (ui[iz  ][ix  ]+ui[iz+1][ix  ])*(-2.0f*b)+
          (um[iz+1][ix  ])*(b-a)
        );
      }
    }
    for (int iz=_izd-1; iz>=_izc+1; --iz) {
      for (int ix=_ixb; ix<_ixc; ++ix) {
        float si = _s[iz][ix];
        float m =  1.0f-_w[iz][ix];
        float a = -0.50f*oxt;
        float b =  0.50f*ott*si;
        float c = -0.25f*oxx/si;
        float d =  1.0f/(a-b)*m;
        up[iz][ix] += -m*up[iz][ix]+d*(
          (um[iz  ][ix-1]+um[iz  ][ix+1]+cp[iz-1][ix+1]+cp[iz-1][ix-1])*c+
          (um[iz  ][ix  ]+cp[iz-1][ix  ])*(a+b-2.0f*c)+
          (ui[iz  ][ix  ]+ui[iz-1][ix  ])*(-2.0f*b)+
          (um[iz-1][ix  ])*(b-a)
        );
      }
    }
    // Corners
    for (int iz=_iza; iz<_izb; ++iz) {
      for (int ix=_ixa; ix<_ixb; ++ix) {
        float m = 1.0f-_w[iz][ix];
        float r = _dt/(SQRT2*_dx*_s[iz][ix]);
        up[iz][ix] += -m*up[iz][ix]+m*(
          ui[iz][ix]+r*(ui[iz+1][ix]+ui[iz][ix+1]-2.0f*ui[iz][ix])
        );
      }
    }
    for (int iz=_izc; iz<_izd; ++iz) {
      for (int ix=_ixa; ix<_ixb; ++ix) {
        float m = 1.0f-_w[iz][ix];
        float r = _dt/(SQRT2*_dx*_s[iz][ix]);
        up[iz][ix] += -m*up[iz][ix]+m*(
          ui[iz][ix]+r*(ui[iz-1][ix]+ui[iz][ix+1]-2.0f*ui[iz][ix])
        );
      }
    }
    for (int iz=_iza; iz<_izb; ++iz) {
      for (int ix=_ixc; ix<_ixd; ++ix) {
        float m = 1.0f-_w[iz][ix];
        float r = _dt/(SQRT2*_dx*_s[iz][ix]);
        up[iz][ix] += -m*up[iz][ix]+m*(
          ui[iz][ix]+r*(ui[iz+1][ix]+ui[iz][ix-1]-2.0f*ui[iz][ix])
        );
      }
    }
    for (int iz=_izc; iz<_izd; ++iz) {
      for (int ix=_ixc; ix<_ixd; ++ix) {
        float m = 1.0f-_w[iz][ix];
        float r = _dt/(SQRT2*_dx*_s[iz][ix]);
        up[iz][ix] += -m*up[iz][ix]+m*(
          ui[iz][ix]+r*(ui[iz-1][ix]+ui[iz][ix-1]-2.0f*ui[iz][ix])
        );
      }
    }
  }

  public void adjointAbsorb(
  final float[][] um, final float[][] ui, final float[][] up) {
    final float oxt = (float)(1.0/(_dx*_dt));
    final float oxx = (float)(1.0/(_dx*_dx));
    final float ott = (float)(1.0/(_dt*_dt));
    final float[][] cp = copy(up); // possible to avoid copy?

    for (int ix=_ixa; ix<_ixb-1; ++ix) {
      for (int iz=_izb; iz<_izc; ++iz) {
        float si = _s[iz][ix];
        float m =  1.0f-_w[iz][ix];
        float a =  0.50f*oxt;
        float b = -0.50f*ott*si;
        float c =  0.25f*oxx/si;
        float d =  1.0f/(a-b)*m;
        up[iz  ][ix  ] -= cp[iz][ix]*m;
        um[iz-1][ix  ] += cp[iz][ix]*c*d;
        um[iz+1][ix  ] += cp[iz][ix]*c*d;
        up[iz+1][ix+1] += cp[iz][ix]*c*d;
        up[iz-1][ix+1] += cp[iz][ix]*c*d;
        um[iz  ][ix  ] += cp[iz][ix]*(a+b-2.0f*c)*d;
        up[iz  ][ix+1] += cp[iz][ix]*(a+b-2.0f*c)*d;
        ui[iz  ][ix  ] += cp[iz][ix]*(-2.0f*b)*d;
        ui[iz  ][ix+1] += cp[iz][ix]*(-2.0f*b)*d;
        um[iz  ][ix+1] += cp[iz][ix]*(b-a)*d;
      }
    }
    for (int ix=_ixd-1; ix>=_ixc+1; --ix) {
      for (int iz=_izb; iz<_izc; ++iz) {
        float si = _s[iz][ix];
        float m =  1.0f-_w[iz][ix];
        float a = -0.50f*oxt;
        float b =  0.50f*ott*si;
        float c = -0.25f*oxx/si;
        float d =  1.0f/(a-b)*m;
        up[iz  ][ix  ] -= cp[iz][ix]*m;
        um[iz-1][ix  ] += cp[iz][ix]*c*d;
        um[iz+1][ix  ] += cp[iz][ix]*c*d;
        up[iz+1][ix-1] += cp[iz][ix]*c*d;
        up[iz-1][ix-1] += cp[iz][ix]*c*d;
        um[iz  ][ix  ] += cp[iz][ix]*(a+b-2.0f*c)*d;
        up[iz  ][ix-1] += cp[iz][ix]*(a+b-2.0f*c)*d;
        ui[iz  ][ix  ] += cp[iz][ix]*(-2.0f*b)*d;
        ui[iz  ][ix-1] += cp[iz][ix]*(-2.0f*b)*d;
        um[iz  ][ix-1] += cp[iz][ix]*(b-a)*d;
      }
    }
    for (int iz=_iza; iz<_izb-1; ++iz) {
      for (int ix=_ixb; ix<_ixc; ++ix) {
        float si = _s[iz][ix];
        float m =  1.0f-_w[iz][ix];
        float a = -0.50f*oxt;
        float b =  0.50f*ott*si;
        float c = -0.25f*oxx/si;
        float d =  1.0f/(a-b)*m;
        up[iz  ][ix  ] -= cp[iz][ix]*m;
        um[iz  ][ix-1] += cp[iz][ix]*c*d;
        um[iz  ][ix+1] += cp[iz][ix]*c*d;
        up[iz+1][ix+1] += cp[iz][ix]*c*d;
        up[iz+1][ix-1] += cp[iz][ix]*c*d;
        um[iz  ][ix  ] += cp[iz][ix]*(a+b-2.0f*c)*d;
        up[iz+1][ix  ] += cp[iz][ix]*(a+b-2.0f*c)*d;
        ui[iz  ][ix  ] += cp[iz][ix]*(-2.0f*b)*d;
        ui[iz+1][ix  ] += cp[iz][ix]*(-2.0f*b)*d;
        um[iz+1][ix  ] += cp[iz][ix]*(b-a)*d;
      }
    }
    for (int iz=_izd-1; iz>=_izc+1; --iz) {
      for (int ix=_ixb; ix<_ixc; ++ix) {
        float si = _s[iz][ix];
        float m =  1.0f-_w[iz][ix];
        float a = -0.50f*oxt;
        float b =  0.50f*ott*si;
        float c = -0.25f*oxx/si;
        float d =  1.0f/(a-b)*m;
        up[iz  ][ix  ] -= cp[iz][ix]*m;
        um[iz  ][ix-1] += cp[iz][ix]*c*d;
        um[iz  ][ix+1] += cp[iz][ix]*c*d;
        up[iz-1][ix+1] += cp[iz][ix]*c*d;
        up[iz-1][ix-1] += cp[iz][ix]*c*d;
        um[iz  ][ix  ] += cp[iz][ix]*(a+b-2.0f*c)*d;
        up[iz-1][ix  ] += cp[iz][ix]*(a+b-2.0f*c)*d;
        ui[iz  ][ix  ] += cp[iz][ix]*(-2.0f*b)*d;
        ui[iz-1][ix  ] += cp[iz][ix]*(-2.0f*b)*d;
        um[iz-1][ix  ] += cp[iz][ix]*(b-a)*d;
      }
    }
    // Corners
    for (int iz=_iza; iz<_izb; ++iz) {
      for (int ix=_ixa; ix<_ixb; ++ix) {
        float m = 1.0f-_w[iz][ix];
        float r = _dt/(SQRT2*_dx*_s[iz][ix]);
        up[iz  ][ix  ] -= cp[iz][ix]*m;
        ui[iz  ][ix  ] += cp[iz][ix]*m;
        ui[iz+1][ix  ] += cp[iz][ix]*m*r;
        ui[iz  ][ix+1] += cp[iz][ix]*m*r;
        ui[iz  ][ix  ] -= cp[iz][ix]*m*r*2.0f;
      }
    }
    for (int iz=_izc; iz<_izd; ++iz) {
      for (int ix=_ixa; ix<_ixb; ++ix) {
        float m = 1.0f-_w[iz][ix];
        float r = _dt/(SQRT2*_dx*_s[iz][ix]);
        up[iz  ][ix  ] -= cp[iz][ix]*m;
        ui[iz  ][ix  ] += cp[iz][ix]*m;
        ui[iz-1][ix  ] += cp[iz][ix]*m*r;
        ui[iz  ][ix+1] += cp[iz][ix]*m*r;
        ui[iz  ][ix  ] -= cp[iz][ix]*m*r*2.0f;
      }
    }
    for (int iz=_iza; iz<_izb; ++iz) {
      for (int ix=_ixc; ix<_ixd; ++ix) {
        float m = 1.0f-_w[iz][ix];
        float r = _dt/(SQRT2*_dx*_s[iz][ix]);
        up[iz  ][ix  ] -= cp[iz][ix]*m;
        ui[iz  ][ix  ] += cp[iz][ix]*m;
        ui[iz+1][ix  ] += cp[iz][ix]*m*r;
        ui[iz  ][ix-1] += cp[iz][ix]*m*r;
        ui[iz  ][ix  ] -= cp[iz][ix]*m*r*2.0f;
      }
    }
    for (int iz=_izc; iz<_izd; ++iz) {
      for (int ix=_ixc; ix<_ixd; ++ix) {
        float m = 1.0f-_w[iz][ix];
        float r = _dt/(SQRT2*_dx*_s[iz][ix]);
        up[iz  ][ix  ] -= cp[iz][ix]*m;
        ui[iz  ][ix  ] += cp[iz][ix]*m;
        ui[iz-1][ix  ] += cp[iz][ix]*m*r;
        ui[iz  ][ix-1] += cp[iz][ix]*m*r;
        ui[iz  ][ix  ] -= cp[iz][ix]*m*r*2.0f;
      }
    }
  }

  private void absorb(
  final float[][] um, final float[][] ui, final float[][] up) {
    final float oxt = (float)(1.0/(_dx*_dt));
    final float oxx = (float)(1.0/(_dx*_dx));
    final float ott = (float)(1.0/(_dt*_dt));
    //final float[][] cp = copy(up);
    final float[][] cp = up;
    for (int ix=_ixa; ix<_ixb-1; ++ix) {
      for (int iz=_izb; iz<_izc; ++iz) {
        float si = _s[iz][ix];
        float a =  0.50f*oxt;
        float b = -0.50f*ott*si;
        float c =  0.25f*oxx/si;
        float vp = 1.0f/(a-b)*(
          (um[iz-1][ix  ]+um[iz+1][ix  ]+cp[iz+1][ix+1]+cp[iz-1][ix+1])*c+
          (um[iz  ][ix  ]+cp[iz  ][ix+1])*(a+b-2.0f*c)+
          (ui[iz  ][ix  ]+ui[iz  ][ix+1])*(-2.0f*b)+
          (um[iz  ][ix+1])*(b-a)
        );
        float w = _w[iz][ix];
        up[iz][ix] = w*up[iz][ix]+(1.0f-w)*vp;
      }
    }
    for (int ix=_ixd-1; ix>=_ixc+1; --ix) {
      for (int iz=_izb; iz<_izc; ++iz) {
        float si = _s[iz][ix];
        float a = -0.50f*oxt;
        float b =  0.50f*ott*si;
        float c = -0.25f*oxx/si;
        float vp = 1.0f/(a-b)*(
          (um[iz-1][ix  ]+um[iz+1][ix  ]+cp[iz+1][ix-1]+cp[iz-1][ix-1])*c+
          (um[iz  ][ix  ]+cp[iz  ][ix-1])*(a+b-2.0f*c)+
          (ui[iz  ][ix  ]+ui[iz  ][ix-1])*(-2.0f*b)+
          (um[iz  ][ix-1])*(b-a)
        );
        float w = _w[iz][ix];
        up[iz][ix] = w*up[iz][ix]+(1.0f-w)*vp;
      }
    }
    for (int iz=_iza; iz<_izb-1; ++iz) {
      for (int ix=_ixb; ix<_ixc; ++ix) {
        float si = _s[iz][ix];
        float a =  0.50f*oxt;
        float b = -0.50f*ott*si;
        float c =  0.25f*oxx/si;
        float vp = 1.0f/(a-b)*(
          (um[iz  ][ix-1]+um[iz  ][ix+1]+cp[iz+1][ix+1]+cp[iz+1][ix-1])*c+
          (um[iz  ][ix  ]+cp[iz+1][ix  ])*(a+b-2.0f*c)+
          (ui[iz  ][ix  ]+ui[iz+1][ix  ])*(-2.0f*b)+
          (um[iz+1][ix  ])*(b-a)
        );
        float w = _w[iz][ix];
        up[iz][ix] = w*up[iz][ix]+(1.0f-w)*vp;
      }
    }
    for (int iz=_izd-1; iz>=_izc+1; --iz) {
      for (int ix=_ixb; ix<_ixc; ++ix) {
        float si = _s[iz][ix];
        float a = -0.50f*oxt;
        float b =  0.50f*ott*si;
        float c = -0.25f*oxx/si;
        float vp = 1.0f/(a-b)*(
          (um[iz  ][ix-1]+um[iz  ][ix+1]+cp[iz-1][ix+1]+cp[iz-1][ix-1])*c+
          (um[iz  ][ix  ]+cp[iz-1][ix  ])*(a+b-2.0f*c)+
          (ui[iz  ][ix  ]+ui[iz-1][ix  ])*(-2.0f*b)+
          (um[iz-1][ix  ])*(b-a)
        );
        float w = _w[iz][ix];
        up[iz][ix] = w*up[iz][ix]+(1.0f-w)*vp;
      }
    }
    // Corners
    for (int iz=_iza; iz<_izb; ++iz) {
      for (int ix=_ixa; ix<_ixb; ++ix) {
        float w = _w[iz][ix];
        float r = _dt/(SQRT2*_dx*_s[iz][ix]);
        float vp = ui[iz][ix]+r*(ui[iz+1][ix]+ui[iz][ix+1]-2.0f*ui[iz][ix]);
        up[iz][ix] = w*up[iz][ix]+(1.0f-w)*vp;
      }
    }
    for (int iz=_izc; iz<_izd; ++iz) {
      for (int ix=_ixa; ix<_ixb; ++ix) {
        float w = _w[iz][ix];
        float r = _dt/(SQRT2*_dx*_s[iz][ix]);
        float vp = ui[iz][ix]+r*(ui[iz-1][ix]+ui[iz][ix+1]-2.0f*ui[iz][ix]);
        up[iz][ix] = w*up[iz][ix]+(1.0f-w)*vp;
      }
    }
    for (int iz=_iza; iz<_izb; ++iz) {
      for (int ix=_ixc; ix<_ixd; ++ix) {
        float w = _w[iz][ix];
        float r = _dt/(SQRT2*_dx*_s[iz][ix]);
        float vp = ui[iz][ix]+r*(ui[iz+1][ix]+ui[iz][ix-1]-2.0f*ui[iz][ix]);
        up[iz][ix] = w*up[iz][ix]+(1.0f-w)*vp;
      }
    }
    for (int iz=_izc; iz<_izd; ++iz) {
      for (int ix=_ixc; ix<_ixd; ++ix) {
        float w = _w[iz][ix];
        float r = _dt/(SQRT2*_dx*_s[iz][ix]);
        float vp = ui[iz][ix]+r*(ui[iz-1][ix]+ui[iz][ix-1]-2.0f*ui[iz][ix]);
        up[iz][ix] = w*up[iz][ix]+(1.0f-w)*vp;
      }
    }
  }

  private float[][] makeWeights() {
    int fd = FD_ORDER/2;
    float ob = 1.0f/(_b-1.0f);
    float[][] w = new float[_nz][_nx];
    for (int ix=_ixa; ix<_ixb; ++ix) {
      for (int iz=ix; iz<_nz-ix; ++iz) {
        w[iz][ix      ] = (ix-fd)*ob;
        w[iz][_nx-1-ix] = (ix-fd)*ob;
      }
    }
    for (int iz=_iza; iz<_izb; ++iz) {
      for (int ix=iz; ix<_nx-iz; ++ix) {
        w[iz      ][ix] = (iz-fd)*ob;
        w[_nz-1-iz][ix] = (iz-fd)*ob;
      }
    }
    for (int iz=_izb; iz<_izc; ++iz) {
      for (int ix=_ixb; ix<_ixc; ++ix) {
        w[iz][ix] = 1.0f;
      }
    }
    //SimplePlot.asPixels(w).addColorBar();
    return w;
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
    //SimplePlot.asPixels(v).addColorBar();
    return v;
  }
}
