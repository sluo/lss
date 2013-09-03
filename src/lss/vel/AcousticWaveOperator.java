package lss.vel;

import java.util.Random;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

// testing
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;

public class AcousticWaveOperator {

  // count number of same digits
  public static int compareDigits(float xa, float xb) {
    int significantDigits = 10;
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

  public static float dot(float[][][] u, float[][][] a) {
    return dot(u,a,0);
  }
  public static float dot(float[][][] u, float[][][] a, int nabsorb) {
    int nz = u[0][0].length;
    int nx = u[0].length;
    int nt = u.length;
    float sum = 0.0f;
    for (int it=0; it<nt; ++it)
      for (int ix=nabsorb; ix<nx-nabsorb; ++ix)
        for (int iz=nabsorb; iz<nz-nabsorb; ++iz)
          sum += u[it][ix][iz]*a[it][ix][iz];
    return sum;
  }

  //////////////////////////////////////////////////////////////////////////

  public AcousticWaveOperator(
  float[][] s, double dx, double dt, int nabsorb) {
    int nz = s[0].length;
    int nx = s.length;
    Check.argument(nabsorb>=FD_ORDER/2,"nabsorb>=FD_ORDER/2");
    _nabsorb = nabsorb;
    _b = nabsorb-FD_ORDER/2;
    _nz = nz+2*nabsorb;
    _nx = nx+2*nabsorb;
    _dx = dx;
    _dt = dt;
    //_fz = -nabsorb*_dx;
    //_fx = -nabsorb*_dx;
    //_sz = new Sampling(_nz,_dx,_fz);
    //_sx = new Sampling(_nx,_dx,_fx);
    _ixa = FD_ORDER/2;
    _ixb = _ixa+_b;
    _ixc = _ixb+nx;
    _ixd = _ixc+_b;
    _iza = FD_ORDER/2;
    _izb = _iza+_b;
    _izc = _izb+nz;
    _izd = _izc+_b;
    _w = makeWeights();
    float scale = (float)(_dt*_dt/(_dx*_dx));
    _s = extendModel(s,nabsorb);
    //addRandomBoundary(0.15f,_s);
    _r = copy(_s);
    mul(_r,_r,_r);
    div(scale,_r,_r);
  }

  /* Indexing example for FD_ABSORB/2=2, nabsorb=3, nx=3:

          |ixa|ixb        |ixc|ixd
  |---|---|---|---|---|---|---|---|---|
  | S | S | B | M | M | M | B | S | S |
  |---|---|---|---|---|---|---|---|---|
          |   |           |   |

    S - stencil cells for finite-difference stencil
    B - boundary cells for absorbing boundary
    M - model cells for slowness model

  */

  public void applyForward(
  Source source, float[][][] u) {
    apply(true,source,null,u);
  }

  public void applyForward(
  Source source, Receiver receiver) {
    apply(true,source,receiver,null);
  }

  public void applyForward(
  Source source, Receiver receiver, float[][][] u) {
    apply(true,source,receiver,u);
  }

  public void applyAdjoint(
  Source source, float[][][] u) {
    apply(false,source,null,u);
  }

  public void applyAdjoint(
  Source source, Receiver receiver) {
    apply(false,source,receiver,null);
  }

  public void applyAdjoint(
  Source source, Receiver receiver, float[][][] u) {
    apply(false,source,receiver,u);
  }

  //////////////////////////////////////////////////////////////////////////
  // static

  public static float[][] crop(float[][] x, int nabsorb) {
    int nz = x[0].length;
    int nx = x.length;
    return copy(nz-2*nabsorb,nx-2*nabsorb,nabsorb,nabsorb,x);
  }

  public static float[][] collapse(
  float[][][] u, float[][][] a, int nabsorb) {
    int nz = u[0][0].length;
    int nx = u[0].length;
    int nt = u.length;
    float[][] y = new float[nx-2*nabsorb][nz-2*nabsorb];
    for (int it=0; it<nt; ++it)
      for (int ix=nabsorb; ix<nx-nabsorb; ++ix)
        for (int iz=nabsorb; iz<nz-nabsorb; ++iz)
          y[ix-nabsorb][iz-nabsorb] += u[it][ix][iz]*a[it][ix][iz];
    return y;
  }

  public static interface Source {
    public void add(float[][] ui, int it, int nabsorb);
  }

  public static class RickerSource implements Source {
    public RickerSource(int xs, int zs, double dt, double fpeak) {
      _tdelay = 1.0/fpeak;
      _fpeak = fpeak;
      _dt = dt;
      _xs = xs;
      _zs = zs;
    }
    public void add(float[][] ui, int it, int nabsorb) {
      ui[_xs+nabsorb][_zs+nabsorb] += ricker(it*_dt-_tdelay);
    }
    private float ricker(double t) {
      double x = PI*_fpeak*t;
      double xx = x*x;
      return (float)((1.0-2.0*xx)*exp(-xx));
    }
    private int _xs,_zs;
    private double _dt,_fpeak,_tdelay;
  }

  public static class ReceiverSource implements Source {
    public ReceiverSource(Receiver receiver) {
      int[][] c = receiver.getIndices();
      _xr = c[0];
      _zr = c[1];
      _nr = receiver.getNr();
      _nt = receiver.getNt();
      _data = receiver.getData();
    }
    public void add(float[][] ui, int it, int nabsorb) {
      if (it>=_nt) 
        return;
      for (int ir=0; ir<_nr; ++ir) {
        int xs = _xr[ir]+nabsorb;
        int zs = _zr[ir]+nabsorb;
        ui[xs][zs] += _data[ir][it];
      }
    }
    private int _nr,_nt;
    private int[] _xr, _zr;
    private float[][] _data;
  }

  public static class WavefieldSource implements Source {
    public WavefieldSource(float[][][] u) {
      _nz = u[0][0].length;
      _nx = u[0].length;
      _u = u;
    }
    public void add(float[][] ui, int it, int nabsorb) {
      for (int ix=0; ix<_nx; ++ix)
        for (int iz=0; iz<_nz; ++iz)
          ui[ix][iz] += _u[it][ix][iz];
    }
    private float[][][] _u;
    private int _nz,_nx;
  }

  public static class Receiver {
    public Receiver(int xr, int zr, int nt) {
      this(new int[]{xr},new int[]{zr},nt);
    }
    public Receiver(int[] xr, int[] zr, int nt) {
      int nxr = xr.length;
      int nzr = zr.length;
      Check.argument(nxr==nzr,"nxr==nzr");
      _nr = nxr;
      _nt = nt;
      _xr = xr;
      _zr = zr;
      _data = new float[_nr][nt];
    }
    public void setData(float[][] ui, int it, int nabsorb) {
      for (int ir=0; ir<_nr; ++ir) {
        int xr = _xr[ir]+nabsorb;
        int zr = _zr[ir]+nabsorb;
        _data[ir][it] = ui[xr][zr];
      }
    }
    public int[][] getIndices() {
      return new int[][]{_xr,_zr};
    }
    public float[][] getData() {
      return _data;
    }
    public int getNt() {
      return _nt;
    }
    public int getNr() {
      return _nr;
    }
    private int _nr,_nt;
    private int[] _xr, _zr;
    private float[][] _data;
  }

  //////////////////////////////////////////////////////////////////////////
  // private

  private int _b,_nabsorb; // absorbing boundary
  private int _nx,_nz;
  private double _dx,_dt;
  private float[][] _s; // slowness
  private float[][] _r; // dt^2/(dx^2*s^2)
  private float[][] _w; // weights for boundary
  private int _ixa,_ixb,_ixc,_ixd;
  private int _iza,_izb,_izc,_izd;

//  private static final int FD_ORDER = 2;
//  private static final float C0 = -4.0f;
//  private static final float C1 =  1.0f;
//  private void apply(
//  boolean forward, Source source, Receiver receiver, float[][][] u) {
//    int nt = (receiver==null)?u.length:receiver.getNt();
//    float[][] um,ui,up;
//    if (forward) {
//      for (int it=0; it<nt; ++it) {
//        um = (it>0)?u[it-1]:new float[_nx][_nz];
//        ui = u[it];
//        up = (it<nt-1)?u[it+1]:new float[_nx][_nz];
//        forwardStep(um,ui,up);
//      }
//    /*
//    } else {
//      for (int it=nt-1; it>=0; --it) {
//        up = (it>0)?u[it-1]:new float[_nx][_nz];
//        ui = u[it];
//        um = (it<nt-1)?u[it+1]:new float[_nx][_nz];
//        forwardStep(um,ui,up);
//      }
//    */
//    } else {
//      for (int it=nt-1; it>=0; --it) {
//        um = (it>0)?u[it-1]:new float[_nx][_nz];
//        ui = u[it];
//        up = (it<nt-1)?u[it+1]:new float[_nx][_nz];
//        adjointStep(um,ui,up);
//      }
//    }
//  }
//  private void forwardStep(
//  final float[][] um, final float[][] ui, final float[][] up) {
//    Parallel.loop(_ixa,_ixd,new Parallel.LoopInt() {
//    public void compute(int ix) {
//      int ixm1 = ix-1;
//      int ixp1 = ix+1;
//      int ixm2 = ix-2;
//      int ixp2 = ix+2;
//      for (int iz=_iza; iz<_izd; ++iz) {
//        int izm1 = iz-1;
//        int izp1 = iz+1;
//        int izm2 = iz-2;
//        int izp2 = iz+2;
//        float r = _r[ix][iz];
//
//        /*
//        up[ix][iz] += (2.0f+C0*r)*ui[ix][iz]-um[ix][iz]+r*(
//          C1*(ui[ix][izm1]+ui[ix][izp1]+ui[ixm1][iz]+ui[ixp1][iz]));
//        */
//        up[ix][iz] -= um[ix  ][iz  ];
//        up[ix][iz] += ui[ix  ][iz  ]*(2.0f+C0*r);
//        up[ix][iz] += ui[ix  ][izm1]*C1*r;
//        up[ix][iz] += ui[ix  ][izp1]*C1*r;
//        up[ix][iz] += ui[ixp1][iz  ]*C1*r;
//        up[ix][iz] += ui[ixm1][iz  ]*C1*r;
//      }
//    }});
//  }
//  private void adjointStep(
//  final float[][] um, final float[][] ui, final float[][] up) {
//    Parallel.loop(_ixa,_ixd,new Parallel.LoopInt() {
//    public void compute(int ix) {
//      int ixm1 = ix-1;
//      int ixp1 = ix+1;
//      int ixm2 = ix-2;
//      int ixp2 = ix+2;
//      for (int iz=_iza; iz<_izd; ++iz) {
//        int izm1 = iz-1;
//        int izp1 = iz+1;
//        int izm2 = iz-2;
//        int izp2 = iz+2;
//        float r = _r[ix][iz];
//
//        /*
//        um[ix][iz] += (2.0f+C0*r)*ui[ix][iz]-up[ix][iz]+r*(
//          C1*(ui[ix][izm1]+ui[ix][izp1]+ui[ixm1][iz]+ui[ixp1][iz]));
//        */
//        um[ix  ][iz  ] -= up[ix  ][iz  ];
//        ui[ix  ][iz  ] += up[ix  ][iz  ]*(2.0f+C0*r);
//        ui[ix  ][iz  ] += up[ix  ][izp1]*C1*_r[ix  ][izp1];
//        ui[ix  ][iz  ] += up[ix  ][izm1]*C1*_r[ix  ][izm1];
//        ui[ix  ][iz  ] += up[ixm1][iz  ]*C1*_r[ixm1][iz  ];
//        ui[ix  ][iz  ] += up[ixp1][iz  ]*C1*_r[ixp1][iz  ];
//      }
//    }});
//  }

  private void apply(
  boolean forward, Source source, Receiver receiver, float[][][] u) {
    int nt = (u==null)?receiver.getNt():u.length;
    float[][] ut = new float[_nx][_nz]; // temp
    float[][] um = new float[_nx][_nz]; // it-1
    float[][] ui = new float[_nx][_nz]; // it
    float[][] up;                       // it+1
    int fit = (forward)?0:nt-1;
    int pit = (forward)?1:-1;
    for (int it=fit, count=0; count<nt; it+=pit, ++count) {
      if (u==null) { // next time
        up = ut;
        zero(ut); // necessary only if using += for injection
      } else {
        up = u[it];
      }
      //up = (u==null)?ut:u[it]; // next time 
      step(forward,um,ui,up); // time step 
      source.add(up,it,_nabsorb); // inject source (off for adjoint test)
      absorb(um,ui,up); // absorbing boundaries (off for adjoint test)
      if (receiver!=null)
        receiver.setData(up,it,_nabsorb); // data
      ut = um; um = ui; ui = up; // rotate arrays
    }
  }

  private void step(
  final boolean forward,
  final float[][] um, final float[][] ui, final float[][] up) {
    if (forward) {
      Parallel.loop(_ixa,_ixd,new Parallel.LoopInt() {
      public void compute(int ix) {
        forwardStepSliceX(ix,um,ui,up);
      }});
    } else {
      Parallel.loop(_ixa,_ixd,new Parallel.LoopInt() {
      public void compute(int ix) {
        adjointStepSliceX(ix,um,ui,up);
        //forwardStepSliceX(ix,um,ui,up);
      }});
    }
  }

//  // 20th order coefficients
//  private static final int FD_ORDER = 20;
//  private static final float C00 = -0.32148051f*10.0f*2.0f;
//  private static final float C01 =  0.19265816f*10.0f;
//  private static final float C02 = -0.43052632f*1.0f;
//  private static final float C03 =  0.15871000f*1.0f;
//  private static final float C04 = -0.68711400f*0.1f;
//  private static final float C05 =  0.31406935f*0.1f;
//  private static final float C06 = -0.14454222f*0.1f;
//  private static final float C07 =  0.65305182f*0.01f;
//  private static final float C08 = -0.28531535f*0.01f;
//  private static final float C09 =  0.11937032f*0.01f;
//  private static final float C10 = -0.47508613f*0.001f;
//  private void forwardStepSliceX(
//  int ix, float[][] um, float[][] ui, float[][] up) {
//    for (int iz=_iza; iz<_izd; ++iz) {
//      float r = _r[ix][iz];
//      up[ix][iz] += (2.0f+C00*r)*ui[ix][iz]-um[ix][iz]+r*(
//        C01*(ui[ix-1][iz]+ui[ix][iz-1]+ui[ix][iz+1]+ui[ix+1][iz])+
//        C02*(ui[ix-2][iz]+ui[ix][iz-2]+ui[ix][iz+2]+ui[ix+2][iz])+
//        C03*(ui[ix-3][iz]+ui[ix][iz-3]+ui[ix][iz+3]+ui[ix+3][iz])+
//        C04*(ui[ix-4][iz]+ui[ix][iz-4]+ui[ix][iz+4]+ui[ix+4][iz])+
//        C05*(ui[ix-5][iz]+ui[ix][iz-5]+ui[ix][iz+5]+ui[ix+5][iz])+
//        C06*(ui[ix-6][iz]+ui[ix][iz-6]+ui[ix][iz+6]+ui[ix+6][iz])+
//        C07*(ui[ix-7][iz]+ui[ix][iz-7]+ui[ix][iz+7]+ui[ix+7][iz])+
//        C08*(ui[ix-8][iz]+ui[ix][iz-8]+ui[ix][iz+8]+ui[ix+8][iz])+
//        C09*(ui[ix-9][iz]+ui[ix][iz-9]+ui[ix][iz+9]+ui[ix+9][iz])+
//        C10*(ui[ix-10][iz]+ui[ix][iz-10]+ui[ix][iz+10]+ui[ix+10][iz]));
//    }
//  }
//  private void adjointStepSliceX(
//  int ix, float[][] um, float[][] ui, float[][] up) {
//    for (int iz=_iza; iz<_izd; ++iz) {
//      float r = _r[ix][iz];
//      up[ix][iz] += (2.0f+C00*r)*ui[ix][iz]-um[ix][iz]+
//        C01*(_r[ix-1][iz  ]*ui[ix-1][iz  ]+
//             _r[ix  ][iz-1]*ui[ix  ][iz-1]+
//             _r[ix  ][iz+1]*ui[ix  ][iz+1]+
//             _r[ix+1][iz  ]*ui[ix+1][iz  ])+
//        C02*(_r[ix-2][iz  ]*ui[ix-2][iz  ]+
//             _r[ix  ][iz-2]*ui[ix  ][iz-2]+
//             _r[ix  ][iz+2]*ui[ix  ][iz+2]+
//             _r[ix+2][iz  ]*ui[ix+2][iz  ])+
//        C03*(_r[ix-3][iz  ]*ui[ix-3][iz  ]+
//             _r[ix  ][iz-3]*ui[ix  ][iz-3]+
//             _r[ix  ][iz+3]*ui[ix  ][iz+3]+
//             _r[ix+3][iz  ]*ui[ix+3][iz  ])+
//        C04*(_r[ix-4][iz  ]*ui[ix-4][iz  ]+
//             _r[ix  ][iz-4]*ui[ix  ][iz-4]+
//             _r[ix  ][iz+4]*ui[ix  ][iz+4]+
//             _r[ix+4][iz  ]*ui[ix+4][iz  ])+
//        C05*(_r[ix-5][iz  ]*ui[ix-5][iz  ]+
//             _r[ix  ][iz-5]*ui[ix  ][iz-5]+
//             _r[ix  ][iz+5]*ui[ix  ][iz+5]+
//             _r[ix+5][iz  ]*ui[ix+5][iz  ])+
//        C06*(_r[ix-6][iz  ]*ui[ix-6][iz  ]+
//             _r[ix  ][iz-6]*ui[ix  ][iz-6]+
//             _r[ix  ][iz+6]*ui[ix  ][iz+6]+
//             _r[ix+6][iz  ]*ui[ix+6][iz  ])+
//        C07*(_r[ix-7][iz  ]*ui[ix-7][iz  ]+
//             _r[ix  ][iz-7]*ui[ix  ][iz-7]+
//             _r[ix  ][iz+7]*ui[ix  ][iz+7]+
//             _r[ix+7][iz  ]*ui[ix+7][iz  ])+
//        C08*(_r[ix-8][iz  ]*ui[ix-8][iz  ]+
//             _r[ix  ][iz-8]*ui[ix  ][iz-8]+
//             _r[ix  ][iz+8]*ui[ix  ][iz+8]+
//             _r[ix+8][iz  ]*ui[ix+8][iz  ])+
//        C09*(_r[ix-9][iz  ]*ui[ix-9][iz  ]+
//             _r[ix  ][iz-9]*ui[ix  ][iz-9]+
//             _r[ix  ][iz+9]*ui[ix  ][iz+9]+
//             _r[ix+9][iz  ]*ui[ix+9][iz  ])+
//        C10*(_r[ix-10][iz  ]*ui[ix-10][iz  ]+
//             _r[ix  ][iz-10]*ui[ix  ][iz-10]+
//             _r[ix  ][iz+10]*ui[ix  ][iz+10]+
//             _r[ix+10][iz  ]*ui[ix+10][iz  ]);
//    }
//  }

  // First 21-point Laplacian stencil from
  // Patra, M. and M. Karttunen, 2005, Stencils
  // with Isotropic Error for Differential Operators.
  private static final int FD_ORDER = 4;
  private static final float C0 = -9.0f/2.0f;
  private static final float C1 =  16.0f/15.0f;
  private static final float C2 =  2.0f/15.0f;
  private static final float C3 = -1.0f/15.0f;
  private static final float C4 = -1.0f/120.0f;
  private void forwardStepSliceX(
  int ix, float[][] um, float[][] ui, float[][] up) {
    int ixm1 = ix-1;
    int ixp1 = ix+1;
    int ixm2 = ix-2;
    int ixp2 = ix+2;
    for (int iz=_iza; iz<_izd; ++iz) {
      int izm1 = iz-1;
      int izp1 = iz+1;
      int izm2 = iz-2;
      int izp2 = iz+2;
      float r = _r[ix][iz];
      // FIXME: Source injection replaces += needed to pass adjoint test.
      up[ix][iz] += (2.0f+C0*r)*ui[ix][iz]-um[ix][iz]+r*(
        C1*(ui[ix  ][izm1]+ui[ix  ][izp1]+ui[ixm1][iz  ]+ui[ixp1][iz  ])+
        C2*(ui[ixm1][izm1]+ui[ixm1][izp1]+ui[ixp1][izm1]+ui[ixp1][izp1])+
        C3*(ui[ix  ][izm2]+ui[ix  ][izp2]+ui[ixm2][iz  ]+ui[ixp2][iz  ])+
        C4*(ui[ixm2][izm2]+ui[ixm2][izp2]+ui[ixp2][izm2]+ui[ixp2][izp2]));
    }
  }
  private void adjointStepSliceX(
  int ix, float[][] um, float[][] ui, float[][] up) {
    int ixm1 = ix-1;
    int ixp1 = ix+1;
    int ixm2 = ix-2;
    int ixp2 = ix+2;
    for (int iz=_iza; iz<_izd; ++iz) {
      int izm1 = iz-1;
      int izp1 = iz+1;
      int izm2 = iz-2;
      int izp2 = iz+2;
      float r = _r[ix][iz];
      // FIXME: Source injection replaces += needed to pass adjoint test.
      up[ix][iz] += (2.0f+C0*r)*ui[ix][iz]-um[ix][iz]+
        C1*(_r[ix  ][izm1]*
            ui[ix  ][izm1]+
            _r[ix  ][izp1]*
            ui[ix  ][izp1]+
            _r[ixm1][iz  ]*
            ui[ixm1][iz  ]+
            _r[ixp1][iz  ]*
            ui[ixp1][iz  ])+
        C2*(_r[ixm1][izm1]*
            ui[ixm1][izm1]+
            _r[ixm1][izp1]*
            ui[ixm1][izp1]+
            _r[ixp1][izm1]*
            ui[ixp1][izm1]+
            _r[ixp1][izp1]*
            ui[ixp1][izp1])+
        C3*(_r[ix  ][izm2]*
            ui[ix  ][izm2]+
            _r[ix  ][izp2]*
            ui[ix  ][izp2]+
            _r[ixm2][iz  ]*
            ui[ixm2][iz  ]+
            _r[ixp2][iz  ]*
            ui[ixp2][iz  ])+
        C4*(_r[ixm2][izm2]*
            ui[ixm2][izm2]+
            _r[ixm2][izp2]*
            ui[ixm2][izp2]+
            _r[ixp2][izm2]*
            ui[ixp2][izm2]+
            _r[ixp2][izp2]*
            ui[ixp2][izp2]);

    }
  }

//  // Second 21-point Laplacian stencil from
//  // Patra, M. and M. Karttunen, 2005, Stencils
//  // with Isotropic Error for Differential Operators.
//  private static final int FD_ORDER = 4;
//  private static final float C0 = -21.0f/5.0f;
//  private static final float C1 =  13.0f/15.0f;
//  private static final float C2 =  4.0f/15.0f;
//  private static final float C3 = -1.0f/60.0f;
//  private static final float C4 = -1.0f/30.0f;
//  private void forwardStepSliceX(
//  int ix, float[][] um, float[][] ui, float[][] up) {
//    int ixm1 = ix-1;
//    int ixp1 = ix+1;
//    int ixm2 = ix-2;
//    int ixp2 = ix+2;
//    for (int iz=_iza; iz<_izd; ++iz) {
//      int izm1 = iz-1;
//      int izp1 = iz+1;
//      int izm2 = iz-2;
//      int izp2 = iz+2;
//      float r = _r[ix][iz];
//      // FIXME: Source injection replaces += needed to pass adjoint test.
//      up[ix][iz] += (2.0f+C0*r)*ui[ix][iz]-um[ix][iz]+r*(
//        C1*(ui[ix  ][izm1]+ui[ix  ][izp1]+ui[ixm1][iz  ]+ui[ixp1][iz  ])+
//        C2*(ui[ixm1][izm1]+ui[ixm1][izp1]+ui[ixp1][izm1]+ui[ixp1][izp1])+
//        C3*(ui[ix  ][izm2]+ui[ix  ][izp2]+ui[ixm2][iz  ]+ui[ixp2][iz  ])+
//        C4*(ui[ixm2][izm1]+ui[ixm2][izp1]+ui[ixm1][izm2]+ui[ixm1][izp2]+
//            ui[ixp1][izm2]+ui[ixp1][izp2]+ui[ixp2][izm1]+ui[ixp2][izp1]));
//    }
//  }
//  private void adjointStepSliceX(
//  int ix, float[][] um, float[][] ui, float[][] up) {
//    int ixm1 = ix-1;
//    int ixp1 = ix+1;
//    int ixm2 = ix-2;
//    int ixp2 = ix+2;
//    for (int iz=_iza; iz<_izd; ++iz) {
//      int izm1 = iz-1;
//      int izp1 = iz+1;
//      int izm2 = iz-2;
//      int izp2 = iz+2;
//      float r = _r[ix][iz];
//      // FIXME: Source injection replaces += needed to pass adjoint test.
//      up[ix][iz] += (2.0f+C0*r)*ui[ix][iz]-um[ix][iz]+
//        C1*(_r[ix  ][izm1]*
//            ui[ix  ][izm1]+
//            _r[ix  ][izp1]*
//            ui[ix  ][izp1]+
//            _r[ixm1][iz  ]*
//            ui[ixm1][iz  ]+
//            _r[ixp1][iz  ]*
//            ui[ixp1][iz  ])+
//        C2*(_r[ixm1][izm1]*
//            ui[ixm1][izm1]+
//            _r[ixm1][izp1]*
//            ui[ixm1][izp1]+
//            _r[ixp1][izm1]*
//            ui[ixp1][izm1]+
//            _r[ixp1][izp1]*
//            ui[ixp1][izp1])+
//        C3*(_r[ix  ][izm2]*
//            ui[ix  ][izm2]+
//            _r[ix  ][izp2]*
//            ui[ix  ][izp2]+
//            _r[ixm2][iz  ]*
//            ui[ixm2][iz  ]+
//            _r[ixp2][iz  ]*
//            ui[ixp2][iz  ])+
//        C4*(_r[ixm2][izm1]*
//            ui[ixm2][izm1]+
//            _r[ixm2][izp1]*
//            ui[ixm2][izp1]+
//            _r[ixm1][izm2]*
//            ui[ixm1][izm2]+
//            _r[ixm1][izp2]*
//            ui[ixm1][izp2]+
//            _r[ixp1][izm2]*
//            ui[ixp1][izm2]+
//            _r[ixp1][izp2]*
//            ui[ixp1][izp2]+
//            _r[ixp2][izm1]*
//            ui[ixp2][izm1]+
//            _r[ixp2][izp1]*
//            ui[ixp2][izp1]);
//    }
//  }

  // Liu, Y. and M. K. Sen, 2010, A hybrid scheme for absorbing
  // edge reflections in numerical modeling of wave propagation.
  private void absorb(
  final float[][] um, final float[][] ui, final float[][] up) {
    final float oxt = (float)(1.0/(_dx*_dt));
    final float oxx = (float)(1.0/(_dx*_dx));
    final float ott = (float)(1.0/(_dt*_dt));
    //final float[][] cp = copy(up);
    final float[][] cp = up;
    for (int ix=_ixa; ix<_ixb-1; ++ix) {
      for (int iz=_izb; iz<_izc; ++iz) {
        float si = _s[ix][iz];
        float a =  0.50f*oxt;
        float b = -0.50f*ott*si;
        float c =  0.25f*oxx/si;
        float vp = 1.0f/(a-b)*(
          (um[ix][iz-1]+um[ix][iz+1]+cp[ix+1][iz+1]+cp[ix+1][iz-1])*c+
          (um[ix][iz]+cp[ix+1][iz])*(a+b-2.0f*c)+
          (ui[ix][iz]+ui[ix+1][iz])*(-2.0f*b)+
          (um[ix+1][iz])*(b-a)
        );
        float w = _w[ix][iz];
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }
    for (int ix=_ixd-1; ix>=_ixc+1; --ix) {
      for (int iz=_izb; iz<_izc; ++iz) {
        float si = _s[ix][iz];
        float a = -0.50f*oxt;
        float b =  0.50f*ott*si;
        float c = -0.25f*oxx/si;
        float vp = 1.0f/(a-b)*(
          (um[ix][iz-1]+um[ix][iz+1]+cp[ix-1][iz+1]+cp[ix-1][iz-1])*c+
          (um[ix][iz]+cp[ix-1][iz])*(a+b-2.0f*c)+
          (ui[ix][iz]+ui[ix-1][iz])*(-2.0f*b)+
          (um[ix-1][iz])*(b-a)
        );
        float w = _w[ix][iz];
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }
    for (int iz=_iza; iz<_izb-1; ++iz) {
      for (int ix=_ixb; ix<_ixc; ++ix) {
        float si = _s[ix][iz];
        float a =  0.50f*oxt;
        float b = -0.50f*ott*si;
        float c =  0.25f*oxx/si;
        float vp = 1.0f/(a-b)*(
          (um[ix-1][iz]+um[ix+1][iz]+cp[ix+1][iz+1]+cp[ix-1][iz+1])*c+
          (um[ix][iz]+cp[ix][iz+1])*(a+b-2.0f*c)+
          (ui[ix][iz]+ui[ix][iz+1])*(-2.0f*b)+
          (um[ix][iz+1])*(b-a)
        );
        float w = _w[ix][iz];
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }
    for (int iz=_izd-1; iz>=_izc+1; --iz) {
      for (int ix=_ixb; ix<_ixc; ++ix) {
        float si = _s[ix][iz];
        float a = -0.50f*oxt;
        float b =  0.50f*ott*si;
        float c = -0.25f*oxx/si;
        float vp = 1.0f/(a-b)*(
          (um[ix-1][iz]+um[ix+1][iz]+cp[ix+1][iz-1]+cp[ix-1][iz-1])*c+
          (um[ix][iz]+cp[ix][iz-1])*(a+b-2.0f*c)+
          (ui[ix][iz]+ui[ix][iz-1])*(-2.0f*b)+
          (um[ix][iz-1])*(b-a)
        );
        float w = _w[ix][iz];
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }
    // Corners
    for (int ix=_ixa; ix<_ixb; ++ix) {
      for (int iz=_iza; iz<_izb; ++iz) {
        float w = _w[ix][iz];
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[ix][iz]));
        float vp = ui[ix][iz]+r*(ui[ix][iz+1]+ui[ix+1][iz]-2.0f*ui[ix][iz]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }
    for (int ix=_ixa; ix<_ixb; ++ix) {
      for (int iz=_izc; iz<_izd; ++iz) {
        float w = _w[ix][iz];
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[ix][iz]));
        float vp = ui[ix][iz]+r*(ui[ix][iz-1]+ui[ix+1][iz]-2.0f*ui[ix][iz]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }
    for (int ix=_ixc; ix<_ixd; ++ix) {
      for (int iz=_iza; iz<_izb; ++iz) {
        float w = _w[ix][iz];
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[ix][iz]));
        float vp = ui[ix][iz]+r*(ui[ix][iz+1]+ui[ix-1][iz]-2.0f*ui[ix][iz]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }
    for (int ix=_ixc; ix<_ixd; ++ix) {
      for (int iz=_izc; iz<_izd; ++iz) {
        float w = _w[ix][iz];
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[ix][iz]));
        float vp = ui[ix][iz]+r*(ui[ix][iz-1]+ui[ix-1][iz]-2.0f*ui[ix][iz]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }
  }

  private static float[][] extendModel(float[][] c, int b) {
    int nz = c[0].length;
    int nx = c.length;
    float[][] v = new float[nx+2*b][nz+2*b];
    copy(nz,nx,0,0,c,b,b,v);
    for (int ix=b; ix<nx+b; ++ix) {
      for (int iz=0, jz=nz+b; iz<b; ++iz, ++jz) {
        v[ix][iz] = v[ix][b];
        v[ix][jz] = v[ix][nz+b-1];
      }
    }
    for (int ix=0, jx=nx+b; ix<b; ++ix, ++jx) {
      copy(v[b],v[ix]);
      copy(v[nx+b-1],v[jx]);
    }
    return v;
  }

  private void addRandomBoundary(float maxPerturb, float[][] s) {
    float[][] r = new float[_nx][_nz];
    Random random = new Random(012345);
    for (int ix=0; ix<_nx; ++ix) {
      for (int iz=0; iz<_nz; ++iz) {
        float wi = 1.0f-_w[ix][iz];
        wi *= wi;
        if (wi>0.0f) {
          s[ix][iz] += wi*maxPerturb*(1.1f*random.nextFloat()-1.0f);
        }
      }
    }
    SimplePlot.asPixels(s);
  }

  private float[][] makeWeights() {
    int fd = FD_ORDER/2;
    float ob = 1.0f/(_b-1.0f);
    float[][] w = new float[_nx][_nz];
    for (int ix=_ixa; ix<_ixb; ++ix) {
      for (int iz=ix; iz<_nz-ix; ++iz) {
        w[ix      ][iz] = (ix-fd)*ob;
        w[_nx-1-ix][iz] = (ix-fd)*ob;
      }
    }
    for (int iz=_iza; iz<_izb; ++iz) {
      for (int ix=iz; ix<_nx-iz; ++ix) {
        w[ix][iz      ] = (iz-fd)*ob;
        w[ix][_nz-1-iz] = (iz-fd)*ob;
      }
    }
    for (int ix=_ixb; ix<_ixc; ++ix) {
      for (int iz=_izb; iz<_izc; ++iz) {
        w[ix][iz] = 1.0f;
      }
    }
    //SimplePlot.asPixels(w).addColorBar();
    return w;
  }
}
