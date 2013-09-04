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

public class AcousticWaveOperator {

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

  public AcousticWaveOperator(
  float[][] s, double dx, double dt, int nabsorb) {
    int nx = s[0].length;
    int nz = s.length;
    Check.argument(nabsorb>=FD_ORDER/2,"nabsorb>=FD_ORDER/2");
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
    float scale = (_dt*_dt/(_dx*_dx));
    _s = extendModel(s,nabsorb);
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
  // static

  public static float[][] crop(float[][] x, int nabsorb) {
    int nx = x[0].length;
    int nz = x.length;
    return copy(nx-2*nabsorb,nz-2*nabsorb,nabsorb,nabsorb,x);
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

  public static interface Source {
    public void add(float[][] ui, int it, int nabsorb);
  }

  public static class RickerSource implements Source {
    public RickerSource(int xs, int zs, float dt, float fpeak) {
      _tdelay = 1.0f/fpeak;
      _fpeak = fpeak;
      _dt = dt;
      _xs = xs;
      _zs = zs;
    }
    public void add(float[][] ui, int it, int nabsorb) {
      ui[_zs+nabsorb][_xs+nabsorb] += ricker(it*_dt-_tdelay);
    }
    private float ricker(float t) {
      float x = FLT_PI*_fpeak*t;
      float xx = x*x;
      return (float)((1.0-2.0*xx)*exp(-xx));
    }
    private int _xs,_zs;
    private float _dt,_fpeak,_tdelay;
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
        ui[zs][xs] += _data[ir][it];
      }
    }
    private int _nr,_nt;
    private int[] _xr, _zr;
    private float[][] _data;
  }

  public static class WavefieldSource implements Source {
    public WavefieldSource(float[][][] u) {
      _nx = u[0][0].length;
      _nz = u[0].length;
      _u = u;
    }
    public WavefieldSource(float[][][] u, float[][] r) {
      _nx = u[0][0].length;
      _nz = u[0].length;
      _u = u;
      _r = r;
    }
    public void add(float[][] ui, int it, int nabsorb) {
      if (_r==null) {
        for (int ix=0; ix<_nx; ++ix)
          for (int iz=0; iz<_nz; ++iz)
            ui[iz][ix] += _u[it][iz][ix];
      } else {
        for (int ix=nabsorb; ix<_nx; ++ix)
          for (int iz=0; iz<_nz; ++iz)
            ui[iz][ix] += _u[it][iz][ix]*_r[iz-nabsorb][ix-nabsorb];
      }
    }
    private float[][] _r = null; //reflectivity
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
        _data[ir][it] = ui[zr][xr];
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
  private float _dx,_dt;
  private float[][] _s; // slowness
  private float[][] _r; // dt^2/(dx^2*s^2)
  private float[][] _w; // weights for boundary
  private int _ixa,_ixb,_ixc,_ixd;
  private int _iza,_izb,_izc,_izd;

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
    for (int it=fit, count=0; count<nt-2; it+=pit, ++count) {

      // Set next time step.
      if (u==null) {
        up = ut;
        zero(ut); // necessary if using += for injection
      } else {
        up = u[it];
      }

      // Step.
      if (forward) {
        forwardStep(um,ui,up);
        if (source!=null)
          source.add(up,it,_nabsorb);
        forwardAbsorb(um,ui,up);
      } else {
        adjointAbsorb(up,ui,um);
        if (source!=null)
          source.add(up,it,_nabsorb);
        adjointStep(up,ui,um);
        //forwardStep(um,ui,up);
        //if (source!=null)
        //  source.add(up,it,_nabsorb);
        //forwardAbsorb(um,ui,up);
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

  // 21-point Laplacian stencil from
  // Patra, M. and M. Karttunen, 2005, Stencils
  // with Isotropic Error for Differential Operators.
  private static final int FD_ORDER = 4;
  private static final float C0 = -9.0f/2.0f;
  private static final float C1 =  16.0f/15.0f;
  private static final float C2 =  2.0f/15.0f;
  private static final float C3 = -1.0f/15.0f;
  private static final float C4 = -1.0f/120.0f;

  private void forwardStep(
  final float[][] um, final float[][] ui, final float[][] up) {
    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
    public void compute(int iz) {
      int izm1 = iz-1;
      int izp1 = iz+1;
      int izm2 = iz-2;
      int izp2 = iz+2;
      for (int ix=_ixa; ix<_ixd; ++ix) {
        int ixm1 = ix-1;
        int ixp1 = ix+1;
        int ixm2 = ix-2;
        int ixp2 = ix+2;
        float r = _r[iz][ix];
        up[iz][ix] += (2.0f+C0*r)*ui[iz][ix]-um[iz][ix]+r*(
          C1*(ui[izm1][ix  ]+ui[izp1][ix  ]+ui[iz  ][ixm1]+ui[iz  ][ixp1])+
          C2*(ui[izm1][ixm1]+ui[izp1][ixm1]+ui[izm1][ixp1]+ui[izp1][ixp1])+
          C3*(ui[izm2][ix  ]+ui[izp2][ix  ]+ui[iz  ][ixm2]+ui[iz  ][ixp2])+
          C4*(ui[izm2][ixm2]+ui[izp2][ixm2]+ui[izm2][ixp2]+ui[izp2][ixp2]));
        //up[iz][ix] -= um[iz  ][ix  ];
        //up[iz][ix] += ui[iz  ][ix  ]*(2.0f+C0*r);
        //up[iz][ix] += ui[izm1][ix  ]*C1*r;
        //up[iz][ix] += ui[izp1][ix  ]*C1*r;
        //up[iz][ix] += ui[iz  ][ixm1]*C1*r;
        //up[iz][ix] += ui[iz  ][ixp1]*C1*r;
        //up[iz][ix] += ui[izm1][ixm1]*C2*r;
        //up[iz][ix] += ui[izp1][ixm1]*C2*r;
        //up[iz][ix] += ui[izm1][ixp1]*C2*r;
        //up[iz][ix] += ui[izp1][ixp1]*C2*r;
        //up[iz][ix] += ui[izm2][ix  ]*C3*r;
        //up[iz][ix] += ui[izp2][ix  ]*C3*r;
        //up[iz][ix] += ui[iz  ][ixm2]*C3*r;
        //up[iz][ix] += ui[iz  ][ixp2]*C3*r;
        //up[iz][ix] += ui[izm2][ixm2]*C4*r;
        //up[iz][ix] += ui[izp2][ixm2]*C4*r;
        //up[iz][ix] += ui[izm2][ixp2]*C4*r;
        //up[iz][ix] += ui[izp2][ixp2]*C4*r;
      }
    }});
  }

  private void adjointStep(
  final float[][] um, final float[][] ui, final float[][] up) {

    //for (int iz=_iza; iz<_izd; ++iz) {
    //  int izm1 = iz-1;
    //  int izp1 = iz+1;
    //  int izm2 = iz-2;
    //  int izp2 = iz+2;
    //  for (int ix=_ixa; ix<_ixd; ++ix) {
    //    int ixm1 = ix-1;
    //    int ixp1 = ix+1;
    //    int ixm2 = ix-2;
    //    int ixp2 = ix+2;
    //    float r = _r[iz][ix];
    //    um[iz  ][ix  ] -= up[iz][ix];
    //    ui[iz  ][ix  ] += up[iz][ix]*(2.0f+C0*r);
    //    ui[izm1][ix  ] += up[iz][ix]*C1*r;
    //    ui[izp1][ix  ] += up[iz][ix]*C1*r;
    //    ui[iz  ][ixm1] += up[iz][ix]*C1*r;
    //    ui[iz  ][ixp1] += up[iz][ix]*C1*r;
    //    ui[izm1][ixm1] += up[iz][ix]*C2*r;
    //    ui[izp1][ixm1] += up[iz][ix]*C2*r;
    //    ui[izm1][ixp1] += up[iz][ix]*C2*r;
    //    ui[izp1][ixp1] += up[iz][ix]*C2*r;
    //    ui[izm2][ix  ] += up[iz][ix]*C3*r;
    //    ui[izp2][ix  ] += up[iz][ix]*C3*r;
    //    ui[iz  ][ixm2] += up[iz][ix]*C3*r;
    //    ui[iz  ][ixp2] += up[iz][ix]*C3*r;
    //    ui[izm2][ixm2] += up[iz][ix]*C4*r;
    //    ui[izp2][ixm2] += up[iz][ix]*C4*r;
    //    ui[izm2][ixp2] += up[iz][ix]*C4*r;
    //    ui[izp2][ixp2] += up[iz][ix]*C4*r;
    //  }
    //}

    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
    public void compute(int iz) {
      for (int ix=_ixa; ix<_ixd; ++ix) {
        float r = _r[iz][ix];
        ui[iz-2][ix  ] += up[iz][ix]*C3*r;
        ui[iz-2][ix-2] += up[iz][ix]*C4*r;
        ui[iz-2][ix+2] += up[iz][ix]*C4*r;
      }
    }});
    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
    public void compute(int iz) {
      for (int ix=_ixa; ix<_ixd; ++ix) {
        float r = _r[iz][ix];
        ui[iz-1][ix  ] += up[iz][ix]*C1*r;
        ui[iz-1][ix-1] += up[iz][ix]*C2*r;
        ui[iz-1][ix+1] += up[iz][ix]*C2*r;
      }
    }});
    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
    public void compute(int iz) {
      for (int ix=_ixa; ix<_ixd; ++ix) {
        float r = _r[iz][ix];
        um[iz  ][ix  ] -= up[iz][ix];
        ui[iz  ][ix  ] += up[iz][ix]*(2.0f+C0*r);
        ui[iz  ][ix-1] += up[iz][ix]*C1*r;
        ui[iz  ][ix+1] += up[iz][ix]*C1*r;
        ui[iz  ][ix-2] += up[iz][ix]*C3*r;
        ui[iz  ][ix+2] += up[iz][ix]*C3*r;
      }
    }});
    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
    public void compute(int iz) {
      for (int ix=_ixa; ix<_ixd; ++ix) {
        float r = _r[iz][ix];
        ui[iz+1][ix  ] += up[iz][ix]*C1*r;
        ui[iz+1][ix-1] += up[iz][ix]*C2*r;
        ui[iz+1][ix+1] += up[iz][ix]*C2*r;
      }
    }});
    Parallel.loop(_iza,_izd,new Parallel.LoopInt() {
    public void compute(int iz) {
      for (int ix=_ixa; ix<_ixd; ++ix) {
        float r = _r[iz][ix];
        ui[iz+2][ix  ] += up[iz][ix]*C3*r;
        ui[iz+2][ix-2] += up[iz][ix]*C4*r;
        ui[iz+2][ix+2] += up[iz][ix]*C4*r;
      }
    }});

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
//        um[iz][ix] -= up[iz  ][ix  ];
//        ui[iz][ix] += up[iz  ][ix  ]*(2.0f+C0*r);
//        ui[iz][ix] += up[izp1][ix  ]*C1*r;
//        ui[iz][ix] += up[izm1][ix  ]*C1*r;
//        ui[iz][ix] += up[iz  ][ixp1]*C1*r;
//        ui[iz][ix] += up[iz  ][ixm1]*C1*r;
//        ui[iz][ix] += up[izp1][ixp1]*C2*r;
//        ui[iz][ix] += up[izm1][ixp1]*C2*r;
//        ui[iz][ix] += up[izp1][ixm1]*C2*r;
//        ui[iz][ix] += up[izm1][ixm1]*C2*r;
//        ui[iz][ix] += up[izp2][ix  ]*C3*r;
//        ui[iz][ix] += up[izm2][ix  ]*C3*r;
//        ui[iz][ix] += up[iz  ][ixp2]*C3*r;
//        ui[iz][ix] += up[iz  ][ixm2]*C3*r;
//        ui[iz][ix] += up[izp2][ixp2]*C4*r;
//        ui[iz][ix] += up[izm2][ixp2]*C4*r;
//        ui[iz][ix] += up[izp2][ixm2]*C4*r;
//        ui[iz][ix] += up[izm2][ixm2]*C4*r;
//    }});
//    for (int iz=_iza; iz<_izd; ++iz) {
////      ui[iz  ][ix+1  ] += up[iz][ix    ]*C1*r;
////      ui[iz  ][_ixd  ] += up[iz][_ixd-1]*C1*_r[iz][_ixd-1];
////      ui[iz  ][_ixa  ] -= up[iz][_ixa-1]*C1*_r[iz][_ixa-1];
//
////      ui[iz  ][ix-1  ] += up[iz][ix  ]*C1*r;
////      ui[iz  ][_ixa-1] += up[iz][_ixa]*C1*_r[iz][_ixa];
////      ui[iz  ][_ixd-1] -= up[iz][_ixd]*C1*_r[iz][_ixd];
//    }

  }

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
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[iz][ix]));
        up[iz][ix] += -m*up[iz][ix]+m*(
          ui[iz][ix]+r*(ui[iz+1][ix]+ui[iz][ix+1]-2.0f*ui[iz][ix])
        );
      }
    }
    for (int iz=_izc; iz<_izd; ++iz) {
      for (int ix=_ixa; ix<_ixb; ++ix) {
        float m = 1.0f-_w[iz][ix];
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[iz][ix]));
        up[iz][ix] += -m*up[iz][ix]+m*(
          ui[iz][ix]+r*(ui[iz-1][ix]+ui[iz][ix+1]-2.0f*ui[iz][ix])
        );
      }
    }
    for (int iz=_iza; iz<_izb; ++iz) {
      for (int ix=_ixc; ix<_ixd; ++ix) {
        float m = 1.0f-_w[iz][ix];
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[iz][ix]));
        up[iz][ix] += -m*up[iz][ix]+m*(
          ui[iz][ix]+r*(ui[iz+1][ix]+ui[iz][ix-1]-2.0f*ui[iz][ix])
        );
      }
    }
    for (int iz=_izc; iz<_izd; ++iz) {
      for (int ix=_ixc; ix<_ixd; ++ix) {
        float m = 1.0f-_w[iz][ix];
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[iz][ix]));
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
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[iz][ix]));
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
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[iz][ix]));
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
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[iz][ix]));
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
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[iz][ix]));
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
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[iz][ix]));
        float vp = ui[iz][ix]+r*(ui[iz+1][ix]+ui[iz][ix+1]-2.0f*ui[iz][ix]);
        up[iz][ix] = w*up[iz][ix]+(1.0f-w)*vp;
      }
    }
    for (int iz=_izc; iz<_izd; ++iz) {
      for (int ix=_ixa; ix<_ixb; ++ix) {
        float w = _w[iz][ix];
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[iz][ix]));
        float vp = ui[iz][ix]+r*(ui[iz-1][ix]+ui[iz][ix+1]-2.0f*ui[iz][ix]);
        up[iz][ix] = w*up[iz][ix]+(1.0f-w)*vp;
      }
    }
    for (int iz=_iza; iz<_izb; ++iz) {
      for (int ix=_ixc; ix<_ixd; ++ix) {
        float w = _w[iz][ix];
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[iz][ix]));
        float vp = ui[iz][ix]+r*(ui[iz+1][ix]+ui[iz][ix-1]-2.0f*ui[iz][ix]);
        up[iz][ix] = w*up[iz][ix]+(1.0f-w)*vp;
      }
    }
    for (int iz=_izc; iz<_izd; ++iz) {
      for (int ix=_ixc; ix<_ixd; ++ix) {
        float w = _w[iz][ix];
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[iz][ix]));
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
