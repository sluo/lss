package lss.vel;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

// testing
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;

public class AcousticWaveOperator {

  public AcousticWaveOperator(
  float[][] s, double dx, double dt) {
    int nabsorb = NABSORB;
    int nz = s[0].length;
    int nx = s.length;
    Check.argument(nabsorb>=FD_ORDER/2,"nabsorb>=FD_ORDER/2");
    _b = nabsorb-FD_ORDER/2;
    _nz = nz+2*nabsorb;
    _nx = nx+2*nabsorb;
    _dx = dx;
    _dt = dt;
    _fz = -_b*_dx;
    _fx = -_b*_dx;
    _sz = new Sampling(_nz,_dx,_fz);
    _sx = new Sampling(_nx,_dx,_fx);
    _ixa = FD_ORDER/2;
    _ixb = _ixa+_b;
    _ixc = _ixb+nx;
    _ixd = _ixc+_b;
    _iza = FD_ORDER/2;
    _izb = _iza+_b;
    _izc = _izb+nz;
    _izd = _izc+_b;
    float scale = (float)(_dt*_dt/(_dx*_dx));
    _s = extendModel(s,nabsorb);
    _r = copy(_s);
    mul(_r,_r,_r);
    div(scale,_r,_r);
  }

  /* Indexing example for FD_ABSORB/2=2, nabsorb=3, sx.count=3:

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

  private void apply(
  boolean forward, Source source, Receiver receiver, float[][][] u) {
    int nz = _nz-2*NABSORB;
    int nx = _nx-2*NABSORB;
    int nt = (receiver==null)?u.length:receiver.getNt();
    float[][] um = new float[_nx][_nz];
    float[][] ui = new float[_nx][_nz];
    //float[][] up = (u==null)?new float[_nx][_nz]:u[0];
    float[][] up = new float[_nx][_nz];
    int itf = (forward)?0:nt-1;
    int itp = (forward)?1:-1;
    for (int it=itf,count=0; count<nt; it+=itp,++count) {
      
      // Time step.
      forwardStep(um,ui,up);

      // Source injection.
      source.add(_dt*it,up,_sx,_sz);

      // Absorbing boundaries.
      absorb(um,ui,up);

      // Copy data and wavefield.
      if (receiver!=null)
        receiver.setData(it,up,_sx,_sz);
      if (u!=null)
        copy(nz,nx,_b,_b,up,0,0,u[it]);

      // Rotate arrays
      float[][] ut = um;
      um = ui;
      ui = up;
      //up = (u!=null&&it<nt-1)?u[it+1]:ut;
      up = ut;
    }
  }

  public static class Receiver {
    public Receiver(float xr, float zr, int nt) {
      this(new float[]{xr},new float[]{zr},nt);
    }
    public Receiver(float[] xr, float[] zr, int nt) {
      int nxr = xr.length;
      int nzr = zr.length;
      Check.argument(nxr==nzr,"nxr==nzr");
      _nr = nxr;
      _nt = nt;
      _xr = xr;
      _zr = zr;
      _data = new float[_nr][nt];
    }
    public int getNt() {
      return _nt;
    }
    public void setData(int it, float[][] ui, Sampling sx, Sampling sz) {
      double dz = sz.getDelta();
      double dx = sx.getDelta();
      double fz = sz.getFirst();
      double fx = sx.getFirst();
      for (int ir=0; ir<_nr; ++ir) {
        int kxr = (int)round((_xr[ir]-fx)/dx);
        int kzr = (int)round((_zr[ir]-fz)/dz);
        _data[ir][it] = ui[kxr][kzr];
      }
    }
    public float[][] getData() {
      return _data;
    }
    public int _nr,_nt;
    public float[] _xr, _zr;
    public float[][] _data;
  }

  public static interface Source {
    public void add(double t, float[][] f, Sampling sx, Sampling sz);
  }
  public static class RickerSource implements Source {
    public RickerSource(double fpeak, float x, float z) {
      _tdelay = 1.0/fpeak;
      _fpeak = fpeak;
      _x = x;
      _z = z;
    }
    public void add(double t, float[][] ui, Sampling sx, Sampling sz) {
      double dz = sz.getDelta();
      double dx = sx.getDelta();
      double fz = sz.getFirst();
      double fx = sx.getFirst();
      int kxs = (int)round((_x-fx)/dx);
      int kzs = (int)round((_z-fz)/dz);
      ui[kxs][kzs] += ricker(t-_tdelay);
    }
    private float ricker(double t) {
      double x = PI*_fpeak*t;
      double xx = x*x;
      return (float)((1.0-2.0*xx)*exp(-xx));
    }
    private int _kzs,_kxs;
    private float _x,_z;
    private double _fpeak,_tdelay;
  }

  private void forwardStep(
  final float[][] um, final float[][] ui, final float[][] up) {
    Parallel.loop(_ixa,_ixd,new Parallel.LoopInt() {
    public void compute(int ix) {
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
        up[ix][iz] = (2.0f+C5*r)*ui[ix][iz]-um[ix][iz]+r*(
          C4*(ui[ix  ][izm1]+ui[ix  ][izp1]+ui[ixm1][iz  ]+ui[ixp1][iz  ])+
          C3*(ui[ixm1][izm1]+ui[ixm1][izp1]+ui[ixp1][izm1]+ui[ixp1][izp1])+
          C2*(ui[ix  ][izm2]+ui[ix  ][izp2]+ui[ixm2][iz  ]+ui[ixp2][iz  ])+
          C1*(ui[ixm2][izm1]+ui[ixm2][izp1]+ui[ixm1][izm2]+ui[ixm1][izp2]+
              ui[ixp1][izm2]+ui[ixp1][izp2]+ui[ixp2][izm1]+ui[ixp2][izp1]));
      }
    }});
  }

  // Coefficients for 21-point Laplacian from
  // Patra, M. and M. Karttunen, 2005, Stencils
  // with Isotropic Error for Differential Operators.
  private static final float C1 = -1.0f/30.0f;
  private static final float C2 = -1.0f/60.0f;
  private static final float C3 =  4.0f/15.0f;
  private static final float C4 =  13.0f/15.0f;
  private static final float C5 = -21.0f/5.0f;

  // Liu, Y. and M. K. Sen, 2010, A hybrid scheme for absorbing
  // edge reflections in numerical modeling of wave propagation.
  private void absorb(
  final float[][] um, final float[][] ui, final float[][] up) {

    final float oxt = (float)(1.0/(_dx*_dt));
    final float oxx = (float)(1.0/(_dx*_dx));
    final float ott = (float)(1.0/(_dt*_dt));
    final float[][] cp = copy(up);

    final float fd = (float)(FD_ORDER/2);
    final float ob = 1.0f/_b;

    Parallel.loop(_ixa,_ixb,new Parallel.LoopInt() {
    public void compute(int ix) {
      float w = (ix-fd)*ob;
      for (int iz=ix+1; iz<_nz-1-ix; ++iz) {
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
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }});

    Parallel.loop(_ixc,_ixd,new Parallel.LoopInt() {
    public void compute(int ix) {
      float w = ((_nx-1)-(ix+fd))*ob;
      for (int iz=_nx-ix; iz<_nz-1-(_nx-1-ix); ++iz) {
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
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }});

    Parallel.loop(_iza,_izb,new Parallel.LoopInt() {
    public void compute(int iz) {
      float w = (iz-fd)*ob;
      for (int ix=iz+1; ix<_nx-1-iz; ++ix) {
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
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }});

    Parallel.loop(_izc,_izd,new Parallel.LoopInt() {
    public void compute(int iz) {
      float w = ((_nz-1)-(iz+fd))*ob;
      for (int ix=_nz-iz; ix<_nx-1-(_nz-1-iz); ++ix) {
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
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }});

    // Corners
    for (int ix=_ixa; ix<_ixb; ++ix) {
      float w = (ix-fd)*ob;
      for (int iz=_iza; iz<_izb; ++iz) {
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[ix][iz]));
        float vp = ui[ix][iz]+r*(ui[ix][iz+1]+ui[ix+1][iz]-2.0f*ui[ix][iz]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    //}
    //for (int ix=_ixa; ix<_ixb; ++ix) {
    //  float w = (float)(ix-_ixa)/(float)_b;
      for (int iz=_izc; iz<_izd; ++iz) {
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[ix][iz]));
        float vp = ui[ix][iz]+r*(ui[ix][iz-1]+ui[ix+1][iz]-2.0f*ui[ix][iz]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }
    for (int ix=_ixc; ix<_ixd; ++ix) {
      float w = ((_nx-1)-(ix+fd))*ob;
      for (int iz=_iza; iz<_izb; ++iz) {
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[ix][iz]));
        float vp = ui[ix][iz]+r*(ui[ix][iz+1]+ui[ix-1][iz]-2.0f*ui[ix][iz]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    //}
    //for (int ix=_ixc; ix<_ixd; ++ix) {
    //  float w = (float)((_nx-1)-ix)/(float)_b;
      for (int iz=_izc; iz<_izd; ++iz) {
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[ix][iz]));
        float vp = ui[ix][iz]+r*(ui[ix][iz-1]+ui[ix-1][iz]-2.0f*ui[ix][iz]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }

  }

  //////////////////////////////////////////////////////////////////////////

  private static final int NABSORB = 12; // size of absorbing boundary
  private static final int FD_ORDER = 4; // finite-difference order

  private int _b; // absorbing boundary size
  //private Source _source;
  //private Receiver _receiver;
  private double _t; // current time
  private int _it; // current time index
  private int _nx,_nz,_nt;
  private double _dx,_dz,_dt;
  private double _fx,_fz,_ft;
  private Sampling _sz,_sx;
  private float[][] _s; // slowness
  private float[][] _r; // dt^2/(dx^2*s^2)
  private int _ixa,_ixb,_ixc,_ixd;
  private int _iza,_izb,_izc,_izd;

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
}
