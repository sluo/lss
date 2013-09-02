package lss.vel;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

// testing
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;

public class Wavefield {

  public Wavefield(Sampling sz, Sampling sx, Sampling st) {
    _dz = sz.getDelta();
    _dx = sx.getDelta();
    _dt = st.getDelta();
    Check.argument(_dx==_dz,"dx==dz");
    _nz = sz.getCount()+2*_b;
    _nx = sx.getCount()+2*_b;
    _nt = st.getCount();
  }

  public void modelAcousticData(
    Source source, Receiver receiver, float[][] s, float[][] d)
  {
    modelAcousticDataAndWavefield(source,receiver,s,d,null);
  }

  public void modelAcousticWavefield(
    Source source, float[][] s, float[][][] u)
  {
    modelAcousticDataAndWavefield(source,null,s,null,u);
  }

  public void modelAcousticDataAndWavefield(
    Source source, Receiver receiver,
    float[][] s, float[][] d, float[][][] u)
  {
    updateParameters(source,receiver,s,null);
    modelAcoustic(d,u);
  }
  
  public static interface Source {
    public void add(double t, float[][] f);
  }

  public static class Receiver {
    public Receiver(int kzr, int kxr) {
      this(new int[]{kzr},new int[]{kxr});
    }
    public Receiver(int[] kzr, int[] kxr) {
      this.kzr = kzr;
      this.kxr = kxr;
    }
    public int[] kzr,kxr;
  }

  //////////////////////////////////////////////////////////////////////////

  private static int _b = 20; // absorbing boundary size

  private Source _source;
  private Receiver _receiver;
  private double _t; // current time
  private int _it; // current time index
  private int _nx,_nz,_nt; // number of samples
  private double _dx,_dz,_dt; // sampling intervals
  private float[][] _s; // slowness
  private float[][] _r; // dt^2/(dx^2*s^2)

  private void updateParameters(
    Source source, Receiver receiver, float[][] s, float[][] s1)
  {
    _source = source;
    _receiver = receiver;
    float r = (float)(_dt*_dt/(_dx*_dx));
    _s = extendModel(s);
    _r = copy(_s);
    mul(_r,_r,_r);
    div(  r,_r,_r);
  }

  private void modelAcoustic(float[][] d, float[][][] u) {
    int nz = _nz-2*_b;
    int nx = _nx-2*_b;
    int nt = _nt;
    float[][] um = new float[_nx][_nz];
    float[][] ui = new float[_nx][_nz];
    float[][] up = new float[_nx][_nz];
    for (int it=0; it<nt; ++it) {
      
      // One time step
      step(um,ui,up);

      // Source injection
      _source.add(_dt*it,up);

      // Absorbing boundaries
      absorb(um,ui,up);

      // Copy data and wavefield
      if (_receiver!=null) {
        int nr = _receiver.kzr.length;
        for (int ir=0; ir<nr; ++ir) {
          int kz = _receiver.kzr[ir]+_b;
          int kx = _receiver.kxr[ir]+_b;
          d[ir][it] = ui[kx][kz];
        }
      }
      if (u!=null) {
        copy(nz,nx,_b,_b,ui,0,0,u[it]);
      }

      // Rotate arrays
      float[][] t = um;
      um = ui;
      ui = up;
      up = t;
    }

    // Time reverse if necessary
    if (_source instanceof AdjointSource) {
      reverse3(u);
    }
  }

  public static class RickerSource implements Source {
    public RickerSource(double fpeak, int kzs, int kxs) {
      _fpeak = fpeak;
      _kzs = kzs+_b;
      _kxs = kxs+_b;
      _tdelay = 1.0/fpeak;
    }
    public void add(double t, float[][] f) {
      f[_kxs][_kzs] += ricker(t-_tdelay);
    }
    private float ricker(double t) {
      double x = PI*_fpeak*t;
      double xx = x*x;
      return (float)((1.0-2.0*xx)*exp(-xx));
    }
    private int _kzs,_kxs;
    private double _fpeak,_tdelay;
  }

  private static abstract class TimeFunctionSource implements Source {
    public TimeFunctionSource(double fpeak, int kzs, int kxs) {
      _fpeak = fpeak;
      _kzs = kzs+_b;
      _kxs = kxs+_b;
      _tdelay = 1.0/fpeak;
    }
    public void add(double t, float[][] f) {
      f[_kxs][_kzs] += function(t-_tdelay);
    }
    public abstract float function(double t);
    double _fpeak,_tdelay;
    int _kzs,_kxs;
  }
  public static class GaussianSource extends TimeFunctionSource {
    public GaussianSource(double fpeak, int kzs, int kxs) {
      super(fpeak,kzs,kxs);
    }
    public float function(double t) {
      double x = PI*_fpeak*t;
      double xx = x*x;
      return (float)exp(-xx);
    }
  }
  public static class Gaussian2Source extends TimeFunctionSource {
    public Gaussian2Source(double fpeak, int kzs, int kxs) {
      super(fpeak,kzs,kxs);
    }
    public float function(double t) {
      double x = PI*_fpeak*t;
      double xx = x*x;
      return (float)((1.0-2.0*xx)*exp(-xx));
    }
  }
  public static class Gaussian4Source extends TimeFunctionSource {
    public Gaussian4Source(double fpeak, int kzs, int kxs) {
      super(fpeak,kzs,kxs);
    }
    public float function(double t) {
      double a = PI*_fpeak;
      double aa = a*a;
      double x = a*t;
      double xx = x*x;
      //return (float)(2.0*aa*(3.0-12.0*xx+4.0*xx*xx)*exp(-xx));
      return (float)((3.0-12.0*xx+4.0*xx*xx)*exp(-xx));
    }
  }

  // 2nd derivative of Ricker wavelet
  public static class Ricker2Source implements Source {
    public Ricker2Source(double fpeak, int kzs, int kxs) {
      _fpeak = fpeak;
      _kzs = kzs+_b;
      _kxs = kxs+_b;
      _tdelay = 1.0/fpeak;
    }
    public void add(double t, float[][] f) {
      f[_kxs][_kzs] += ricker(t-_tdelay);
    }
    private float ricker(double t) {
      double a = PI*_fpeak;
      double aa = a*a;
      double x = a*t;
      double xx = x*x;
      //return (float)(2.0*aa*(3.0-12.0*xx+4.0*xx*xx)*exp(-xx));
      return (float)((3.0-12.0*xx+4.0*xx*xx)*exp(-xx));
    }
    private int _kzs,_kxs;
    private double _fpeak,_tdelay;
  }

  public static class PlaneWaveSource implements Source {
    public PlaneWaveSource(
      double angle, double tdelay, double fpeak,
      double dx, int kzs, int kxs, float[][] v)
    {
      _fpeak = fpeak;
      _angle = angle;
      _kzs = kzs;
      _kxs = kxs; // where is plane-wave centered?
      _tdelay = tdelay;
      _nx = v.length;
      _dx = dx;
      _ps = sin(DBL_PI*angle/180.0)/v[kxs][kzs]; // TODO: function of x?
      _w = fillfloat(1.0f,_nx);
      //new RecursiveGaussianFilter(_nx/4).apply0(_w,_w); // taper
      div(_w,max(_w),_w);
    }
    public void add(double t, float[][] f) {
      //double p = sin(_angle)/_v[_kzs][_nx/2]; // TODO: function of x?
      for (int ix=0; ix<_nx; ++ix) {
        double t0 = _ps*(ix-_kxs)*_dx;
        // TODO + or - t0?
        f[ix+_b][_kzs+_b] += _w[ix]*ricker(t-_tdelay+t0);
      }
    }
    private float ricker(double t) {
      double x = PI*_fpeak*t;
      double xx = x*x;
      return (float)((1.0-2.0*xx)*exp(-xx));
    }
    private int _nx,_kzs,_kxs;
    private double _fpeak,_angle,_tdelay,_ps,_dx;
    private float[] _w; // weight for tapering sources
  }

  public static class DefaultSource implements Source {
    public DefaultSource(double dt, int kzs, int kxs, float[] s) {
      this(dt,new int[]{kzs},new int[]{kxs},new float[][]{s});
    }
    public DefaultSource(double dt, int[] kzs, int[] kxs, float[][] s) {
      Check.argument(kzs.length==kxs.length,"kzs.length=kxs.length");
      _ns = kzs.length;
      _dt = dt;
      _kzs = kzs;
      _kxs = kxs;
      _s = s;
    }
    public void add(double t, float[][] f) {
      int it = (int)(t/_dt);
      for (int is=0; is<_ns; is++) {
        int kzs = _b+_kzs[is];
        int kxs = _b+_kxs[is];
        f[kxs][kzs] += _s[is][it];
      }
    }
    private double _dt;
    private int _ns;
    private int[] _kzs = null;
    private int[] _kxs = null;
    private float[][] _s = null;
  }

  public static class WavefieldSource implements Source {
    public WavefieldSource(double dt, float[][] r, float[][][] u) {
      int nz = u[0][0].length;
      int nx = u[0].length;
      int nt = u.length;
      _dt = dt;
      _r = r;
      _u = u;
    }
    public void add(final double t, final float[][] f) {
      final float[][] uf = _u[(int)(t/_dt)];
      final float[][] rf = _r;
      final int nz = uf[0].length;
      final int nx = uf.length;
      Parallel.loop(nx,new Parallel.LoopInt() {
      public void compute(int ix) {
        for (int iz=0; iz<nz; ++iz) {
          int kzs = _b+iz;
          int kxs = _b+ix;
          f[kxs][kzs] += uf[ix][iz]*rf[ix][iz];
        }
      }});
    }
    private double _dt;
    private float[][] _r;
    private float[][][] _u;
  }

  public static class AdjointSource implements Source {
    public AdjointSource(double dt, int kzs, int kxs, float[] s) {
      this(dt,new int[]{kzs},new int[]{kxs},new float[][]{s});
    }
    public AdjointSource(double dt, int[] kzs, int[] kxs, float[][] s) {
      _source = new DefaultSource(dt,kzs,kxs,reverse(s));
    }
    public void add(double t, float[][] f) {
      _source.add(t,f);
    }
    private static float[][] reverse(float[][] x) {
      int n1 = x[0].length;
      int n2 = x.length;
      float[][] y = new float[n2][n1];
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0,j1=n1-1; i1<n1; ++i1,--j1)
          y[i2][j1] = x[i2][i1];
      return y;
    }
    private DefaultSource _source;
  }

  // Coefficients for 21-point Laplacian
  private static final float C1 =  0.0f;
  private static final float C2 = -1.0f/30.0f;
  private static final float C3 = -1.0f/60.0f;
  private static final float C4 =  4.0f/15.0f;
  private static final float C5 =  13.0f/15.0f;
  private static final float C6 = -21.0f/5.0f;

  // Coefficients for 5-point Laplacian
  private static final float D1 =  1.0f/6.0f;
  private static final float D2 =  2.0f/3.0f;
  private static final float D3 = -10.0f/3.0f;

  private void step(
    final float[][] um, final float[][] ui, final float[][] up)
  {

    Parallel.loop(1,_nx-1,new Parallel.LoopInt() {
    public void compute(int ix) {

      // ix = 1 or ix = nx-2
      if (ix==1 || ix==_nx-2) {
        int ixm1 = ix-1;
        int ixp1 = ix+1;
        for (int iz=1; iz<_nz-1; ++iz) {
          int izm1 = iz-1;
          int izp1 = iz+1;
          float r = _r[ix][iz];
          up[ix][iz] = (2.0f+D3*r)*ui[ix][iz]-um[ix][iz]+r*(
            D2*(ui[ixm1][iz  ]+ui[ix  ][izm1]+ui[ix  ][izp1]+ui[ixp1][iz  ])+
            D1*(ui[ixm1][izm1]+ui[ixm1][izp1]+ui[ixp1][izm1]+ui[ixp1][izp1]));
        }

      // 1 < ix < nx-2
      } else {
        int ixm1 = ix-1;
        int ixp1 = ix+1;
        int ixm2 = ix-2;
        int ixp2 = ix+2;

        // iz = 1
        {
          int iz = 1;
          int izm1 = iz-1;
          int izp1 = iz+1;
          float r = _r[ix][iz];
          up[ix][iz] = (2.0f+D3*r)*ui[ix][iz]-um[ix][iz]+r*(
            D2*(ui[ixm1][iz  ]+ui[ix  ][izm1]+ui[ix  ][izp1]+ui[ixp1][iz  ])+
            D1*(ui[ixm1][izm1]+ui[ixm1][izp1]+ui[ixp1][izm1]+ui[ixp1][izp1]));
        }
        
        // 1 < iz < nz-2
        for (int iz=2; iz<_nz-2; ++iz) {
          int izm1 = iz-1;
          int izp1 = iz+1;
          int izm2 = iz-2;
          int izp2 = iz+2;
          float r = _r[ix][iz];
          up[ix][iz] = (2.0f+C6*r)*ui[ix][iz]-um[ix][iz]+r*(
            C5*(ui[ix  ][izm1]+ui[ix  ][izp1]+ui[ixm1][iz  ]+ui[ixp1][iz  ])+
            C4*(ui[ixm1][izm1]+ui[ixm1][izp1]+ui[ixp1][izm1]+ui[ixp1][izp1])+
            C3*(ui[ix  ][izm2]+ui[ix  ][izp2]+ui[ixm2][iz  ]+ui[ixp2][iz  ])+
            C2*(ui[ixm2][izm1]+ui[ixm2][izp1]+ui[ixm1][izm2]+ui[ixm1][izp2]+
                ui[ixp1][izm2]+ui[ixp1][izp2]+ui[ixp2][izm1]+ui[ixp2][izp1]));
        }

        // iz = nz-2
        {
          int iz = _nz-2;
          int izm1 = iz-1;
          int izp1 = iz+1;
          float r = _r[ix][iz];
          up[ix][iz] = (2.0f+D3*r)*ui[ix][iz]-um[ix][iz]+r*(
            D2*(ui[ixm1][iz  ]+ui[ix  ][izm1]+ui[ix  ][izp1]+ui[ixp1][iz  ])+
            D1*(ui[ixm1][izm1]+ui[ixm1][izp1]+ui[ixp1][izm1]+ui[ixp1][izp1]));
        }
      }
    }});
  }

  private void absorb(
    final float[][] um, final float[][] ui, final float[][] up)
  {
    //absorbReynolds(um,ui,up);
    absorbClayton(um,ui,up);
    //absorbHale(um,ui,up);
    //Parallel.loop(0,_nx,new Parallel.LoopInt() {
    //  public void compute(int ix) {
    //    absorbSliceX(ix,_b,um);
    //    absorbSliceX(ix,_b,ui);
    //    absorbSliceX(ix,_b,up);
    //  }
    //});
  }

  private void absorbReynolds(float[][] um, float[][] ui, float[][] up) {
    float a = (float)(_dt/_dx);

    for (int ix=0; ix<_b; ++ix) {
      float w = (float)ix/(float)_b;
      for (int iz=ix+1; iz<_nz-1-ix; ++iz) {
        float vp = ui[ix][iz]+ui[ix+1][iz]-um[ix+1][iz]+
          (a/_s[ix][iz])*(ui[ix+1][iz]-ui[ix][iz]-um[ix+2][iz]+um[ix+1][iz]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }

    for (int ix=_nx-_b; ix<_nx; ++ix) {
      float w = (float)((_nx-1)-ix)/(float)_b;
      for (int iz=_nx-ix; iz<_nz-1-(_nx-1-ix); ++iz) {
        float vp = ui[ix][iz]+ui[ix-1][iz]-um[ix-1][iz]+
          (a/_s[ix][iz])*(ui[ix-1][iz]-ui[ix][iz]-um[ix-2][iz]+um[ix-1][iz]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }

    for (int iz=0; iz<_b; ++iz) {
      float w = (float)iz/(float)_b;
      for (int ix=iz+1; ix<_nx-1-iz; ++ix) {
        float vp = ui[ix][iz]+ui[ix][iz+1]-um[ix][iz+1]+
          (a/_s[ix][iz])*(ui[ix][iz+1]-ui[ix][iz]-um[ix][iz+2]+um[ix][iz+1]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }

    for (int iz=_nz-_b; iz<_nz; ++iz) {
      float w = (float)((_nz-1)-iz)/(float)_b;
      for (int ix=_nz-iz; ix<_nx-1-(_nz-1-iz); ++ix) {
        float vp = ui[ix][iz]+ui[ix][iz-1]-um[ix][iz-1]+
          (a/_s[ix][iz])*(ui[ix][iz-1]-ui[ix][iz]-um[ix][iz-2]+um[ix][iz-1]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }

//    // Corners
////    iz = 0;
////    ix = 0;
//    for (ix=0,iz=0; ix<_b; ++ix,++iz) {
//      vp[ix][iz] = ui[ix][iz]+(a/_s[ix][iz])*(ui[ix+1][iz+1]-ui[ix][iz]);
//    }
////
////    iz = 0;
////    ix = _nx-1;
//    for (ix=_nx-1,iz=0; iz<_b; --ix,++iz) {
//      vp[ix][iz] = ui[ix][iz]+(a/_s[ix][iz])*(ui[ix-1][iz+1]-ui[ix][iz]);
//    }
////
////    iz = _nz-1;
////    ix = 0;
//    for (ix=0,iz=_nz-1; ix<_b; ++ix,--iz) {
//      vp[ix][iz] = ui[ix][iz]+(a/_s[ix][iz])*(ui[ix+1][iz-1]-ui[ix][iz]);
//    }
////
////    iz = _nz-1;
////    ix = _nx-1;
//    for (ix=_nx-1,iz=_nx-1; ix>=_nx-_b; --ix,--iz) {
//      vp[ix][iz] = ui[ix][iz]+(a/_s[ix][iz])*(ui[ix-1][iz-1]-ui[ix][iz]);
//    }



  }

  private void absorbClayton(
    final float[][] um, final float[][] ui, final float[][] up)
  {
    final float oxt = (float)(1.0/(_dx*_dt));
    final float oxx = (float)(1.0/(_dx*_dx));
    final float ott = (float)(1.0/(_dt*_dt));
    final float[][] cp = copy(up);

    Parallel.loop(_b,new Parallel.LoopInt() {
    public void compute(int ix) {
    //for (int ix=0; ix<_b; ++ix) {
      float w = (float)ix/(float)_b;
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

    Parallel.loop(_nx-_b,_nx,new Parallel.LoopInt() {
    public void compute(int ix) {
    //for (int ix=_nx-_b; ix<_nx; ++ix) {
      float w = (float)((_nx-1)-ix)/(float)_b;
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

    Parallel.loop(_b,new Parallel.LoopInt() {
    public void compute(int iz) {
    //for (int iz=0; iz<_b; ++iz) {
      float w = (float)iz/(float)_b;
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

    Parallel.loop(_nz-_b,_nz,new Parallel.LoopInt() {
    public void compute(int iz) {
    //for (int iz=_nz-_b; iz<_nz; ++iz) {
      float w = (float)((_nz-1)-iz)/(float)_b;
      for (int ix=_nz-iz; ix<_nx-1-(_nz-1-iz); ++ix) {
      //for (int ix=_nz+1-iz; ix<_nx-2-(_nz-1-iz); ++ix) {
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
    for (int ix=0; ix<_b; ++ix) {
      float w = (float)ix/(float)_b;
      for (int iz=0; iz<_b; ++iz) {
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[ix][iz]));
        float vp = ui[ix][iz]+r*(ui[ix][iz+1]+ui[ix+1][iz]-2.0f*ui[ix][iz]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }
    for (int ix=0; ix<_b; ++ix) {
      float w = (float)ix/(float)_b;
      for (int iz=_nz-_b; iz<_nz; ++iz) {
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[ix][iz]));
        float vp = ui[ix][iz]+r*(ui[ix][iz-1]+ui[ix+1][iz]-2.0f*ui[ix][iz]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }
    for (int ix=_nx-_b; ix<_nx; ++ix) {
      float w = (float)((_nx-1)-ix)/(float)_b;
      for (int iz=0; iz<_b; ++iz) {
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[ix][iz]));
        float vp = ui[ix][iz]+r*(ui[ix][iz+1]+ui[ix-1][iz]-2.0f*ui[ix][iz]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }
    for (int ix=_nx-_b; ix<_nx; ++ix) {
      float w = (float)((_nx-1)-ix)/(float)_b;
      for (int iz=_nz-_b; iz<_nz; ++iz) {
        float r = (float)(_dt/(sqrt(2.0)*_dx*_s[ix][iz]));
        float vp = ui[ix][iz]+r*(ui[ix][iz-1]+ui[ix-1][iz]-2.0f*ui[ix][iz]);
        up[ix][iz] = w*up[ix][iz]+(1.0f-w)*vp;
      }
    }

  }

  private void absorbSliceX(int ix ,int b, float[][] u) {
    int nzmb = _nz-b;
    //float w = 0.005f;
    //float w = 0.0010f;
    float w = 0.0015f;
    //float w = 0.0020f;
    //float w = 0.0025f;
    //float w = 0.0050f;
    float ws = w*w;
    float[] uix = u[ix];
    if (ix<b) {
      float x = (float)(b-ix);
      float xx = x*x;
      for (int iz=0; iz<b; ++iz) {
        float z = (float)(b-iz);
        float rs = (xx+z*z)*ws;
        float e = exp(-rs);
        uix[iz      ] *= e;
        uix[_nz-1-iz] *= e;
      }
      for (int iz=b; iz<nzmb; ++iz) {
        float r = x*w;
        uix[iz] *= exp(-r*r);
      } 
    } else if (ix<_nx-b) {
      for (int iz=0; iz<b; ++iz) {
        float r = (float)(b-iz)*w;
        float e = exp(-r*r);
        uix[iz      ] *= e;
        uix[_nz-1-iz] *= e;
      }
    } else {
      float x = 1.0f+(float)(b+ix-_nx);
      float xx = x*x;
      for (int iz=0; iz<b; ++iz) {
        float z = (float)(b-iz);
        float rs = (xx+z*z)*ws;
        float e = exp(-rs);
        uix[iz      ] *= e;
        uix[_nz-1-iz] *= e;
      }
      for (int iz=b; iz<nzmb; ++iz) {
        float r = x*w;
        uix[iz] *= exp(-r*r);
      }
    }
  }

  private void absorbHale(float[][] um, float[][] ui, float[][] up) {
    int ix,iz;
    float d = (float)(_dx/_dt);
    float odx = (float)(0.5/_dx);
    float odt = (float)(0.5/_dt);

    ix = 0;
    for (iz=1; iz<_nz-1; ++iz) {
      float si = _s[ix][iz];
      float pz = odx*(ui[ix+1][iz+1]-ui[ix+1][iz-1]);
      float pt = odt*(up[ix+1][iz  ]-um[ix+1][iz  ]);
      float a = 1.0f-pz*pz/(si*si*pt*pt);
      float c = (a>0)?sqrt(a):0.0f;
      float r = si*c*d;
      float g = (1.0f-r)/(1.0f+r);
      up[ix][iz] = ui[ix+1][iz]+g*(up[ix+1][iz]-ui[ix][iz]);
    }

    ix = _nx-1;
    for (iz=1; iz<_nz-1; ++iz) {
      float si = _s[ix][iz];
      float pz = odx*(ui[ix-1][iz+1]-ui[ix-1][iz-1]);
      float pt = odt*(up[ix-1][iz  ]-um[ix-1][iz  ]);
      float a = 1.0f-pz*pz/(si*si*pt*pt);
      float c = (a>0)?sqrt(a):0.0f;
      float r = si*c*d;
      float g = (1.0f-r)/(1.0f+r);
      up[ix][iz] = ui[ix-1][iz]+g*(up[ix-1][iz]-ui[ix][iz]);
    }

    iz = 0;
    for (ix=1; ix<_nx-1; ++ix) {
      float si = _s[ix][iz];
      float px = odx*(ui[ix+1][iz+1]-ui[ix-1][iz+1]);
      float pt = odt*(up[ix  ][iz+1]-um[ix  ][iz+1]);
      float a = 1.0f-px*px/(si*si*pt*pt);
      float c = (a>0)?sqrt(a):0.0f;
      float r = si*c*d;
      float g = (1.0f-r)/(1.0f+r);
      up[ix][iz] = ui[ix][iz+1]+g*(up[ix][iz+1]-ui[ix][iz]);
    }

    iz = _nz-1;
    for (ix=1; ix<_nx-1; ++ix) {
      float si = _s[ix][iz];
      float px = odx*(ui[ix+1][iz-1]-ui[ix-1][iz-1]);
      float pt = odt*(up[ix  ][iz-1]-um[ix  ][iz-1]);
      float a = 1.0f-px*px/(si*si*pt*pt);
      float c = (a>0)?sqrt(a):0.0f;
      float r = si*c*d;
      float g = (1.0f-r)/(1.0f+r);
      up[ix][iz] = ui[ix][iz-1]+g*(up[ix][iz-1]-ui[ix][iz]);
    }
  }

  private float[][] extendModel(float[][] c) {
    int nz = c[0].length;
    int nx = c.length;
    float[][] v = new float[nx+2*_b][nz+2*_b];
    copy(nz,nx,0,0,c,_b,_b,v);
    for (int ix=_b; ix<nx+_b; ++ix) {
      for (int iz=0, jz=nz+_b; iz<_b; ++iz, ++jz) {
        v[ix][iz] = v[ix][_b];
        v[ix][jz] = v[ix][nz+_b-1];
      }
    }
    for (int ix=0, jx=nx+_b; ix<_b; ++ix, ++jx) {
      copy(v[_b],v[ix]);
      copy(v[nx+_b-1],v[jx]);
    }
    return v;
  }

  private void reverse3(float[][][] f) {
    int n3 = f.length;
    float[][] t;
    for (int i3=0, j3=n3-1; i3<n3/2; ++i3, --j3) {
      t = f[i3];
      f[i3] = f[j3];
      f[j3] = t;
    }
  }

}
