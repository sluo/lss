package lss.vel;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

// testing
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
    _u0m = new float[_nx][_nz];
    _u0i = new float[_nx][_nz];
    _u0p = new float[_nx][_nz];
    _u1m = new float[_nx][_nz];
    _u1i = new float[_nx][_nz];
    _u1p = new float[_nx][_nz];
  }

  public void modelAcousticData(
    Source source, Receiver receiver, float[][] s0, float[][] d)
  {
    modelAcousticDataAndWavefield(source,receiver,s0,d,null);
  }

  public void modelAcousticWavefield(
    Source source, float[][] s0, float[][][] u0)
  {
    modelAcousticDataAndWavefield(source,null,s0,null,u0);
  }

  public void modelAcousticDataAndWavefield(
    Source source, Receiver receiver,
    float[][] s0, float[][] d, float[][][] u0)
  {
    updateParameters(source,receiver,s0,null);
    modelAcoustic(d,u0);
  }

  public void modelBornData(
    Source source, Receiver receiver, float[][] s0, float[][] s1, float[][] d)
  {
    modelBornDataAndWavefield(source,receiver,s0,s1,d,null,null);
  }

  public void modelBornWavefield(
    Source source, float[][] s0, float[][] s1, float[][][] u0, float[][][] u1)
  {
    modelBornDataAndWavefield(source,null,s0,s1,null,u0,u1);
  }

  public void modelBornDataAndWavefield(
    Source source, Receiver receiver, float[][] s0, float[][] s1, 
    float[][] d, float[][][] u0, float[][][] u1)
  {
    updateParameters(source,receiver,s0,s1);
    modelBorn(d,u0,u1);
  }

  public static abstract class Source {
    public void add(double t, float[][] f) {
      add(1.0f,t,f);
    }
    public void sub(double t, float[][] f) {
      add(-1.0f,t,f);
    }
    abstract void add(float sign, double t, float[][] f);
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

  private static int _b = 50; // absorbing boundary size

  private Source _source;
  private Receiver _receiver;
  private double _t; // current time
  private int _it; // current time index
  private int _nx,_nz,_nt; // number of samples
  private double _dx,_dz,_dt; // sampling intervals
  private float[][] _s0; // slowness
  private float[][] _s1; // perturbation slowness
  private float[][] _r0; // dt^2/(6*dx^2*s0^2)
  private float[][] _u0m,_u0i,_u0p; // u0(x;t-1), u0(x;t), u0(x;t+1)
  private float[][] _u1m,_u1i,_u1p; // u1(x;t-1), u1(x;t), u1(x;t+1)

  private void updateParameters(
    Source source, Receiver receiver, float[][] s0, float[][] s1)
  {
    _source = source;
    _receiver = receiver;
    float r = (float)(_dt*_dt/(6.0*_dx*_dx));
    _s0 = extendModel(s0);
    _r0 = copy(_s0);
    mul(_r0,_r0,_r0);
    div(  r,_r0,_r0);
    if (s1!=null) {
      _s1 = extendModel(s1);
    }
  }

  private void modelAcoustic(float[][] d, float[][][] u0) {
    int nz = _nz-2*_b;
    int nx = _nx-2*_b;
    int nt = _nt;

    zero(_u0m);
    zero(_u0i);
    zero(_u0p);
    for (int it=0; it<nt; ++it) {
      
      // One time step
      step(_u0m,_u0i,_u0p);

      // Source injection
      _source.add(_dt*it,_u0p);
      if (it>=2) {
        _source.sub(_dt*(it-2),_u0m); // SEP trick?
      }

      // Absorbing boundaries
      absorb(_u0m,_u0i,_u0p);

      // Copy data and wavefield
      if (_receiver!=null) {
        int nr = _receiver.kzr.length;
        for (int ir=0; ir<nr; ++ir) {
          int kz = _receiver.kzr[ir]+_b;
          int kx = _receiver.kxr[ir]+_b;
          d[ir][it] = _u0p[kx][kz];
        }
      }
      if (u0!=null) {
        copy(nz,nx,_b,_b,_u0p,0,0,u0[it]);
      }

      // Rotate arrays
      float[][] t = _u0m;
      _u0m = _u0i;
      _u0i = _u0p;
      _u0p = t;
    }

    // Time reverse if necessary
    if (_source instanceof AdjointSource) {
      reverse3(u0);
    }
  }

  private void modelBorn(float[][] d, float[][][] u0, float[][][] u1) {
    int nz = _nz-2*_b;
    int nx = _nx-2*_b;
    int nt = _nt;

    zero(_u0m);
    zero(_u0i);
    zero(_u0p);
    zero(_u1m);
    zero(_u1i);
    zero(_u1p);
    for (int it=0; it<nt; ++it) {
      
      // One time step
      step(_u0m,_u0i,_u0p);
      step(_u1m,_u1i,_u1p);

      // Source injection for u0
      _source.add(_dt*it,_u0p);
      _source.sub(_dt*(it-2),_u0m); // SEP trick?

      // Virtual source injection for u1
      //for (int ix=0; ix<_nx; ++ix) {
      Parallel.loop(_nx,new Parallel.LoopInt() {
      public void compute(int ix) {
        for (int iz=0; iz<_nz; ++iz) {
          _u1p[ix][iz] -= 2.0f*(_s1[ix][iz]/_s0[ix][iz])*(
            _u0p[ix][iz]-2.0f*_u0i[ix][iz]+_u0m[ix][iz]);
        }
      }});
      //}

      // Absorbing boundaries
      absorb(_u0m,_u0i,_u0p);
      absorb(_u1m,_u1i,_u1p);

      // Copy data and wavefield
      if (_receiver!=null) {
        int nr = _receiver.kzr.length;
        for (int ir=0; ir<nr; ++ir) {
          int kz = _receiver.kzr[ir]+_b;
          int kx = _receiver.kxr[ir]+_b;
          d[ir][it] = _u1p[kx][kz];
        }
      }
      if (u0!=null) {
        copy(nz,nx,_b,_b,_u0p,0,0,u0[it]);
        copy(nz,nx,_b,_b,_u1p,0,0,u1[it]);
      }

      // Rotate arrays
      float[][] t0 = _u0m;
      _u0m = _u0i;
      _u0i = _u0p;
      _u0p = t0;
      float[][] t1 = _u1m;
      _u1m = _u1i;
      _u1i = _u1p;
      _u1p = t1;
    }

    // Time reverse if necessary
    if (_source instanceof AdjointSource) {
      reverse3(u0);
      reverse3(u1);
    }
  }

  public static class RickerSource extends Source {
    public RickerSource(double fpeak, int kzs, int kxs) {
      _fpeak = fpeak;
      _kzs = kzs+_b;
      _kxs = kxs+_b;
      _tdelay = 1.0/fpeak;
    }
    public void add(float sign, double t, float[][] f) {
      f[_kxs][_kzs] += sign*ricker(t-_tdelay);
    }
    private float ricker(double t) {
      double x = PI*_fpeak*t;
      double xx = x*x;
      return (float)((1.0-2.0*xx)*exp(-xx));
    }
    private int _kzs,_kxs;
    private double _fpeak,_tdelay;
  }

  public static class PlaneWaveSource extends Source {
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
    public void add(float sign, double t, float[][] f) {
      //double p = sin(_angle)/_v[_kzs][_nx/2]; // TODO: function of x?
      for (int ix=0; ix<_nx; ++ix) {
        double t0 = _ps*(ix-_kxs)*_dx;
        // TODO + or - t0?
        f[ix+_b][_kzs+_b] += sign*_w[ix]*ricker(t-_tdelay+t0);
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

  public static class DefaultSource extends Source {
    public DefaultSource(double dt, int kzs, int kxs, float[] s) {
      this(dt,new int[]{kzs},new int[]{kxs},new float[][]{s});
    }
    public DefaultSource(double dt, int[] kzs, int[] kxs, float[][] s) {
      Check.argument(kzs.length==kxs.length,"kzs.length=kxs.length");
      _ns = kzs.length;
      _nt = s[0].length;
      _dt = dt;
      _kzs = kzs;
      _kxs = kxs;
      _s = s;
    }
    public void add(float sign, double t, float[][] f) {
      int it = (int)(t/_dt);
      for (int is=0; is<_ns; is++) {
        int kzs = _b+_kzs[is];
        int kxs = _b+_kxs[is];
        f[kxs][kzs] += sign*_s[is][it];
      }
    }
    private double _dt;
    private int _ns,_nt;
    private int[] _kzs = null;
    private int[] _kxs = null;
    private float[][] _s = null;
  }

  public static class AdjointSource extends Source {
    public AdjointSource(double dt, int kzs, int kxs, float[] s) {
      this(dt,new int[]{kzs},new int[]{kxs},new float[][]{s});
    }
    public AdjointSource(double dt, int[] kzs, int[] kxs, float[][] s) {
      _source = new DefaultSource(dt,kzs,kxs,reverse(s));
    }
    public void add(float sign, double t, float[][] f) {
      _source.add(sign,t,f);
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

  private void step(
    final float[][] um, final float[][] ui, final float[][] up)
  {
    Parallel.loop(1,_nx-1,new Parallel.LoopInt() {
    public void compute(int ix) {
      int ixm = ix-1;
      int ixp = ix+1;
      for (int iz=1; iz<_nz-1; ++iz) {
        int izm = iz-1;
        int izp = iz+1;
        float r = _r0[ix][iz];
        up[ix][iz] = (2.0f-r*20.0f)*ui[ix][iz]-um[ix][iz]+r*(
          1.0f*(ui[ixm][izm]+ui[ixm][izp]+ui[ixp][izm]+ui[ixp][izp])+
          4.0f*(ui[ixm][iz ]+ui[ix ][izm]+ui[ix ][izp]+ui[ixp][iz ])
        );
      }
    }});
  }

  private void absorb(
    final float[][] um, final float[][] ui, final float[][] up)
  {
    absorbHale(um,ui,up);
    Parallel.loop(0,_nx,new Parallel.LoopInt() {
      public void compute(int ix) {
        absorbSliceX(ix,_b,um);
        absorbSliceX(ix,_b,ui);
        absorbSliceX(ix,_b,up);
      }
    });
  }

  private void absorbSliceX(int ix ,int b, float[][] u) {
    int nzmb = _nz-b;
    //float w = 0.005f;
    float w = 0.0015f;
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
      float si = _s0[ix][iz];
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
      float si = _s0[ix][iz];
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
      float si = _s0[ix][iz];
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
      float si = _s0[ix][iz];
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
    int nxp = nx+2*_b;
    int nzp = nz+2*_b;
    float[][] v = new float[nxp][nzp];
    copy(nz,nx,0,0,c,_b,_b,v);
    for (int ix=_b; ix<nx+_b; ++ix) {
      float dv1 = v[ix][_b+1]-v[ix][_b];
      float dv2 = v[ix][nz+_b-1]-v[ix][nz+_b-2];
      for (int iz=0, jz=nz+_b; iz<_b; ++iz, ++jz) {
        //v[ix][iz] = v[ix][_b]-(_b-iz)*dv1;
        //v[ix][jz] = v[ix][nz+_b-1]+(jz-(nz+_b-1))*dv2;
        v[ix][iz] = v[ix][_b];
        v[ix][jz] = v[ix][nz+_b-1];
      }
    }
    for (int ix=0, jx=nx+_b; ix<_b; ++ix, ++jx) {
      copy(v[_b],v[ix]);
      copy(v[nx+_b-1],v[jx]);
    }

    float[][] w = new float[nxp][nzp];
    for (int ix=0; ix<=nxp/2; ++ix) {
      for (int iz=0; iz<=nzp/2; ++iz) {
        float xf = (ix<_b)?(float)(_b-ix):0.0f;
        float zf = (iz<_b)?(float)(_b-iz):0.0f;
        float wi = 1.0f-min(sqrt(xf*xf+zf*zf)/_b,1.0f);
        w[ix      ][iz      ] = wi;
        w[nxp-1-ix][iz      ] = wi;
        w[ix      ][nzp-1-iz] = wi;
        w[nxp-1-ix][nzp-1-iz] = wi;
      }
    }
    mul(w,w,w);
    float[][] t = copy(v);
    new RecursiveExponentialFilter(_b).apply(t,t);
    v = add(mul(w,v),mul(sub(1.0f,w),t));
    //SimplePlot.asPixels(t);
    //SimplePlot sp = SimplePlot.asPixels(v); sp.addColorBar();
    //SimplePlot.asPixels(w);
    return v;
  }
  private float[][] xextendModel(float[][] c) {
    int nz = c[0].length;
    int nx = c.length;
    float[][] v = new float[nx+2*_b][nz+2*_b];
    copy(nz,nx,0,0,c,_b,_b,v);
    for (int ix=_b; ix<nx+_b; ++ix) {
      float dv1 = v[ix][_b+1]-v[ix][_b];
      float dv2 = v[ix][nz+_b-1]-v[ix][nz+_b-2];
      for (int iz=0, jz=nz+_b; iz<_b; ++iz, ++jz) {
        //v[ix][iz] = v[ix][_b]-(_b-iz)*dv1;
        //v[ix][jz] = v[ix][nz+_b-1]+(jz-(nz+_b-1))*dv2;
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
