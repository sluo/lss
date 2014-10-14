package lss.vel;

import java.awt.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class AcousticWavefieldY {

  public enum SourceWavelet {
    RICKER,
    PLANE_WAVE
  }

  /**
   * Source function.
   */
  public static abstract class Source {
    public void add(double t, float[][] f) {
      add(1.0f,t,f);
    }
    public void sub(double t, float[][] f) {
      add(-1.0f,t,f);
    }
    abstract void add(float sign, double t, float[][] f);
  }

  public AcousticWavefieldY(
    Sampling sz, Sampling sx, Sampling st)
  {
    _dz = sz.getDelta();
    _dx = sx.getDelta();
    _dt = st.getDelta();
    _nz = sz.getCount()+2*_b;
    _nx = sx.getCount()+2*_b;
    _nt = st.getCount();
    Check.argument(_dx==_dz,"dx==dz");
    _u = new float[_nt][_nx][_nz];
  }

  public void propagate(Source source, float[][] v) {
    Check.argument(v[0].length==_nz-2*_b,"dimensions");
    Check.argument(v.length==_nx-2*_b,"dimensions");
    updateAndPropagate(source,v);
    if (source instanceof AdjointSource) {
      reverse3(_u);
    }
  }

  public void forwardPropagate(Source source, float[][] v) {
    propagate(source,v);
  }

  public void backPropagate(Source source, float[][] v) {
    propagate(source,v);
    reverse3(_u);
  }

/*
  public void forwardPropagate(
    double freq, int kzs, int kxs, float[][] c)
  {
    Source source = new RickerSource(freq,kzs,kxs);
    propagate(source,c);
  }

  public void forwardPropagate(
    SourceWavelet wavelet, double freq, int kzs, int kxs, float[][] c)
  {
    if (wavelet==SourceWavelet.RICKER) {
      Source source = new RickerSource(freq,kzs,kxs);
      propagate(source,c);
    } else {
      // TODO
    }
  }

  public void backPropagate(
    float[] s, int kzs, int kxs, float[][] c)
  {
    Source source = new AdjointSource(_dt,kzs,kxs,s);
    propagate(source,c);
    reverse3(_u);
  }

  public void backPropagate(
    float[][] s, int kzs[], int[] kxs, float[][] c)
  {
    Source source = new AdjointSource(_dt,kzs,kxs,s);
    propagate(source,c);
    reverse3(_u);
  }
*/

  public float[] getWavefield(int kzr, int kxr) {
    float[] y = new float[_nt];
    kzr += _b; kxr += _b;
    for (int it=0; it<_nt; it++)
      y[it] = _u[it][kxr][kzr];
    return y;
  }

  public float[][] getWavefield(int[] kzr, int[] kxr) {
    Check.argument(kzr.length==kxr.length,"kzr.length=kxr.length");
    int nr = kzr.length;
    float[][] y = new float[nr][_nt];
    for (int it=0; it<_nt; it++) {
      for (int ir=0; ir<nr; ir++) {
        int kzri = kzr[ir]+_b;
        int kxri = kxr[ir]+_b;
        y[ir][it] = _u[it][kxri][kzri];
      }
    }
    return y;
  }

  public float[][][] getWavefield() {
    //return _u;
    return copy(_nz-2*_b,_nx-2*_b,_nt,_b,_b,0,_u);
  }

  public float[][] correlate(AcousticWavefield wavefield) {
    float[][][] u = wavefield.getWavefield();
    float[][] y = new float[_nx-2*_b][_nz-2*_b];
    for (int it=0; it<_nt; ++it) {
      for (int ix=_b, jx=0; ix<_nx-_b; ++ix, ++jx) {
        for (int iz=_b, jz=0; iz<_nz-_b; ++iz, ++jz) {
          y[jx][jz] += _u[it][ix][iz]*u[it][ix][iz];
        }
      }
    }
    return y;
  }

  /////////////////////////////////////////////////////////////////////////
  // private

  private static int _b = 50; // absorbing boundary size
  private Source _source;
  private double _t; // current time
  private int _it; // current time index
  private int _nx,_nz,_nt; // number of samples
  private double _dx,_dz,_dt; // sampling intervals
  private float[][] _v; // velocity
  private float[][] _s; // v^2*dt^2/(6*dx^2)
  private float[][] _um,_ui,_up; // u(x;t-1), u(x;t), u(x;t+1)
  private float[][][] _u; // wavefield

  private void updateAndPropagate(Source source, float[][] v) {
    //double odx = 1.0/_dx;
    //double r = 0.5*_dt*_dt*odx*odx;
    //_rvs = extendModel(v);
    //_v = copy(_rvs);
    //mul(_rvs,_rvs,_rvs);
    //mul((float)r,_rvs,_rvs);
    //_source = source;
    //propagate();


    float r = (float)(_dt*_dt/(6.0*_dx*_dx));
    _v = extendModel(v);
    _s = copy(_v);
    mul(_s,_s,_s);
    mul( r,_s,_s);
    _source = source;
    propagate();
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
    //public AdjointSource(double dt, int kzs, int kxs, float[] s) {
    //  this(dt,new int[]{kzs},new int[]{kxs},new float[][]{s});
    //}
    //public AdjointSource(double dt, int[] kzs, int[] kxs, float[][] s) {
    //  Check.argument(kzs.length==kxs.length,"kzs.length=kxs.length");
    //  _ns = kzs.length;
    //  _nt = s[0].length;
    //  _dt = dt;
    //  _kzs = kzs;
    //  _kxs = kxs;
    //  _s = s;
    //}
    //public void add(float sign, double t, float[][] f) {
    //  int it = (int)((t)/_dt);
    //  for (int is=0; is<_ns; is++) {
    //    int kzs = _b+_kzs[is];
    //    int kxs = _b+_kxs[is];
    //    if (it>=0 && it<_nt-1)
    //      f[kxs][kzs] += sign*_s[is][_nt-1-it];
    //  }
    //}
    //private double _dt;
    //private int _ns,_nt;
    //private int[] _kzs;
    //private int[] _kxs;
    //private float[][] _s;

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

  private void propagate() {
    _um = _u[0];
    _ui = _u[1];
    _up = _u[2];
    zero(_um);
    zero(_ui);
    zero(_up);
    for (_it=0; _it<_nt-3; ++_it) {
      System.out.format("\r%d",_it+3);
      step();
    }
    System.out.print("\n");
  }
  
  private void step() {
    Parallel.loop(1,_nx-1,new Parallel.LoopInt() {
    public void compute(int ix) {
      int ixm = ix-1;
      int ixp = ix+1;
      for (int iz=1; iz<_nz-1; ++iz) {
        int izm = iz-1;
        int izp = iz+1;
        float si = _s[ix][iz];
        _up[ix][iz] = (2.0f-si*20.0f)*_ui[ix][iz]-_um[ix][iz]+si*(
          1.0f*(_ui[ixm][izm]+_ui[ixm][izp]+_ui[ixp][izm]+_ui[ixp][izp])+
          4.0f*(_ui[ixm][iz ]+_ui[ix ][izm]+_ui[ix ][izp]+_ui[ixp][iz ])
        );
      }
    }});
    _source.sub(_dt*_it,_up);
    if (_it>=2) _source.add(_dt*(_it-2),_um); // SEP trick?
    absorb();
    _um = _u[_it+1];
    _ui = _u[_it+2];
    _up = _u[_it+3];
  }

  private void absorb() {
    //absorbReynolds(_um,_ui,_up);
    absorbHale();
    /*
    */
    Parallel.loop(0,_nx,new Parallel.LoopInt() {
      public void compute(int ix) {
        absorbSliceX(ix,_b,_um);
        absorbSliceX(ix,_b,_ui);
        absorbSliceX(ix,_b,_up);
      }
    });
    //absorbClayton();
  }

  private static void absorbSliceX(int ix ,int b, float[][] u) {
    int nz = u[0].length;
    int nx = u.length;
    int nzmb = nz-b;
    //float w = 0.005f;
    float w = 0.0025f;
    //float w = 0.0050f;
    //float w = 0.30f/b;
    float ws = w*w;
    float[] uix = u[ix];
    if (ix<b) {
      float x = (float)(b-ix);
      float xx = x*x;
      for (int iz=0; iz<b; ++iz) {
        float z = (float)(b-iz);
        float rs = (xx+z*z)*ws;
        float e = exp(-rs);
        uix[iz     ] *= e;
        uix[nz-1-iz] *= e;
      }
      for (int iz=b; iz<nzmb; ++iz) {
        float r = x*w;
        uix[iz] *= exp(-r*r);
      } 
    } else if (ix<nx-b) {
      for (int iz=0; iz<b; ++iz) {
        float r = (float)(b-iz)*w;
        float e = exp(-r*r);
        uix[iz     ] *= e;
        uix[nz-1-iz] *= e;
      }
    } else {
      float x = 1.0f+(float)(b+ix-nx);
      float xx = x*x;
      for (int iz=0; iz<b; ++iz) {
        float z = (float)(b-iz);
        float rs = (xx+z*z)*ws;
        float e = exp(-rs);
        uix[iz     ] *= e;
        uix[nz-1-iz] *= e;
      }
      for (int iz=b; iz<nzmb; ++iz) {
        float r = x*w;
        uix[iz] *= exp(-r*r);
      }
    }
  }

  private void absorbClayton() {
    // A1: ix = 0
    for (int iz=0, kx=0; iz<_nz; ++iz) {
      float a = _v[kx][iz]*(float)(_dt/_dx);
      _up[kx][iz] = _ui[kx][iz]+a*(_ui[kx+1][iz]-_ui[kx][iz]);
    }

    // A1: ix = nx-1
    for (int iz=0, kx=_nx-1; iz<_nz; ++iz) {
      float a = _v[kx][iz]*(float)(_dt/_dx);
      _up[kx][iz] = _ui[kx][iz]-a*(_ui[kx][iz]-_ui[kx-1][iz]);
    }

    // A1: iz = 0
    for (int ix=0, kz=0; ix<_nx; ++ix) {
      float a = _v[ix][kz]*(float)(_dt/_dz);
      _up[ix][kz] = _ui[ix][kz]+a*(_ui[ix][kz+1]-_ui[ix][kz]);
    }

    // A1: iz = nz-1
    for (int ix=0, kz=_nz-1; ix<_nx; ++ix) {
      float a = _v[ix][kz]*(float)(_dt/_dz);
      _up[ix][kz] = _ui[ix][kz]-a*(_ui[ix][kz]-_ui[ix][kz-1]);
    }
  }

  private void absorbHale() {
    int ix,iz;
    float d = (float)(_dx/_dt);
    float odx = (float)(0.5/_dx);
    float odt = (float)(0.5/_dt);

    ix = 0;
    for (iz=1; iz<_nz-1; ++iz) {
      float vi = _v[ix][iz];
      float pz = odx*(_ui[ix+1][iz+1]-_ui[ix+1][iz-1]);
      float pt = odt*(_up[ix+1][iz  ]-_um[ix+1][iz  ]);
      float a = 1.0f-vi*vi*pz*pz/(pt*pt);
      float c = (a>0)?sqrt(a):0.0f;
      float g = (vi-c*d)/(vi+c*d);
      _up[ix][iz] = _ui[ix+1][iz]+g*(_up[ix+1][iz]-_ui[ix][iz]);
    }

    ix = _nx-1;
    for (iz=1; iz<_nz-1; ++iz) {
      float vi = _v[ix][iz];
      float pz = odx*(_ui[ix-1][iz+1]-_ui[ix-1][iz-1]);
      float pt = odt*(_up[ix-1][iz  ]-_um[ix-1][iz  ]);
      float a = 1.0f-vi*vi*pz*pz/(pt*pt);
      float c = (a>0)?sqrt(a):0.0f;
      float g = (vi-c*d)/(vi+c*d);
      _up[ix][iz] = _ui[ix-1][iz]+g*(_up[ix-1][iz]-_ui[ix][iz]);
    }

    iz = 0;
    for (ix=1; ix<_nx-1; ++ix) {
      float vi = _v[ix][iz];
      float px = odx*(_ui[ix+1][iz+1]-_ui[ix-1][iz+1]);
      float pt = odt*(_up[ix  ][iz+1]-_um[ix  ][iz+1]);
      float a = 1.0f-vi*vi*px*px/(pt*pt);
      float c = (a>0)?sqrt(a):0.0f;
      float g = (vi-c*d)/(vi+c*d);
      _up[ix][iz] = _ui[ix][iz+1]+g*(_up[ix][iz+1]-_ui[ix][iz]);
    }

    iz = _nz-1;
    for (ix=1; ix<_nx-1; ++ix) {
      float vi = _v[ix][iz];
      float px = odx*(_ui[ix+1][iz-1]-_ui[ix-1][iz-1]);
      float pt = odt*(_up[ix  ][iz-1]-_um[ix  ][iz-1]);
      float a = 1.0f-vi*vi*px*px/(pt*pt);
      float c = (a>0)?sqrt(a):0.0f;
      float g = (vi-c*d)/(vi+c*d);
      _up[ix][iz] = _ui[ix][iz-1]+g*(_up[ix][iz-1]-_ui[ix][iz]);
    }
  }

  private void absorbReynolds(
    float[][] um, float[][] ui, float[][] up)
  {
    int km,ki,kp;

    // ix = 0
    kp = 0;
    ki = 1;
    km = 2;
    for (int iz=0; iz<_nz; ++iz) {
      float t = _v[kp][iz]*(float)(_dt/_dx);
      float r = 1.0f+t;
      float s = 1.0f-t;
      up[kp][iz] = s*(ui[kp][iz]-um[ki][iz])+r*ui[ki][iz]-t*um[km][iz];
    }

    // ix = nx-1
    kp = _nx-1;
    ki = _nx-2;
    km = _nx-3;
    for (int iz=0; iz<_nz; ++iz) {
      float t = _v[kp][iz]*(float)(_dt/_dx);
      float r = 1.0f+t;
      float s = 1.0f-t;
      up[kp][iz] = s*(ui[kp][iz]-um[ki][iz])+r*ui[ki][iz]-t*um[km][iz];
    }

    // iz = 0
    kp = 0;
    ki = 1;
    km = 2;
    for (int ix=0; ix<_nx; ++ix) {
      float t = _v[ix][kp]*(float)(_dt/_dx);
      float r = 1.0f+t;
      float s = 1.0f-t;
      up[ix][kp] = s*(ui[ix][kp]-um[ix][ki])+r*ui[ix][ki]-t*um[ix][km];
    }

    // iz = nz-1
    kp = _nz-1;
    ki = _nz-2;
    km = _nz-3;
    for (int ix=0; ix<_nx; ++ix) {
      float t = _v[ix][kp]*(float)(_dt/_dx);
      float r = 1.0f+t;
      float s = 1.0f-t;
      up[ix][kp] = s*(ui[ix][kp]-um[ix][ki])+r*ui[ix][ki]-t*um[ix][km];
    }

  }

  private float[][] extendModel(float[][] c) {
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
