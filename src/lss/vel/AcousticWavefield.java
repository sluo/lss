package lss.vel;

import java.awt.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class AcousticWavefield {

  public enum SourceWavelet {
    RICKER
  }

  /**
   * Source function.
   */
  public interface Source {
    public void add(double t, float[][] f);
  }

  public AcousticWavefield(
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
    Source source = new AdjointSource(s,kzs,kxs);
    propagate(source,c);
    reverseIndex3(_u);
  }

  public void backPropagate(
    float[][] s, int kzs[], int[] kxs, float[][] c)
  {
    Source source = new AdjointSource(s,kzs,kxs);
    propagate(source,c);
    reverseIndex3(_u);
  }

  public float[] getWavefield(int kzr, int kxr) {
    float[] y = new float[_nt];
    kzr += _b; kxr += _b;
    for (int it=0; it<_nt; it++)
      y[it] = _u[it][kxr][kzr];
    return y;
  }

  public float[][] getWavefield(int[] kzr, int[] kxr) {
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

  private static int _b = 42; // absorbing boundary size
  private Source _source;
  private double _t; // current time
  private int _it; // current time index
  private int _nx,_nz,_nt; // number of samples
  private double _dx,_dz,_dt; // sampling intervals
  private float[][] _v; // velocity
  private float[][] _rvs; // v * v * (0.5 * dt * dt) / (dx * dx)
  private float[][] _um,_ui,_up; // u(x;t-1), u(x;t), u(x;t+1)
  private float[][][] _u; // wavefield

  private void propagate(Source source, float[][] v) {
    _source = source;
    setModel(v);
    propagate();
  }

  private void setModel(float[][] v) {
    double odx = 1.0/_dx;
    double r = 0.5*_dt*_dt*odx*odx;
    _rvs = extendModel(v);
    _v = copy(_rvs);
    mul(_rvs,_rvs,_rvs);
    mul((float)r,_rvs,_rvs);
  }

  private class RickerSource implements Source {
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

  private class AdjointSource implements Source {
    public AdjointSource(float[] s1, int kzs, int kxs) {
      _kzs = kzs+_b;
      _kxs = kxs+_b;
      _s1 = s1;
      _ns = 1;
    }
    public AdjointSource(float[][] s2, int[] kzsa, int[] kxsa) {
      Check.argument(kzsa.length==kxsa.length,"kzs.length=kxs.length");
      _ns = kzsa.length;
      _kzsa = copy(kzsa);
      _kxsa = copy(kxsa);
      for (int is=0; is<_ns; is++) {
        _kzsa[is] += _b;
        _kxsa[is] += _b;
      }
      _s2 = s2;
    }
    public void add(double t, float[][] f) {
      int it = (int)(t/_dt);
      if (_kzsa==null && _kxsa==null) {
        float a = (it<0)?0.0f:_s1[_nt-1-it];
        f[_kxs][_kzs] += a;
      } else if (_kzsa!=null && _kxsa!=null) {
        for (int is=0; is<_ns; is++) {
          int kzs = _kzsa[is];
          int kxs = _kxsa[is];
          float a = (it<0)?0.0f:_s2[is][_nt-1-it];
          f[kxs][kzs] += a;
        }
      }
    }
    private int _ns,_kzs,_kxs;
    private int[] _kzsa = null;
    private int[] _kxsa = null;
    private float[] _s1 = null;
    private float[][] _s2 = null;
  }

  private void propagate() {
    _um = _u[0];
    _ui = _u[1];
    _up = _u[2];
    zero(_um);
    zero(_ui);
    zero(_up);
    for (_it=0; _it<_nt-3; ++_it) {
      System.out.format("\r%d",_it);
      step();
    }
    System.out.print("\n");
  }
  
  private void step() {
    Parallel.loop(1,_nx-1,new Parallel.LoopInt() {
    public void compute(int ix) {
      stepSliceX(ix);
    }});
    _source.add(_it*_dt,_up);
    absorb();
    _um = _u[_it+1];
    _ui = _u[_it+2];
    _up = _u[_it+3];
  }

  // 4th-order stencil coefficients for 2nd-order derivative
  private static final float C0 = -2.5000000000000000f;
  private static final float C1 =  1.3333333730697632f;
  private static final float C2 = -0.0833333358168602f;

  private void stepSliceX(int ix) {

    // At second and second-to-last sample, use 2nd-order stencil.
    if (ix==1 || ix==_nx-2) {
      int ixm1 = ix-1, ixp1 = ix+1;
      for (int iz=1; iz<_nz-1; ++iz) {
        int izm1 = iz-1, izp1 = iz+1;
        float rvsi = _rvs[ix][iz];
        _up[ix][iz] = (2.0f-6.0f*rvsi)*_ui[ix][iz]-_um[ix][iz]+
          rvsi*(((_ui[ixp1][iz  ]+_ui[ix  ][izp1]+
                  _ui[ixm1][iz  ]+_ui[ix  ][izm1])+0.5f*
                 (_ui[ixp1][izp1]+_ui[ixm1][izp1]+
                  _ui[ixm1][izm1]+_ui[ixp1][izm1])));
      }
    }

    // In the interior, use (mostly) 4th-order stencil.
    else {
      int ixm1 = ix-1, ixp1 = ix+1;
      int ixm2 = ix-2, ixp2 = ix+2;

      //  iz = 1
      {
        int iz = 1;
        int izm1 = iz-1, izp1 = iz+1;
        float rvsi = _rvs[ix][iz];
        _up[ix][iz] = (2.0f-6.0f*rvsi)*_ui[ix][iz]-_um[ix][iz]+
          rvsi*(((_ui[ixp1][iz  ]+_ui[ix  ][izp1]+
                  _ui[ixm1][iz  ]+_ui[ix  ][izm1])+0.5f*
                 (_ui[ixp1][izp1]+_ui[ixm1][izp1]+
                  _ui[ixm1][izm1]+_ui[ixp1][izm1])));
      }

      // 1 < iz < nz-2
      for (int iz=2; iz<_nz-2; ++iz) {
        int izm1 = iz-1, izp1 = iz+1;
        int izm2 = iz-2, izp2 = iz+2;
        float rvsi = _rvs[ix][iz];
        _up[ix][iz] = (2.0f+3.0f*C0*rvsi)*_ui[ix][iz]-_um[ix][iz]+
          rvsi*(C1*((_ui[ixp1][iz  ]+_ui[ix  ][izp1]+
                    _ui[ixm1][iz  ]+_ui[ix  ][izm1])+0.5f*
                   (_ui[ixp1][izp1]+_ui[ixm1][izp1]+
                    _ui[ixm1][izm1]+_ui[ixp1][izm1]))+
               C2*((_ui[ixp2][iz  ]+_ui[ix  ][izp2]+
                    _ui[ixm2][iz  ]+_ui[ix  ][izm2])+0.5f*
                   (_ui[ixp2][izp2]+_ui[ixm2][izp2]+
                    _ui[ixm2][izm2]+_ui[ixp2][izm2])));
      }

      // iz = nz-2
      {
        int iz = _nz-2;
        int izm1 = iz-1, izp1 = iz+1;
        float rvsi = _rvs[ix][iz];
        _up[ix][iz] = (2.0f-6.0f*rvsi)*_ui[ix][iz]-_um[ix][iz]+
          rvsi*(((_ui[ixp1][iz  ]+_ui[ix  ][izp1]+
                  _ui[ixm1][iz  ]+_ui[ix  ][izm1])+0.5f*
                 (_ui[ixp1][izp1]+_ui[ixm1][izp1]+
                  _ui[ixm1][izm1]+_ui[ixp1][izm1])));
      }
    }
  }

  private void absorb() {
    //absorbReynolds(_um,_ui,_up);
    Parallel.loop(0,_nx,new Parallel.LoopInt() {
      public void compute(int ix) {
        absorbSliceX(ix,_b,_ui);
        absorbSliceX(ix,_b,_up);
      }
    });
    absorbClayton();
  }

  private static void absorbSliceX(int ix ,int b, float[][] u) {
    int nz = u[0].length;
    int nx = u.length;
    int nzmb = nz-b;
    float w = 0.005f;
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

  private void reverseIndex3(float[][][] f) {
    int n3 = f.length;
    float[][][] ft = new float[n3][][];
    for (int i3=0; i3<n3; i3++)
      ft[i3] = f[n3-1-i3];
    for (int i3=0; i3<n3; i3++)
      f[i3] = ft[i3];
  }

}
