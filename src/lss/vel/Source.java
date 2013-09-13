package lss.vel;

import static edu.mines.jtk.util.ArrayMath.*;

public interface Source {

  public void add(float[][] ui, int it, int nabsorb);

  ////////////////////////////////////////////////////////////////////////////

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

  public static class Gaussian4Source implements Source {
    public Gaussian4Source(int xs, int zs, float dt, float fpeak) {
      _tdelay = 1.0f/fpeak;
      _fpeak = fpeak;
      _dt = dt;
      _xs = xs;
      _zs = zs;
    }
    public void add(float[][] ui, int it, int nabsorb) {
      ui[_zs+nabsorb][_xs+nabsorb] += gaussian(it*_dt-_tdelay);
    }
    private float gaussian(float t) {
      double x = PI*_fpeak*t;
      double xx = x*x;
      //return (float)((3.0-12.0*xx+4.0*xx*xx)*exp(-xx));
      return (float)((0.075-0.3*xx+0.1*xx*xx)*exp(-xx)); // 0.025 scaling
    }
    private int _xs,_zs;
    private float _dt,_fpeak,_tdelay;
  }

  public static class WaveletSource implements Source {
    public WaveletSource(int xs, int zs, float[] w) {
      _xs = xs;
      _zs = zs;
      _nt = w.length;
      _w = w;
    }
    public void add(float[][] ui, int it, int nabsorb) {
      if (it>=_nt) 
        return;
      ui[_zs+nabsorb][_xs+nabsorb] += _w[it];
    }
    private float[] _w;
    private int _xs,_zs,_nt;
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
        for (int ix=nabsorb; ix<_nx-nabsorb; ++ix)
          for (int iz=nabsorb; iz<_nz-nabsorb; ++iz)
            ui[iz][ix] += _u[it][iz][ix]*_r[iz-nabsorb][ix-nabsorb];
      }
    }
    private float[][] _r = null; //reflectivity
    private float[][][] _u;
    private int _nz,_nx;
  }


}
