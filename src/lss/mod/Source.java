package lss.mod;

import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

public interface Source {

  public void add(float[][] ui, int it, int nabsorb);

  ////////////////////////////////////////////////////////////////////////////
  // implementation

  static abstract class PointSource implements Source {
    public PointSource(int xs, int zs, float dt, float fpeak) {
      _tdelay = 1.5f/fpeak;
      _fpeak = fpeak;
      _dt = dt;
      _xs = xs;
      _zs = zs;
    }
    public void add(float[][] ui, int it, int nabsorb) {
      ui[_zs+nabsorb][_xs+nabsorb] += function(it*_dt-_tdelay);
    }
    abstract float function(float t);
    float _dt,_fpeak,_tdelay;
    int _xs,_zs;
  }

  public static class RickerSource extends PointSource {
    public RickerSource(int xs, int zs, float dt, float fpeak) {
      super(xs,zs,dt,fpeak);
    }
    float function(float t) {
      float x = FLT_PI*_fpeak*t;
      float xx = x*x;
      return (float)((1.0-2.0*xx)*exp(-xx));
    }
  }

  public static class Gaussian4Source extends PointSource {
    public Gaussian4Source(int xs, int zs, float dt, float fpeak) {
      super(xs,zs,dt,fpeak);
    }
    float function(float t) {
      double x = PI*_fpeak*t;
      double xx = x*x;
      //return (float)((3.0-12.0*xx+4.0*xx*xx)*exp(-xx));
      return (float)((0.075-0.3*xx+0.1*xx*xx)*exp(-xx)); // 0.025 scaling
    }
  }

  public static class SimultaneousSource implements Source {
    public SimultaneousSource(Source[] source) {
      _ns = source.length;
      _source = source;
    }
    public void add(final float[][] ui, final int it, final int nabsorb) {
      for (int isou=0; isou<_ns; ++isou) {
        _source[isou].add(ui,it,nabsorb);
      }
    }
    private int _ns;
    private Source[] _source;
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
    public void add(final float[][] ui, final int it, final int nabsorb) {
      if (_r==null) {
        Parallel.loop(_nz,new Parallel.LoopInt() {
        public void compute(int iz) {
          for (int ix=0; ix<_nx; ++ix) {
            ui[iz][ix] += _u[it][iz][ix];
          }
        }});
      } else {
        Parallel.loop(nabsorb,_nz-nabsorb,new Parallel.LoopInt() {
        public void compute(int iz) {
          for (int ix=nabsorb; ix<_nx-nabsorb; ++ix) {
            ui[iz][ix] += _u[it][iz][ix]*_r[iz-nabsorb][ix-nabsorb];
          }
        }});
      }
    }
    private float[][] _r = null; //reflectivity
    private float[][][] _u;
    private int _nz,_nx;
  }

}
