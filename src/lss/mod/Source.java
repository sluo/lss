package lss.mod;

import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.dsp.SincInterpolator;
import edu.mines.jtk.dsp.Sampling;

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



  /* this wavelet source injects the source in a 2D sinc
     window around the true location (which may not lie on 
     a grid point
     Graham. J. Hicks (2002). 
     ”Arbitrary source and receiver positioning in finite‐difference schemes 
     using Kaiser windowed sinc functions.” 67, SPECIAL
     SECTION—SEISMIC SIGNATURES OF FLUID TRANSPORT, 156-165.
     doi: 10.1190/1.1451454
  */
  public static class WaveletSourceSinc implements Source {
    public WaveletSourceSinc(float xs, float zs, float[] w) {
      _xs = xs;
      _ixs = (int) (xs+0.5f); // round to nearest grid point
      _zs = zs; 
      _izs = (int) (zs+0.5f); // round to nearest grid point
      _nt = w.length;
      _w = w;
      buildTable();
    }


    public void add(float[][] ui, int it, int nabsorb) {
      if (it>=_nt) 
        return;

      int n = (int) _n/2;

      /* Scatter wavelet in the surounding locations */
      for (int i2=0; i2<_n; i2++){
        int i2u = -n +i2;        
        for(int i1=0; i1<_n; i1++){
          int i1u = -n +i1;
          ui[_izs+nabsorb+i2u][_ixs+nabsorb+i1u] += _sinc_weights[i2][i1]*_w[it];
        }
      }
    }

    private void buildTable(){
      _n = 9;
      float o =(int) (-_n/2)+1;
      SincInterpolator sinc = new SincInterpolator();

      /* Here I build the input spike */
      Sampling is1 = new Sampling(_n,1,o);
      Sampling is2 = is1;
      float[][] ispike = new float[_n][_n];
      ispike[(int)_n/2][(int)_n/2] = 1.0f;


      /* Here I build the shifted spike */
      float d1 = (float) _xs - _ixs;
      float d2 = (float) _zs - _izs; 

      Sampling os1 = new Sampling(_n,1,o+d1);
      Sampling os2 = new Sampling(_n,1,o+d2);

      _sinc_weights = new float[_n][_n];

      for (int i2=0; i2< os2.getCount(); i2++){
        float z = (float) os2.getValue(i2);
        for (int i1=0; i1< os1.getCount(); i1++){
          float x = (float) os1.getValue(i1);
          _sinc_weights[i2][i1] = sinc.interpolate(is1,is2,ispike,x,z);
        }
      }
    }
    
    private float[] _w;
    private float[][] _sinc_weights;
    private float _xs,_zs;
    private int _ixs,_izs,_nt;
    private int _n;
  }








  public static class ReceiverSource implements Source {
    public ReceiverSource(Receiver receiver) {
      int[][] c = receiver.getIndices();
      _xr = c[0];
      _zr = c[1];
      _nr = receiver.getNr();
      _nt = receiver.getNt();
      _data = receiver.getData();
      _rec = receiver;
    }

    public void add(float[][] ui, int it, int nabsorb) {
      if(_rec.isSincRec()) {
        _rec.setDataRadj(ui, it, nabsorb);
      }else{ 
        addI(ui, it, nabsorb);
      }
    }

    public void addI(float[][] ui, int it, int nabsorb) {
      if (it>=_nt) 
        return;
      for (int ir=0; ir<_nr; ++ir) {
        int xs = _xr[ir]+nabsorb;
        int zs = _zr[ir]+nabsorb;
        ui[zs][xs] += _data[ir][it];
      }
    }
    private Receiver _rec;
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
