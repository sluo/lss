package lss.mod;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

public abstract class SincSource implements Source {

  /**
   * Implementation of the Source interface.
   */
  public void add(float[][] ui, int it, int nabsorb) {
    add(it,ui);
  }

  /**
   * Injects the source.
   * @param it current time step.
   * @param ui current time slice of the wavefield.
   */
  abstract void add(int it, float[][] ui);

  /**
   * Sets the parameters describing the modeling domain.
   * @param dx spatial sampling interval in the x-direction.
   * @param dz spatial sampling interval in the z-direction.
   * @param dt time sampling interval.
   * @param nabsorb size of absorbing boundary.
   */
  abstract void setupForDomain(
    int nx, double dx, int nz, double dz, int nt, double dt, int nabsorb);

  ///////////////////////////////////////////////////////////////////////////
  // Point sources

  /**
   * A Ricker wavelet.
   */
  public static class RickerWavelet extends PointSource {
    public RickerWavelet(double xs, double zs, double fpeak) {
      super(xs,zs);
      _tdelay = 1.5f/(float)fpeak;
      _fpeak = (float)fpeak;
    }
    public void add(int it, float[][] ui) {
      float x = FLT_PI*_fpeak*((float)_dt*it-_tdelay);
      float xx = x*x;
      float r = (float)((1.0-2.0*xx)*exp(-xx));
      for (int iz=0; iz<_lsinc; ++iz) {
        for (int ix=0; ix<_lsinc; ++ix) {
          float ss = _sincx[ix]*_sincz[iz];
          ui[_ifz+iz][_ifx+ix] += ss*r;
        }
      }
    }
    private float _fpeak,_tdelay;
  }

  /**
   * An arbitrary wavelet.
   */
  public class Wavelet extends PointSource {
    public Wavelet(double xs, double zs, double dt, double ft, float[] w) {
      super(xs,zs);
      _w = w;
      _sw = new Sampling(w.length,_dt,ft);
      _si = new SincInterpolator();
    }
    public void add(int it, float[][] ui) {
      float r = _si.interpolate(_sw,_w,_dt*it);
      for (int iz=0; iz<_lsinc; ++iz) {
        for (int ix=0; ix<_lsinc; ++ix) {
          float ss = _sincx[ix]*_sincz[iz];
          ui[_ifz+iz][_ifx+ix] += ss*r;
        }
      }
    }
    private float[] _w; // array containint the wavelet
    private Sampling _sw; // time sampling for wavelet
    private SincInterpolator _si;
  }

  /**
   * A point source is injected at a single location.
   */
  static abstract class PointSource extends SincSource {
    public PointSource(double xs, double zs) {
      _xs = xs;
      _zs = zs;
      _asinc = new SincInterpolator().getTable();
      _lsinc = _asinc[0].length;
      _nsinc = _asinc.length;
    }
    public void setupForDomain(
      int nx, double dx, int nz, double dz, int nt, double dt, int nabsorb)
    {
      _dt = (float)dt;
      _ifx = (int)(_xs/dx)-3+nabsorb; // first samples that have
      _ifz = (int)(_zs/dz)-3+nabsorb; // nonzero sinc coefficients
      _sincx = getSinc(_xs,dx,_asinc);
      _sincz = getSinc(_zs,dz,_asinc);
      //edu.mines.jtk.mosaic.SimplePlot.asPoints(_sincx);
      //edu.mines.jtk.mosaic.SimplePlot.asPoints(_sincz);
    }
    double _xs,_zs,_dt;
    int _ifx,_ifz;
    int _lsinc; // length of sinc approximations
    int _nsinc; // number of sinc approximations
    float[][] _asinc; // array of sinc approximations
    float[] _sincx; // sinc coefficients for x direction
    float[] _sincz; // sinc coefficients for z direction
  }

  ///////////////////////////////////////////////////////////////////////////

  /**
   * Gets the sinc approximation for the given position
   * relative to the nearest sample.
   * @param x the position.
   * @param dx the sampling interval.
   * @param asinc precomputed table of sinc coefficients.
   */
  private static float[] getSinc(double x, double dx, float[][] asinc) {
    int lsinc = asinc[0].length;
    int nsinc = asinc.length;
    int nsincm1 = nsinc-1;
    double xscale = 1.0/dx;
    double xshift = lsinc;
    double xn = xshift+x*xscale;
    int ixn = (int)xn;
    double frac = xn-ixn;
    if (frac<0.0) frac += 1.0;
    int ksinc = (int)(frac*nsincm1+0.5);
    //System.out.println("xscale="+xscale);
    //System.out.println("xshift="+xshift);
    //System.out.println("xn="+xn);
    //System.out.println("ixn="+ixn);
    //System.out.println("frac="+frac);
    //System.out.println("ksinc="+ksinc);
    return asinc[ksinc];
  }

}
