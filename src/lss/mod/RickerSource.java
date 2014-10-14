package lss.mod;

import edu.mines.jtk.dsp.SincInterpolator;
import static edu.mines.jtk.util.ArrayMath.*;

public class RickerSource implements Source {

  public RickerSource(double x, double z, double fpeak) {
    _tdelay = 1.5f/(float)fpeak;
    _fpeak = (float)fpeak;
    _x = x;
    _z = z;

    SincInterpolator si = new SincInterpolator();
    _asinc = si.getTable();
    _lsinc = _asinc[0].length;
    _nsinc = _asinc.length;
    _nsincm1 = _nsinc-1;
  }

  public void setupForDomain(
    double dx, double dz, double dt, int nabsorb)
  {
    _dt = (float)dt;
    _ifx = (int)(_x/dx)-3+nabsorb;
    _ifz = (int)(_z/dz)-3+nabsorb;
    _sincx = getSinc(_x,dx);
    _sincz = getSinc(_z,dz);
    edu.mines.jtk.mosaic.SimplePlot.asPoints(_sincx);
    edu.mines.jtk.mosaic.SimplePlot.asPoints(_sincz);
  }

  public void add(float[][] ui, int it, int nabsorb) {
    float x = FLT_PI*_fpeak*((float)_dt*it-_tdelay);
    float xx = x*x;
    float r = (float)((1.0-2.0*xx)*exp(-xx));
    for (int iz=0; iz<_lsinc; ++iz) {
      for (int ix=0; ix<_lsinc; ++ix) {
        float w = _sincx[ix]*_sincz[iz];
        ui[_ifz+iz][_ifx+ix] += w*r;
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////

  private float _fpeak,_tdelay;
  private double _x,_z,_dt;
  private int _ifx,_ifz;

  private int _lsinc; // length of sinc approximations
  private int _nsinc; // number of sinc approximations
  private int _nsincm1; // _nsinc-1
  private float[][] _asinc; // array of sinc approximations
  private float[] _sincx; // sinc coefficients for x direction
  private float[] _sincz; // sinc coefficients for z direction

  // See SincInterpolator
  private float[] getSinc(double x, double dx) {
    double xscale = 1.0/dx;
    double xshift = _lsinc;
    double xn = xshift+x*xscale;
    int ixn = (int)xn;
    double fracx = xn-ixn;
    if (fracx<0.0) fracx += 1.0;
    int ksincx = (int)(fracx*_nsincm1+0.5);
    //System.out.println("xscale="+xscale);
    //System.out.println("xshift="+xshift);
    //System.out.println("xn="+xn);
    //System.out.println("ixn="+ixn);
    //System.out.println("fracx="+fracx);
    //System.out.println("ksincx="+ksincx);
    return _asinc[ksincx];
  }
}
