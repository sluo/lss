package lss.mod;

import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.dsp.SincInterpolator;
import edu.mines.jtk.dsp.Sampling;

public abstract class ArbitrarySource implements Source {

  public ArbitrarySource() {
    SincInterpolator si = new SincInterpolator();
    _asinc = si.getTable();
    _lsinc = _asinc[0].length;
    _nsinc = _asinc.length;
    _nsincm1 = _nsinc-1;
  }

  // See SincInterpolator
  public float[] getSinc(double x, double dx) {
    double xn = x*dx;
    int ixn = (int)xn;
    double fracx = xn-(int)ixn;
    if (fracx<0.0) fracx += 1.0;
    int ksincx = (int)(fracx*_nsincm1+0.5);
    return _asinc[ksincx];
  }

  abstract void setupForDomain(double dx, double dz, double dt, int nabsorb);

  private int _lsinc; // length of sinc approximations
  private int _nsinc; // number of sinc approximations
  private int _nsincm1; // _nsinc-1
  private float[][] _asinc; // array of sinc approximations
  double _dx,_dz,_dt;

}
