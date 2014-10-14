package lss.mod;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

public class SincReceiver extends Receiver {

  public SincReceiver(double xr, double zr, int nt) {
    this(new double[]{xr},new double[]{zr},nt);
  }

  public SincReceiver(double[] xr, double[] zr, int nt) {
    super(0,0,nt); // TODO: get rid of this eventually
    Check.argument(xr.length==zr.length,"xr.length==zr.length");
    _xr = xr;
    _zr = zr;
    _nr = xr.length;
    _nt = nt;
    _data = new float[_nr][_nt];
    _si = new SincInterpolator();
  }

  public void setupForDomain(
    int nx, double dx, int nz, double dz, int nt, double dt, int nabsorb)
  {
    //_nx = nx;
    //_dx = dx;
    //_nz = nz;
    //_dz = dz;
    //_fx = -_dx*nabsorb;
    //_fz = -_dz*nabsorb;
    _nt = nt;
    _dt = dt;
    _sx = new Sampling(nx+2*nabsorb,dx,-dx*nabsorb);
    _sz = new Sampling(nz+2*nabsorb,dz,-dz*nabsorb);
  }

  @Override
  public void setData(float[][] ui, int it, int nabsorb) {
    setData(it,ui);
  }

  public void setData(int it, float[][] ui) {
    for (int ir=0; ir<_nr; ++ir) {
      double xr = _xr[ir];
      double zr = _zr[ir];
      this._data[ir][it] = _si.interpolate(_sx,_sz,ui,xr,zr);
    }
  }

  public float[][] getData() {
    return this._data;
  }

  ////////////////////////////////////////////////////////////////////////////

  //private int _nx,_nz,_nt,_nr;
  //private double _dx,_dz;
  //private double _fx,_fz;

  private int _nr;
  private int _nt;
  private double _dt;
  private Sampling _sx,_sz;

  private double[] _xr, _zr;
  private float[][] _data;
  private SincInterpolator _si;
}
