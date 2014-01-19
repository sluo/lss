package lss.mod;

import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

public class Receiver {

  /**
   * Residual.
   */
  public static interface Residual {

    /**
     * Computes the residual.
     * @param recp input receiver containing predicted data.
     * @param reco input receiver containing observed data.
     * @return receiver containing residuals.
     */
    public Receiver compute(Receiver rcp, Receiver rco);
  }

  ////////////////////////////////////////////////////////////////////////////
  // Constructors

  public Receiver(int xr, int zr, int nt) {
    this(new int[]{xr},new int[]{zr},nt,null);
  }

  public Receiver(int[] xr, int[] zr, int nt) {
    this(xr,zr,nt,null);
  }

  public Receiver(int[] xr, int[] zr, float[][] data) {
    this(xr,zr,-1,data);
  }

  public Receiver(Receiver rc) {
    this(copy(rc.getXIndices()),copy(rc.getZIndices()),copy(rc.getData()));
  }

  private Receiver(int[] xr, int[] zr, int nt, float[][] data) {
    Check.argument(xr.length==zr.length,"xr.length==zr.length");
    _nr = xr.length;
    _xr = xr;
    _zr = zr;
    for (int ir=0; ir<_nr; ++ir) {
      Check.argument(_xr[ir]>=0,"xr[ir]>=0");
      Check.argument(_zr[ir]>=0,"zr[ir]>=0");
    }
    if (data!=null) {
      Check.argument(data.length==_nr,"data.length==nr");
      _nt = data[0].length;
      _data = data;
    } else {
      _nt = nt;
      _data = new float[_nr][_nt];
    }
  }

  ////////////////////////////////////////////////////////////////////////////

  public void setData(float[][] ui, int it, int nabsorb) {
    for (int ir=0; ir<_nr; ++ir) {
      int xr = _xr[ir]+nabsorb;
      int zr = _zr[ir]+nabsorb;
      _data[ir][it] = ui[zr][xr];
    }
  }

  public void setData(float[][] d) {
    Check.argument(d[0].length==_nt,"d[0].length==_nt");
    Check.argument(d.length==_nr,"d.length==_nt");
    copy(d,_data);
  }

  public int[][] getIndices() {
    return new int[][]{_xr,_zr};
  }

  public int[] getXIndices() {
    return _xr;
  }

  public int[] getZIndices() {
    return _zr;
  }

  public float[][] getData() {
    return _data;
  }

  public int getNt() {
    return _nt;
  }

  public int getNr() {
    return _nr;
  }

  public Receiver clone() {
    return new Receiver(copy(_xr),copy(_zr),copy(_data));
  }

  private int _nr,_nt;
  private int[] _xr, _zr;
  private float[][] _data;
}