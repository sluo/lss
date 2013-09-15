package lss.vel;

import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.ArrayMath;

public class Receiver {

  public Receiver(int xr, int zr, int nt) {
    this(new int[]{xr},new int[]{zr},nt,null);
  }

  public Receiver(int[] xr, int[] zr, int nt) {
    this(xr,zr,nt,null);
  }

  public Receiver(int[] xr, int[] zr, float[][] data) {
    this(xr,zr,-1,data);
  }

  private Receiver(int[] xr, int[] zr, int nt, float[][] data) {
    int nxr = xr.length;
    int nzr = zr.length;
    Check.argument(nxr==nzr,"nxr==nzr");
    _nr = nxr;
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
    ArrayMath.copy(d,_data);
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

  public Receiver copy() {
    return new Receiver(
      ArrayMath.copy(_xr),ArrayMath.copy(_zr),ArrayMath.copy(_data));
  }

  private int _nr,_nt;
  private int[] _xr, _zr;
  private float[][] _data;
}
