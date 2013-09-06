package lss.vel;

import edu.mines.jtk.util.Check;

public class Receiver {

  public Receiver(int xr, int zr, int nt) {
    this(new int[]{xr},new int[]{zr},nt);
  }

  public Receiver(int[] xr, int[] zr, int nt) {
    int nxr = xr.length;
    int nzr = zr.length;
    Check.argument(nxr==nzr,"nxr==nzr");
    _nr = nxr;
    _nt = nt;
    _xr = xr;
    _zr = zr;
    _data = new float[_nr][nt];
  }

  public Receiver(int[] xr, int[] zr, float[][] data) {
    int nxr = xr.length;
    int nzr = zr.length;
    Check.argument(nxr==nzr,"nxr==nzr");
    _nr = nxr;
    Check.argument(data.length==_nr,"data.length==nr");
    _nt = data[0].length;
    _xr = xr;
    _zr = zr;
    _data = data;
  }

  public void setData(float[][] ui, int it, int nabsorb) {
    for (int ir=0; ir<_nr; ++ir) {
      int xr = _xr[ir]+nabsorb;
      int zr = _zr[ir]+nabsorb;
      _data[ir][it] = ui[zr][xr];
    }
  }

  public int[][] getIndices() {
    return new int[][]{_xr,_zr};
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

  private int _nr,_nt;
  private int[] _xr, _zr;
  private float[][] _data;
}
