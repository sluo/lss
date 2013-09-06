package lss.util;

import edu.mines.jtk.util.Check;

/**
 * A 4-D array of floats in which the 3-D arrays of floats
 * indexed by the fourth dimension are shared.
 * Values of index i4 that are greated than n4 will return a 3-D array
 * corresponding to a lower index value. For example, for n4=10, i4=11
 * will return the array corresponding to i4=1.
 * @author Simon Luo
 * @version 2013.09.06
 */
public class SharedFloat4 {

  public SharedFloat4(int n1, int n2, int n3, int n4) {
    this(new float[n4][n3][n2][n1]);
  }

  public SharedFloat4(float[][][][] f) {
    _n4 = f.length;
    _f = f;
  }

  public int getN4() {
    return _n4;
  }

  public float[][][] get(int i4) {
    while (i4>=_n4) 
      i4 -= _n4;
    return _f[i4];

  }

  public float[][][][] get() {
    return _f;
  }

  private int _n4;
  private float[][][][] _f;
}
