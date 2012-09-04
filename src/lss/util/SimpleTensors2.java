package lss.util;

import edu.mines.jtk.dsp.Eigen;
import edu.mines.jtk.dsp.EigenTensors2;
import edu.mines.jtk.dsp.Tensors2;
import edu.mines.jtk.util.Parallel;

/**
 * An array of tensors.
 * Each tensor is a symmetric positive-semidefinite 2-by-2 matrix:
 * <pre><code>
 * A = |a11 a12|
 *     |a12 a22|
 * </code></pre>
 * @author Simon Luo
 * @version 2012.02.01
 */
public class SimpleTensors2 implements Tensors2 {

  /**
   * Constructs tensors.
   * @param a array of {a11,a12,a22} of tensor elements.
   */
  public SimpleTensors2(float[][][] a) {
    _a = a;
  }

  /**
   * Gets tensor elements for specified indices.
   * @param i1 index for 1st dimesion.
   * @param i2 index for 2nd dimesion.
   * @param a array {a11,a12,a22} of tensor elements.
   */
  public void getTensor(int i1, int i2, float[] a) {
    a[0] = _a[0][i2][i1];
    a[1] = _a[1][i2][i1];
    a[2] = _a[2][i2][i1];
  }

  /**
   * Gets all tensor elements.
   * @param a11 tensor elements a11.
   * @param a12 tensor elements a12.
   * @param a22 tensor elements a22.
   */
  public void getTensors(float[][] a11, float[][] a12, float[][] a22) {
    edu.mines.jtk.util.ArrayMath.copy(_a[0],a11);
    edu.mines.jtk.util.ArrayMath.copy(_a[1],a12);
    edu.mines.jtk.util.ArrayMath.copy(_a[2],a22);
  }

  /**
   * Gets the eigen-decomposition of tensors.
   * @return eigen-decomposition.
   */
  public EigenTensors2 asEigenTensors() {
    if (_et==null)
      _et = getEigenTensors(_a);
    return _et;
  }

  public static EigenTensors2 getEigenTensors(float[][][] a) {
    int n1 = a[0][0].length;
    int n2 = a[0].length;
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    float[][] eu = new float[n2][n1];
    float[][] ev = new float[n2][n1];
    solveEigenproblems(a,u1,u2,eu,ev);
    return new EigenTensors2(u1,u2,eu,ev);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  
  private float[][][] _a;
  private EigenTensors2 _et = null;

  private static void solveEigenproblems(
    float[][][] as,
    float[][] u1, final float[][] u2,
    float[][] eu, final float[][] ev)
  {
    int n1 = as[0][0].length;
    int n2 = as[0].length;
    double[][] a = new double[2][2];
    double[][] z = new double[2][2];
    double[] e = new double[2];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        a[0][0] = as[0][i2][i1];
        a[0][1] = as[1][i2][i1];
        a[1][0] = as[1][i2][i1];
        a[1][1] = as[2][i2][i1];
        Eigen.solveSymmetric22(a,z,e);
        float u1i = (float)z[0][0];
        float u2i = (float)z[0][1];
        if (u1i<0.0f) {
          u1i = -u1i;
          u2i = -u2i;
        }
        float eui = (float)e[0];
        float evi = (float)e[1];
        if (evi<0.0f) evi = 0.0f;
        if (eui<evi) eui = evi;
        u1[i2][i1] = u1i;
        u2[i2][i1] = u2i;
        eu[i2][i1] = eui;
        ev[i2][i1] = evi;
      }
    }
  }
  
}
