package lss.util;

import edu.mines.jtk.dsp.Eigen;
import edu.mines.jtk.dsp.EigenTensors3;
import edu.mines.jtk.dsp.Tensors3;
import edu.mines.jtk.util.ArrayMath;
import edu.mines.jtk.util.Parallel;

/**
 * An array of tensors.
 * Each tensor is a symmetric positive-semidefinite 3-by-3 matrix:
 * <pre><code>
 *     |a11 a12 a13|
 * A = |a12 a22 a23|
 *     |a13 a23 a33|
 * </code></pre>
 * @author Simon Luo
 * @version 2012.02.01
 */
public class SimpleTensors3 implements Tensors3 {

  /**
   * Constructs tensors.
   * @param a array of {a11,a12,a13,a22,a23,a33} of tensor elements.
   */
  public SimpleTensors3(float[][][][] a) {
    _a = a;
  }

  /**
   * Gets tensor elements for specified indices.
   * @param i1 index for 1st dimesion.
   * @param i2 index for 2nd dimesion.
   * @param i3 index for 3rd dimesion.
   * @param a array {a11,a12,a13,a22,a23,a33} of tensor elements.
   */
  public void getTensor(int i1, int i2, int i3, float[] a) {
    a[0] = _a[0][i3][i2][i1];
    a[1] = _a[1][i3][i2][i1];
    a[2] = _a[2][i3][i2][i1];
    a[3] = _a[3][i3][i2][i1];
    a[4] = _a[4][i3][i2][i1];
    a[5] = _a[5][i3][i2][i1];
  }

  /**
   * Gets all tensor elements.
   * @param a11 tensor elements a11.
   * @param a12 tensor elements a12.
   * @param a13 tensor elements a22.
   * @param a22 tensor elements a22.
   * @param a23 tensor elements a22.
   * @param a33 tensor elements a22.
   */
  public void getTensors(
    float[][][] a11, float[][][] a12, float[][][] a13,
    float[][][] a22, float[][][] a23, float[][][] a33)
  {
    if (a11!=null) ArrayMath.copy(_a[0],a11);
    if (a12!=null) ArrayMath.copy(_a[1],a12);
    if (a13!=null) ArrayMath.copy(_a[2],a13);
    if (a22!=null) ArrayMath.copy(_a[3],a22);
    if (a23!=null) ArrayMath.copy(_a[4],a23);
    if (a33!=null) ArrayMath.copy(_a[5],a33);
  }

  /**
   * Gets the eigen-decomposition of tensors.
   * @return eigen-decomposition.
   */
  public EigenTensors3 asEigenTensors() {
    if (_et==null)
      _et = getEigenTensors(_a);
    return _et;
  }

  public static EigenTensors3 getEigenTensors(float[][][][] a) {
    int n1 = a[0][0][0].length;
    int n2 = a[0][0].length;
    int n3 = a[0].length;
    float[][][] u1 = new float[n3][n2][n1];
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] w1 = new float[n3][n2][n1];
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] eu = new float[n3][n2][n1];
    float[][][] ev = new float[n3][n2][n1];
    float[][][] ew = new float[n3][n2][n1];
    solveEigenproblems(a,u1,u2,w1,w2,eu,ev,ew);
    return new EigenTensors3(u1,u2,w1,w2,eu,ev,ew,false);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  
  private float[][][][] _a;
  private EigenTensors3 _et = null;

  private static void solveEigenproblems(
    final float[][][][] as,
    final float[][][] u1, final float[][][] u2,
    final float[][][] w1, final float[][][] w2,
    final float[][][] eu, final float[][][] ev, final float[][][] ew)
  {
    final int n1 = as[0][0][0].length;
    final int n2 = as[0][0].length;
    final int n3 = as[0].length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      double[][] a = new double[3][3];
      double[][] z = new double[3][3];
      double[] e = new double[3];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          a[0][0] = as[0][i3][i2][i1];
          a[0][1] = as[1][i3][i2][i1];
          a[0][2] = as[2][i3][i2][i1];
          a[1][0] = as[1][i3][i2][i1];
          a[1][1] = as[3][i3][i2][i1];
          a[1][2] = as[4][i3][i2][i1];
          a[2][0] = as[2][i3][i2][i1];
          a[2][1] = as[4][i3][i2][i1];
          a[2][2] = as[5][i3][i2][i1];
          Eigen.solveSymmetric33(a,z,e);
          float u1i = (float)z[0][0];
          float u2i = (float)z[0][1];
          float u3i = (float)z[0][2];
          float w1i = (float)z[2][0];
          float w2i = (float)z[2][1];
          float w3i = (float)z[2][2];
          if (u1i<0.0f) {
            u1i = -u1i;
            u2i = -u2i;
            u3i = -u3i;
          }
          if (w3i<0.0f) {
            w1i = -w1i;
            w2i = -w2i;
            w3i = -w3i;
          }
          float eui = (float)e[0];
          float evi = (float)e[1];
          float ewi = (float)e[2];
          if (ewi<0.0f) ewi = 0.0f;
          if (evi<ewi) evi = ewi;
          if (eui<evi) eui = evi;
          u1[i3][i2][i1] = u1i;
          u2[i3][i2][i1] = u2i;
          w1[i3][i2][i1] = w1i;
          w2[i3][i2][i1] = w2i;
          eu[i3][i2][i1] = eui;
          ev[i3][i2][i1] = evi;
          ew[i3][i2][i1] = ewi;
        }
      }
    }});
  }
  
}
