package lss.vel;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class Util {

  public static void sdiv(
  final float num, final float[][] x, final float[][] y) {
    Check.argument(num!=0.0f,"num!=0");
    final int n1 = x[0].length;
    final int n2 = x.length;
    final Almost almost = Almost.FLOAT;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        double den = x[i2][i1];
        double sign = 1.0;
        if ((num>0.0 && den<0.0) || (num<0.0 && den>0.0))
          sign = -1.0;
        double q = almost.divide(num,den,0.0);
        if (q==sign*0.01*Float.MAX_VALUE)
          y[i2][i1] = 0.0f;
        else
          y[i2][i1] = (float)q;
      }
    }});
  }

  public static float[][][][] zeros(int n1, int n2, int n3, int n4) {
    return new float[n4][n3][n2][n1];
  }

  public static float[] like(float[] f) {
    return new float[f.length];
  }
  public static float[][] like(float[][] f) {
    return new float[f.length][f[0].length];
  }
  public static float[][][] like(float[][][] f) {
    return new float[f.length][f[0].length][f[0][0].length];
  }

}
