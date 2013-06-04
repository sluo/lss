package lss.vel;

public class Util {

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
