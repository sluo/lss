package lss.util;

import edu.mines.jtk.dsp.Sampling;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;

// FOR TESTING
import java.util.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;

/**
 * Gridding by nearest-neighbor interpolation or discrete Sibson 
 * interpolation of scattered samples f(x1,x2).
 * @author Simon Luo
 * @version 2011.12.08
 */
public class Grid2 {

  /**
   * Constructs a gridder.
   * @param f sample values f(x1,x2).
   * @param x1 sample coordinates x1
   * @param x2 sample coordinates x2
   */
  public Grid2(
    float[] f, float[] x1, float[] x2)
  {
    _f = f;
    _x = new float[][]{x1,x2};
    _tree = new KdTree(_x);
  }

  /**
   * Computes gridded sample values using 
   * discrete Sibson interpolation.
   * @param s1 sampling of x1.
   * @param s2 sampling of x2.
   * @return array of gridded sample values.
   */
  public float[][] gridSibson(Sampling s1, Sampling s2) {
    return _parallel?gridSibsonP(s1,s2):gridSibsonS(s1,s2);
  }

  /**
   * Computes gridded sample values using 
   * nearest-neighbor interpolation.
   * @param s1 sampling of x1.
   * @param s2 sampling of x2.
   * @return array of gridded sample values.
   */
  public float[][] gridNearest(Sampling s1, Sampling s2) {
    return _parallel?gridNearestP(false,s1,s2):gridNearestS(false,s1,s2);
  }

  /**
   * Computes gridded sample values using a
   * naive nearest-neighbor interpolation.
   * @param s1 sampling of x1.
   * @param s2 sampling of x2.
   * @return array of gridded sample values.
   */
  public float[][] gridNearestSlow(Sampling s1, Sampling s2) {
    return gridNearestP(true,s1,s2);
  }

  /**
   * Enables or disables or parallel processing.
   * @param parallel true for parallel processing
   */
  public static void setParallel(boolean parallel) {
    _parallel = parallel;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final int NREDUCE = 8;
  private static boolean _parallel = true;

  private float[] _f;   // values
  private float[][] _x; // coordinates
  private KdTree _tree;

  private float[][] gridNearestS(
    boolean slow, Sampling s1, Sampling s2)
  {
    int n1 = s1.getCount(); 
    int n2 = s2.getCount();
    final float[][] y = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      gridNearestSlice2(i2,y,slow,s1,s2);
    }
    return y;
  }

  private float[][] gridNearestP(
    final boolean slow, final Sampling s1, final Sampling s2)
  {
    final int n1 = s1.getCount(); 
    final int n2 = s2.getCount();
    final float[][] y = new float[n2][n1];
    loop(n2,new LoopInt() {
    public void compute(int i2) {
      gridNearestSlice2(i2,y,slow,s1,s2);
    }});
    return y;
  }

  private void gridNearestSlice2(
    int i2, float[][] y,
    final boolean slow, Sampling s1, Sampling s2)
  {
    int n1 = s1.getCount();
    for (int i1=0; i1<n1; ++i1) {
      float v1 = (float)s1.getValue(i1);
      float v2 = (float)s2.getValue(i2);
      float[] xq = new float[]{v1,v2};
      int nearest = (slow)?_tree.findNearestSlow(xq):
                           _tree.findNearest(xq);
      y[i2][i1] = _f[nearest];
    }
  }

  private float[][] gridSibsonS(Sampling s1, Sampling s2) {
    int n2 = s2.getCount();
    float[][][] yc = accumulateSibson(0,n2,s1,s2);
    return div(yc[0],yc[1]);
  }

  public float[][] gridSibsonP(final Sampling s1, final Sampling s2) {
    final int n2 = s2.getCount();
    final int k2 = n2/NREDUCE;
    float[][][] yc = reduce(NREDUCE,new ReduceInt<float[][][]>() {
      public float[][][] compute(int i) {
        int min2 = i*k2;
        int max2 = (i<NREDUCE-1)?min2+k2:n2;
        return accumulateSibson(min2,max2,s1,s2);
      }
      public float[][][] combine(float[][][] yc1, float[][][] yc2) {
        add(yc1[0],yc2[0],yc2[0]);
        add(yc1[1],yc2[1],yc2[1]);
        return yc2;
      }
    });
    return div(yc[0],yc[1]);
  }

  private float[][][] accumulateSibson(
    int p2, int q2, Sampling s1, Sampling s2)
  {
    int n1 = s1.getCount(); 
    int n2 = s2.getCount();
    double d1 = s1.getDelta();
    double d2 = s2.getDelta();
    double f1 = s1.getFirst();
    double f2 = s2.getFirst();
    float[] xq = new float[2];
    float[][] c = new float[n2][n1];
    float[][] y = new float[n2][n1];
    for (int i2=p2; i2<q2; ++i2) {
      float v2 = (float)s2.getValue(i2);
      xq[1] = v2;
      for (int i1=0; i1<n1; ++i1) {
        float v1 = (float)s1.getValue(i1);
        xq[0] = v1;

        // Find distance to nearest neighbor
        int nearest = _tree.findNearest(xq);
        float rs = distance(nearest,_x,xq);
        float f = _f[nearest];

        // Fill circle and increment counter
        float r = sqrt(rs);
        int min1 = s1.indexOfNearest(v1-r);
        int max1 = s1.indexOfNearest(v1+r);
        int min2 = s2.indexOfNearest(v2-r);
        int max2 = s2.indexOfNearest(v2+r);
        for (int j2=min2; j2<=max2; ++j2) {
          float w2 = (float)s2.getValue(j2);
          float r2 = w2-v2;
          float r2s = r2*r2;
          for (int j1=min1; j1<=max1; ++j1) {
            float w1 = (float)s1.getValue(j1);
            float r1 = w1-v1;
            float r1s = r1*r1;
            if (r1s+r2s<rs) {
              y[j2][j1] += f;
              c[j2][j1] += 1.0f;
            }
          }
        }
      }
    }
    return new float[][][]{y,c};
  }

  private static float distance(int i, float[][] x, float[] xq) {
    int k = x.length;
    double ds = 0.0f;
    for (int ik=0; ik<k; ++ik) {
      double d = x[ik][i]-xq[ik];
      ds += d*d;
    }
    return (float)ds;
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  private static boolean _png = false;

  public static void main(String[] args) {
    goGrid();
  }

  private static void goGrid() {
    int ns = 1000;
    int n1 = 1001;
    int n2 = 1001;
    Sampling s1 = new Sampling(n1,1.0,0.0);
    Sampling s2 = new Sampling(n2,1.0,0.0);
    Stopwatch s = new Stopwatch();

    Gaussian gaussian = new Gaussian(n1/5.0,n1,n2);
    float[][] g = gaussian.toArray();
    float[] f = new float[ns];
    float[] x1 = new float[ns];
    float[] x2 = new float[ns];
    Random random = new Random();
    for (int i=0; i<ns; ++i) {
      float f1 = (n1-1)*random.nextFloat();
      float f2 = (n2-1)*random.nextFloat();
      f[i] = gaussian.evaluate(f1,f2);
      x1[i] = f2;
      x2[i] = f1;
    }
    plot(g,x1,x2,"samples");

    Grid2 gridder = new Grid2(f,x1,x2);

    /*
    // naive nearest neighbor
    System.out.print("nearest neighbor (naive): ");
    s.restart();
    float[][] naive = gridder.gridNearestSlow(s1,s2);
    s.stop();
    System.out.format("%.4fs\n",s.time());
    plot(naive,"naivenn");
    */

    /*
    // serial nearest neighbor
    System.out.print("nearest neighbor (serial): ");
    gridder.setParallel(false);
    s.restart();
    float[][] nears = gridder.gridNearest(s1,s2);
    s.stop();
    System.out.format("%.4fs\n",s.time());
    plot(nears,"serialnn");
    */

    /*
    */
    // parallel nearest neighbor
    System.out.print("nearest neighbor (parallel): ");
    gridder.setParallel(true);
    s.restart();
    float[][] nearp = gridder.gridNearest(s1,s2);
    s.stop();
    System.out.format("%.4fs\n",s.time());
    plot(nearp,"parallelnn");

    //assertEqual(naive,nears);
    //assertEqual(naive,nearp);

    /*
    // serial discrete Sibson
    System.out.print("discrete sibson (serial): ");
    gridder.setParallel(false);
    s.restart();
    float[][] sibs = gridder.gridSibson(s1,s2);
    s.stop();
    System.out.format("%.4fs\n",s.time());
    plot(sibs,"serialsib");
    */

    /*
    // parallel discrete Sibson
    System.out.print("discrete sibson (parallel): ");
    gridder.setParallel(true);
    s.restart();
    float[][] sibp = gridder.gridSibson(s1,s2);
    s.stop();
    System.out.format("%.4fs\n",s.time());
    plot(sibp,"parallelsib");
    */

    //assertEqual(sibs,sibp);
  }

  private static void assertEqual(float[][] x, float[][] y) {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float xi = x[i2][i1];
        float yi = y[i2][i1];
        if (xi!=yi) System.out.println("(i1,i2) = ("+i1+","+i2+")");
        assert x[i2][i1]==y[i2][i1];
      }
    }
  }

  private static class Gaussian {
    public Gaussian(double sigma, int n1, int n2) {
      this.s = -0.5/(sigma*sigma);
      this.n1 = n1;
      this.n2 = n2;
      this.c1 = n1/2;
      this.c2 = n2/2;
    }
    public float evaluate(float x1, float x2) {
      double k1 = x1-c1;
      double k2 = x2-c2;
      double en = k1*k1+k2*k2;
      return (float)exp(s*en);
    }
    public float[][] toArray() {
      float[][] g = new float[n2][n1];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          double k1 = i1-c1;
          double k2 = i2-c2;
          double en = k1*k1+k2*k2;
          g[i2][i1] = (float)exp(s*en);
        }
      }
      return g;
    }
    private double s;
    private int n1,n2,c1,c2;
  }

  private static void plot(float[][] f, String png) {
    frame(panel(f),png);
  }

  private static void plot(
    float[][] f, float[] x1, float[] x2, String png)
  {
    PlotPanel panel = panel(f);
    PointsView point = panel.addPoints(x1,x2);
    point.setLineStyle(PointsView.Line.NONE);
    point.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    point.setMarkSize(4.0f);
    frame(panel,png);
  }

  private static PlotPanel panel(float[][] f) {
    PlotPanel panel = new PlotPanel();
    panel.setLimits(0,0,f[0].length,f.length);
    panel.addColorBar().setWidthMinimum(70);
    PixelsView pixel = panel.addPixels(f);
    pixel.setInterpolation(PixelsView.Interpolation.NEAREST);
    pixel.setColorModel(edu.mines.jtk.awt.ColorMap.JET);
    return panel;
  }

  private static void frame(PlotPanel p, String png) {
    PlotFrame frame = new PlotFrame(p);
    frame.setSize(800,700);
    frame.setVisible(true);
    if (_png)
      frame.paintToPng(360.0,3.0,"/Users/sluo/Desktop/"+png+".png");
  }
}

