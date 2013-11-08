package lss.eni;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class FirstBreaks {

  /**
   * Constructs a first-break picker.
   * @param ewindow length of energy window.
   */
  public FirstBreaks(int ewindow) {
    _ewindow = ewindow;
  }

  public FirstBreaks() {
    this(EWINDOW_DEFAULT);
  }

  /**
   * Picks first breaks.
   * See Wong, Han, Bancroft, and Stewart, 2009, Automatic
   * time-picking of first arrivals on noisy microseismic data.
   * NB: Equation 4 is wrong; use reciprocal.
   */
  public float pick(float[] f) {
    final int nt = f.length;
    final int ne = _ewindow;
    final float eps = 1.0E-1f*max(abs(f)); // epsilon to stabilize division
    final float[] er = new float[nt]; // energy ratio
    // it=0
    float num = f[0];
    float den = 0.0f;
    for (int it=0; it<min(ne+1,nt); ++it) {
      float fi = f[it];
      den += fi*fi;
    }
    er[0] = abs(f[0])*den/(num+eps);
    // 0<it<nt
    for (int it=1; it<nt; ++it) {
      float fa = (it-ne-1>=0)?f[it-ne-1]:0.0f;
      float fb = f[it-1];
      float fc = f[it];
      float fd = (it+ne<nt)?f[it+ne]:0.0f;
      num += fc*fc-fa*fa;
      den += fd*fd-fb*fb;
      er[it] = abs(f[it])*den/(num+eps);
    }
    int[] i = new int[1];
    max(er,i);
    return i[0];
  }

  public float[] pick(final float[][] f) {
    final int nt = f[0].length;
    final int nr = f.length;
    final float[] pp = new float[nr]; // picks
    final float[][] er = new float[nr][nt]; // energy ratio
    Parallel.loop(nr,new Parallel.LoopInt() {
    public void compute(int ir) {
      pp[ir] = pick(f[ir]);
    }});
    return pp;
  }

  private int _ewindow; // energy window
  private static final int EWINDOW_DEFAULT = 1000; // default window
}
