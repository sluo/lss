package lss.vel;

import edu.mines.jtk.util.Parallel;
import lss.util.*;

public class WaveOperatorS {

  public WaveOperatorS(
  float[][] s, double dx, double dt, 
  int nabsorb, SharedFloat4 u, SharedFloat4 a) {
    _wave = new WaveOperator(s,dx,dt,nabsorb);
    _u = u;
    _a = a;
    _s = s;
    _nx = s[0].length;
    _nz = s.length;
    _nabsorb = nabsorb;
  }

  public void applyForward(
  final Source[] source, final Receiver[] receiver) {
    final int ns = source.length;
    final int np = _u.getN4(); // number of parallel shots
    PartialParallel parallel = new PartialParallel(np);
    parallel.loop(ns,new Parallel.LoopInt() {
      public void compute(int isou) {
        _wave.applyForward(source[isou],receiver[isou]);
      }
    });
  }

  ////////////////////////////////////////////////////////////////////////////  
  // private

  private final int _nabsorb;
  private final int _nx,_nz;
  private final float[][] _s;
  private final SharedFloat4 _u;
  private final SharedFloat4 _a;
  private final WaveOperator _wave;

}
