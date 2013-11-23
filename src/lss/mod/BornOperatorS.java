package lss.mod;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import lss.util.*;

// testing
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;

/**
 * Acoustic constant-density Born modeling with an optional 
 * time-shift operator, parallelized over sources.
 * @see lss.vel.BornOperator
 * @author Simon Luo, Colorado SChool
 * @version 2013.11.20
 */
public class BornOperatorS {

  public static Source[] getSourceArray(int n) {
    return new Source[n];
  }

  public static Receiver[] getReceiverArray(int n) {
    return new Receiver[n];
  }

  ////////////////////////////////////////////////////////////////////////////  

  public BornOperatorS(
  float[][] s, double dx, double dt, 
  int nabsorb, SharedFloat4 u, SharedFloat4 a) {
    Check.argument(s[0].length==u.getN1()-2*nabsorb,"consistent nx");
    Check.argument(s.length==u.getN2()-2*nabsorb,"consistent nz");
    _born = new BornOperator(s,dx,dt,nabsorb);
    _u = u;
    _a = a;
    _s = s;
    _nx = s[0].length;
    _nz = s.length;
    _nabsorb = nabsorb;
  }

  public void setAdjoint(boolean adjoint) {
    _born.setAdjoint(adjoint);
  }

  public int[] getNxNz() {
    return new int[]{_nx,_nz};
  }

  public int[] getDimensions() {
    return new int[]{_nx,_nz};
  }

  public void setSlowness(float[][] s) {
    copy(s,_s);
    _born.setSlowness(s);
  }

  /**
   * Applies the forward operator.
   * @param source input sources for computing background wavefields. 
   * @param rx input reflectivity image.
   * @param ts input time shifts.
   * @param receiver output receivers.
   */
  public void applyForward(
    final Source[] source, final float[][] rx, final float[][][] ts,
    final Receiver[] receiver)
  {
    Check.argument(rx[0].length==_nx,"consistent nx");
    Check.argument(rx.length==_nz,"consistent nz");
    final int ns = source.length;
    final int np = _u.getN4(); // number of parallel shots
    PartialParallel parallel = new PartialParallel(np);
    parallel.loop(ns,new Parallel.LoopInt() {
      public void compute(int isou) {
        float[][][] u = _u.get(isou);
        if (ts==null) {
          _born.applyForward(source[isou],u,rx,receiver[isou]);
        } else {
          _born.applyForward(source[isou],u,rx,ts[isou],receiver[isou]);
        }
      }
    });
  }

  public void applyForward(
    Source[] source, float[][] rx,
    Receiver[] receiver)
  {
    applyForward(source,rx,null,receiver);
  }

  /**
   * Applies the adjoint operator.
   * @param source input sources for computing background wavefields. 
   * @param receiver input receivers containing data to be migrated.
   * @param ts input time shifts.
   * @param ry output reflectivity image.
   */
  public void applyAdjoint(
    final Source[] source, final Receiver[] receiver, final float[][][] ts,
    final float[][] ry)
  {
    Check.argument(ry[0].length==_nx,"consistent nx");
    Check.argument(ry.length==_nz,"consistent nz");
    final int ns = source.length;
    final int nx = ry[0].length;
    final int nz = ry.length;
    final int np = _a.getN4(); // number of parallel shots
    PartialParallel parallel = new PartialParallel(np);
    float[][] rz = parallel.reduce(ns,new Parallel.ReduceInt<float[][]>() {
      public float[][] compute(int isou) {
        float[][][] u = _u.get(isou);
        float[][][] a = _a.get(isou);
        float[][] rt = new float[nz][nx];
        if (ts==null) {
          _born.applyAdjoint(source[isou],u,a,receiver[isou],rt);
        } else {
          _born.applyAdjoint(source[isou],u,a,receiver[isou],ts[isou],rt);
        }
        return rt;
      }
      public float[][] combine(float[][] ra, float[][] rb) {
        return add(ra,rb);
      }
    });
    copy(rz,ry);
  }

  public void applyAdjoint(
    Source[] source, Receiver[] receiver,
    float[][] ry)
  {
    applyAdjoint(source,receiver,null,ry);
  }

  /**
   * Applies the Hessian operator.
   * @param source input sources for computing background wavefields. 
   * @param receiver input receivers containing data to be migrated.
   * @param rx input reflectivity image.
   * @param ts input time shifts.
   * @param ry output reflectivity image.
   */
  public void applyHessian(
    final Source[] source, final Receiver[] receiver,
    final float[][] rx, final float[][][] ts,
    final float[][] ry)
  {
    Check.argument(rx[0].length==_nx,"consistent nx");
    Check.argument(rx.length==_nz,"consistent nz");
    Check.argument(ry[0].length==_nx,"consistent nx");
    Check.argument(ry.length==_nz,"consistent nz");
    final int ns = receiver.length;
    final int nx = ry[0].length;
    final int nz = ry.length;
    final int np = _a.getN4(); // number of parallel shots
    PartialParallel parallel = new PartialParallel(np);
    float[][] rz = parallel.reduce(ns,new Parallel.ReduceInt<float[][]>() {
      public float[][] compute(int isou) {
        float[][] rt = new float[nz][nx];
        float[][][] u = _u.get(isou);
        float[][][] a = _a.get(isou);
        if (ts==null) {
          _born.applyHessian(source[isou],receiver[isou],u,a,rx,rt);
        } else {
          _born.applyHessian(source[isou],receiver[isou],u,a,rx,ts[isou],rt);
        }
        return rt;
      }
      public float[][] combine(float[][] ra, float[][] rb) {
        return add(ra,rb);
      }
    });
    copy(rz,ry);
  }

  public void applyHessian(
    Source[] source, Receiver[] receiver, float[][] rx,
    float[][] ry)
  {
    applyHessian(source,receiver,rx,null,ry);
  }

  // Mostly for adjoint test.
  public void applyForward(
  final float[][] rx, final Receiver[] receiver) {
    final int ns = receiver.length;
    final int np = _u.getN4(); // number of parallel shots
    PartialParallel parallel = new PartialParallel(np);
    parallel.loop(ns,new Parallel.LoopInt() {
      public void compute(int isou) {
        float[][][] u = _u.get(isou);
        _born.applyForward(u,rx,receiver[isou]);
      }
    });
  }
  public void applyAdjoint(
  final Receiver[] receiver, final float[][] ry) {
    final int ns = receiver.length;
    final int nx = ry[0].length;
    final int nz = ry.length;
    final int np = _a.getN4(); // number of parallel shots
    PartialParallel parallel = new PartialParallel(np);
    float[][] rz = parallel.reduce(ns,new Parallel.ReduceInt<float[][]>() {
      public float[][] compute(int isou) {
        float[][][] u = _u.get(isou);
        float[][][] a = _a.get(isou);
        float[][] rt = new float[nz][nx];
        _born.applyAdjoint(u,a,receiver[isou],rt);
        return rt;
      }
      public float[][] combine(float[][] ra, float[][] rb) {
        return add(ra,rb);
      }
    });
    copy(rz,ry);
  }

  ////////////////////////////////////////////////////////////////////////////
  // preconditioning

  public void applyForIllumination(
  final Source[] source, final float[][] m) {
    final int ns = source.length;
    final int nx = m[0].length;
    final int nz = m.length;
    final int np = _u.getN4(); // number of parallel shots
    PartialParallel parallel = new PartialParallel(np);
    float[][] mm = parallel.reduce(ns,new Parallel.ReduceInt<float[][]>() {
      public float[][] compute(int isou) {
        float[][][] u = _u.get(isou);
        float[][] mt = new float[nz][nx];
        _born.applyForIllumination(source[isou],u,mt);
        return mt;
      }
      public float[][] combine(float[][] ma, float[][] mb) {
        return add(ma,mb);
      }
    });
    copy(mm,m);
  }
  public void xapplyForIllumination(Source[] source, float[][] m) {
    float[][][] b = _u.get(0);
    Check.argument(b[0][0].length-m[0].length==2*_nabsorb,"consistent nx");
    Check.argument(b[0].length-m.length==2*_nabsorb,"consistent nz");
    Source simultaneousSource = new Source.SimultaneousSource(source);
    _born.getWaveOperator().applyForward(simultaneousSource,b);
    WaveOperator.collapse(b,b,_nabsorb,m);
  }

  ////////////////////////////////////////////////////////////////////////////  
  // private

  private final int _nabsorb;
  private final int _nx,_nz;
  private final float[][] _s;
  private final SharedFloat4 _u;
  private final SharedFloat4 _a;
  private final BornOperator _born;

}
