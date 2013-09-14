package lss.vel;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import lss.util.*;
import lss.vel.BornWaveOperator;

// testing
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;


public class BornWaveOperatorS {

  public static Source[] getSourceArray(int n) {
    return new Source[n];
  }

  public static Receiver[] getReceiverArray(int n) {
    return new Receiver[n];
  }

  ////////////////////////////////////////////////////////////////////////////  

  public BornWaveOperatorS(
  float[][] s, double dx, double dt, 
  int nabsorb, SharedFloat4 u, SharedFloat4 a) {
    _born = new BornWaveOperator(s,dx,dt,nabsorb);
    _born.setIlluminationCompensation(true);
    _born.setReflectivityRoughening(0.0);
    _u = u;
    _a = a;
    _s = s;
    _dx = dx;
    _dt = dt;
    _nabsorb = nabsorb;
  }

  public void setReflectivityRoughening(double sigma) {
    if (sigma>0.0) {
      _ref = new RecursiveExponentialFilter(sigma/_dx);
      _ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_VALUE);
    } else {
      _ref = null;
    }
  }

  public void applyForward(
  final Source[] source, final float[][] rx, final Receiver[] receiver) {
    final float[][] qx = copy(rx);
    applyForwardRoughen(rx,qx);
    final int ns = source.length;
    final int np = _u.getN4(); // number of parallel shots
    PartialParallel parallel = new PartialParallel(np);
    parallel.loop(ns,new Parallel.LoopInt() {
      public void compute(int isou) {
        float[][][] u = _u.get(isou);
        _born.applyForward(source[isou],u,qx,receiver[isou]);
      }
    });
  }

  public void applyAdjoint(
  final Source[] source, final Receiver[] receiver, final float[][] ry) {
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
        _born.applyAdjoint(source[isou],u,a,receiver[isou],rt);
        return rt;
      }
      public float[][] combine(float[][] ra, float[][] rb) {
        return add(ra,rb);
      }
    });
    //copy(rz,ry);
    applyAdjointRoughen(rz,ry);
  }

  public void applyHessian(
  final Source[] source, final Receiver[] receiver,
  final float[][] rx, final float[][] ry) {
    final float[][] qx = copy(rx);
    applyForwardRoughen(rx,qx);
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
        _born.applyHessian(source[isou],u,a,receiver[isou],qx,rt);
        return rt;
      }
      public float[][] combine(float[][] ra, float[][] rb) {
        return add(ra,rb);
      }
    });
    //copy(rz,ry);
    applyAdjointRoughen(rz,ry);
  }

  // Mostly for adjoint test.
  public void applyForward(
  final float[][] rx, final Receiver[] receiver) {
    final float[][] qx = copy(rx);
    applyForwardRoughen(rx,qx);
    final int ns = receiver.length;
    final int np = _u.getN4(); // number of parallel shots
    PartialParallel parallel = new PartialParallel(np);
    parallel.loop(ns,new Parallel.LoopInt() {
      public void compute(int isou) {
        float[][][] u = _u.get(isou);
        _born.applyForward(u,qx,receiver[isou]);
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
    //copy(rz,ry);
    applyAdjointRoughen(rz,ry);
  }

  ////////////////////////////////////////////////////////////////////////////  
  // private

  private RecursiveExponentialFilter _ref = null;
  private final int _nabsorb;
  private final double _dx,_dt;
  private final float[][] _s;
  private final SharedFloat4 _u;
  private final SharedFloat4 _a;
  private final BornWaveOperator _born;

  private void applyForwardRoughen(float[][] rx, float[][] ry) {
    if (_ref!=null) {
      Check.argument(rx!=ry,"rx!=ry");
      _ref.apply1(rx,ry);
      _ref.apply2(ry,ry);
      sub(rx,ry,ry);
    } else {
      copy(rx,ry);
    }
  } 

  private void applyAdjointRoughen(float[][] rx, float[][] ry) {
    if (_ref!=null) {
      Check.argument(rx!=ry,"rx!=ry");
      _ref.apply2(rx,ry);
      _ref.apply1(ry,ry);
      sub(rx,ry,ry);
    } else {
      copy(rx,ry);
    }
  }
}
