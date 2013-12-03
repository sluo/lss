package lss.mig;

import edu.mines.jtk.opt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import lss.mod.*;
import lss.opt.*;
import lss.util.SharedFloat4;

/**
 * Linearized waveform inversion, with the option to
 * include time-shifts within the forward modeling operator.
 * @author Simon Luo, Colorado School of Mines
 * @version 2013.11.22
 */
public class BornSolver {

  /**
   * Constructs a solver.
   * @param born Born modeling operator.
   * @param src sources.
   * @param rcp receivers for predicted data.
   * @param rco receivers containing observed data.
   * @param ref roughening filter.
   * @param mp model preconditioning, e.g., water-layer mask.
   * @param ts time shifts.
   */
  public BornSolver(
    BornOperatorS born,
    Source[] src, Receiver[] rcp, Receiver[] rco,
    RecursiveExponentialFilter ref, float[][] mp, float[][][] ts)
  {
    Check.argument(src.length==rco.length,"src.length==rco.length");
    int[] nxz = born.getNxNz();
    _nx = nxz[0];
    _nz = nxz[1];
    _born = born;
    _src = src;
    _rcp = rcp;
    _rco = rco;
    _ref = ref;
    _mp = mp;
    _ts = ts;
  }

  public BornSolver(
    BornOperatorS born,
    Source[] src, Receiver[] rcp, Receiver[] rco,
    RecursiveExponentialFilter ref, float[][] mp)
  {
    this(born,src,rcp,rco,ref,mp,null);
  }

  /**
   * Solves for the reflectivity image.
   * @param niter number of CG iterations to perform.
   * @param rx input/output reflectivity image.
   */
  public void solve(int niter, float[][] rx) {
    CgSolver cg = new CgSolver(0.0,niter);
    VecArrayFloat2 vx = new VecArrayFloat2(rx);
    VecArrayFloat2 vb = new VecArrayFloat2(makeB());
    CgSolver.A ma = new A();
    CgSolver.A mm = new M();
    cg.solve(ma,mm,vb,vx);
  }

  public float[][] solve(int niter) {
    float[][] rx = new float[_nz][_nx];
    solve(niter,rx);
    return rx;
  }

  // Old method using QuadraticSolver.
  public float[][] xsolve(int niter) {
    QuadraticSolver qs = new QuadraticSolver(new Q());
    ArrayVect2f v2y = (ArrayVect2f)qs.solve(niter,null);
    return v2y.getData();
  }

  public void setTimeShifts(float[][][] ts) {
    _ts = ts;
  }

  public void setObservedData(Receiver[] rco) {
    _rco = rco;
  }

  public void setTrueReflectivity(float[][] r) {
    _r = r;
  }

  public void setMask(float[][] m) {
    _mp = m;
  }

  public void setTrueBornOperator(BornOperatorS bornt) {
    _bornt = bornt;
  }

  ////////////////////////////////////////////////////////////////////////////
  // private

  // required parameters
  private final int _nz,_nx;
  private final RecursiveExponentialFilter _ref; // roughening filter
  private final BornOperatorS _born;
  private final Source[] _src;
  private final Receiver[] _rcp;
  private Receiver[] _rco;

  // optional parameters
  private float[][] _mp = null; // model preconditioner
  private float[][] _r = null; // true reflectivity
  private float[][][] _ts = null; // time shifts
  private BornOperatorS _bornt = null; // Born operator with true slowness

  // LHS Hessian operator.
  private void applyA(float[][] rx, float[][] ry) {
    //System.out.println("Applying LHS...");
    _born.applyHessian(_src,_rcp,rx,_ts,ry);
  }

  // Approximate inverse Hessian or preconditioner.
  private void applyM(float[][] rx, float[][] ry) {
    float[][] rz = new float[_nz][_nx];
    if (_ref!=null) {
      applyAdjointRoughen(_ref,rx,rz);
    }
    if (_mp!=null) {
      mul(_mp,rz,rz); // model precondition (mask) if non-null
    }
    if (_ref!=null) {
      applyForwardRoughen(_ref,rz,ry);
    }
  }

  // RHS vector.
  private float[][] makeB() {
    float[][] rb = new float[_nz][_nx];
    if (_r!=null) {
      if (_bornt!=null) {
        _bornt.applyForward(_src,_r,_rco);
      } else {
        _born.applyForward(_src,_r,_rco);
      }
    }
    _born.applyAdjoint(_src,_rco,_ts,rb);
    return rb;
  }

  private static void applyForwardRoughen(
  RecursiveExponentialFilter ref, float[][] rx, float[][] ry) {
    float[][] cx = (rx==ry)?copy(rx):rx;
    ref.apply1(cx,ry);
    ref.apply2(ry,ry);
    sub(cx,ry,ry);
  } 

  private static void applyAdjointRoughen(
  RecursiveExponentialFilter ref, float[][] rx, float[][] ry) {
    float[][] cx = (rx==ry)?copy(rx):rx;
    ref.apply2(cx,ry);
    ref.apply1(ry,ry);
    sub(cx,ry,ry);
  }

  private void computeInverseIllumination(float[][] ii) {
    System.out.println("computing illumination");
    _born.applyForIllumination(_src,ii); // illumination
    mul(ii,ii,ii); // squared
    add(1.0e-8f*max(ii),ii,ii); // stabilized
    div(1.0f,ii,ii); // inverted
    mul(1.0f/max(ii),ii,ii); // normalized
  }

  ////////////////////////////////////////////////////////////////////////////
  // QuadraticSolver

  private class Q implements Quadratic {
    public void multiplyHessian(Vect vx) {
      float[][] rx = ((ArrayVect2f)vx).getData();
      applyA(rx,rx);
    }
    public void inverseHessian(Vect vx) {
      float[][] rx = ((ArrayVect2f)vx).getData();
      applyM(rx,rx);
    }
    public Vect getB() {
      float[][] rb = makeB();
      mul(-1.0f,rb,rb); // Harlan's B defined as negative of RHS of Ax=b.
      return new ArrayVect2f(rb,1.0);
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // CgSolver

  private class A implements CgSolver.A {
    public void apply(Vec vx, Vec vy) {
      float[][] rx = ((VecArrayFloat2)vx).getArray();
      float[][] ry = ((VecArrayFloat2)vy).getArray();
      applyA(rx,ry);
    }
  }

  private class M implements CgSolver.A {
    public void apply(Vec vx, Vec vy) {
      float[][] rx = ((VecArrayFloat2)vx).getArray();
      float[][] ry = ((VecArrayFloat2)vy).getArray();
      applyM(rx,ry);
    }
  }
}
