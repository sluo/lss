package lss.vel;

import edu.mines.jtk.opt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import lss.util.SharedFloat4;

// TESTING
import edu.mines.jtk.mosaic.*;

public class BornSolver {

  public BornSolver(
    BornOperatorS born, Source[] src,
    Receiver[] rcp, Receiver[] rco, RecursiveExponentialFilter ref)
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
    _ii = new float[_nz][_nx];
    computeInverseIllumination(_ii);
    _qs = new QuadraticSolver(new Q());
  }

  public void setObservedData(final Receiver[] rco) {
    Check.argument(_rco.length==rco.length,"_rco.length==rco.length");
    final int ns = rco.length;
    Parallel.loop(ns,new Parallel.LoopInt() {
    public void compute(int is) {
      copy(rco[is].getData(),_rco[is].getData());
    }});
  }

  public float[][] solve(int niter) {
    ArrayVect2f v2y = (ArrayVect2f)_qs.solve(niter,null);
    return v2y.getData();
  }

  public void setTrueReflectivity(float[][] r) {
    _r = r;
  }

  public void setMask(float[][] m) {
    _m = m;
  }

  public void setTrueBornOperator(BornOperatorS bornt) {
    _bornt = bornt;
  }

  ////////////////////////////////////////////////////////////////////////////
  // private

  // required parameters
  private final int _nz,_nx;
  private final float[][] _ii; // inverse illumination
  private final BornOperatorS _born;
  private final Source[] _src;
  private final Receiver[] _rco;
  private final Receiver[] _rcp;
  private final RecursiveExponentialFilter _ref; // roughening filter
  private final QuadraticSolver _qs;

  // optional parameters
  private float[][] _m = null; // mask
  private float[][] _r = null; // true reflectivity
  private BornOperatorS _bornt = null; // Born operator with true slowness

  private class Q implements Quadratic {
    public void multiplyHessian(Vect vx) {
      ArrayVect2f v2x = (ArrayVect2f)vx;
      float[][] rx = v2x.getData();
      _born.applyHessian(_src,_rcp,rx,rx);
    }
    public void inverseHessian(Vect vx) {
      ArrayVect2f v2x = (ArrayVect2f)vx;
      float[][] rx = v2x.getData();
      float[][] ry = new float[_nz][_nx];
      _born.applyAdjointRoughen(_ref,rx,ry);
      mul(_ii,ry,ry); // inverse illumination
      if (_m!=null)
        mul(_m,ry,ry); // mask if non-null
      _born.applyForwardRoughen(_ref,ry,rx);
    }
    public Vect getB() {
      float[][] rb = new float[_nz][_nx];
      ArrayVect2f vb = new ArrayVect2f(rb,1.0);
      if (_r!=null) {
        if (_bornt!=null) {
          _bornt.applyForward(_src,_r,_rcp);
          _born.applyAdjoint(_src,_rcp,rb);
        } else {
          _born.applyHessian(_src,_rco,_r,rb);
        }
      } else {
        Check.argument(_bornt==null,"bornt==null");
        _born.applyAdjoint(_src,_rco,rb);
      }
      mul(-1.0f,rb,rb); // need negative gradient (?)
      SimplePlot.asPixels(transpose(rb));
      return vb;
    }
  }

  private void computeInverseIllumination(float[][] ii) {
    System.out.println("computing illumination");
    _born.applyForIllumination(_src,ii); // illumination
    mul(ii,ii,ii); // squared
    add(1.0e-8f*max(ii),ii,ii); // stabilized
    div(1.0f,ii,ii); // inverted
    mul(1.0f/max(ii),ii,ii); // normalized
  }
}
