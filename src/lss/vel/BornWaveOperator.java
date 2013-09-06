package lss.vel;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

// testing
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;

public class BornWaveOperator {

  public BornWaveOperator(
  float[][] s, double dx, double dt, int nabsorb) {
    _wave = new AcousticWaveOperator(s,dx,dt,nabsorb);
    int nx = s[0].length;
    int nz = s.length;
    _nabsorb = nabsorb;
    _nz = nz+2*nabsorb;
    _nx = nx+2*nabsorb;
    _dx = (float)dx;
    _dt = (float)dt;
    _s = s;
  }

  public void applyForward(
  Source source, float[][][] u, float[][] rx, Receiver receiver) {
    _wave.applyForward(source,u);
    applyForward(u,rx,receiver);
  }

  public void applyForward(
  float[][][] u, float[][] rx, Receiver receiver) {
    Check.argument(u[0][0].length-rx[0].length==2*_nabsorb,
      "x-dimension inconsistent");
    Check.argument(u[0].length-rx.length==2*_nabsorb,
      "z-dimension inconsistent");
    _wave.applyForward(new Source.WavefieldSource(u,rx),receiver);
  }

  public void applyAdjoint(
  Source source, float[][][] u,
  float[][][] a, Receiver receiver, float[][] ry) {
    _wave.applyForward(source,u);
    applyAdjoint(u,a,receiver,ry);
  }

  public void applyAdjoint(
  float[][][] u, float[][][] a, Receiver receiver, float[][] ry) {
    Check.argument(u[0][0].length-ry[0].length==2*_nabsorb,
      "x-dimension inconsistent");
    Check.argument(u[0].length-ry.length==2*_nabsorb,
      "z-dimension inconsistent");
    _wave.applyAdjoint(new Source.ReceiverSource(receiver),a);
    AcousticWaveOperator.collapse(u,a,_nabsorb,ry);
  }

  public void applyHessian(
  Source source, float[][][] u, float[][][] a,
  Receiver receiver, float[][] rx, float[][] ry) {
    applyForward(source,u,rx,receiver);
    applyAdjoint(u,a,receiver,ry);
  }

  public void applyHessian(
  float[][][] u, float[][][] a,
  Receiver receiver, float[][] rx, float[][] ry) {
    applyForward(u,rx,receiver);
    applyAdjoint(u,a,receiver,ry);
  }

  //////////////////////////////////////////////////////////////////////////
  // private

  private AcousticWaveOperator _wave;
  private int _nabsorb;
  private int _nx,_nz;
  private float _dx,_dt;
  private float[][] _s; // background slowness

}
