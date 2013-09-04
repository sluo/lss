package lss.vel;

import java.util.Random;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

// testing
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;

import lss.vel.AcousticWaveOperator.*;

public class BornWaveOperator {

  public BornWaveOperator(
  Source source, float[][] s, double dx, double dt, int nabsorb) {
    this(source,null,s,dx,dt,nabsorb);
  }

  public BornWaveOperator(
  float[][][] u, float[][] s, double dx, double dt, int nabsorb) {
    this(null,u,s,dx,dt,nabsorb);
  }

  private BornWaveOperator(
  Source source, float[][][] u,
  float[][] s, double dx, double dt, int nabsorb) {
    _wave = new AcousticWaveOperator(s,dx,dt,nabsorb);
    int nx = s[0].length;
    int nz = s.length;
    if (u!=null)
      Check.argument(u[0][0].length-nx==2*nabsorb,
        "x-dimension inconsistent");
      Check.argument(u[0].length-nz==2*nabsorb,
        "z-dimension inconsistent");
    _nabsorb = nabsorb;
    _nz = nz+2*nabsorb;
    _nx = nx+2*nabsorb;
    _dx = (float)dx;
    _dt = (float)dt;
    _s = s;
    _u = u;
    _source = source;
  }

  private float[][][] getBackgroundWavefield(int nt) {
    if (_u==null)
      _u = new float[nt][_nz][_nx];
      _wave.applyForward(_source,_u);
    return _u;
  }

  public void applyForward(float[][] r, Receiver receiver) {
    int nt = receiver.getNt();
    float[][][] u = getBackgroundWavefield(nt);
    Source source = new WavefieldSource(u,r);
    _wave.applyForward(source,receiver);
  }

  public void applyAdjoint(
  float[][][] a, Receiver receiver, float[][] r) {
    int nt = a.length;
    float[][][] u = getBackgroundWavefield(nt);
    Source source = new ReceiverSource(receiver);
    _wave.applyAdjoint(source,a);
    AcousticWaveOperator.collapse(u,a,_nabsorb,r);
  }

  private void rmul(final float[][] r, final float[][][] u) {
    final int nx = r[0].length;
    final int nz = r.length;
    final int nt = u.length;
    Parallel.loop(nt,new Parallel.LoopInt() {
    public void compute(int it) { 
      for (int iz=0; iz<nz; ++iz)
        for (int ix=0; ix<nx; ++ix)
          u[it][iz+_nabsorb][ix+_nabsorb] *= r[iz][ix];
    }});
  }

  //////////////////////////////////////////////////////////////////////////
  // private

  private AcousticWaveOperator _wave;
  private AcousticWaveOperator.Source _source;
  private int _nabsorb;
  private int _nx,_nz;
  private float _dx,_dt;
  private float[][] _s; // background slowness
  private float[][][] _u = null; // background wavefield

}
