package lss.dev;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

import lss.dev.WaveOperator;

/**
 * Forward, backward, and adjoint migration operator.
 * <p>
 * The forward operator models data from a model, while the backward
 * and adjoint operators migrate data to produce an image of the model.
 * Here, the model is a reflectivity model, or a model of the perturbation
 * to the background slowness model.
 * <p>
 * @author Simon Luo
 * @version 2013.08.29
 */
public class WaveMigrationOperator {

  public WaveMigrationOperator(Sampling sz, Sampling st, float[] s) {
    _nz = sz.getCount();
    _nt = st.getCount();
    _dz = sz.getDelta();
    _dt = st.getDelta();
    _a = div((float)(_dt*_dt/(_dz*_dz)),mul(s,s));
    _wave = new WaveOperator(sz,st,s);
  }


  public float[] applyForward(float[] m) {
    float[][] u = new float[_nt][];
    float[][] b = getBackgroundWavefield();

    SimplePlot.asPixels(b);
    for (int it=0; it<_nt; ++it)
      u[it] = mul(m,b[it]);
    SimplePlot.asPixels(u);

    float[] d = _wave.applyForward(u);
    SimplePlot.asPixels(u);

    return d;

  }

  public float[] applyAdjoint(float[] d) {
    float[][] u = new float[_nt][_nz];
    _wave.applyAdjoint(d,u);
    SimplePlot.asPixels(u);
    float[] m = new float[_nz];
    for (int it=0; it<_nt; ++it)
      for (int iz=0; iz<_nz; ++iz)
        m[iz] += _b[it][iz]*u[it][iz];
    return m;
  }

  //////////////////////////////////////////////////////////////////////////
  // private

  private static final double FPEAK = 10.0; // Ricker peak frequency
  private float[][] _b = null; // background wavefield
  private float[] _a;
  private double _dz;
  private double _dt;
  private int _nz;
  private int _nt;
  private WaveOperator _wave;

  public float[][] getBackgroundWavefield() {
    if (_b==null) {
      _b = new float[_nt][_nz];
      //float[] src = new float[_nt];
      for (int it=0; it<_nt; ++it) {
        double t = _dt*it;
        //src[it] = ricker(t);
        _b[it][0] = ricker(t);
      }
      //plotDepthSlice(0,_b);
      //SimplePlot.asPixels(_b);
      _wave.applyForward(_b);
    }
    return _b;
  }

  private static float ricker(double t) {
    double tdelay = 1.0/FPEAK;
    double x = PI*FPEAK*(t-tdelay);
    double xx = x*x;
    return (float)((1.0-2.0*xx)*exp(-xx));
  }

  //////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) { 
    //test();
    adjointTest();
  }

  public static void test() {
    Sampling sz = new Sampling(501,0.0120,0.0);
    Sampling st = new Sampling(1201,0.0015,0.0);
    //Sampling sz = new Sampling(1001,0.0060,0.0);
    //Sampling st = new Sampling(3001,0.0005,0.0);
    int nz = sz.getCount();
    int nt = st.getCount();
    double dz = sz.getDelta();
    double dt = st.getDelta();
    float[] s = fillfloat(0.25f,nz); // slowness

    WaveMigrationOperator wmo = new WaveMigrationOperator(sz,st,s);
    float[][] b = wmo.getBackgroundWavefield();
    //SimplePlot.asPixels(b);

    //plotDepthSlice(0,b);

    // Reflectivity model.
    float[] m = new float[nz];
    m[nz/2] = 1.0f;

    // Forward operator (modeling).
    float[] d = wmo.applyForward(m);
    SimplePlot.asPoints(d);

    // Adjoint operator (migration).
    float[] r = wmo.applyAdjoint(d);
    SimplePlot.asPoints(r);


  }

  public static void adjointTest() {
    Sampling sz = new Sampling(1001,0.0120,0.0);
    Sampling st = new Sampling(1201,0.0015,0.0);
    int nz = sz.getCount();
    int nt = st.getCount();
    double dz = sz.getDelta();
    double dt = st.getDelta();
    float[] s = fillfloat(0.25f,nz); // slowness
    WaveMigrationOperator wave = new WaveMigrationOperator(sz,st,s);
    Random random = new Random(0123);
    float[] m = sub(randfloat(random,nz),0.5f);
    float[] d = sub(randfloat(random,nt),0.5f);
    float sum1 = dot(wave.applyForward(m),d);
    float sum2 = dot(wave.applyAdjoint(d),m);
    System.out.println("adjoint test:");
    System.out.println(sum1);
    System.out.println(sum2);
  }

  private static void plotDepthSlice(int ix, float[][] u) {
    int nt = u.length;
    float[] d = new float[nt];
    for (int it=0; it<nt; ++it)
      d[it] = u[it][ix];
    SimplePlot.asPoints(d);
  }

  private static float dot(float[] x, float[] y) {
    int n = x.length;
    float sum = 0.0f;
    for (int i=0; i<n; ++i)
      sum += x[i]*y[i];
    return sum;
  }

}
