package lss.dev;

import java.util.Random;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

// TESTING
import edu.mines.jtk.mosaic.*;

/**
 * Semblance spectra for velocity analysis.
 * @author Simon Luo
 * @version 2012.02.16
 */
public class SemblanceSpectrum {

  public SemblanceSpectrum(
    Sampling st, Sampling sx, Sampling sv)
  {
    _st = st;
    _sx = sx;
    _sv = sv;
  }

  public float[][] apply(float[][] p) {
    int nv = _sv.getCount();
    float[][] s = new float[nv][];
    for (int iv=0; iv<nv; ++iv) {
      double v = _sv.getValue(iv);
      float[][] q = nmo(v,_st,_sx,p);
      s[iv] = semblance(_sigma,q);
    }
    return s;
  }

  ///////////////////////////////////////////////////////////////////////////

  private static double _sigma = 2.0; // smoother half-width
  private Sampling _st,_sx,_sv;

  private static float[][] nmo(
    double vnmo, Sampling st, Sampling sx, float[][] p)
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    int nx = sx.getCount();
    double dx = sx.getDelta();
    double fx = sx.getFirst();
    float[] t = new float[nt];
    float[][] q = new float[nx][nt];
    SincInterpolator si = new SincInterpolator();
    for (int ix=0; ix<nx; ++ix) {
      double x = sx.getValue(ix);
      double xxg = (x*x)/(vnmo*vnmo);
      for (int it=0; it<nt; ++it) {
        double t0 = st.getValue(it);
        t[it] = (float)sqrt(t0*t0+xxg); 
      }
      si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
    }
    return q;
  }

  private static float[] semblance(double tsigma, float[][] q) {
    int nx = q.length;
    int nt = q[0].length;
    float[] sn = new float[nt];
    float[] sd = new float[nt];
    for (int ix=0; ix<nx; ++ix) {
      for (int it=0; it<nt; ++it) {
        float qi = q[ix][it];
        if (qi==0.0f) continue;
        sn[it] += qi;
        sd[it] += qi*qi;
      }
    }
    mul(sn,sn,sn);
    mul(nx,sd,sd);
    esmooth(tsigma,sn,sn);
    esmooth(tsigma,sd,sd);
    float[] s = sn;
    for (int it=0; it<nt; ++it)
      s[it] = (sd[it]>0.0)?sn[it]/sd[it]:0.0f;
    return s;
  }

  public static void esmooth(double sigma, float[] p, float[] q) {
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE);
    //ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_VALUE);
    //ref.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE);
    ref.apply(p,q);
  }

  /////////////////////////////////////////////////////////////////////////
  // testing

  public static float[][] makeRickerGather(
    double fpeak, double snr, float[] vnmo, Sampling st, Sampling sx) 
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    double lt = st.getLast();
    int nx = sx.getCount();
    double dx = sx.getDelta();
    double fx = sx.getFirst();
    float[][] p = new float[nx][nt];
    double thalf =  1.0/fpeak;
    Random random = new Random(314);
    for (int jt=0; jt<nt; ++jt) {
      if (random.nextDouble()<0.1) // if reflection, ...
        continue;
      double t0 = st.getValue(jt);
      double a0 = 2.0*random.nextDouble()-1.0;
      double v0 = vnmo[st.indexOfNearest(t0)];
      double gamma = 1.0/(v0*v0);
      for (int ix=0; ix<nx; ++ix) {
        double x = (float)sx.getValue(ix);
        double t = sqrt(t0*t0+x*x*gamma);
        int itlo = max(0,(int)((t-thalf-ft)/dt));
        int ithi = min(nt-1,(int)((t+thalf-ft)/dt));
        for (int it=itlo; it<=ithi; ++it) {
          double twave = st.getValue(it)-t;
          p[ix][it] += (float)(a0*ricker(fpeak,twave));
        }
      }
    }
    return addRandomNoise((float)snr,p); // noise
  }

  public static float[] makeLinearVelocity(
    double vmin, double vmax, Sampling st) 
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    double lt = st.getLast();
    float[] v = new float[nt];
    for (int it=0; it<nt; ++it) {
      v[it] = (float)((it*vmax+(nt-1-it)*vmin)/(nt-1));
    }
    return v;
  }

  private static float[][] addRandomNoise(float r, float[][] p) {
    int nt = p[0].length;
    int nx = p.length;
    //float pmax = max(abs(p)); // peak signal
    float prms = sqrt(sum(mul(p,p))/nt/nx); // rms of signal
    Random random = new Random(3);
    float[][] s = sub(randfloat(random,nt,nx),0.5f);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply0X(s,s); // noise, bandlimited in time only
    float srms = sqrt(sum(mul(s,s))/nt/nx); // rms of noise
    return add(mul(prms/(srms*r),s),p); // r = rms-signal / rms-noise
  }

  private static float ricker(double fpeak, double time) {
    double x = PI*fpeak*time;
    double xx = x*x;
    return (float)((1.0-2.0*xx)*exp(-xx));
  }
}
