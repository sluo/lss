/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lss.sem;

import java.awt.*;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Tests use of weighted semblance in velocity analysis.
 * @author Dave Hale, Colorado School of Mines
 * @version 2011.04.05
 */
public class Velan {

  public static float[][] pickPeaks(
    final Sampling st,
    final float[] vmin, final float[] vmax, final float[][][] s)
  {
    final int nt = s[0][0].length;
    final int ngather = s.length;
    final int nthread = Runtime.getRuntime().availableProcessors();
    final float[][] picks = new float[ngather][];
    final AtomicInteger ai = new AtomicInteger();
    Thread[] threads = new Thread[nthread];
    for (int ithread=0; ithread<nthread; ++ithread) {
      threads[ithread] = new Thread() {
        public void run() {
          for (int igather=ai.getAndIncrement(); igather<ngather; 
                   igather=ai.getAndIncrement()) {
            picks[igather] = pickPeaks(st,vmin,vmax,s[igather]);
          }
        }
      };
    }
    Threads.startAndJoin(threads);
    return picks;

    //float[] p = new float[nt];
    //for (int igather=0; igather<ngather; ++igather)
    //  add(p,picks[igather],p);
    //div(p,ngather,p);
    //return p;
  }

  public static float[] pickPeaks(
    Sampling sv, float[] vmin, float[] vmax, float[][] s)
  {
    int nt = s[0].length;
    int nv = s.length;
    double dv = sv.getDelta();
    double fv = sv.getFirst();
    float[] p = new float[nt];
    for (int it=0; it<nt; ++it) {
      int va = (int)((vmin[it]-fv)/dv);
      int vb = (int)((vmax[it]-fv)/dv);
      float smax = 0.0f;
      for (int iv=va; iv<vb; ++iv) {
        double v = fv+iv*dv;
        float si = s[iv][it];
        if (si>smax) {
          smax = si;
          p[it] = (float)v;
        }
      }
    }
    return p;
  }

  // XXX: for flattening
  public static float[][] makeLinearGather(
    double fpeak, float[] vnmo, Sampling st, Sampling sx)
  {
    return makeLinearGather(false,fpeak,vnmo,st,sx);
  }
  public static float[][] makeLinearGather(
    boolean fault, double fpeak, float[] vnmo, Sampling st, Sampling sx) 
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    double lt = st.getLast();
    int nx = sx.getCount();
    double dx = sx.getDelta();
    double fx = sx.getFirst();
    Random random = new Random(314159);
    //Random random = new Random(12111);
    //Random random = new Random();
    float[][] p = new float[nx][nt];
    double thalf =  1.0/fpeak;
    for (int jt=0; jt<nt; ++jt) {
      if (random.nextDouble()<0.1) // if reflection, ...
        continue;
      double t0 = st.getValue(jt);
      double a0 = 2.0*random.nextDouble()-1.0;
      double v0 = vnmo[st.indexOfNearest(t0)];
      double gamma = 1.0/(v0*v0);
      for (int ix=0; ix<nx; ++ix) {
        double x = (float)sx.getValue(ix);
        double t = t0+x*gamma;

        if (ix<nx/2 && fault) // fault
          t -= 0.2;

        int itlo = max(0,(int)((t-thalf-ft)/dt));
        int ithi = min(nt-1,(int)((t+thalf-ft)/dt));
        for (int it=itlo; it<=ithi; ++it) {
          double twave = st.getValue(it)-t;
          p[ix][it] += (float)(a0*ricker(fpeak,twave));
        }
      }
    }
    return p;
  }

  private static final float SCL = 0.1f;

  public static float[][][] velocitySpectrum(
    final Sampling st, final Sampling sx, final float[][][] p,
    final Sampling sv, final double tsigma, final boolean weighted)
  {
    final int ngather = p.length;
    final int nthread = Runtime.getRuntime().availableProcessors();
    final float[][][] s = new float[ngather][][];
    final AtomicInteger ai = new AtomicInteger();
    Thread[] threads = new Thread[nthread];
    for (int ithread=0; ithread<nthread; ++ithread) {
      threads[ithread] = new Thread() {
        public void run() {
          for (int igather=ai.getAndIncrement(); igather<ngather; 
                   igather=ai.getAndIncrement()) {
            s[igather] = velocitySpectrum(st,sx,p[igather],sv,tsigma,weighted);
          }
        }
      };
    }
    Threads.startAndJoin(threads);
    return s;
  }

  public static float[][] velocitySpectrum(
    Sampling st, Sampling sx, float[][] p,
    Sampling sv, double tsigma, boolean weighted)
  {
    return velocitySpectrum(st,sx,p,sv,tsigma,weighted,null);
  }
  public static float[][] velocitySpectrum(
    Sampling st, Sampling sx, float[][] p, 
    Sampling sv, double tsigma, boolean weighted,
    float[][] b)
  {
    int nv = sv.getCount();
    int nt = st.getCount();
    float[] scl = new float[nt];

    float[][] s = new float[nv][];
    float[][] sn = new float[nv][];
    float[][] sd = new float[nv][];

    for (int iv=0; iv<nv; ++iv) {
      double v = sv.getValue(iv);
      float[][] q = nmo(v,st,sx,p);
      if (weighted) {
        float[] biv = (b!=null)?b[iv]:null;
        s[iv] = semblance(st,sx,v,tsigma,q,biv,scl);
      } else {
        s[iv] = semblance(tsigma,q);
      }
    }
    if (weighted) {
      for (int iv=0; iv<nv; ++iv) {
        for (int it=0; it<nt; ++it) {
          float scli = scl[it];
          float oscl = (scli!=0.0f)?1.0f/scli:1.0f;
          s[iv][it] *= oscl;
        }
      }
    }
    return s;
  }

  private static void splot(float[][] x) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    sp.setSize(500,800);
    PixelsView pv = sp.addPixels(x);
    pv.setColorModel(ColorMap.JET);
    pv.setClips(0.0f,1.0f);
    sp.addColorBar();
  }

  // Weighted semblance
  public static float[] semblance(
    Sampling st, Sampling sx, double vnmo, double tsigma, 
    float[][] q, float[] biv, float[] scl)
  {
    int nx = q.length;
    int nt = q[0].length;
    float[] r = new float[nt];
    float xxsum = 0.0f;
    float xsum = 0.0f;
    for (int ix=0; ix<nx; ++ix) {
      float x = (float)sx.getValue(ix);
      xxsum += x*x;
      xsum += x;
      for (int it=0; it<nt; ++it) {
        r[it] += q[ix][it];
      }
    }
    float xxscl = SCL*xxsum/nx;
    float[] arr = new float[nt];
    float[] arq = new float[nt];
    float[] aqq = new float[nt];
    float[] brr = new float[nt];
    float[] brq = new float[nt];
    float[] bqq = new float[nt];
    float gamma = (float)(1.0/(vnmo*vnmo));
    for (int ix=0; ix<nx; ++ix) {
      float x = (float)sx.getValue(ix);
      float xx = x*x;
      float xxg = xx*gamma;
      for (int it=0; it<nt; ++it) {
        float t0 = (float)st.getValue(it);
        float ti = sqrt(t0*t0+xxg);
        float qi = q[ix][it];
        float ri = r[it];
        if (qi==0.0f) continue;
        float ui = (xx*t0)/(ti*xxscl);
        float rr = ri*ri;
        float rq = ri*qi;
        float qq = qi*qi;
        arr[it] += rr;
        arq[it] += rq;
        aqq[it] += qq;
        brr[it] += ui*rr;
        brq[it] += ui*rq;
        bqq[it] += ui*qq;
      }
    }
    esmooth(tsigma,arr,arr);
    esmooth(tsigma,arq,arq);
    esmooth(tsigma,aqq,aqq);
    esmooth(tsigma,brr,brr);
    esmooth(tsigma,brq,brq);
    esmooth(tsigma,bqq,bqq);
    float[] s = new float[nt];
    for (int it=0; it<nt; ++it) {
      double arri = arr[it];
      double arqi = arq[it];
      double aqqi = aqq[it];
      double brri = brr[it];
      double brqi = brq[it];
      double bqqi = bqq[it];
      double saaa = arri*arqi*aqqi;
      double saab = arri*arqi*bqqi;
      double saba = arri*brqi*aqqi;
      double sbaa = brri*arqi*aqqi;
      double sabb = arri*brqi*bqqi;
      double sbab = brri*arqi*bqqi;
      double sbba = brri*brqi*aqqi;
      double t0 = saaa-saba-saab+sabb;
      double t1 = saaa-sbaa-saab+sbab;
      double t2 = saaa-sbaa-saba+sbba;
      double b = 0.0;

      //if (t0<t1 && t1<t2 || t2<t1 && t1<t0) {
      //  double bnum = sbaa-2.0*saba+saab;
      //  double bden = sabb-2.0*sbab+sbba+bnum;
      //  b = (bden!=0.0)?bnum/bden:1.0;
      //} else {
      //  double bnum = arqi;
      //  double bden = -brqi+bnum;
      //  b = (bden!=0.0)?bnum/bden:0.0;
      //}

      double d = arqi*(sbab-sabb)+brqi*(saba-sbaa);
      //double d1 = (arqi*bqqi-aqqi*brqi)*(-arri*brqi+arqi*brri);
      //System.out.println("d1="+d1+", d2="+d);
      if (d<0.0) {
      //if (t0<t1 && t1<t2 || t2<t1 && t1<t0) {

      //double rrq = arqi/(arqi-brqi);
      //double rrr = arri/(arri-brri);
      //double rqq = aqqi/(aqqi-bqqi);
      //if (rrr<rrq && rrq<rqq || rqq<rrq && rrq<rrr) {

        double bnum = sbaa-2.0*saba+saab;
        double bden = sabb-2.0*sbab+sbba+bnum;
        b = (bden!=0.0)?bnum/bden:1.0;
        //b = 1.0; // XXX
      } else {
        double bnum = arqi;
        double bden = -brqi+bnum;
        b = (bden!=0.0)?bnum/bden:0.0;
        //b = 0.0; // XXX
      }

      if (b<0.0 || b>1.0){
        double snuma = arqi*arqi;
        double sdena = arri*aqqi;
        double sa = (sdena>0.0)?(float)(snuma/sdena):0.0f;
        double snumb = brqi*brqi;
        double sdenb = brri*bqqi;
        double sb = (sdenb>0.0)?(float)(snumb/sdenb):0.0f;
        s[it] = (float)min(sa,sb);
        if (sa<=sb) {
          s[it] = (float)sa;
          if (biv!=null) biv[it] = 0.0f;
        } else {
          s[it] = (float)sb;
          if (biv!=null) biv[it] = 1.0f;
        }
        //s[it] = 1.0f; // XXX
      } else {
        double a = 1.0-b;
        double srri = a*arri+b*brri;
        double srqi = a*arqi+b*brqi;
        double sqqi = a*aqqi+b*bqqi;
        double snum = srqi*srqi;
        double sden = srri*sqqi;
        s[it] = (sden>0.0)?(float)(snum/sden):0.0f;
        if (biv!=null) biv[it] = (float)b;
        //s[it] = 0.0f; // XXX
      }

      // Scaling
      double snum = arqi*arqi;
      double sden = arri*aqqi;
      float sc = (sden>0.0)?(float)(snum/sden):0.0f;
      float sw = s[it];
      if (sc<0.0f) sc = 0.0f;
      if (sc>1.0f) sc = 1.0f;
      float scli = scl[it];
      float sclt = (sc!=0.0f)?sw/sc:0.0f;
      if (sclt>scli) {
        scl[it] = sclt;
      }
      //scl[it] = 1.0f; // XXX

    }
    return s;
  }

  // Conventional semblance
  public static float[] semblance(double tsigma, float[][] q) {
    int nx = q.length;
    int nt = q[0].length;
    float[] sn = new float[nt];
    float[] sd = new float[nt];
    //float[] sx = new float[nt];
    for (int ix=0; ix<nx; ++ix) {
      for (int it=0; it<nt; ++it) {
        float qi = q[ix][it];
        //sx[it] += 1.0f;
        if (qi==0.0f) continue;
        sn[it] += qi;
        sd[it] += qi*qi;
      }
    }
    mul(sn,sn,sn);
    //mul(sd,sx,sd);
    mul(nx,sd,sd);
    esmooth(tsigma,sn,sn);
    esmooth(tsigma,sd,sd);
    //rsmooth(tsigma,sn,sn);
    //rsmooth(tsigma,sd,sd);
    float[] s = sn;
    for (int it=0; it<nt; ++it) {
      s[it] = (sd[it]>0.0)?sn[it]/sd[it]:0.0f;
      if (s[it]<0.0) s[it] = 0.0f;
      if (s[it]>1.0) s[it] = 1.0f;
    }
    //if (num)
    //  return sn;
    //else
    //  return sd;
    return s;
  }

  private static void rsmooth(double tsigma, float[] x, float[] y) {
    double ts = 0.5*(sqrt(1.0+12*tsigma*tsigma)-1.0);
    //System.out.println("ts="+ts);
    RecursiveRectangleFilter rrf = 
      new RecursiveRectangleFilter((int)-ts,(int)ts);
    rrf.apply(x,y);
  }

  private static void xrsmooth(double tsigma, float[] x, float[] y) {
    int nt = x.length;
    if (x==y) x = copy(x);
    double ts = 0.5*(sqrt(1.0+12*tsigma*tsigma)-1.0);
    for (int it=0; it<nt; ++it) {
      double s = 0.0;
      int count = 0;
      for (int jt=max(0,(int)(it-ts/2)); jt<min((int)(it+ts/2+1),nt); ++jt) {
        s += x[jt];
        ++count;
      }
      y[it] = (float)(s/count);
    }
  }

  public static float[][] nmo(
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

  public static float[][] nmo(
    float[] vnmo, Sampling st, Sampling sx, float[][] p)
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
      for (int it=0; it<nt; ++it) {
        double t0 = st.getValue(it);
        double v0 = vnmo[it];
        double xxg = (x*x)/(v0*v0);
        t[it] = (float)sqrt(t0*t0+xxg); 
      }
      si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
    }

    return q;
  }

  public static float[][] addRandomNoise(float r, float[][] p) {
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

  public static float[][] makeRickerGather(
    double fpeak, float[] vnmo, Sampling st, Sampling sx) 
  {
    return makeRickerGather(new Random(314159),fpeak,vnmo,st,sx);
  }

  public static float[][] makeRickerGather(
    Random random, double fpeak, float[] vnmo, Sampling st, Sampling sx) 
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
    return p;
  }

  public static float[][] xmakeRickerGather(
    double fpeak, float[] vnmo, Sampling st, Sampling sx) 
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    double lt = st.getLast();
    int nx = sx.getCount();
    double dx = sx.getDelta();
    double fx = sx.getFirst();
    Random random = new Random(314159);
    float[][] p = new float[nx][nt];
    double thalf =  1.0/fpeak;
    //for (int jt=0; jt<nt; ++jt) {
    for (int jt=nt/2; jt<nt; jt+=nt) {
      double t0 = st.getValue(jt);
      //double a0 = 2.0*random.nextDouble()-1.0;
      double a0 = 1.0;
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
    return p;
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

      //XXX
      //double t = ft+it*dt;
      //if (t==3.5) 
      //  System.out.println("it="+it+", t="+t+", v="+v[it]);

    }
    return v;
  }

  public static void esmooth(double sigma, float[] p, float[] q) {
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE);
    //ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_VALUE);
    //ref.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE);
    ref.apply(p,q);
  }

  public static void xesmooth(double sigma, float[] p, float[] q) {
    if (p==q) p = copy(p);
    float sigmas = (float)(sigma*sigma);
    float a = (sigmas>0.0)?(1.0f+sigmas-sqrt(1.0f+2.0f*sigmas))/sigmas:0.0f;
    float b = (1.0f-a)/(1.0f+a);
    int n = p.length;
    float qi = b/(1.0f-a)*p[0];
    q[0] = qi;
    for (int i=1; i<n; ++i) {
      qi = a*qi+b*p[i];
      q[i] = qi;
    }
    qi = a*b/(1.0f-a)*p[n-1];
    q[n-1] += qi;
    for (int i=n-2; i>=0; --i) {
      qi = a*(qi+b*p[i+1]);
      q[i] += qi;
    }
  }

  private static float ricker(double fpeak, double time) {
    double x = PI*fpeak*time;
    double xx = x*x;
    return (float)((1.0-2.0*xx)*exp(-xx));
  }

  /////////////////////////////////////////////////////////////////////////

  private static void testGather() {
    Sampling st = new Sampling(1001,0.004,0.0);
    Sampling sx = new Sampling(60,0.050,0.050);
    Sampling sv = new Sampling(101,0.020,1.5);
    int nt = st.getCount();
    int nv = sv.getCount();
    float[][] b = new float[nv][nt];
    double[] vp = {2.00,3.00}; // velocity of primaries
    double[] vm = {1.98,2.70}; // velocity of multiples
    double fpeak = 25.0;
    double tsigma = 8.0;
    float snr = 1.0e6f;
    float[] vps = Velan.makeLinearVelocity(vp[0],vp[1],st);
    float[] vms = Velan.makeLinearVelocity(vm[0],vm[1],st);
    float[][] p = Velan.makeRickerGather(fpeak,vps,st,sx);
    p  = add(p,Velan.makeRickerGather(fpeak,vms,st,sx)); // multiples
    p = Velan.addRandomNoise(snr,p); // noise
    float[][] sc = velocitySpectrum(st,sx,p,sv,tsigma,false,b);
    float[][] sw = velocitySpectrum(st,sx,p,sv,tsigma,true,b);
    cplot(p,st,sx);
    splot(sc,st,sv);
    splot(sw,st,sv);
    bplot(b,st,sv);
  }

  private static void splot(float[][] f, Sampling s1, Sampling s2) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    sp.setHLimits(s2.getFirst(),s2.getLast());
    sp.setVLimits(s1.getFirst(),s1.getLast());
    sp.setHLabel("Velocity (km/s)");
    sp.setVLabel("Time (s)");
    sp.addColorBar("Semblance");
    sp.setSize(600,800);
    PixelsView pv = sp.addPixels(s1,s2,f);
    pv.setColorModel(ColorMap.JET);
  }

  private static void cplot(float[][] f, Sampling s1, Sampling s2) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    sp.setHLimits(s2.getFirst(),s2.getLast());
    sp.setVLimits(s1.getFirst(),s1.getLast());
    sp.setHLabel("Offset (km)");
    sp.setVLabel("Time (s)");
    sp.addColorBar("Amplitude");
    sp.setSize(600,800);
    PixelsView pv = sp.addPixels(s1,s2,f);
    pv.setColorModel(ColorMap.GRAY);
  }

  private static void bplot(float[][] f, Sampling s1, Sampling s2) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    sp.setHLimits(s2.getFirst(),s2.getLast());
    sp.setVLimits(s1.getFirst(),s1.getLast());
    sp.setHLabel("Velocity (km/s)");
    sp.setVLabel("Time (s)");
    sp.addColorBar("b value");
    sp.setSize(600,800);
    PixelsView pv = sp.addPixels(s1,s2,f);
    pv.setColorModel(ColorMap.JET);
  }

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        testGather();
      }
    });
  }

}
