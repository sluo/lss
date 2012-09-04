package lss.vel;

import java.awt.*;
import java.util.*;
import javax.swing.*;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.io.ArrayOutputStream;
import static edu.mines.jtk.util.ArrayMath.*;

public class VelocityModel { 

  public static float[][] marmousi(
    Sampling s1, Sampling s2, double sigma)
  {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    double d1 = s1.getDelta();
    double d2 = s2.getDelta();
    double f1 = s1.getFirst();
    double f2 = s2.getFirst();
    float[][] m = marmousi(sigma);
    float[][] y = fillfloat(1500f,n1,n2);
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    si.setUniform(m[0].length,0.004,0.0,m.length,0.004,0.0,m);
    int k1 = (int)(0.1/d1); // water
    int i2=0; double t2=0.0; while (t2<2301*0.004 && i2<n2) {
      t2 = f2+i2*d2;
      int i1=0; double t1=0.0; while(t1<751*0.004 && k1+i1<n1) {
        t1 = f1+i1*d1;
        y[i2][k1+i1] = si.interpolate(t1,t2);
        ++i1;
      }
      ++i2;
    }
    return div(y,1000.f);
  }

  private static float[][] marmousi(double sigma) {
    float[][] m = new float[2301][751];
    float[][] t = new float[2301][743];
    readImage("/data/seis/marmousi/marmousi.dat",m);
    copy(743,2301,8,0,m,0,0,t);
    if (sigma>=1.0) {
      int e = 3*(int)sigma;
      m = extendVelocityModel(e,t);
      reciprocal(m);
      new RecursiveGaussianFilter(sigma).apply00(m,m); // smooth slowness
      reciprocal(m);
      copy(743,2301,e,e,m,0,0,t);
    }
    return t;
  }

  //////////////////////////////////////////////////////////////////////////
  // private

  private static void reciprocal(float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float xi = x[i2][i1];
        x[i2][i1] = (xi!=0.0f)?1.0f/xi:0.0f;
      }
    }
  }

  private static float[][] extendVelocityModel(
    int e, float[][] c)
  {
    int n1 = c[0].length;
    int n2 = c.length;
    int n1p = n1+2*e;
    int n2p = n2+2*e;
    float[][] v = new float[n2p][n1p];
    copy(n1,n2,0,0,c,e,e,v);
    for (int i2=e; i2<n2+e; ++i2) {
      for (int i1=0, j1=n1+e; i1<e; ++i1, ++j1) {
        v[i2][i1] = v[i2][e];
        v[i2][j1] = v[i2][n1-1+e];
      }
    }
    for (int i2=0, j2=n2+e; i2<e; ++i2, ++j2) {
      copy(v[e],v[i2]);
      copy(v[n2-1+e],v[j2]);
    }
    return v;
  }

  private static void readImage(String name, float[][] y) {
    try {
      ArrayInputStream ais = new ArrayInputStream(name);
      ais.readFloats(y);
      ais.close();
    } catch (java.io.IOException e) {
      e.printStackTrace();
    }
  }

  private static void writeImage(String name, float[][] y) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(name);
      aos.writeFloats(y);
      aos.close();
    } catch (java.io.IOException e) {
      e.printStackTrace();
    }
  }
}
