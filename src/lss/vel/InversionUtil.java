package lss.vel;

import edu.mines.jtk.io.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class InversionUtil {

  // u: state variable
  // a: adjoint-state variable
  public static float[][] imagingCondition(
    double dt, float[][][] u, float[][][] a)
  {
    int n1 = u[0][0].length;
    int n2 = u[0].length;
    int n3 = u.length;
    float[][] y = new float[n2][n1];
    float oden = 1.0f/(float)(dt*dt);

    // i3=0
    for (int i2=0; i2<n2; i2++) {
      for (int i1=0; i1<n1; i1++) {
        float d2ui = oden*
          (u[0][i2][i1]-2.0f*u[1][i2][i1]+u[2][i2][i1]);
        y[i2][i1] += d2ui*a[0][i2][i1];
      }
    }
    
    // i3=1:n3-1
    for (int i3=1; i3<n3-1; i3++) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1; i1++) {
          float d2ui = oden*
            (u[i3+1][i2][i1]-2.0f*u[i3][i2][i1]+u[i3-1][i2][i1]);
          y[i2][i1] += d2ui*a[i3][i2][i1];
        }
      }
    }

    // i3=n3-1
    for (int i2=0; i2<n2; i2++) {
      for (int i1=0; i1<n1; i1++) {
        float d2ui = oden*
          (u[n3-1][i2][i1]-2.0f*u[n3-2][i2][i1]+u[n3-3][i2][i1]);
        y[i2][i1] += d2ui*a[n3-1][i2][i1];
      }
    }

    return y;
  }
  public static float[][] imagingCondition(float[][][] u, float[][][] a) {
    int n1 = u[0][0].length;
    int n2 = u[0].length;
    int n3 = u.length;
    float[][] y = new float[n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1; i1++) {
          y[i2][i1] -= u[i3][i2][i1]*a[i3][i2][i1];
        }
      }
    }
    return y;
  }

  public static float dot(float[][] a, float[][] b) {
    float y = 0.0f;
    for (int i2=0, n2=a.length; i2<n2; i2++)
      y += dot(a[i2],b[i2]);
    return y;
  }

  public static float dot(float[] a, float[] b) {
    return sum(mul(a,b));
  }

  public static float[][] rotateGradient(final int na, final float[][] x) {
    int n1 = x[0].length;
    int n2 = x.length;
    final double dtheta = 2.0*PI/na;
    final SincInterpolator si = new SincInterpolator();
    float[][] g = Parallel.reduce(na-1,new Parallel.ReduceInt<float[][]>() {
      public float[][] compute(int ia) {
        double theta = dtheta*(ia+1);
        return rotateGradientOnce(theta,si,x);
      }
      public float[][] combine(float[][] f1, float[][] f2) {
        return add(f1,f2);
      }
    });

    return add(x,g);
  }
  public static float[][] rotateGradientOnce(
    double theta, SincInterpolator si, float[][] x)
  {
    int n1 = x[0].length;
    int n2 = x.length;
    float[][] y = new float[n2][n1];
    double r1 = (n1-1)/2.0;
    double r2 = (n2-1)/2.0;
    double s = sin(theta);
    double c = cos(theta);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        double q1 = i1-r1;
        double q2 = i2-r2;
        double t1 = c*q1+s*q2;
        double t2 = c*q2-s*q1;
        y[i2][i1] = si.interpolate(n1,1.0,0.0,n2,1.0,0.0,x,t1+r1,t2+r2);
      }
    }
    return y;
  }

  public static void muteTopAndBottom(float[][] f) {
    muteTop(f);
    muteBottom(f);
  }
  public static void muteTopAndBottom(float[][] f, int zero) {
    muteTop(f,zero);
    muteBottom(f,zero);
  }
  public static void muteTopAndBottom(float[][] f, int zero, int taper) {
    muteTop(f,zero,taper);
    muteBottom(f,zero,taper);
  }

  public static void muteTop(float[][] f) {
    muteTop(f,5);
  }
  public static void muteTop(float[][] f, int zero) {
    muteTop(f,zero,5);
  }
  public static void muteTop(float[][] f, int zero, int taper) {
    int nx = f.length;
    for (int ix=0; ix<nx; ++ix) {
      for (int iz=0; iz<zero; ++iz)
        f[ix][iz] = 0.0f;
      for (int iz=zero; iz<zero+taper; ++iz)
        f[ix][iz] *= (float)(iz+1-zero)/(1.0f+taper);
    }
  }

  public static void muteBottom(float[][] f) {
    muteBottom(f,5);
  }
  public static void muteBottom(float[][] f, int zero) {
    muteBottom(f,zero,5);
  }
  public static void muteBottom(float[][] f, int zero, int taper) {
    int nx = f.length;
    int nz = f[0].length;
    for (int ix=0; ix<nx; ++ix) {
      for (int iz=nz-zero; iz<nz; ++iz)
        f[ix][iz] = 0.0f;
      for (int iz=nz-zero-taper, nzmz=nz-zero; iz<nzmz; ++iz)
        f[ix][iz] *= 1.0f-(float)(iz+1-nz+zero+taper)/(1.0f+taper);
    }
  }

  public static void muteLeftAndRight(float[][] f) {
    muteLeft(f);
    muteRight(f);
  }
  public static void muteLeftAndRight(float[][] f, int zero) {
    muteLeft(f,zero);
    muteRight(f,zero);
  }
  public static void muteLeftAndRight(float[][] f, int zero, int taper) {
    muteLeft(f,zero,taper);
    muteRight(f,zero,taper);
  }

  public static void muteLeft(float[][] f) {
    muteLeft(f,5);
  }
  public static void muteLeft(float[][] f, int zero) {
    muteLeft(f,zero,5);
  }
  public static void muteLeft(float[][] f, int zero, int taper) {
    int nx = f.length;
    int nz = f[0].length;
    for (int ix=0; ix<zero; ++ix) {
      mul(0.0f,f[ix],f[ix]);
    }
    for (int ix=zero, zpt=zero+taper; ix<zpt; ++ix) {
      float r = (float)(ix+1-zero)/(1.0f+taper);
      mul(r,f[ix],f[ix]);
    }
  }

  public static void muteRight(float[][] f) {
    muteRight(f,5);
  }
  public static void muteRight(float[][] f, int zero) {
    muteRight(f,zero,5);
  }
  public static void muteRight(float[][] f, int zero, int taper) {
    int nx = f.length;
    int nz = f[0].length;
    for (int ix=nx-zero; ix<nx; ++ix) {
      mul(0.0f,f[ix],f[ix]);
    }
    for (int ix=nx-zero-taper, nxmz=nx-zero; ix<nxmz; ++ix) {
      float r = 1.0f-(float)(ix+1-nx+zero+taper)/(1.0f+taper);
      mul(r,f[ix],f[ix]);
    }
  }

  public static void mirror(float[][] f) {
    int n2 = f.length;
    float[][] temp = copy(f);
    for (int i2=0, j2=n2-1; i2<n2; i2++, j2--)
      add(f[i2],temp[j2],f[i2]);
    mul(0.5f,f,f);
  }

  public static void window(float[][] f) {
    int nx = f.length;
    for (int ix=0; ix<nx; ++ix) {
      window(f[ix]);
    }
  }
//  public static void window(float[] f) {
//    int nt = f.length;
//    int kt = 0;
//    float max = 0.0f;
//    for (int it=0; it<nt; ++it) {
//      float fi = f[it];
//      if (fi>max) {
//        max = fi;
//        kt = it;
//      }
//    }
//    int delay = 150;
//    float taper = 50.0f;
//    kt += delay;
//    for (int it=kt; it<nt; ++it) {
//      float r = 1.0f-(it-kt)/taper;
//      if (r<0.0f) r=0.0f;
//      f[it] *= r;
//    }
//  }

  // Picking using Wong (2009)
  public static void window(float[] f) {
    window(300,f);
  }
  public static void window(int ne, float[] f) {
    int n1 = f.length;
    int a1 = pickFirstArrival(ne,f);
    int delay = 300;
    float taper = 50.0f;
    a1 += delay;
    for (int i1=a1; i1<n1; ++i1) {
      float r = 1.0f-(i1-a1)/taper;
      if (r<0.0f) r=0.0f;
        f[i1] *= r;
    }
  }
  private static int pickFirstArrival(int ne, float[] x) {
    int n1 = x.length;
    float[] er = new float[n1]; // energy ratio
    for (int i1=0; i1<n1; ++i1) {
      float num = 0.0f;
      float den = 0.0f;
      for (int ie=0; ie<ne; ++ie) {
        int im = i1-ie;
        int ip = i1+ie;
        float xm = (im<0   )?0.5f*(x[0   ]+x[1   ]):x[im];
        float xp = (ip>n1-1)?0.5f*(x[n1-1]+x[n1-2]):x[ip];
        num += xp*xp;
        den += xm*xm;
      }
      float eri = (den>0.0f)?abs(x[i1])*num/den:0.0f;
      er[i1] = eri*eri*eri;
    }
    int[] fa = new int[1];
    max(er,fa);
    return fa[0];
  }

  public static void readFloats(String name, float[] image) {
    try {
      ArrayInputStream ais = new ArrayInputStream(name);
      ais.readFloats(image);
      ais.close();
    } catch (java.io.IOException e) {
      System.out.println(e.toString());
    }
  }
  public static void readFloats(String name, float[][] image) {
    try {
      ArrayInputStream ais = new ArrayInputStream(name);
      ais.readFloats(image);
      ais.close();
    } catch (java.io.IOException e) {
      System.out.println(e.toString());
    }
  }

  public static void writeFloats(String name, float[] image) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(name);
      aos.writeFloats(image);
      aos.close();
    } catch (java.io.IOException e) {
      System.out.println(e.toString());
    }
  }
  public static void writeFloats(String name, float[][] image) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(name);
      aos.writeFloats(image);
      aos.close();
    } catch (java.io.IOException e) {
      System.out.println(e.toString());
    }
  }

  // Normalize
  public static float[] normal(float[] x) {
    return div(x,max(abs(x)));
  }
  public static void normal(float[] x, float[] y) {
    div(x,max(abs(x)),y);
  }
  public static float[][] normal(float[][] x) {
    return div(x,max(abs(x)));
  }
  public static void normal(float[][] x, float[][] y) {
    div(x,max(abs(x)),y);
  }
  public static float[][][] normal(float[][][] x) {
    return div(x,max(abs(x)));
  }
  public static void normal(float[][][] x, float[][][] y) {
    div(x,max(abs(x)),y);
  }

}
