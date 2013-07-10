package lss.vel;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Computes a measure of correlation between two 1-D sequences.
 * The similarity between two sequences can be quantified with their
 * crosscorrelation, while their dissimilarity can be quantified with
 * either the sum of absolute differences or the sum of squared differences
 * between the two sequences as a function of lag. Adapted from 
 * {@link http://dhale.github.io/jtk/api/edu/mines/jtk/dsp/Conv.html}.
 * @version 2013.07.01
 * @author Simon Luo
 */
public class Correlation {

  public enum Measure {
    PRODUCT,
    ABSOLUTE_DIFFERENCE,
    SQUARED_DIFFERENCE
  }

  public Correlation() {
    this(Measure.PRODUCT);
  }

  public Correlation(Measure measure) {
    _measure = measure;
  }

  public void correlate(
    int lx, int kx, float[] x,
    int ly, int ky, float[] y,
    int lz, int kz, float[] z)
  {
    boolean copy = x==y;
    x = reverse(lx,x,copy);
    kx = 1-kx-lx;
    convFast(lx,kx,x,ly,ky,y,lz,kz,z);
    if (!copy)
      reverse(lx,x,false);
  }

  //////////////////////////////////////////////////////////////////////////
  // private

  private Measure _measure;

  private float function(float a, float b) {
    if (_measure==Measure.ABSOLUTE_DIFFERENCE) {
      return abs(a-b);
    } else if (_measure==Measure.SQUARED_DIFFERENCE) {
      return (a-b)*(a-b);
    } else {
      return a*b;
    }
  }

  private void convFast(
    int lx, int kx, float[] x,
    int ly, int ky, float[] y,
    int lz, int kz, float[] z)
  {
    // If necessary, swap x and y so that x is the shorter sequence.
    // This simplifies the logic below.
    if (lx>ly) {
      int lt = lx;  lx = ly;  ly = lt;
      int kt = kx;  kx = ky;  ky = kt;
      float[] t = x;  x = y;  y = t;
    }

    // Bounds for index i.
    int imin = kz-kx-ky;
    int imax = imin+lz-1;

    // Variables that we expect to reside in registers.
    int i,ilo,ihi,j,jlo,jhi,iz;
    float sa,sb,xa,xb,ya,yb;

    // Off left: imin <= i <= -1
    ilo = imin;
    ihi = min(-1,imax);
    for (i=ilo,iz=i-imin; i<=ihi; ++i,++iz)
      z[iz] = 0.0f;

    // Rolling on: 0 <= i <= lx-2 and 0 <= j <= i
    ilo = max(0,imin);
    ihi = min(lx-2,imax);
    jlo = 0;
    jhi = ilo;
    for (i=ilo,iz=i-imin; i<ihi; i+=2,iz+=2,jhi+=2) {
      sa = 0.0f;
      sb = 0.0f;
      yb = y[i-jlo+1];
      for (j=jlo; j<jhi; j+=2) {
        xa = x[j];
        sb += function(xa,yb);
        ya = y[i-j];
        sa += function(xa,ya);
        xb = x[j+1];
        sb += function(xb,ya);
        yb = y[i-j-1];
        sa += function(xb,yb);
      }
      xa = x[j];
      sb += function(xa,yb);
      if (j==jhi) {
        ya = y[i-j];
        sa += function(xa,ya);
        xb = x[j+1];
        sb += function(xb,ya);
      }
      z[iz  ] = sa;
      z[iz+1] = sb;
    }
    if (i==ihi) {
      jlo = 0;
      jhi = i;
      sa = 0.0f;
      for (j=jlo; j<=jhi; ++j)
        sa += function(x[j],y[i-j]);
      z[iz] = sa;
    }

    // Middle: lx-1 <= i <= ly-1 and 0 <= j <= lx-1
    ilo = max(lx-1,imin);
    ihi = min(ly-1,imax);
    jlo = 0;
    jhi = lx-1;
    for (i=ilo,iz=i-imin; i<ihi; i+=2,iz+=2) {
      sa = 0.0f;
      sb = 0.0f;
      yb = y[i-jlo+1];
      for (j=jlo; j<jhi; j+=2) {
        xa = x[j];
        sb += function(xa,yb);
        ya = y[i-j];
        sa += function(xa,ya);
        xb = x[j+1];
        sb += function(xb,ya);
        yb = y[i-j-1];
        sa += function(xb,yb);
      }
      if (j==jhi) {
        xa = x[j];
        sb += function(xa,yb);
        ya = y[i-j];
        sa += function(xa,ya);
      }
      z[iz  ] = sa;
      z[iz+1] = sb;
    }
    if (i==ihi) {
      sa = 0.0f;
      for (j=jlo; j<=jhi; ++j)
        sa += function(x[j],y[i-j]);
      z[iz] = sa;
    }

    // Rolling off: ly <= i <= lx+ly-2 and i-ly+1 <= j <= lx-1
    ilo = max(ly,imin);
    ihi = min(lx+ly-2,imax);
    jlo = ihi-ly+1;
    jhi = lx-1;
    for (i=ihi,iz=i-imin; i>ilo; i-=2,iz-=2,jlo-=2) {
      sa = 0.0f;
      sb = 0.0f;
      yb = y[i-jhi-1];
      for (j=jhi; j>jlo; j-=2) {
        xa = x[j];
        sb += function(xa,yb);
        ya = y[i-j];
        sa += function(xa,ya);
        xb = x[j-1];
        sb += function(xb,ya);
        yb = y[i-j+1];
        sa += function(xb,yb);
      }
      xa = x[j];
      sb += function(xa,yb);

      if (j==jlo) {
        ya = y[i-j];
        sa += function(xa,ya);
        xb = x[j-1];
        sb += function(xb,ya);
      }
      z[iz  ] = sa;
      z[iz-1] = sb;
    }
    if (i==ilo) {
      jlo = i-ly+1;
      jhi = lx-1;
      sa = 0.0f;
      for (j=jhi; j>=jlo; --j)
    	sa += function(x[j],y[i-j]);
      z[iz] = sa;
    }
	
    // Off right: lx+ly-1 <= i <= imax
    ilo = max(lx+ly-1,imin);
    ihi = imax;
    for (i=ilo,iz=i-imin; i<=ihi; ++i,++iz)
      z[iz] = 0.0f;

  }

  private static float[] reverse(int n1, float[] z, boolean copy) {
    if (copy) 
      z = copy(n1,z);
    for (int i1=0,j1=n1-1; i1<j1; ++i1,--j1) {
      float zt = z[i1];
      z[i1] = z[j1];
      z[j1] = zt;
    }
    return z;
  }

}
