package lss.vel;
import java.lang.StringBuilder;
import edu.mines.jtk.lapack.DMatrix;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Coefficients for finite-difference stencils.
 */
public class Coefficients {

  public static void main(String[] args) {
    int derivativeOrder = 2;
    int halfStencilOrder = 1; // half stencil order
    getCoefficients(derivativeOrder,halfStencilOrder);
  }

  public static void getCoefficients(
    int derivativeOrder, int halfStencilOrder)
  {
    double[] c = computeCoefficients(derivativeOrder,halfStencilOrder);
    printCoefficients(c);
  }

  ////////////////////////////////////////////////////////////////////////////////////
  // private

  private static double[] computeCoefficients(
    int derivativeOrder, int halfStencilOrder)
  {
    int stencilOrder = 2*halfStencilOrder;
    int n = stencilOrder+1;
    DMatrix mata = new DMatrix(n,n);
    double r = halfStencilOrder;
    for (int j=0; j<n; ++j) {
      for (int i=0; i<n; ++i) {
        mata.set(i,j,pow(r,i)/factorial(i));
      }
      r -= 1.0;
    }
    DMatrix mate = new DMatrix(n,1);
    if (derivativeOrder==1) {
      mate.set(1,0,1.0);
    } else if (derivativeOrder==2) {
      mate.set(2,0,1.0);
    }
    DMatrix matc = mata.solve(mate);
    double[] c = matc.getArray();
    return reverse(c);
  }

  private static int factorial(int k) {
    if (k==0) {
      return 1;
    } else {
      return k*factorial(k-1);
    }
  }

  private static void printCoefficients(double[] c) {
    int n = c.length;
    int factor = 1000;
    for (int i=0; i<n; ++i) {
      System.out.println(toFraction(c[i],factor));
    }
  }

  // http://stackoverflow.com/questions/5968636/converting-a-float-
  // into-a-string-fraction-representation/5968920#5968920
  private static String toFraction(double d, int factor) {
    StringBuilder sb = new StringBuilder();
    if (d < 0) {
      sb.append('-');
      d = -d;
    } else {
      sb.append(' ');
    }
    long l = (long) d;
    //if (l != 0) sb.append(l); d -= l; // prints whole numbers separately
    double error = Math.abs(d);
    int bestDenominator = 1;
    for(int i=2;i<=factor;i++) {
      double error2 = Math.abs(d - (double) Math.round(d * i) / i);
      if (error2 < error) {
        error = error2;
        bestDenominator = i;
      }
    }
    if (bestDenominator > 1) {
      sb.append(Math.round((d)*bestDenominator)).append('/').append(bestDenominator);
    } else {
      sb.append('0');
    }
    return sb.toString();
  }
}
