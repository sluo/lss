package lss.util;

import java.util.Arrays;
import java.util.Comparator;
import java.util.ArrayList;

import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.MedianFinder;
import static edu.mines.jtk.util.ArrayMath.*;

// TESTING
import java.util.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;

/**
 * A k-d tree.
 * A k-d tree is a binary tree in which every node is a k-dimensional point.
 * This implementation is adapted from Friedman, J.H., J.L. Bentley, and 
 * R.A. Finkel 1977, An Algorithm for Finding Best Matches in Logarithmic 
 * Expected Time: ACM Transactions on Mathematical Software, 3, 209--226.
 * @author Simon Luo
 * @version 2011.12.08
 */
public class KdTree {

  /**
   * Constructs a k-d tree.
   * @param x input array of sample coordinates.
   */
  public KdTree(float[][] x) {
    int n = x[0].length;
    _k = x.length;
    _x = x;
    _i = rampint(0,1,n);
    _root = new Node(0,n-1);
  }

  /**
   * Finds the nearest sample.
   * @param xq coordinates of query point.
   * @return index of the nearest sample.
   */
  public int findNearest(float[] xq) {
    return new Search(xq).findNearest();
  }

  /**
   * Finds the nearest sample by checking 
   * the distance to all samples.
   * @param xq coordinates of query point.
   * @return index of the nearest sample.
   */
  public int findNearestSlow(float[] xq) {
    int n = _x[0].length;
    int nearest = 0;
    float rmin = distance(nearest,_x,xq);
    for (int i=1; i<n; ++i) {
      float r = distance(i,_x,xq);
      if (r<rmin) {
        nearest = i;
        rmin = r;
      }
    }
    return nearest;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _k;
  private int[] _i;
  private float[][] _x;
  private Node _root;

  private class Search {
    public Search(float[] x) {
      xq = x;
      xmin = fillfloat(-FLT_MAX,_k);
      xmax = fillfloat(FLT_MAX,_k);
    }
    public int findNearest() {
      r = FLT_MAX;
      search(_root);
      return nearest;
    }
    private boolean search(Node node) {
      // if leaf node, check distances
      if (node.isLeaf) {
        int p = node.p;
        int q = node.q;
        for (int i=p; i<=q; ++i)
          checkDistance(_i[i]);
        // return true if nearest point has been found
        return (ballWithinBounds())?true:false;
      } 
      // else continue searching
      else {
        Node left = node.left;
        Node right = node.right;
        int d = node.d; // split index (discriminator)
        float v = node.v; // split value (partition)
        float xqd = xq[d];
        // if query point is left of split value, search left
        if (xqd<=v) {
          float t = xmax[d];
          xmax[d] = v;
          if (search(left))
            return true;
          xmax[d] = t;
        }
        // else search right
        else {
          float t = xmin[d];
          xmin[d] = v;
          if (search(right))
            return true;
          xmin[d] = t;
        }
        // if bounds overlap ball, search other child
        if (xqd<=v) {
          float t = xmin[d];
          xmin[d] = v;
          if (boundsOverlapBall())
            if (search(right))
              return true;
          xmin[d] = t;
        } else {
          float t = xmax[d];
          xmax[d] = v;
          if (boundsOverlapBall())
            if (search(left))
              return true;
          xmax[d] = t;
        }
        // return true if nearest point has been found
        return (ballWithinBounds())?true:false;
      }
    }
    private void checkDistance(int i) {
      float t = distance(i,_x,xq);
      if (t<r) {
        r = t;
        nearest = i;
      }
    }
    private boolean ballWithinBounds() {
      for (int d=0; d<_k; ++d) {
        if (coordinateDistance(d,xq,xmin)<=r ||
            coordinateDistance(d,xq,xmax)<=r) {
          return false;
        }
      }
      return true;
    }
    private boolean boundsOverlapBall() {
      float s = 0.0f;
      for (int d=0; d<_k; ++d) {
        float xqd = xq[d];
        // if lower than low boundary
        if (xqd<xmin[d]) {
          s += coordinateDistance(d,xq,xmin);
          if (s>r) return false;
        }
        // if higher than high boundary
        if (xqd>xmax[d]) {
          s += coordinateDistance(d,xq,xmax);
          if (s>r) return false;
        }

      }
      return true;
    }
    private int nearest; // index of nearest sample
    private float r; // distance to nearest sample
    private float[] xq,xmin,xmax;
  }

  private class Node {
    public Node(int p, int q) {
      // This node is a leaf.
      if (q-p<NLEAF) {
        this.p = p;
        this.q = q;
        isLeaf = true;
        left = null;
        right = null;
      }
      // Too many samples, so split.
      else {
        d = findSpreadestCoordinate(p,q,_i,_x);
        int m = medianSplit(p,q,_i,_x[d]);
        v = _x[d][_i[m]];
        isLeaf = false;
        left = new Node(p,m);
        right = new Node(m+1,q);
      }
    }
    public Node left,right; // if leaf, null
    public boolean isLeaf; // if leaf, true
    public int p; // if leaf, lower index of array of samples i[]
    public int q; // if leaf, upper index of array of samples i[]
    public int d; // if not leaf, dimension that is split (discriminator)
    public float v; // if not leaf, value of the split (partition)
    private static final int NLEAF = 8; // TODO: optimal?
  }

  private static float distance(int i, float[][] x, float[] xq) {
    int k = x.length;
    double ds = 0.0f;
    for (int ik=0; ik<k; ++ik) {
      double d = x[ik][i]-xq[ik];
      ds += d*d;
    }
    return (float)ds;
  }

  private static float coordinateDistance(
    int i, float[] x1, float[] x2)
  {
    double d = x1[i]-x2[i];
    return (float)(d*d);
  }

  // Median search algorithm from Wirth, N., Algorithms + Data 
  // Structures = Programs. Implementation adapted from Devillard, N.,
  // 1998, Fast median search: an ANSI C implementation.
  public static int medianSplit(int p, int q, int[] i, float[] x) {
    int m = (p+q)/2;
    while (p<q) {
      float xmid = x[i[m]];
      int r = p;
      int s = q;
      do {
        while (x[i[r]]<xmid) ++r;
        while (xmid<x[i[s]]) --s;
        if (r<=s) {
          int t = i[r];
          i[r++] = i[s];
          i[s--] = t;
        }
      } while (r<=s);
      if (s<m) p=r;
      if (m<r) q=s;
    }
    return m;
  }
  
  private static int findSpreadestCoordinate(
    int p, int q, int[] i, float[][] x)
  {
    int nk = x.length;
    int spreadest = 0;
    float maxspread = 0.0f;
    for (int k=0; k<nk; ++k) {
      float max = x[k][i[p]];
      float min = max;
      for (int ii=p+1; ii<=q; ++ii) {
        float xi = x[k][i[ii]];
        if (xi>max) max = xi;
        if (xi<min) min = xi;
      }
      float spread = max-min;
      if (spread>maxspread) {
        maxspread = spread;
        spreadest = k;
      }
    }
    return spreadest;
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    goMedian();
  }

  private static void goMedian() {
    int n = 100;
    int[] i = rampint(0,1,n);
    Random r = new Random(3);
    Random s = new Random(4);
    float[] x = randfloat(r,n); mul(10.0f,x,x);
    float[] y = randfloat(s,n); mul(10.0f,y,y);
    for (int j=0; j<4; ++j) {
      x[j] = 0.0f;
      y[j] = 0.0f;
    }

    int p = 0;
    //int p = 2;
    int q = n-1;
    int m = medianSplit(p,q,i,x);
    System.out.println("median="+x[i[m]]);
    for (int j=0; j<n; ++j)
      System.out.println(x[i[j]]);
  }

}
