package lss.vel;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import dnp.*;

public class Fletcher {

  public static void main(String[] args) {
    //Sampling sz = new Sampling(501,0.010,0.0);
    //Sampling sx = new Sampling(501,0.010,0.0);
    int nz = 501;
    int nx = 501;
    float[][] s = new float[nx][nz];
    for (int ix=0; ix<nx; ++ix) {
      for (int iz=0; iz<nz/2; ++iz) {
        s[ix][iz] = 1.0f;
      }
      for (int iz=nz/2; iz<nz; ++iz) {
        s[ix][iz] = 0.5f;
      }
    } 

    // RHS
    float[][] b = makeRhs(s);
    VecArrayFloat2 vb = new VecArrayFloat2(b);

    // Unknown vector (model)
    float[][] x = new float[nx][nz];
    //copy(s,x);
    VecArrayFloat2 vx = new VecArrayFloat2(x);

    // LHS operator
    A2 a2 = new A2();
    CgSolver cg = new CgSolver(0.001f,1000);
    
    // Solver
    cg.solve(a2,vb,vx);

    SimplePlot.asPixels(s);
    SimplePlot.asPixels(b);
    SimplePlot.asPixels(x);

  }

  private static float[][] makeRhs(float[][] f) {
    int nz = f[0].length;
    int nx = f.length;
    float[][] g = new float[nx][nz];

    float r = 0.5f;
    float s = 0.5f;
    for (int i2=1; i2<nx; ++i2) {
      for (int i1=1; i1<nz; ++i1) {
        float f1r = f[i2  ][i1  ]-f[i2  ][i1-1];
        float f1s = f[i2-1][i1  ]-f[i2-1][i1-1];
        //float f2r = f[i2  ][i1  ]-f[i2-1][i1  ];
        //float f2s = f[i2  ][i1-1]-f[i2-1][i1-1];
        float f2r = 0.0f;
        float f2s = 0.0f;
        float g1 = r*f1r+s*f1s;
        float g2 = r*f2r+s*f2s; // gather ends
        float g1r = g1*r;
        float g1s = g1*s;
        float g2r = g2*r;
        float g2s = g2*s;
        g[i2  ][i1  ]  = g1r+g2r;
        g[i2  ][i1-1] -= g1r-g2s;
        g[i2-1][i1  ] += g1s-g2r;
        g[i2-1][i1-1] -= g1s+g2s;
      }
    }
    return g;
  }

  private static class A2 implements CgSolver.A {
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      VecArrayFloat2 v2z = v2x.clone();
      v2y.zero();
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = v2z.getArray();
      applyLhs(z,y);
    }
  }
  
  private static void applyLhs(float[][] f, float[][] g) {
    int nz = f[0].length;
    int nx = f.length;

    float r = 0.5f;
    float s = 0.5f;

    for (int i2=1; i2<nx; ++i2) {
      for (int i1=1; i1<nz; ++i1) {
        float f1r = f[i2  ][i1  ]-f[i2  ][i1-1];
        float f1s = f[i2-1][i1  ]-f[i2-1][i1-1];
        //float f2r = f[i2  ][i1  ]-f[i2-1][i1  ];
        //float f2s = f[i2  ][i1-1]-f[i2-1][i1-1];
        float f2r = 0.0f;
        float f2s = 0.0f;
        float g1 = r*f1r+s*f1s;
        float g2 = r*f2r+s*f2s; // gather ends
        float g1r = g1*r;
        float g1s = g1*s;
        float g2r = g2*r;
        float g2s = g2*s;
        g[i2  ][i1  ]  = g1r+g2r;
        g[i2  ][i1-1] -= g1r-g2s;
        g[i2-1][i1  ] += g1s-g2r;
        g[i2-1][i1-1] -= g1s+g2s;
      }
    }
  }

}
