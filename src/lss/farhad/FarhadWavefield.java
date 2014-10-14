package lss.farhad;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.awt.ColorMap;

public class FarhadWavefield {

  public static int compareDigits(double xa, double xb) {
    int significantDigits = 10;
    boolean matches = false;
    while (!matches && significantDigits>0) {
      Almost almost = new Almost(significantDigits);
      matches = almost.equal(xa, xb);
      if (!matches) {
        --significantDigits;
      }
    }
    return significantDigits;
  }

  public static double dot(double[][][] u, double[][][] a) {
    int nz = u[0][0].length;
    int nx = u[0].length;
    int nt = u.length;
    double sum = 0.0f;
    for (int it=0; it<nt; ++it)
      for (int ix=0; ix<nx; ++ix)
        for (int iz=0; iz<nz; ++iz)
          sum += u[it][ix][iz]*a[it][ix][iz];
    return sum;
  }

  //////////////////////////////////////////////////////////////////////////

  public FarhadWavefield(Sampling sz, Sampling sx, Sampling st) {
    _dz = sz.getDelta();
    _dx = sx.getDelta();
    _dt = st.getDelta();
    _dtdx2 = (_dt*_dt)/(_dx*_dx);
    _dtdz2 = (_dt*_dt)/(_dz*_dz);
    _nz = sz.getCount(); 
    _nx = sx.getCount(); 
    _nt = st.getCount();
    _adjoint = false;
  }

  public void modelAcousticWavefield(
    Source source, double[][] v, double[][][] u)
  {
    setSourceAndVelocity(source,v);
    modelAcoustic(u);
  }

  public void setAdjoint(boolean adjointOption) {
    _adjoint = adjointOption;  
    if (adjointOption) {
      System.out.println("ADJOINT PROPAGATOR");
    }
    else{
      System.out.println("FORWARD PROPAGATOR");
    }
  }
  
  private void modelAcoustic(double[][][] u) {
    int nz = _nz; 
    int nx = _nx;   
    int nt = _nt;
    
    double[][][] psi = zerodouble(_nz,_nx,2);
    
    //Stencil_3 stencil = new Stencil_3(); 
    //double[] sten = stencil.getStencil();
    
    //Stencil_15 stencil = new Stencil_15(); 
    //double[] sten = stencil.getStencil();
    
    //Stencil_Anu stencil = new Stencil_Anu(); 
    //double[] sten = stencil.getStencil();
    
    Stencil_21 stencil = new Stencil_21(); 
    double[] sten = stencil.getStencil();
    
    double [][] absorber = getAbsorber(factor,_b,_b,true);
    //plot(absorber,"A");

    for (int it=0; it<nt; ++it) {
      //System.out.format("\r%d",it);
      if (!_adjoint) {
        progressForward(u,psi,_source,sten,absorber,it);
      } else {
        progressAdjoint(u,psi,_source,sten,absorber,it);
      }
    }
    System.out.println();
    
    if (_source instanceof AdjointSource) {
      //System.out.println("\nin time reverse of adjoint source"); 
      reverse3(u);
    }
  }

  private void progressForward
    (double[][][] u, final double[][][] psi, Source source,
    final double[] sten, final double[][] A, int step)  {
    double dt = _dt;
    int nz = _nz;
    final int nx = _nx;
    final int prev = (step  )%2;    
    final int curr = (step+1)%2;    
    final int next = prev;
    final double dtdx2 = _dtdx2;
    final double dtdz2 = _dtdz2;
    final double[][] vp2 = mul(_v,_v); // velocity squared
    final int move = sten.length/2;

    
    Parallel.loop(move,nz-move,new Parallel.LoopInt() {
    public void compute (int iz) {  
      for (int ix=move;ix<nx-move;++ix) {
        double accum = 0;
        for (int sti=0;sti<sten.length;++sti) {
           double stTemp = sten[sti];
           int bump = sti-move;
           //Laplacian
           accum = accum + dtdx2 * stTemp * psi[curr][ix+bump][iz     ];
           accum = accum + dtdz2 * stTemp * psi[curr][ix     ][iz+bump];
        }
        accum = accum * vp2[ix][iz];
        //accum = accum + 2.0f * psi[curr][ix][iz]-A[ix][iz] * psi[prev][ix][iz]; // XXX
        //psi[next][ix][iz] = accum * A[ix][iz]; 
        accum = accum + 2.0f * psi[curr][ix][iz]-psi[prev][ix][iz];
        psi[next][ix][iz] = accum; 
      }
    }}); // end of parallel loop strucure
    source.add(step,psi[next]);
    copy(psi[curr],u[step]);
  }



  private void progressAdjoint
    (double[][][] u, final double[][][] psi, Source source,
    final double[] sten, final double[][] A, int step)  {
    double dt = _dt;
    int nz = _nz;
    final int nx = _nx;
    final int prev = (step  )%2;    
    final int curr = (step+1)%2;    
    final int next = prev;
    final double dtdx2 = _dtdx2;
    final double dtdz2 = _dtdz2;
    final double[][] vp2 = mul(_v,_v); 

    final int move = sten.length/2;

    Parallel.loop(move,nz-move,new Parallel.LoopInt() {
    public void compute (int iz) {  
      for (int ix=move;ix<nx-move;++ix) {
        double accum = 0;
        for (int sti=0;sti<sten.length;++sti) {
           double stTemp = sten[sti];
           int bump = sti-move;
           accum = accum+dtdx2*stTemp*vp2[ix+bump][iz     ]*psi[curr][ix+bump][iz     ];
           accum = accum+dtdz2*stTemp*vp2[ix     ][iz+bump]*psi[curr][ix     ][iz+bump];
        }
        //accum = accum + 2.0f * psi[curr][ix][iz]-A[ix][iz] * psi[prev][ix][iz]; // XXX
        //psi[next][ix][iz] = accum * A[ix][iz]; 
        accum = accum + 2.0f * psi[curr][ix][iz]-psi[prev][ix][iz];
        psi[next][ix][iz] = accum; 
      }
    }});
    source.add(step,psi[next]);
    copy(psi[curr],u[step]);
  }

  private interface Stencil {
    double[] getStencil();
  }
  

  public static class Stencil_Anu implements Stencil {
    public Stencil_Anu() {
      _stencil = new double[15]; 
      makeStencil(); 
    }
    public double[] getStencil() {
      return _stencil;
    }
   
    private void makeStencil() {
      _stencil[7 ] = (-3.16677852 * 1);
      _stencil[8 ] = ( 1.88007854 * 1);
      _stencil[9 ] = (-3.89038422 * 0.1);
      _stencil[10] = ( 1.24526825 * 0.1);
      _stencil[11] = (-4.28440709 * 0.01);
      _stencil[12] = ( 1.35799582 * 0.01);
      _stencil[13] = (-3.39721227 * 0.001);
      _stencil[14] = ( 4.83639697 * 0.0001);
      for (int i=0;i<7;++i) 
        _stencil[i] = _stencil[14-i];
    }
    private double[] _stencil; 
  }



  public static class Stencil_15 implements Stencil {
    public Stencil_15() {
      _stencil = new double[15]; 
      makeStencil(); 
    }
    public double[] getStencil() {
      return _stencil;
    }
    
    private void makeStencil() {
      _stencil[7 ] = (-0.31038384 * 10);
      _stencil[8 ] = ( 0.18222336 * 10);
      _stencil[9 ] = (-0.34456101 * 1 );
      _stencil[10] = ( 0.96150994 * 0.1);
      _stencil[11] = (-0.28189965 * 0.1);
      _stencil[12] = ( 0.78062201 * 0.01);
      _stencil[13] = (-0.19471198 * 0.01);
      _stencil[14] = ( 0.42654175 * 0.001);
      for (int i=0;i<7;++i) 
        _stencil[i] = _stencil[14-i];
    }
    private double[] _stencil; 
  }
 
  public static class Stencil_21 implements Stencil {
    public Stencil_21() {
      _stencil = new double[21]; 
      makeStencil(); 
    }
    public double[] getStencil() {
      return _stencil;
    }
    
    private void makeStencil() {
      _stencil[10] = (-0.32148051 * 10);   
      _stencil[11] = ( 0.19265816 * 10);
      _stencil[12] = (-0.43052632 * 1); 
      _stencil[13] = ( 0.15871000 * 1);
      _stencil[14] = (-0.68711400 * 0.1);
      _stencil[15] = ( 0.31406935 * 0.1);
      _stencil[16] = (-0.14454222 * 0.1);
      _stencil[17] = ( 0.65305182 * 0.01);
      _stencil[18] = (-0.28531535 * 0.01);
      _stencil[19] = ( 0.11937032 * 0.01);
      _stencil[20] = (-0.47508613 * 0.001);
      for (int i=0;i<10;++i) 
        _stencil[i] = _stencil[20-i];
    }
    private double[] _stencil; 
  }
  
  public static class Stencil_3 implements Stencil {
    public Stencil_3() {
      _stencil = new double[3]; 
      makeStencil(); 
    }
    public double[] getStencil() {
      return _stencil;
    }
    
    private void makeStencil() {
      _stencil[0] = ( 1.0);   
      _stencil[1] = (-2.0);
      _stencil[2] = ( 1.0); 
    }
    private double[] _stencil; 
  }


  //ABSORBER
  private double[][] getAbsorber (double factor, int thickx, int thickz,
                    boolean top_absorber) {
      int nz = _nz;
      int nx = _nx;
      double dt = _dt;
      double dz = _dz;
      double dx = _dx;
      double vmax = max(_v);
      //System.out.println("vmax = " + vmax); 
      double[][] A = filldouble(1.0,nz,nx);
      
      double ttrav = 2.0f * thickx*dx / (vmax*dt);
      //System.out.println("ttrav = " + ttrav); 
      double power = 0.0; 
      double tmp = 0.0;

      for (int iz=0; iz<nz; ++iz) {
        for (int ix=0; ix<thickx; ++ix) {
          tmp = cos((double)(ix) * 3.14159/thickx); 
          power = (tmp+1.0) / ttrav; 
          A[ix][iz] = A[ix][iz]*Math.pow(factor,power); 
          A[nx-ix-1][iz]=A[ix][iz]; 
        }
      }
  
      ttrav = 2.0f * thickz*dz / (vmax*dt);
      for (int ix=0; ix<nx; ++ix) {
        for (int iz=0; iz<thickz; ++iz) {
          tmp = cos((double)(iz) * 3.14159f/thickz); 
          power = (tmp+1.0) / ttrav;  
          A[ix][nz-iz-1] = A[ix][nz-iz-1]* Math.pow(factor,power); 
          if (top_absorber) 
            A[ix][iz] = A[ix][nz-iz-1];
        }
      }
    return A;
  }

  ////////////////////////////////////SOURCES/////////////////////////////
  public static interface Source {
    public void add(int it, double[][] f);
  }
  
  public static class RickerSource implements Source {
    public RickerSource(double fpeak, double dt, int kzs, int kxs) {
      _fpeak = fpeak;
      _kzs = kzs;
      _kxs = kxs;
      _tdelay = 1.0/fpeak;
      _dt = dt;
    }
    public void add(int it,double[][] f) {
      double t = it * _dt;
      f[_kxs][_kzs] += ricker(it*_dt-_tdelay);//*_dt*_dt;
    }
   

    private double ricker(double t) {
      double x = PI*_fpeak*t;
      double xx = x*x;
      return (1.0-2.0*xx)*exp(-xx);
    }
    
    private int _kzs,_kxs;
    private double _fpeak,_tdelay;
    private double _dt;
  }
  
  public static class DefaultSource implements Source {
    public DefaultSource(double[][] v, double dt, int kzs, int kxs, double[] s) {
      this(v,dt,new int[]{kzs},new int[]{kxs},new double[][]{s});
    }
    public DefaultSource(double[][] v, double dt, int[] kzs, int[] kxs, double[][] s) {
      Check.argument(kzs.length==kxs.length,"kzs.length=kxs.length");
       _ns = kzs.length;
       _dt = dt;
       _kzs = kzs;
       _kxs = kxs;
       _s = s;
       _v = v;
    }
    
    public void add(int it, double[][] f) {
      for (int is=0; is<_ns; is++) {
        int ix = _kxs[is];
        int iz = _kzs[is];
        f[ix][iz] += _s[is][it];//*_dt*_dt;// *_v[ix][iz]*_v[ix][iz];
      }
    }
    
    private double _dt;
    private int _ns;
    private double[][] _v;
    private int[] _kzs;
    private int[] _kxs; 
    private double[][] _s;
  }
  
  public static class AdjointSource implements Source {
    public AdjointSource(double[][] v, double dt, int kzs, int kxs, double[] s) {
      this(v, dt,new int[]{kzs},new int[]{kxs},new double[][]{s});
    }
    public AdjointSource(double[][] v, double dt, int[] kzs, int[] kxs, double[][] s) {
      _source = new DefaultSource(v,dt,kzs,kxs,reverse(s));
    }
    public void add(int it, double[][] f) {
      _source.add(it,f);
    }
    
    // reverse the time axis
    private static double[][] reverse(double[][] x) {
      int n1 = x[0].length;
      int n2 = x.length;
      double[][] y = new double[n2][n1];
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0,j1=n1-1; i1<n1; ++i1,--j1)
          y[i2][j1] = x[i2][i1];
      return y;
    }
    private DefaultSource _source;
  }

  public static class WavefieldAdjointSource implements Source {
    public WavefieldAdjointSource(double[][][] w) {
      _source = new WavefieldSource(reverse3(w));
    }
    public void add(int it, double[][] f) {
      _source.add(it,f);
    }
    private static double[][][] reverse3(double[][][] x) {
      int n1 = x[0][0].length;
      int n2 = x[0].length;
      int n3 = x.length;
      double[][][] y = new double[n3][n2][n1];
      for (int i3=0; i3<n3; ++i3)
        for (int i2=0; i2<n2; ++i2)
          for (int i1=0,j1=n1-1; i1<n1; ++i1,--j1)
            y[i3][i2][j1] = x[i3][i2][i1];
      return y;
    }
    private WavefieldSource _source;
  }

  public static class WavefieldSource implements Source {
    public WavefieldSource(double[][][] w) {
      _w = w;
    }
    public void add(final int it, final double[][] f) {
      final double[][] wf = _w[it]; 
      final int nz = wf[0].length;
      final int nx = wf.length;
      Parallel.loop(nx,new Parallel.LoopInt() {
      public void compute(int ix) {
        for (int iz=0; iz<nz; ++iz) 
          f[ix][iz] += wf[ix][iz];
      }});
    }
    
    private double _dt;
    private double[][] _r;
    private double[][][] _w;
    private double[][] _v;
  }
  
  public static class WavefieldSourceX implements Source {
    public WavefieldSourceX(double[][] v, double dt, double[][] r, double[][][] w) {
      int nz = w[0][0].length;
      int nx = w[0].length;
      int nt = w.length;
      _dt = dt;
      _r = r;
      _w = w;
      _v= v;
    }
    
    public void add(final int it, final double[][] f) {
      final double[][] wf = _w[it]; 
      final double[][] rf = _r;
      final int nz = wf[0].length;
      final int nx = wf.length;
      Parallel.loop(nx,new Parallel.LoopInt() {
      public void compute(int ix) {
        for (int iz=0; iz<nz; ++iz) 
          f[ix][iz] += wf[ix][iz]*rf[ix][iz];//+_dt*_dt;// *_v[ix][iz]*_v[ix][iz];
      }});
    }
    
    private double _dt;
    private double[][] _r;
    private double[][][] _w;
    private double[][] _v;
  }
  
  private void setSourceAndVelocity(
    Source source, double[][] v)
  {
    _source = source;
    _v = v;
  }
  
  private void reverse3(double[][][] f) {
    int n3 = f.length;
    double[][] t;
    for (int i3=0, j3=n3-1; i3<n3/2; ++i3, --j3) {
      t = f[i3];
      f[i3] = f[j3];
      f[j3] = t;
    }
  }

  private void plot(double[][] x, String t) {
    SimplePlot plot = SimplePlot.asPixels(x);
    PixelsView pv =plot.addPixels(x);
    pv.setColorModel(ColorMap.JET);
    plot.addColorBar();
    plot.setTitle(t);
  }
  
  private static int _b = 40;// absorbing boundary size
  //private static double factor = 0.0001; 
  private static double factor = 0.001; 
  private Source _source;
  private Stencil stencil;
  private double _t; // current time
  private int _it; // current time index
  private int _nx,_nz,_nt; // number of samples
  private double _dx,_dz,_dt; // sampling intervals
  private double _dtdx2, _dtdz2;
  private double[][] _v; // velocity 
  private boolean _adjoint;
} 
