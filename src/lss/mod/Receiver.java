package lss.mod;

import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.dsp.SincInterpolator;
import edu.mines.jtk.dsp.Sampling;

public class Receiver {

  /**
   * Residual.
   */
  public static interface Residual {

    /**
     * Computes the residual.
     * @param recp input receiver containing predicted data.
     * @param reco input receiver containing observed data.
     * @return receiver containing residuals.
     */
    public Receiver compute(Receiver rcp, Receiver rco);
  }

  ////////////////////////////////////////////////////////////////////////////
  // Constructors

  public Receiver(int xr, int zr, int nt) {
    this(new int[]{xr},new int[]{zr},nt,null);
  }

  public Receiver(int[] xr, int[] zr, int nt) {
    this(xr,zr,nt,null);
  }

  public Receiver(float[] xr, float[] zr, int nt) {
    this(xr,zr,nt,null);
  }
  public Receiver(float[] xr, float[] zr, float[][] data) {
    this(xr,zr,-1,data);
  }
  public Receiver(int[] xr, int[] zr, float[][] data) {
    this(xr,zr,-1,data);
  }

  public Receiver(Receiver rc) {
    this(copy(rc.getXIndices()),copy(rc.getZIndices()),copy(rc.getData()));
  }

  private Receiver(int[] xr, int[] zr, int nt, float[][] data) {
    Check.argument(xr.length==zr.length,"xr.length==zr.length");
    _nr = xr.length;
    _xr = xr;
    _zr = zr;
    for (int ir=0; ir<_nr; ++ir) {
      Check.argument(_xr[ir]>=0,"xr[ir]>=0");
      Check.argument(_zr[ir]>=0,"zr[ir]>=0");
    }
    if (data!=null) {
      Check.argument(data.length==_nr,"data.length==nr");
      _nt = data[0].length;
      _data = data;
    } else {
      _nt = nt;
      _data = new float[_nr][_nt];
    }
    _round_coord = false;
  }

  private Receiver(float[] xr, float[] zr, int nt, float[][] data) {
    Check.argument(xr.length==zr.length,"xr.length==zr.length");
    _nr = xr.length;
    _xr = new int[_nr];
    _zr = new int[_nr];
    
    for (int i=0; i<_nr; ++i){
      _xr[i] = (int) (xr[i]+0.5f);
      _zr[i] = (int) (zr[i]+0.5f);
    }
    _xrf = xr;  
    _zrf = zr;

    for (int ir=0; ir<_nr; ++ir) {
      Check.argument(_xr[ir]>=0,"xr[ir]>=0");
      Check.argument(_zr[ir]>=0,"zr[ir]>=0");
    }
    
    _round_coord = true;
    buildTable();

    if (data!=null) {
      Check.argument(data.length==_nr,"data.length==nr");
      _nt = data[0].length;
      _data = data;
    } else {
      _nt = nt;
      _data = new float[_nr][_nt];
    }
  }


  ////////////////////////////////////////////////////////////////////////////

  
  public void setData(float[][] ui, int it, int nabsorb) {
    if(_round_coord){
      setDataR(ui, it, nabsorb);
    }else{
      setDataI(ui, it, nabsorb);
    }
  }

  public void setDataR(float[][] ui, int it, int nabsorb) {
    int n = (int) _n/2;
 
    for (int ir=0; ir<_nr; ++ir) {
      int xr = _xr[ir]+nabsorb;
      int zr = _zr[ir]+nabsorb;
      /* Gather wavefield around the receiver locations */
      for(int i2=0; i2<_n; i2++){
        int i2u = -n +i2;        
        for(int i1=0; i1<_n; i1++){
          int i1u = -n +i1;
          _data[ir][it] += _sinc_weights[ir][i2][i1]*ui[zr+i2u][xr+i1u]; 
        }
      }
    }
  }

  public boolean isSincRec(){
    return _round_coord;
  }

  public void setDataI(float[][] ui, int it, int nabsorb) {
    for (int ir=0; ir<_nr; ++ir) {
      int xr = _xr[ir]+nabsorb;
      int zr = _zr[ir]+nabsorb;
      _data[ir][it] = ui[zr][xr];
    }
  }

  public void setData(float[][] d) {
    Check.argument(d[0].length==_nt,"d[0].length==_nt");
    Check.argument(d.length==_nr,"d.length==_nt");
    copy(d,_data);
  }

  public int[][] getIndices() {
    return new int[][]{_xr,_zr};
  }

  public int[] getXIndices() {
    return _xr;
  }

  public int[] getZIndices() {
    return _zr;
  }

  public float[][] getData() {
    return _data;
  }

  public int getNt() {
    return _nt;
  }

  public int getNr() {
    return _nr;
  }

  public Receiver clone() {
    return new Receiver(copy(_xr),copy(_zr),copy(_data));
  }

  private void buildTable(){
    _n = 9;
    float o =(int) (-_n/2)+1;
    SincInterpolator sinc = new SincInterpolator();

    /* Here I build the input spike */
    Sampling is1 = new Sampling(_n,1,o);
    Sampling is2 = is1;
    float[][] ispike = new float[_n][_n];
    ispike[(int)_n/2][(int)_n/2] = 1.0f;
    _sinc_weights = new float[_nr][_n][_n];

    for (int i3=0; i3<_nr; ++i3){
      /* Here I build the shifted spike */
      float d1 = (float) _xrf[i3] - _xr[i3];
      float d2 = (float) _zrf[i3] - _zr[i3];
  
      Sampling os1 = new Sampling(_n,1,o+d1);
      Sampling os2 = new Sampling(_n,1,o+d2);
  
      for (int i2=0; i2< os2.getCount(); i2++){
        float z = (float) os2.getValue(i2);
        for (int i1=0; i1< os1.getCount(); i1++){
          float x = (float) os1.getValue(i1);
          _sinc_weights[i3][i2][i1] = sinc.interpolate(is1,is2,ispike,x,z);
        }
      }
    }
  }



  private int _nr,_nt,_n;
  private int[] _xr, _zr;
  private float[][] _data;
  private float[][][] _sinc_weights; 
  private float[] _xrf, _zrf;
  private boolean _round_coord;
}
