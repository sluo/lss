package lss.util;

import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

public class PartialParallel {

  private int _np; // max number of indices to compute in parallel

  public PartialParallel(int nparallel) {
    _np = nparallel;
  }

  ////////////////////////////////////////////////////////////////////////////
  // loop

  public void loop(int end, Parallel.LoopInt body) {
    /*
    int nchunk = (end+_np-1)/_np;
    for (int i=0; i<nchunk; ++i) {
      int fi = i*_np;
      int li = min(fi+_np,end);
      Parallel.loop(fi,li,body);
    }
    */
    loop(0,end,1,body);
  }

  public void loop(int begin, int end, Parallel.LoopInt body) {
    loop(begin,end,1,body);
  }

  public void loop(int begin, int end, int step, Parallel.LoopInt body) {
    int nchunk = ((end-begin+step-1)/step+_np-1)/_np;
    int stride = _np*step;
    for (int i=0; i<nchunk; ++i) {
      int fi = begin+i*stride;
      int li = min(fi+stride,end);
      Parallel.loop(fi,li,step,body);
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // reduce


  public <V> V reduce(int end, Parallel.ReduceInt<V> body) {
    /*
    int nchunk = (end+_np-1)/_np;
    V v = Parallel.reduce(0,min(_np,end),body);
    for (int i=1; i<nchunk; ++i) {
      int fi = i*_np;
      int li = min(fi+_np,end);
      v = body.combine(Parallel.reduce(fi,li,body),v);
    }
    return v;
    */
    return reduce(0,end,1,body);
  }

  public <V> V reduce(int begin, int end, Parallel.ReduceInt<V> body) {
    return reduce(begin,end,1,body);
  }

  public <V> V reduce(
    int begin, int end, int step, Parallel.ReduceInt<V> body)
  {
    int nchunk = ((end-begin+step-1)/step+_np-1)/_np;
    int stride = _np*step;
    V v = Parallel.reduce(begin,min(begin+stride,end),step,body);
    for (int i=1; i<nchunk; ++i) {
      int fi = begin+i*stride;
      int li = min(fi+stride,end);
      v = body.combine(Parallel.reduce(fi,li,step,body),v);
    }
    return v;
  }
  

}
