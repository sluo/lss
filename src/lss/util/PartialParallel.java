package lss.util;

import static edu.mines.jtk.util.ArrayMath.*;

public class PartialParallel {

  public PartialParallel(int nparallel) {
    _np = nparallel;
  }

  public void loop(
  int end, edu.mines.jtk.util.Parallel.LoopInt body) {
    int nchunk = (end+_np-1)/_np;
    for (int i=0; i<nchunk; ++i) {
      int fi = i*_np;
      int li = min(fi+_np,end);
      edu.mines.jtk.util.Parallel.loop(fi,li,body);
    }
  }

  public <V> V reduce(
  int end, edu.mines.jtk.util.Parallel.ReduceInt<V> body) {
    int nchunk = (end+_np-1)/_np;
    V v = edu.mines.jtk.util.Parallel.reduce(0,min(_np,end),body);
    for (int i=1; i<nchunk; ++i) {
      int fi = i*_np;
      int li = min(fi+_np,end);
      v = body.combine(edu.mines.jtk.util.Parallel.reduce(fi,li,body),v);
    }
    return v;
  }
  
  int _np; // max number of indices to compute in parallel

}
