package lss.vel;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import lss.vel.AcousticWaveOperator;
import lss.vel.BornWaveOperator;

// testing
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;


public class BornWaveOperatorS {

//  public BornWaveOperatorS(
//  Source[] source, Receiver[] receiver,
//  float[][] s, double dx, double dt, int nabsorb) {
//    _ns = source.length;
//    _born = new BornWaveOperator[_ns];
//
//  }
//  private float[][][][] _u;
//  private float[][][][] _a;
//
//  private void initialize() {
//
//  }
//
//  private static void Loop(int nsou, int psou, ComputeMethod method) {
//  }
//
//  private static interface ComputeMethod {
//    public void compute(int isou, int fsou, int lsou);
//  }
//
//  private static class ChunkLoopInt {
//    public ChunkLoopInt(ComputeMethod method) {
//      _method = method;
//    }
//    public void loop(
//    private ComputeMethod _method;
//
//  }
//  def ChunkLoopInt(computeMethod,nsou=None):
//    stopwatch = Stopwatch()
//    if nsou is None:
//      nsou = ns
//    nchunk = (nsou+psou-1)/psou
//    for ichunk in range(nchunk):
//      stopwatch.restart()
//      fsou = ichunk*psou # first source
//      lsou = min(fsou+psou,nsou) # last source
//      class Loop(Parallel.LoopInt):
//        def compute(self,isou):
//          computeMethod(isou,fsou,lsou)
//      Parallel.loop(fsou,lsou,Loop())
//      report('lsou=%d'%lsou,stopwatch)
//
//  BornWaveOperator[] _born;
//  int ns;
}
