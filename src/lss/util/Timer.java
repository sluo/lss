package lss.util;

import java.util.HashMap;
import java.util.logging.Logger;
import edu.mines.jtk.util.Stopwatch;

public class Timer {

  public Timer() {
    _timers = new HashMap<String,Stopwatch>();
  }

  public void start(String key) {
    Stopwatch stopwatch = new Stopwatch();
    _timers.put(key,stopwatch);
    log("begin: "+key);
    stopwatch.start();
  }

  public void stop(String key) {
    Stopwatch stopwatch = _timers.get(key);
    stopwatch.stop();
    int s = (int)stopwatch.time();
    int h = s/3600; s -= h*3600;
    int m = s/60; s -= m*60;
    String time = String.format(" %02d:%02d:%02d",h,m,s);
    log("  end: "+key+time);
    _timers.remove(key);
  }

  ///////////////////////////////////////////////////////////////////////////

  private HashMap<String,Stopwatch> _timers;
  private static Logger _log = Logger.getLogger(Timer.class.getName());

  private static void log(String s) {
    _log.finer(s);
  }

}
