/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lss.eni;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Amplitude and phase spectra.
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.10.29
 */
public class Spectrum {

  // Computes the amplitude and phase spectra for the specified sequence.
  public static float[][] computeAmplitudeAndPhase(Sampling st, float[] x) {

    // Time sampling.
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();

    // Frequency sampling.
    int nfft = FftReal.nfftSmall(2*nt);
    int nf = nfft/2+1;
    double df = 1.0/(nfft*dt);
    double ff = 0.0;
    Sampling fs = new Sampling(nf,df,ff);

    // Real-to-complex fast Fourier transform.
    FftReal fft = new FftReal(nfft);
    float[] cf = new float[2*nf];
    copy(nt,x,cf);
    fft.realToComplex(-1,cf,cf);

    // Adjust phase for possibly non-zero time of first sample.
    float[] wft = rampfloat(0.0f,-2.0f*FLT_PI*(float)(df*ft),nf);
    cf = cmul(cf,cmplx(cos(wft),sin(wft)));

    // Amplitude spectrum, normalized.
    float[] af = cabs(cf);
    float amax = max(max(af),FLT_EPSILON);
    af = mul(1.0f/amax,af);

    // Phase spectrum, in cycles.
    float[] pf = carg(cf);
    //pf = mul(0.5f/FLT_PI,pf);

    return new float[][]{af,pf};
  }
}
