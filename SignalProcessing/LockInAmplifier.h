/*
  LockInAmplifier.h

  Class representing a lock-in amplifier
*/

#ifndef LOCK_IN_AMPLIFIER_H
#define LOCK_IN_AMPLIFIER_H

#include "SignalProcessing/Signal.h"
#include "SignalProcessing/LocalOscillator.h"
#include "BasicFunctions/BasicFunctions.h"

#include "TGraph.h"

namespace rad
{
  class LockInAmplifier
  {
  private:
    /// Class members are local oscillator at the reference frequency
    /// TO DO: This should be able to take more advanced and arbitrary signals
    /// although this is perhaps entering the realm of matched filtering
    LocalOscillator refFreq;
    
    /// Should also specify the low pass filter here.
    /// Would require me to write the appropriate class in the first place
    
  public:
    /// \param referenceFreq The processed signal at a given point
    LockInAmplifier(LocalOscillator referenceFreq);

    /// Produces the amplitude vs. time output of the amplifier
    /// \param sig The processed signal
    TGraph* GetAmplitudeTimeDomain(Signal sig);

    /// Produces the phase vs. time output of the amplifier
    /// \param sig The processed signal
    TGraph* GetPhaseTimeDomain(Signal sig);
  };
}

#endif
