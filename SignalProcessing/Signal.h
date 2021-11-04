/*
  Signal.h
  Class for a signal detected on an antenna or antennas
*/

#ifndef SIGNAL_H
#define SIGNAL_H

#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/NoiseFunc.h"
#include "FieldClasses/FieldClasses.h"
#include "BasicFunctions/BasicFunctions.h"

#include <vector>

#include "TGraph.h"

namespace rad
{
  /// Signal base class
  /// Class members are a voltage signal in the time domain and an out of phase component 

  
  class Signal
  {
  private:
    TGraph* grVITime; // In phase component
    TGraph* grVQTime; // Quadrature component
    double sampleRate;
    
  public:
    Signal();
    ~Signal();

    /// \param fp The field point where the antenna is located
    /// \param lo The local oscillator used to do the downmixing
    /// \param srate The sample rate of a hypothetical ADC
    /// \param noiseTerms A vector of noise terms to be added to the signal
    /// \param kUseRetardedTime Boolean on whether or not to use the retarded time
    Signal(FieldPoint fp, LocalOscillator lo, double srate, std::vector<GaussianNoise> noiseTerms={}, const bool kUseRetardedTime=false);
  };
}

#endif
