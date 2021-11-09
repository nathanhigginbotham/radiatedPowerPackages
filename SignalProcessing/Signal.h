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
  /// Class members are:
  /// grInputVoltage A vector of unprocessed voltage signals in the time domain
  /// grVITime In phase processed voltage component
  /// grVQTime Quadrature processed voltage component
  /// sampleRate Hypothetical sample rate of a digitiser
  
  class Signal
  {
  private:
    std::vector<TGraph*> grInputVoltage; 
    TGraph* grVITime; // In phase component
    TGraph* grVQTime; // Quadrature component
    double sampleRate;

    TGraph* DownmixInPhase(TGraph* grInput, LocalOscillator lo);
    TGraph* DownmixQuadrature(TGraph* grInput, LocalOscillator lo);
    
  public:
    ~Signal();
    
    /// \param fp The field points where the antenna is located
    /// \param lo The local oscillator used to do the downmixing
    /// \param srate The sample rate of a hypothetical ADC
    /// \param noiseTerms A vector of noise terms to be added to the signal
    /// \param kUseRetardedTime Boolean on whether or not to use the retarded time
    Signal(std::vector<FieldPoint> fp, LocalOscillator lo, double srate, std::vector<GaussianNoise> noiseTerms={}, const bool kUseRetardedTime=false);

    /// Copy constructor
    Signal(const Signal &s1);

    TGraph* GetInputVoltage(const unsigned int field=0);
    
    TGraph* GetVIUnfilteredTimeDomain(LocalOscillator lo, const unsigned int field=0);
    TGraph* GetVQUnfilteredTimeDomain(LocalOscillator lo, const unsigned int field=0);
    TGraph* GetVIUnsampledTimeDomain(LocalOscillator lo, const unsigned int field=0);
    TGraph* GetVQUnsampledTimeDomain(LocalOscillator lo, const unsigned int field=0);
    
    TGraph* GetVITimeDomain();
    TGraph* GetVQTimeDomain();
  };
}

#endif
