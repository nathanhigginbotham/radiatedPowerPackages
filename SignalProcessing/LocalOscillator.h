/*
  LocalOscillator.h
  Simple model of a local oscillator used for downmixing of signals
*/

#ifndef LOCAL_OSCILLATOR_H
#define LOCAL_OSCILLATOR_H

namespace rad
{
  class LocalOscillator
  {
  private:
    double angularFreq;

  public:
    LocalOscillator();
    LocalOscillator(double angFreq);
    ~LocalOscillator();
    double GetAngularFrequency() { return angularFreq; }
    void SetAngularFrequency(double newFreq) { angularFreq = newFreq; }

    
  };
}

#endif
