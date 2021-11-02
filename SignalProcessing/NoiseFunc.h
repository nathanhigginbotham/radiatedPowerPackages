/*
  NoiseFunc.h
  Contains functions for implementing noise in signal processing
*/
#ifndef NOISE_FUNC_H
#define NOISE_FUNC_H

#include "TRandom3.h"

namespace rad
{

  /// Gaussian random noise which is constant in time
  class GaussianNoise
  {
  private:
    TRandom3 *numGen; 
    double noiseTemp;
    double sampleFreq;
    double resistance;
    double sigma;
    
    void SetSampleFreq(const double fs);
    void SetResistance(const double r);

    void SetSigma();

    double GetNoiseVoltage();
    
  public:
    /// \param T is the noise temperature in kelvin
    /// \param setSeed is the seed for the random number generator
    GaussianNoise(const double T, const double R, const int setSeed=1234); 
    ~GaussianNoise();
  };
}

#endif
