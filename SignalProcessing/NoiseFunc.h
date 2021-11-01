/*
  NoiseFunc.h
  Contains functions for implementing noise in signal processing
*/
#ifndef NOISE_FUNC_H
#define NOISE_FUNC_H

namespace rad
{

  /// Gaussian random noise which is constant in time
  class GaussianNoise
  {
  private:
    double noiseTemp;
    double sampleFreq;
    double resistance;

    void SetSampleFreq(const double fs);
    void SetResistance(const double r);
    
  public:
    /// \param T is the noise temperature in kelvin
    GaussianNoise(const double T); 
    ~GaussianNoise();
  };
}

#endif
