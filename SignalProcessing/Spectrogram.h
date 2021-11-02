/*
  Spectrogram.h
  Class producing a 2 dimensional time vs frequency histogram where the power is the z axis
*/

#ifndef SPECTROGRAM_H
#define SPECTROGRAM_H

#include <vector>

#include "TH2.h"

#include "SignalProcessing/NoiseFunc.h"
#include "BasicFunctions/BasicFunctions.h"
#include "FieldClasses/FieldClasses.h"

namespace rad
{
  class Spectrogram
  {
  private:
    FieldPoint fieldPoint;
    std::vector<GaussianNoise> fNoise;
    
  public:
    Spectrogram(FieldPoint fp, std::vector<GaussianNoise> noiseTerms = {}) : fieldPoint(fp), fNoise(noiseTerms) {}
    ~Spectrogram();

    TH2D* MakeSpectrogram(const int NSamplesPerTimeBin, const bool kUseRetardedTime=false);
  };
}

#endif
