// NoiseFunc.cxx

#include <cmath>

#include "SignalProcessing/NoiseFunc.h"

void rad::GaussianNoise::SetSampleFreq(const double fs) {
  sampleFreq = fs;
}

void rad::GaussianNoise::SetResistance(const double r) {
  resistance = r;
}

rad::GaussianNoise::GaussianNoise(const double T) {
  noiseTemp = T;
}




