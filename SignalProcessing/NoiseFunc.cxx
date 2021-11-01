// NoiseFunc.cxx

#include <cmath>
#include <time.h>

#include "TMath.h"
#include "TRandom3.h"

#include "SignalProcessing/NoiseFunc.h"

rad::GaussianNoise::GaussianNoise(const double T, const int setSeed) {
  noiseTemp = T;
  numGen = (setSeed == 1234) ? new TRandom3(setSeed) : new TRandom3(time(NULL));
}

rad::GaussianNoise::~GaussianNoise() {
  delete numGen;
}

void rad::GaussianNoise::SetSampleFreq(const double fs) {
  sampleFreq = fs;
}

void rad::GaussianNoise::SetResistance(const double r) {
  resistance = r;
}

void rad::GaussianNoise::SetSigma() {
  sigma = TMath::Sqrt( TMath::K() * noiseTemp * sampleFreq );
}

double rad::GaussianNoise::GetNoiseVoltage() {
  return ( sqrt(resistance*0.5) * numGen->Gaus(0, sigma) );
}



