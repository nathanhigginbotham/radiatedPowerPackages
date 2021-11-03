// NoiseFunc.cxx

#include <cmath>
#include <time.h>

#include "TMath.h"
#include "TRandom3.h"

#include "SignalProcessing/NoiseFunc.h"

rad::GaussianNoise::GaussianNoise(double T, double R, int setSeed) {
  noiseTemp = T;
  resistance = R;
  numGen = new TRandom3(1234);//(setSeed == 1234) ? new TRandom3(setSeed) : new TRandom3(time(NULL));
}

rad::GaussianNoise::~GaussianNoise() {
  delete numGen;
}

void rad::GaussianNoise::SetSampleFreq(double fs) {
  sampleFreq = fs;
}

void rad::GaussianNoise::SetResistance(double r) {
  resistance = r;
}

void rad::GaussianNoise::SetSigma() {
  sigma = TMath::Sqrt( TMath::K() * noiseTemp * sampleFreq );
}

double rad::GaussianNoise::GetNoiseVoltage() {
  double premult = TMath::Sqrt(resistance*0.5);
  double volt = premult * numGen->Gaus(0.0, sigma);
  return volt;
}



