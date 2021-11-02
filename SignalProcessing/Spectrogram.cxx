// Spectrogram.cxx

#include <cmath>

#include "TH2.h"

#include "SignalProcessing/Spectrogram.h"
#include "BasicFunctions/BasicFunctions.h"

rad::Spectrogram::~Spectrogram() {
  fNoise.clear();
}

TH2D* rad::Spectrogram::MakeSpectrogram(const double NSamplesPerTimeBin, const bool kUseRetardedTime) {
  const double sampleRate = fieldPoint.GetSampleRate();
  // First check the length of the field histograms
  // (this will be different if the retarded time is used)
  TGraph* grE = fieldPoint.GetEFieldTimeDomain(FieldPoint::Coord_t::kX, kUseRetardedTime);
  // Need to get the total number of bins given the sample rate
  const int nTimeBins = int(floor(grE->GetN() / NSamplesPerTimeBin)); 
  const double firstTime = grE->GetPointX(0);
  const double lastTime  = grE->GetPointX(nTimeBins * NSamplesPerTimeBin - 1);
  const int nFreqBins = (NSamplesPerTimeBin/2) + 1;
  const double deltaF = 1.0 / ((1.0/sampleRate) * NSamplesPerTimeBin);
  const double lastFreq = nFreqBins * deltaF;
  delete grE;
  
  TH2D* h2 = new TH2D("h2", "Spectrogram; Time [s]; Frequency [Hz]; Power [fW]", nTimeBins, firstTime, lastTime, nFreqBins, 0.0, lastFreq);
  
  return h2;
}
