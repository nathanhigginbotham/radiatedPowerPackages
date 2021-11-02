// Spectrogram.cxx

#include <cmath>

#include "TH2.h"

#include "SignalProcessing/Spectrogram.h"
#include "BasicFunctions/BasicFunctions.h"

rad::Spectrogram::~Spectrogram() {
  fNoise.clear();
}

TH2D* rad::Spectrogram::MakeSpectrogram(const int NSamplesPerTimeBin, const bool kUseRetardedTime) {
  const double sampleRate = fieldPoint.GetSampleRate();
  // First check the length of the field histograms
  // (this will be different if the retarded time is used)
  TGraph* grE = fieldPoint.GetEFieldTimeDomain(FieldPoint::Coord_t::kX, kUseRetardedTime);
  // Need to get the total number of bins given the sample rate
  const int nTimeBins = int(floor(double(grE->GetN()) / double(NSamplesPerTimeBin))); 
  const double firstTime = grE->GetPointX(0);
  const double lastTime  = grE->GetPointX(nTimeBins * NSamplesPerTimeBin - 1);
  const int nFreqBins = (NSamplesPerTimeBin/2) + 1;
  const double deltaF = 1.0 / ((1.0/sampleRate) * NSamplesPerTimeBin);
  const double lastFreq = nFreqBins * deltaF;
  delete grE;
  
  TH2D* h2 = new TH2D("h2", "Spectrogram; Time [s]; Frequency [Hz]; Power [W]", nTimeBins, firstTime, lastTime, nFreqBins, 0.0, lastFreq);
  // Loop through time bins and generate a power spectrum at each point
  for (int t = 1; t <= nTimeBins; t++) {
    TGraph* grPower = fieldPoint.GetDipolePowerSpectrumNorm(kUseRetardedTime, (t-1)*NSamplesPerTimeBin, t*NSamplesPerTimeBin - 1);
    // Now add data from the TGraph to the histogram
    for (int i = 0; i < grPower->GetN(); i++) {
      h2->SetBinContent(t, i+1, grPower->GetPointY(i));
    }
    delete grPower;
  }
  
  return h2;
}
