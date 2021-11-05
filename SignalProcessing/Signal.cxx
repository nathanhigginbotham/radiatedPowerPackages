// Signal.cxx

#include "SignalProcessing/Signal.h"
#include "SignalProcessing/LocalOscillator.h"
#include "FieldClasses/FieldClasses.h"
#include "BasicFunctions/BasicFunctions.h"

#include <vector>
#include <iostream>

#include "TGraph.h"

#include "FFTtools.h"

rad::Signal::~Signal() {
  delete grVITime;
  delete grVQTime;
}

rad::Signal::Signal(FieldPoint fp, LocalOscillator lo, double srate,
		    std::vector<GaussianNoise> noiseTerms, const bool kUseRetardedTime) {
  sampleRate = srate;

  // Make sure the noise terms are all set up correctly
  for (int iNoise = 0; iNoise < noiseTerms.size(); iNoise++) {
    (noiseTerms.at(iNoise)).SetSampleFreq(sampleRate);
    (noiseTerms.at(iNoise)).SetSigma();
  }
  
  // Create the signal tgraph in the time domain
  TGraph* grVITimeUnfiltered = new TGraph();
  TGraph* grVQTimeUnfiltered = new TGraph();
  TGraph* grInputVoltage = fp.GetDipoleLoadVoltageTimeDomain(kUseRetardedTime, -1, -1);
  std::cout<<"Performing the downmixing..."<<std::endl;
  for (int i = 0; i < grInputVoltage->GetN(); i++) {
    grVITimeUnfiltered->SetPoint(i, grInputVoltage->GetPointX(i), grInputVoltage->GetPointY(i)*lo.GetInPhaseComponent(grInputVoltage->GetPointX(i)));
    grVQTimeUnfiltered->SetPoint(i, grInputVoltage->GetPointX(i), grInputVoltage->GetPointY(i)*lo.GetQuadratureComponent(grInputVoltage->GetPointX(i)));
  }

  // Now we need to filter and then sample these signals
  std::cout<<"Filtering.."<<std::endl;
  TGraph* grVITimeUnsampled = FFTtools::simplePassBandFilter(grVITimeUnfiltered, 0.0, sampleRate/2.0);
  TGraph* grVQTimeUnsampled = FFTtools::simplePassBandFilter(grVQTimeUnfiltered, 0.0, sampleRate/2.0);

  // Now do sampling
  // Use simple linear interpolation for the job
  grVITime = new TGraph();
  grVQTime = new TGraph();  
  setGraphAttr(grVITime);
  setGraphAttr(grVQTime);

  std::cout<<"Sampling..."<<std::endl;
  double sampleSpacing = 1.0 / sampleRate;
  double sampleTime = grVITimeUnsampled->GetPointX(0);
  for (int i = 0; i < grVITimeUnsampled->GetN(); i++) {
    double time = grVITimeUnsampled->GetPointX(i);
    if (time < sampleTime) continue;
    else if (i == 0) {
      double calcVI = grVITimeUnsampled->GetPointY(0);
      double calcVQ = grVQTimeUnsampled->GetPointY(0);
      for (int iNoise = 0; iNoise < noiseTerms.size(); iNoise++) {
	calcVI += (noiseTerms.at(iNoise)).GetNoiseVoltage();
	calcVQ += (noiseTerms.at(iNoise)).GetNoiseVoltage();
      }
      grVITime->SetPoint(grVITime->GetN(), sampleTime, calcVI);
      grVQTime->SetPoint(grVQTime->GetN(), sampleTime, calcVQ);
      sampleTime += sampleSpacing;
    }
    else {
      // Sample the distribution
      double calcVI = grVITimeUnsampled->GetPointY(i-1) + (sampleTime - grVITimeUnsampled->GetPointX(i-1)) * (grVITimeUnsampled->GetPointY(i) - grVITimeUnsampled->GetPointY(i-1)) / (time - grVITimeUnsampled->GetPointX(i-1));
      double calcVQ = grVQTimeUnsampled->GetPointY(i-1) + (sampleTime - grVQTimeUnsampled->GetPointX(i-1)) * (grVQTimeUnsampled->GetPointY(i) - grVQTimeUnsampled->GetPointY(i-1)) / (time - grVQTimeUnsampled->GetPointX(i-1));
      for (int iNoise = 0; iNoise < noiseTerms.size(); iNoise++) {
	calcVI += (noiseTerms.at(iNoise)).GetNoiseVoltage();
	calcVQ += (noiseTerms.at(iNoise)).GetNoiseVoltage();
      }
      grVITime->SetPoint(grVITime->GetN(), sampleTime, calcVI);
      grVQTime->SetPoint(grVQTime->GetN(), sampleTime, calcVQ);
      sampleTime += sampleSpacing;      
    }
  }
  
  delete grInputVoltage;
  delete grVITimeUnfiltered;
  delete grVQTimeUnfiltered;
  delete grVITimeUnsampled;
  delete grVQTimeUnsampled;
}

rad::Signal::Signal(const Signal &s1) {
  sampleRate = s1.sampleRate;
  grVITime = (TGraph*)s1.grVITime->Clone();
  grVQTime = (TGraph*)s1.grVQTime->Clone();
}
