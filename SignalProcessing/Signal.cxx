// Signal.cxx

#include "SignalProcessing/Signal.h"
#include "SignalProcessing/LocalOscillator.h"
#include "FieldClasses/FieldClasses.h"
#include "BasicFunctions/BasicFunctions.h"

#include <vector>
#include <iostream>
#include <cassert>

#include "TGraph.h"

#include "FFTtools.h"

rad::Signal::~Signal() {
  for (int i = 0; i < grInputVoltage.size(); i++) {
    delete grInputVoltage[i];
  }
  delete grVITime;
  delete grVQTime;
}

rad::Signal::Signal(std::vector<FieldPoint> fp, LocalOscillator lo, double srate,
		    std::vector<GaussianNoise> noiseTerms, const bool kUseRetardedTime) {
  assert(fp.size() > 0);
  
  sampleRate = srate;
  
  // Make sure the noise terms are all set up correctly
  for (int iNoise = 0; iNoise < noiseTerms.size(); iNoise++) {
    (noiseTerms.at(iNoise)).SetSampleFreq(sampleRate);
    (noiseTerms.at(iNoise)).SetSigma();
  }

  // The actual output graphs
  grVITime = new TGraph();
  grVQTime = new TGraph();  
  setGraphAttr(grVITime);
  setGraphAttr(grVQTime);
  grVITime->GetXaxis()->SetTitle("Time [s]");
  grVQTime->GetXaxis()->SetTitle("Time [s]");
  
  std::vector<TGraph*> vecVITime;
  std::vector<TGraph*> vecVQTime;
  
  // For each antenna point generate the field and do the signal processing
  for (int point = 0; point < fp.size(); point++) {
    std::cout<<"Performing the downmixing..."<<std::endl;

    TGraph* grVITimeUnfiltered = DownmixInPhase(grInputVoltage[point], lo);
    TGraph* grVQTimeUnfiltered = DownmixQuadrature(grInputVoltage[point], lo);

    // Now we need to filter and then sample these signals
    std::cout<<"Filtering.."<<std::endl;
    TGraph* grVITimeUnsampled = BandPassFilter(grVITimeUnfiltered, 0.0, sampleRate/2.0);
    TGraph* grVQTimeUnsampled = BandPassFilter(grVQTimeUnfiltered, 0.0, sampleRate/2.0);

    // Now do sampling
    // Use simple linear interpolation for the job
    TGraph* grVITimeTemp = new TGraph();
    TGraph* grVQTimeTemp = new TGraph();  
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
	grVITimeTemp->SetPoint(grVITimeTemp->GetN(), sampleTime, calcVI);
	grVQTimeTemp->SetPoint(grVQTimeTemp->GetN(), sampleTime, calcVQ);
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
	grVITimeTemp->SetPoint(grVITimeTemp->GetN(), sampleTime, calcVI);
	grVQTimeTemp->SetPoint(grVQTimeTemp->GetN(), sampleTime, calcVQ);
	sampleTime += sampleSpacing;      
      }

      // Add the temporary graphs to the vector
      vecVITime.push_back(grVITimeTemp);
      vecVQTime.push_back(grVQTimeTemp);
    }
  
    delete grVITimeUnfiltered;
    delete grVQTimeUnfiltered;
    delete grVITimeUnsampled;
    delete grVQTimeUnsampled;
  } // Generate individual waveforms

  std::cout<<"Summing fields"<<std::endl;
  const int nSampledTimePoints = vecVITime[0]->GetN();
  for (int t = 0; t < nSampledTimePoints; t++) {
    // Loop over the signals from each individual antenna and sum
    double sumVI = 0;
    double sumVQ = 0;
    for (int field = 0; field < vecVITime.size(); field++) {
      sumVI += vecVITime[field]->GetPointY(t);
      sumVQ += vecVQTime[field]->GetPointY(t);     
    }
    grVITime->SetPoint(t, vecVITime[0]->GetPointX(t), sumVI);
    grVQTime->SetPoint(t, vecVQTime[0]->GetPointX(t), sumVQ);
  }
  
  // Remaining cleanup
  // for (auto p : vecVITime) {
  //    delete p;
  // }
  // vecVITime.clear();
  
  // for (auto p : vecVQTime) {
  //    delete p;
  // }
  // vecVQTime.clear();
  
  std::cout<<"Constructor completed"<<std::endl;
  
}

rad::Signal::Signal(const Signal &s1) {
  sampleRate = s1.sampleRate;
  grVITime = (TGraph*)s1.grVITime->Clone();
  grVQTime = (TGraph*)s1.grVQTime->Clone();
  for (int point = 0; point < (s1.grInputVoltage).size(); point++) {
    TGraph* grTemp = (TGraph*)s1.grInputVoltage[point]->Clone(); 
    grInputVoltage.push_back(grTemp);
  }
}

TGraph* rad::Signal::GetInputVoltage(const unsigned int field) {
  TGraph* gr = (TGraph*)grInputVoltage[field]->Clone();
  setGraphAttr(gr);
  return gr;
}

TGraph* rad::Signal::GetVITimeDomain() {
  TGraph* gr = (TGraph*)grVITime->Clone();
  setGraphAttr(gr);
  gr->GetXaxis()->SetTitle("Time [s]");
  gr->GetYaxis()->SetTitle("V_{I} [V]");
  return gr;
}

TGraph* rad::Signal::GetVQTimeDomain() {
  TGraph* gr = (TGraph*)grVQTime->Clone();
  setGraphAttr(gr);
  gr->GetXaxis()->SetTitle("Time [s]");
  gr->GetYaxis()->SetTitle("V_{Q} [V]");
  return gr;
}

TGraph* rad::Signal::DownmixInPhase(TGraph* grInput, LocalOscillator lo) {
  TGraph* grOut = new TGraph();
  for (int i = 0; i < grInput->GetN(); i++) {
    grOut->SetPoint(i, grInput->GetPointX(i), grInput->GetPointY(i)*lo.GetInPhaseComponent(grInput->GetPointX(i)));
  }
  return grOut;
}

TGraph* rad::Signal::DownmixQuadrature(TGraph* grInput, LocalOscillator lo) {
  TGraph* grOut = new TGraph();
  for (int i = 0; i < grInput->GetN(); i++) {
    grOut->SetPoint(i, grInput->GetPointX(i), grInput->GetPointY(i)*lo.GetQuadratureComponent(grInput->GetPointX(i)));
  }
  return grOut;
}

TGraph* rad::Signal::GetVIUnfilteredTimeDomain(LocalOscillator lo,
					       const unsigned int field) {
  return (DownmixInPhase(grInputVoltage[field], lo));
}

TGraph* rad::Signal::GetVQUnfilteredTimeDomain(LocalOscillator lo,
					       const unsigned int field) {
  return (DownmixQuadrature(grInputVoltage[field], lo));
}

TGraph* rad::Signal::GetVIUnsampledTimeDomain(LocalOscillator lo,
					      const unsigned int field) {
  TGraph* grVITimeUnfiltered = GetVIUnfilteredTimeDomain(lo, field);
  TGraph* grVITimeUnsampled = BandPassFilter(grVITimeUnfiltered, 0.0, sampleRate/2.0);
  delete grVITimeUnfiltered;
  return grVITimeUnsampled;
}

TGraph* rad::Signal::GetVQUnsampledTimeDomain(LocalOscillator lo,
					      const unsigned int field) {
  TGraph* grVQTimeUnfiltered = GetVQUnfilteredTimeDomain(lo, field);
  TGraph* grVQTimeUnsampled = BandPassFilter(grVQTimeUnfiltered, 0.0, sampleRate/2.0);
  delete grVQTimeUnfiltered;
  return grVQTimeUnsampled;
}
