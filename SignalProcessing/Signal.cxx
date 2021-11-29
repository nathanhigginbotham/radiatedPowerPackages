// Signal.cxx

#include "SignalProcessing/Signal.h"
#include "SignalProcessing/InducedVoltage.h"
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
  grInputVoltage.clear();
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
    TGraph* grInputVoltageTemp = fp[point].GetAntennaLoadVoltageTimeDomain(kUseRetardedTime, -1, -1);
    grInputVoltage.push_back(grInputVoltageTemp);
    
    std::cout<<"Performing the downmixing..."<<std::endl;
    TGraph* grVITimeUnfiltered = DownmixInPhase(grInputVoltageTemp, lo);
    TGraph* grVQTimeUnfiltered = DownmixQuadrature(grInputVoltageTemp, lo);

    std::cout<<"First downsampling..."<<std::endl;
    TGraph* grVITimeFirstSample = SampleWaveform(grVITimeUnfiltered, 10*sampleRate);
    TGraph* grVQTimeFirstSample = SampleWaveform(grVQTimeUnfiltered, 10*sampleRate);

    delete grVITimeUnfiltered;
    delete grVQTimeUnfiltered;
    
    // Now we need to filter and then sample these signals
    std::cout<<"Filtering.."<<std::endl;
    TGraph* grVITimeUnsampled = BandPassFilter(grVITimeFirstSample, 0.0, sampleRate/2.0);
    TGraph* grVQTimeUnsampled = BandPassFilter(grVQTimeFirstSample, 0.0, sampleRate/2.0);

    delete grVITimeFirstSample;
    delete grVQTimeFirstSample;
    
    // Now do sampling
    // Use simple linear interpolation for the job
    std::cout<<"Sampling..."<<std::endl;
    TGraph* grVITimeTemp = SampleWaveform(grVITimeUnsampled);
    TGraph* grVQTimeTemp = SampleWaveform(grVQTimeUnsampled);

    delete grVITimeUnsampled;
    delete grVQTimeUnsampled;
    
    std::cout<<"Adding noise..."<<std::endl;
    AddGaussianNoise(grVITimeTemp, noiseTerms);
    AddGaussianNoise(grVQTimeTemp, noiseTerms);

    vecVITime.push_back(grVITimeTemp);
    vecVQTime.push_back(grVQTimeTemp);
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

rad::Signal::Signal(FieldPoint fp, LocalOscillator lo, double srate,
		    std::vector<GaussianNoise> noiseTerms, const bool kUseRetardedTime) {
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

  TGraph* grInputVoltageTemp = fp.GetAntennaLoadVoltageTimeDomain(kUseRetardedTime, -1, -1);
  grInputVoltage.push_back(grInputVoltageTemp);
    
  std::cout<<"Performing the downmixing..."<<std::endl;
  TGraph* grVITimeUnfiltered = DownmixInPhase(grInputVoltageTemp, lo);
  TGraph* grVQTimeUnfiltered = DownmixQuadrature(grInputVoltageTemp, lo);

  std::cout<<"First downsampling..."<<std::endl;
  TGraph* grVITimeFirstSample = SampleWaveform(grVITimeUnfiltered, 10*sampleRate);
  TGraph* grVQTimeFirstSample = SampleWaveform(grVQTimeUnfiltered, 10*sampleRate);

  delete grVITimeUnfiltered;
  delete grVQTimeUnfiltered;
  
  // Now we need to filter and then sample these signals
  std::cout<<"Filtering.."<<std::endl;
  TGraph* grVITimeUnsampled = BandPassFilter(grVITimeFirstSample, 0.0, sampleRate/2.0);
  TGraph* grVQTimeUnsampled = BandPassFilter(grVQTimeFirstSample, 0.0, sampleRate/2.0);

  delete grVITimeFirstSample;
  delete grVQTimeFirstSample;
  
  // Now do sampling
  // Use simple linear interpolation for the job
  std::cout<<"Sampling..."<<std::endl;
  grVITime = SampleWaveform(grVITimeUnsampled);
  grVQTime = SampleWaveform(grVQTimeUnsampled);

  delete grVITimeUnsampled;
  delete grVQTimeUnsampled; 
  
  std::cout<<"Adding noise..."<<std::endl;
  AddGaussianNoise(grVITime, noiseTerms);
  AddGaussianNoise(grVQTime, noiseTerms); 
}

rad::Signal::Signal(InducedVoltage iv, LocalOscillator lo, double srate,
		    std::vector<GaussianNoise> noiseTerms) {
  sampleRate = srate;

  // Make sure the noise terms are all set up correctly
  for (int iNoise = 0; iNoise < noiseTerms.size(); iNoise++) {
    (noiseTerms.at(iNoise)).SetSampleFreq(sampleRate);
    (noiseTerms.at(iNoise)).SetSigma();
  }

  TGraph* grInputVoltageTemp = (TGraph*)(iv.GetVoltageGraph())->Clone();
  grInputVoltage.push_back(grInputVoltageTemp);

  std::cout<<"Performing the downmixing..."<<std::endl;
  TGraph* grVITimeUnfiltered = DownmixInPhase(grInputVoltageTemp, lo);
  TGraph* grVQTimeUnfiltered = DownmixQuadrature(grInputVoltageTemp, lo);

  std::cout<<"First downsampling"<<std::endl;
  TGraph* grVITimeFirstSample = SampleWaveform(grVITimeUnfiltered, 10*sampleRate);
  TGraph* grVQTimeFirstSample = SampleWaveform(grVQTimeUnfiltered, 10*sampleRate);  

  delete grVITimeUnfiltered;
  delete grVQTimeUnfiltered;
  
  // Now we need to filter and then sample these signals
  std::cout<<"Filtering.."<<std::endl;
  TGraph* grVITimeUnsampled = BandPassFilter(grVITimeFirstSample, 0.0, sampleRate/2.0);
  TGraph* grVQTimeUnsampled = BandPassFilter(grVQTimeFirstSample, 0.0, sampleRate/2.0);

  delete grVITimeFirstSample;
  delete grVQTimeFirstSample;
  
  // Now do sampling
  // The actual output graphs
  grVITime = new TGraph();
  grVQTime = new TGraph();  
  setGraphAttr(grVITime);
  setGraphAttr(grVQTime);
  grVITime->GetXaxis()->SetTitle("Time [s]");
  grVQTime->GetXaxis()->SetTitle("Time [s]");
  
  // Use simple linear interpolation for the job
  std::cout<<"Sampling..."<<std::endl;
  grVITime = SampleWaveform(grVITimeUnsampled);
  grVQTime = SampleWaveform(grVQTimeUnsampled);

  delete grVITimeUnsampled;
  delete grVQTimeUnsampled; 
  
  std::cout<<"Adding noise..."<<std::endl;
  AddGaussianNoise(grVITime, noiseTerms);
  AddGaussianNoise(grVQTime, noiseTerms); 
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

TGraph* rad::Signal::SampleWaveform(TGraph* grInput) {
  TGraph* grOut = new TGraph();
  double sampleSpacing = 1.0 / sampleRate;
  double sampleTime = grInput->GetPointX(0);

  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    if (time < sampleTime) continue;
    else if (i == 0) {
      double calcV = grInput->GetPointY(0);
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;
    }
    else {
      // Sample the distribution using linear interpolation
      double calcV = grInput->GetPointY(i-1) + (sampleTime - grInput->GetPointX(i-1)) * (grInput->GetPointY(i) - grInput->GetPointY(i-1)) / (time - grInput->GetPointX(i-1));
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;      
    }
  }
  
  return grOut;
}

TGraph* rad::Signal::SampleWaveform(TGraph* grInput, const double sRate) {
  TGraph* grOut = new TGraph();
  double sampleSpacing = 1.0 / sRate;
  double sampleTime = grInput->GetPointX(0);
  
  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    if (time < sampleTime) continue;
    else if (i == 0) {
      double calcV = grInput->GetPointY(0);
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;
    }
    else {
      // Sample the distribution using linear interpolation
      double calcV = grInput->GetPointY(i-1) + (sampleTime - grInput->GetPointX(i-1)) * (grInput->GetPointY(i) - grInput->GetPointY(i-1)) / (time - grInput->GetPointX(i-1));
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;      
    }
  }

  return grOut;
}

void rad::Signal::AddGaussianNoise(TGraph* grInput, std::vector<GaussianNoise> noiseTerms) {
  for (int i = 0; i < grInput->GetN(); i++) {
    double voltage = grInput->GetPointY(i);
    for (int noise = 0; noise < noiseTerms.size(); noise++) {
      voltage += (noiseTerms.at(noise)).GetNoiseVoltage();  
    }
    grInput->SetPointY(i, voltage);
  }
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

TGraph* rad::Signal::GetVIPowerNorm(const double loadResistance) {
  TGraph* grOut = MakePowerSpectrumNorm(grVITime);
  for (int i = 0; i < grOut->GetN(); i++) {
    grOut->SetPointY(i, grOut->GetPointY(i) / loadResistance);
  }
  setGraphAttr(grOut);
  grOut->GetXaxis()->SetTitle("Frequency [Hz]");
  grOut->GetYaxis()->SetTitle("#frac{V_{I}^{2}}{R} #times (#Deltat)^{2} [W s^{2}]");
  return grOut;
}

TGraph* rad::Signal::GetVQPowerNorm(const double loadResistance) {
  TGraph* grOut = MakePowerSpectrumNorm(grVQTime);
  for (int i = 0; i < grOut->GetN(); i++) {
    grOut->SetPointY(i, grOut->GetPointY(i) / loadResistance);
  }
  setGraphAttr(grOut);
  grOut->GetXaxis()->SetTitle("Frequency [Hz]");
  grOut->GetYaxis()->SetTitle("#frac{V_{Q}^{2}}{R} #times (#Deltat)^{2} [W s^{2}]");
  return grOut;
}
