// sineTest.cxx

// STL
#include <iostream>
#include <cmath>

// ROOT includes
#include "TFile.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TMath.h"
#include "TRandom3.h"

// ROOT FFTW Wrapper
#include "FFTtools.h"

#include "FieldClasses/FieldClasses.h"
#include "BasicFunctions/BasicFunctions.h"

int main()
{
  double totalTime1 = 500.0;
  double stepSize1 = 0.01;
  int nSteps1 = totalTime1 / stepSize1;
  double totalTime2 = 300.0;  
  double stepSize2 = 0.005;
  int nSteps2 = totalTime2 / stepSize2;

  TRandom3* randomNum = new TRandom3();
  
  TGraph* gr1 = new TGraph();
  for (int n = 0; n < nSteps1; n++) {
    double time = n * stepSize1;
    double amp = TMath::Sin(2*TMath::Pi()*1.0*time) + 0.2*TMath::Sin(2*TMath::Pi()*5.0*time+1.0) + 0.3*TMath::Sin(2*TMath::Pi()*7.0*time+1.5)/* + randomNum->Gaus(0, 0.05)*/;
    gr1->SetPoint(gr1->GetN(), time, amp);
  }

  TGraph* gr2 = new TGraph();
  for (int n = 0; n < nSteps2; n++) {
    double time = n * stepSize2;
    double amp = TMath::Sin(2*TMath::Pi()*1.0*time) + 0.2*TMath::Sin(2*TMath::Pi()*5.0*time+1.0) + 0.3*TMath::Sin(2*TMath::Pi()*7.0*time+1.5)/* + randomNum->Gaus(0, 0.05)*/;
    gr2->SetPoint(gr2->GetN(), time, amp);
  }

  TGraph* gr1FFT = FFTtools::makePowerSpectrumPeriodogram(gr1);
  TGraph* gr2FFT = FFTtools::makePowerSpectrumPeriodogram(gr2);
  TGraph* gr1FFTVs = rad::MakePowerSpectrumNorm(gr1);
  TGraph* gr2FFTVs = rad::MakePowerSpectrumNorm(gr2);
  gr1FFT->SetLineWidth(2);
  gr2FFT->SetLineWidth(2);
  gr1FFTVs->SetLineWidth(2);
  gr2FFTVs->SetLineWidth(2);
  
  std::cout<<"Sum voltage squared (1, 2)= "<<FFTtools::sumVoltageSquared(gr1, -1, -1)/gr1->GetN()<<", "<<FFTtools::sumVoltageSquared(gr2, -1, -1)/gr2->GetN()<<std::endl;
  std::cout<<"Integral voltage squared (1, 2)= "<<FFTtools::integrateVoltageSquared(gr1, -1, -1)<<", "<<FFTtools::integrateVoltageSquared(gr2, -1, -1)<<std::endl;
  std::cout<<"FFT sum periodogram (1, 2) = "<<FFTtools::sumPower(gr1FFT)<<", "<<FFTtools::sumPower(gr2FFT)<<std::endl;
  std::cout<<"FFT integral (1, 2) = "<<FFTtools::integratePower(gr1FFTVs)/(gr1FFT->GetN()*stepSize1)<<", "<<FFTtools::integratePower(gr2FFTVs)/(gr2FFT->GetN()*stepSize2)<<std::endl;
  
  TFile* fout = new TFile("outputSine.root", "recreate");
  gr1->Write("gr1");
  gr2->Write("gr2");
  gr1FFT->Write("gr1FFT");
  gr2FFT->Write("gr2FFT");
  gr1FFTVs->Write("gr1FFTVs");
  gr2FFTVs->Write("gr2FFTVs");
  
  fout->Close();
  delete fout;
  
  return 0;
}
