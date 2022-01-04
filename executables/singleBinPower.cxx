// singleBinPower.cxx

#include "BasicFunctions/BasicFunctions.h"
#include "SignalProcessing/InducedVoltage.h"
#include "Antennas/HalfWaveDipole.h"
#include "SignalProcessing/Signal.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/NoiseFunc.h"

#include "TFile.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TMath.h"
#include "TString.h"

#include "FFTtools.h"

#include <vector>

using namespace rad;

int main(int argc, char** argv)
{
  std::cout<<"Have "<<argc<<" arguments"<<std::endl;
  for (int i = 0; i < argc; ++i) {
    std::cout<<argv[i]<<std::endl;
  }
  TString outputFile = argv[1];

  const double antennaLowerBandwidth = 26.95e9;
  const double antennaUpperBandwidth = 27.05e9;
  
  TVector3 antennaPoint1(0.02, 0, 0);
  TVector3 antennaDirZ1(0, 1, 0);
  TVector3 antennaDirX1(1, 0, 0);
  HalfWaveDipole* antenna1 = new HalfWaveDipole(antennaPoint1, antennaDirX1, antennaDirZ1, 27.01e9);
  antenna1->SetBandwidth(antennaLowerBandwidth, antennaUpperBandwidth);
  InducedVoltage fp1("/home/sjones/work/qtnm/trajectories/electronTraj60us90Deg.root", antenna1, 54e-6);

  LocalOscillator myLO(26.75e9 * 2 * TMath::Pi());
  const double loadResistance = 70.0;
  const double sampleRate = 750e6; // Hz
  const double noiseTemp = 4.0;
  GaussianNoise noise1(noiseTemp, loadResistance);
  std::vector<GaussianNoise> noiseTerms;
  noiseTerms.push_back(noise1);
  
  Signal signal1(fp1, myLO, sampleRate, noiseTerms);
  Signal signalNoNoise(fp1, myLO, sampleRate);
  
  TFile* fout = new TFile(outputFile, "RECREATE");
  fout->cd();

  TGraph* grVIPowerSpectrum = signal1.GetVIPowerNorm(loadResistance);
  TGraph* grVIPowerSpectrumNoNoise = signalNoNoise.GetVIPowerNorm(loadResistance);
  grVIPowerSpectrum->Write("grVIPowerSpectrum");
  grVIPowerSpectrumNoNoise->Write("grVIPowerSpectrumNoNoise");
  std::cout<<"Power integral "<<IntegratePowerNorm(grVIPowerSpectrum)<<std::endl;
  std::cout<<"Power integral no noise "<<IntegratePowerNorm(grVIPowerSpectrumNoNoise)<<std::endl;

  TGraph* grVIPeriodogram = signal1.GetVIPowerPeriodogram(loadResistance);
  TGraph* grVIPeriodogramNoNoise = signalNoNoise.GetVIPowerPeriodogram(loadResistance);
  std::cout<<"Power integral "<<FFTtools::sumPower(grVIPeriodogram)<<std::endl;
  std::cout<<"Power integral no noise"<<FFTtools::sumPower(grVIPeriodogramNoNoise)<<std::endl;
  grVIPeriodogram->Write("grVIPeriodogram");
  grVIPeriodogramNoNoise->Write("grVIPeriodogramNoNoise");
  
  fout->Close();
  delete fout;
  
  return 0;
}
