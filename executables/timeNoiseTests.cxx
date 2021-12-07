// timeNoiseTests.cxx

#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TAxis.h"
#include "TGraph.h"

#include "SignalProcessing/Signal.h"
#include "SignalProcessing/NoiseFunc.h"
#include "SignalProcessing/LocalOscillator.h"
#include "BasicFunctions/BasicFunctions.h"
#include "FieldClasses/FieldClasses.h"
#include "Antennas/HertzianDipole.h"

using namespace rad;

int main()
{
  TVector3 antennaPoint1(0.02, 0.0, 0.0);
  TVector3 dipoleDirZ1(0.0, 1.0, 0.0);
  TVector3 dipoleDirX1(1.0, 0.0, 0.0);
  HertzianDipole* antenna1 = new HertzianDipole(antennaPoint1, dipoleDirX1, dipoleDirZ1, 27.01e9); 

  const double loadResistance = 70.0;
  const double noiseTemp = 0.01;
  LocalOscillator myLO(26.75e9 * 2 * TMath::Pi());
  GaussianNoise noise1(noiseTemp, loadResistance);
  std::vector<GaussianNoise> noiseTerms;
  noiseTerms.push_back(noise1);
  const double sampleRate = 0.75e9; // Hz

  FieldPoint fpShort("/home/sjones/work/qtnm/trajectories/electronTraj_88deg_800us.root", antenna1);
  fpShort.GenerateFields(0, 5e-7);
  Signal signalShort({fpShort}, myLO, sampleRate, noiseTerms, false);
  
  TGraph* grInputShort = signalShort.GetInputVoltage();
  TGraph* grVIUnfShort = signalShort.GetVIUnfilteredTimeDomain(myLO);
  TGraph* grVIUnsShort = signalShort.GetVIUnsampledTimeDomain(myLO);
  TGraph* grVIShort = signalShort.GetVITimeDomain();
  TGraph* grVQShort = signalShort.GetVQTimeDomain();
  TGraph* grVIShortSpec = MakePowerSpectrumNorm(grVIShort);
  TGraph* grVQShortSpec = MakePowerSpectrumNorm(grVQShort);
  
  FieldPoint fpLong("/home/sjones/work/qtnm/trajectories/electronTraj_88deg_800us.root", antenna1);
  fpLong.GenerateFields(0, 20e-6);
  Signal signalLong({fpLong}, myLO, sampleRate, noiseTerms, false);

  TGraph* grInputLong = signalLong.GetInputVoltage();
  TGraph* grVILong = signalLong.GetVITimeDomain();
  TGraph* grVQLong = signalLong.GetVQTimeDomain();
  TGraph* grVILongSpec = MakePowerSpectrumNorm(grVILong);
  TGraph* grVQLongSpec = MakePowerSpectrumNorm(grVQLong);
  
  TFile* fout = new TFile("timeNoiseOutput.root", "RECREATE");
  fout->cd();

  grInputShort->Write("grInputShort");
  grVIUnfShort->Write("grVIUnfShort");
  grVIUnsShort->Write("grVIUnsShort");
  grVIShort->Write("grVIShort");
  grVQShort->Write("grVQShort");
  grVIShortSpec->Write("grVIShortSpec");
  grVQShortSpec->Write("grVQShortSpec");

  grInputLong->Write("grInputLong");
  grVILong->Write("grVILong");
  grVQLong->Write("grVQLong");
  grVILongSpec->Write("grVILongSpec");
  grVQLongSpec->Write("grVQLongSpec");
  
  fout->Close();
  delete fout;
  
  return 0;
}
