// noiseLevels.cxx

#include <iostream>

#include "TFile.h"
#include "TGraph.h"

#include "FieldClasses/FieldClasses.h"
#include "SignalProcessing/Signal.h"
#include "SignalProcessing/NoiseFunc.h"
#include "Antennas/HertzianDipole.h"
#include "Antennas/HalfWaveDipole.h"

using namespace rad;

int main(int argc, char** argv)
{
  std::cout<<"Have "<<argc<<" arguments"<<std::endl;
  for (int i = 0; i < argc; ++i) {
    std::cout<<argv[i]<<std::endl;
  }
  TString outputFile = argv[1];
  
  const double antennaRadius = 0.02;
  const TVector3 antennaPos(antennaRadius, 0, 0);
  const TVector3 antennaZAxis(0, 1, 0);
  const TVector3 antennaXAxis(-1, 0, 0);

  HertzianDipole* hertzianAntenna = new HertzianDipole(antennaPos, antennaXAxis, antennaZAxis, 27.01e9);
  HalfWaveDipole* halfwaveAntenna = new HalfWaveDipole(antennaPos, antennaXAxis, antennaZAxis, 27.01e9);
  
  FieldPoint fpHertzian("/home/sjones/work/qtnm/trajectories/90DegOnAxis.root", hertzianAntenna);
  fpHertzian.GenerateFields(10e-6);
  std::cout<<"Generated field 1"<<std::endl;
  
  FieldPoint fpHalfwave("/home/sjones/work/qtnm/trajectories/90DegOnAxis.root", halfwaveAntenna);
  fpHalfwave.GenerateFields(10e-6);
  std::cout<<"Generated field 2"<<std::endl;

  TGraph* grVoltageHertzian = fpHertzian.GetAntennaLoadVoltageTimeDomain();
  grVoltageHertzian->GetXaxis()->SetRangeUser(0, 1.5e-10);
  TGraph* grVoltageHalfwave = fpHalfwave.GetAntennaLoadVoltageTimeDomain();
  grVoltageHalfwave->GetXaxis()->SetRangeUser(0, 1.5e-10);
  grVoltageHalfwave->SetLineColor(kRed);

  const double loadResistance = 70.0;
  const double noiseTemp50MHz = 0.26666666;
  LocalOscillator myLO(26.75e9 * 2 * TMath::Pi());
  GaussianNoise noise50MHz(noiseTemp50MHz, loadResistance);

  const double noiseTemp50kHz = 0.00026666666;
  GaussianNoise noise50kHz(noiseTemp50kHz, loadResistance);

  const double sampleRate = 750e6; // Hz
  
  Signal signal50MHz(fpHalfwave, myLO, sampleRate, {noise50MHz}, false);
  Signal signal50kHz(fpHalfwave, myLO, sampleRate, {noise50kHz}, false);

  TGraph* grVITime50MHz = signal50MHz.GetVITimeDomain();
  TGraph* grVITime50kHz = signal50kHz.GetVITimeDomain();
  TGraph* grVIPower50MHz = signal50MHz.GetVIPowerNorm(loadResistance);
  TGraph* grVIPower50kHz = signal50kHz.GetVIPowerNorm(loadResistance);
  
  TFile* fout = new TFile(outputFile, "RECREATE");
  fout->cd();
  grVoltageHertzian->Write("grVoltageHertzian");
  grVoltageHalfwave->Write("grVoltageHalfwave");

  grVITime50MHz->Write("grVITime50MHz");
  grVITime50kHz->Write("grVITime50kHz");
  grVIPower50MHz->Write("grVIPower50MHz");
  grVIPower50kHz->Write("grVIPower50kHz");  
  
  fout->Close();
  delete fout;
  
  return 0;
}
