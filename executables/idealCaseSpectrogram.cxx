/*
  idealCaseSpectrogram.cxx

  Make the spectrogram for the ideal case electron 
    - emitted in the middle of the trap
    - emitted with a pitch angle of 90 degrees
    - 4K thermal noise
*/

#include "SignalProcessing/Signal.h"
#include "Antennas/HalfWaveDipole.h"
#include "SignalProcessing/InducedVoltage.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/NoiseFunc.h"

#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"

using namespace rad;

int main(int argc, char* argv[])
{
  std::cout<<"Have "<<argc<<" arguments"<<std::endl;
  for (int i = 0; i < argc; ++i) {
    std::cout<<argv[i]<<std::endl;
  }
  TString outputFile = argv[1];

  const double antennaRadius = 0.02;
  const double antennaAngle1 = 0*TMath::Pi()/180;
  const double antennaAngle2 = 90*TMath::Pi()/180;

  TVector3 antennaPoint1(antennaRadius*TMath::Cos(antennaAngle1), antennaRadius*TMath::Sin(antennaAngle1), 0.0);
  TVector3 antennaDirZ1(-1*TMath::Sin(antennaAngle1), TMath::Cos(antennaAngle1), 0.0);
  TVector3 antennaDirX1(TMath::Cos(antennaAngle1), TMath::Sin(antennaAngle1), 0.0);
  HalfWaveDipole* antenna1 = new HalfWaveDipole(antennaPoint1, antennaDirX1, antennaDirZ1, 27.01e9);

  const double loadResistance = 70.0;
  const double noiseTemp = 0.0;
  const double optimalTimeBin = 54e-6;
  LocalOscillator myLO(26.75e9 * 2 * TMath::Pi());
  GaussianNoise noise1(noiseTemp, loadResistance);
  std::vector<GaussianNoise> noiseTerms;
  noiseTerms.push_back(noise1);
  const double sampleRate = 0.75e9; // Hz
  const int samplesPerBin = sampleRate*optimalTimeBin;
  
  InducedVoltage iv1("/home/sjones/work/qtnm/trajectories/electronTraj1ms90Deg.root", antenna1);
  
  Signal signal1(iv1, myLO, sampleRate, noiseTerms);
  
  TFile* fout = new TFile(outputFile, "RECREATE");
  fout->cd();

  TH2D* h2Spec = signal1.GetVISpectrogram(loadResistance, samplesPerBin);
  h2Spec->Write();
  
  fout->Close();
  delete fout;
  
  return 0;
}
