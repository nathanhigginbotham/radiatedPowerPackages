// noiseTests.cxx

#include "TFile.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TMath.h"

#include "FieldClasses/FieldClasses.h"
#include "SignalProcessing/NoiseFunc.h"
#include "SignalProcessing/Signal.h"
#include "SignalProcessing/LocalOscillator.h"
#include "BasicFunctions/BasicFunctions.h"

using namespace rad;

int main() {
  TVector3 antennaPoint(0.02, 0.0, 0.0);
  TVector3 dipoleDir(0.0, 1.0, 0.0);
  const double loadResistance = 70.0;
  const double noiseTemp = 0.001;
  LocalOscillator myLO(26.75e9 * 2 * TMath::Pi());
  GaussianNoise noise1(noiseTemp, loadResistance);
  std::vector<GaussianNoise> noiseTerms;
  noiseTerms.push_back(noise1);
  const double sampleRate = 0.75e9; // Hz
  
  FieldPoint fp(antennaPoint, dipoleDir, "/home/sjones/work/qtnm/trajectories/90DegOnAxis.root");
  fp.GenerateFields(4.6e-7);
  std::cout<<"Generated the fields for the field point"<<std::endl;

  TFile *fout = new TFile("noiseTestOutput.root", "RECREATE");
  TGraph* grInputVoltage = fp.GetDipoleLoadVoltageTimeDomain();
  Signal mySignal({fp}, myLO, sampleRate, noiseTerms);
  
  fout->cd();
  grInputVoltage->Write("grInputVoltage");
  TGraph* grVIUnfiltered = mySignal.GetVIUnfilteredTimeDomain(myLO);
  TGraph* grVQUnfiltered = mySignal.GetVQUnfilteredTimeDomain(myLO);
  TGraph* grVIUnsampled = mySignal.GetVIUnsampledTimeDomain(myLO);
  TGraph* grVQUnsampled = mySignal.GetVQUnsampledTimeDomain(myLO);
  TGraph* grVI = mySignal.GetVITimeDomain();
  TGraph* grVQ = mySignal.GetVQTimeDomain();

  TGraph* grVIUnfSpec = MakePowerSpectrumNorm(grVIUnfiltered);
  TGraph* grVQUnfSpec = MakePowerSpectrumNorm(grVQUnfiltered);
  TGraph* grVIUnsSpec = MakePowerSpectrumNorm(grVIUnsampled);
  TGraph* grVQUnsSpec = MakePowerSpectrumNorm(grVQUnsampled);
  TGraph* grVISpec = MakePowerSpectrumNorm(grVI);
  TGraph* grVQSpec = MakePowerSpectrumNorm(grVQ);
  
  grVIUnfiltered->Write("grVIUnf");
  grVQUnfiltered->Write("grVQUnf");
  grVIUnfSpec->Write("grVIUnfSpec");
  grVQUnfSpec->Write("grVQUnfSpec");
  grVIUnsampled->Write("grVIUns");
  grVQUnsampled->Write("grVQUns");
  grVIUnsSpec->Write("grVIUnsSpec");
  grVQUnsSpec->Write("grVQUnsSpec");
  grVI->Write("grVI");
  grVQ->Write("grVQ");
  grVISpec->Write("grVISpec");
  grVQSpec->Write("grVQSpec");
  
  fout->Close();
  delete fout;
  
  return 0;
}
