// signalSum.cxx

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

int main()
{
  TVector3 antennaPoint1(0.05, 0.0, 0.0);
  TVector3 dipoleDir1(0.0, 1.0, 0.0);
  
  const double antennaAngle2 = 10.0 * TMath::Pi()/180;
  TVector3 antennaPoint2(0.05*TMath::Cos(antennaAngle2), 0.05*TMath::Sin(antennaAngle2), 0.0);
  TVector3 dipoleDir2(TMath::Sin(antennaAngle2), TMath::Cos(antennaAngle2), 0.0);
  
  const double antennaAngle3 = -10.0 * TMath::Pi()/180;
  TVector3 dipoleDir3(TMath::Sin(antennaAngle3), TMath::Cos(antennaAngle3), 0.0); 
  TVector3 antennaPoint3(0.05*TMath::Cos(antennaAngle3), 0.05*TMath::Sin(antennaAngle3), 0.0);
  
  const double loadResistance = 70.0;
  const double noiseTemp = 0.0005;
  LocalOscillator myLO(26.75e9 * 2 * TMath::Pi());
  GaussianNoise noise1(noiseTemp, loadResistance);
  std::vector<GaussianNoise> noiseTerms;
  noiseTerms.push_back(noise1);
  const double sampleRate = 0.75e9; // Hz

  FieldPoint fp1(antennaPoint1, dipoleDir1, "/home/sjones/work/qtnm/trajectories/90DegOnAxis.root");
  fp1.GenerateFields(4.6e-7);

  FieldPoint fp2(antennaPoint2, dipoleDir2, "/home/sjones/work/qtnm/trajectories/90DegOnAxis.root");
  fp2.GenerateFields(4.6e-7);
  
  FieldPoint fp3(antennaPoint3, dipoleDir3, "/home/sjones/work/qtnm/trajectories/90DegOnAxis.root");
  fp3.GenerateFields(4.6e-7);
  
  std::cout<<"Generated the fields for the field point"<<std::endl;
  
  TFile* fout = new TFile("signalSumOutput.root", "recreate");

  Signal mySignal({fp1}, myLO, sampleRate, noiseTerms, false);
  fout->cd();
  TGraph* grVI = mySignal.GetVITimeDomain();
  TGraph* grVQ = mySignal.GetVQTimeDomain();
  TGraph* grVISpec = MakePowerSpectrumNorm(grVI);
  TGraph* grVQSpec = MakePowerSpectrumNorm(grVQ);

  Signal mySignalBoth({fp1, fp2}, myLO, sampleRate, noiseTerms, false);
  fout->cd();
  TGraph* grVIBoth = mySignalBoth.GetVITimeDomain();
  TGraph* grVQBoth = mySignalBoth.GetVQTimeDomain();
  TGraph* grVIBothSpec = MakePowerSpectrumNorm(grVIBoth);
  TGraph* grVQBothSpec = MakePowerSpectrumNorm(grVQBoth);

  Signal mySignalAll({fp1, fp2, fp3}, myLO, sampleRate, noiseTerms, false);
  fout->cd();
  TGraph* grVIAll = mySignalAll.GetVITimeDomain();
  TGraph* grVQAll = mySignalAll.GetVQTimeDomain();
  TGraph* grVIAllSpec = MakePowerSpectrumNorm(grVIAll);
  TGraph* grVQAllSpec = MakePowerSpectrumNorm(grVQAll);
  
  grVI->Write("grVI");
  grVQ->Write("grVQ");
  grVISpec->Write("grVISpec");
  grVQSpec->Write("grVQSpec");

  grVIBoth->Write("grVIBoth");
  grVQBoth->Write("grVQBoth");
  grVIBothSpec->Write("grVIBothSpec");
  grVQBothSpec->Write("grVQBothSpec");

  grVIAll->Write("grVIAll");
  grVQAll->Write("grVQAll");
  grVIAllSpec->Write("grVIAllSpec");
  grVQAllSpec->Write("grVQAllSpec");
  
  fout->Close();
  delete fout;
  return 0;
}
