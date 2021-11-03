// spectrogramTest.cxx

#include "TFile.h"
#include "TGraph.h"
#include "TVector3.h"

#include "FieldClasses/FieldClasses.h"
#include "SignalProcessing/Spectrogram.h"
#include "SignalProcessing/NoiseFunc.h"

using namespace rad;

int main() {
  TVector3 antennaPoint(0.02, 0.0, 0.0);
  GaussianNoise* noise1 = new GaussianNoise(1.0, 70.0);
  
  FieldPoint *fp = new FieldPoint(antennaPoint, "/home/sjones/work/qtnm/trajectories/electronTraj600us_89deg.root");
  fp->GenerateFields(50e-6);
  std::cout<<"Generated the fields for the field point"<<std::endl;
  
  TFile* fout = new TFile("spectrogramOutput.root", "RECREATE");
  std::cout<<"Opened the output file"<<std::endl;

  std::vector<GaussianNoise*> noiseTerms;
  noiseTerms.push_back(noise1);
  Spectrogram spec(fp, noiseTerms);
  std::cout<<"Created spectrogram class"<<std::endl;
  TH2D* h2spec = spec.MakeSpectrogram(64768, false);
  std::cout<<"Created spectrogram"<<std::endl;

  fout->cd();
  h2spec->Write("h2spec");
  
  fout->Close();
  delete fout;
  
  return 0;
}
