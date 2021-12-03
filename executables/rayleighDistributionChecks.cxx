// rayleighDistributionChecks.cxx

#include "BasicFunctions/BasicFunctions.h"

#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TH1.h"

using namespace rad;

double rayleighNorm(double *x, double *par)
{
  double xx = sqrt(x[0] * 18.5e3 * 750e6 / (0.5 * 20270));
  double retVal = RayleighPDF(xx, par[0]);
  return retVal;
}

int main(int argc, char** argv)
{
  std::cout<<"Have "<<argc<<" arguments"<<std::endl;
  for (int i = 0; i < argc; ++i) {
    std::cout<<argv[i]<<std::endl;
  }
  TString outputFile = argv[1];

  const double TEff = 4.0;
  const double tAcq = 54e-6;
  const double bandwidth = 1 / tAcq; // Hz
  std::cout<<"Bandwidth "<<bandwidth<<std::endl;
  const double sampleRate = 750e6; // Hz
  const double deltaT = 1 / sampleRate;
  const int N = sampleRate / bandwidth;
  int newLength = (N/2) + 1;
  double* x1 = new double[newLength];
  double* y1 = new double[newLength];
  double tempF = 0;
  
  for (int i = 0; i < newLength; i++) {
    x1[i] = tempF;
    tempF += bandwidth;
    y1[i] = 0;
  }
  
  TGraph* grNoise = new TGraph(newLength, x1, y1);
  AddWhiteNoiseFrequencyDomainPowerNorm(grNoise, TEff);
  setGraphAttr(grNoise);
  TH1D* hNoiseVals = new TH1D("hNoiseVals", "", 100, 0, 5e-9);
  SetHistAttr(hNoiseVals);
  for (int i = 0; i < grNoise->GetN(); i++) {
    double noiseRoot = sqrt(grNoise->GetPointY(i)*bandwidth/(deltaT*0.5*N));
    hNoiseVals->Fill(noiseRoot);
  }
  //hNoiseVals->Scale(1 / hNoiseVals->Integral(1, hNoiseVals->GetNbinsX()));
		    
  const double sigma = TMath::Sqrt(TMath::K() * TEff * bandwidth);
  TF1* fCDF = new TF1("fCDF", RayleighCDFFunc, 0, 4*sigma, 1);
  fCDF->SetParameter(0, sigma);
  fCDF->SetTitle("CDF: Rayleigh distribution with T = 4 K, B = 18.5 kHz; x; f(x; #sigma)");
  fCDF->GetHistogram()->GetXaxis()->SetTitle("#sqrt{k_{B}TB} [W^{0.5}]");

  TF1* fPDF = new TF1("fPDF", RayleighPDFFunc, 0, 4*sigma, 1);//rayleighNorm, 0, 25e-27, 1);
  fPDF->SetParameter(0, sigma);
  fPDF->SetTitle("PDF: Rayleigh distribution with T = 4 K, B = 18.5 kHz; x; f(x; #sigma)");
  fPDF->GetHistogram()->GetXaxis()->SetTitle("#sqrt{k_{B}TB} [W^{0.5}]");
  
  TFile* fout = new TFile(outputFile, "RECREATE");
  fout->cd();
  fCDF->Write();
  fPDF->Write();
  grNoise->Write("grNoise");
  hNoiseVals->Write();
  
  fout->Close();
  delete fout;
  
  return 0;
}
