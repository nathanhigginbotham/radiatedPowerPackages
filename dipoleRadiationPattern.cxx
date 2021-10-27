// dipoleRadiationPattern.cxx
// Basic script for generating components of the electric field
 
// STL
#include <iostream>
#include <cmath>

// ROOT includes
#include "TFile.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TMath.h"

// ROOT FFTW Wrapper
#include "FFTtools.h"

// My files
#include "FieldClasses/FieldClasses.h"
#include "BasicFunctions/BasicFunctions.h"

using namespace rad;

int main()
{
  TVector3 antennaPoint(0.02, 0.0, 0.0);
  // Declare the field point
  FieldPoint fp(antennaPoint);
  fp.GenerateFields("/home/sjones/work/qtnm/trajectories/90DegOnAxis.root", 1e-7);  

  TFile *fout = new TFile("outputFile.root", "RECREATE");
  fout->cd();
  TGraph *grEx = fp.GetEFieldTimeDomain(FieldPoint::Coord_t::kX);
  TGraph *grEy = fp.GetEFieldTimeDomain(FieldPoint::Coord_t::kY);
  TGraph *grEz = fp.GetEFieldTimeDomain(FieldPoint::Coord_t::kZ);
  TGraph *grBx = fp.GetBFieldTimeDomain(FieldPoint::Coord_t::kX);
  TGraph *grBy = fp.GetBFieldTimeDomain(FieldPoint::Coord_t::kY);
  TGraph *grBz = fp.GetBFieldTimeDomain(FieldPoint::Coord_t::kZ);
  grEx->Write("grEx");
  grEy->Write("grEy");
  grEz->Write("grEz");
  grBx->Write("grBx");
  grBy->Write("grBy");
  grBz->Write("grBz");
  TGraph *grExRet = fp.GetEFieldTimeDomain(FieldPoint::Coord_t::kX, true);
  TGraph *grEyRet = fp.GetEFieldTimeDomain(FieldPoint::Coord_t::kY, true);
  TGraph *grEzRet = fp.GetEFieldTimeDomain(FieldPoint::Coord_t::kZ, true);
  TGraph *grBxRet = fp.GetBFieldTimeDomain(FieldPoint::Coord_t::kX, true);
  TGraph *grByRet = fp.GetBFieldTimeDomain(FieldPoint::Coord_t::kY, true);
  TGraph *grBzRet = fp.GetBFieldTimeDomain(FieldPoint::Coord_t::kZ, true);
  grExRet->Write("grExRet");
  grEyRet->Write("grEyRet");
  grEzRet->Write("grEzRet");
  grBxRet->Write("grBxRet");
  grByRet->Write("grByRet");
  grBzRet->Write("grBzRet");

  TGraph *grExPower = FFTtools::makePowerSpectrumPeriodogram(grEx);//fp.GetEFieldPeriodogram(FieldPoint::Coord_t::kX);
  TGraph *grEyPower = FFTtools::makePowerSpectrumPeriodogram(grEy);//fp.GetEFieldPeriodogram(FieldPoint::Coord_t::kY);
  TGraph *grEzPower = fp.GetEFieldPeriodogram(FieldPoint::Coord_t::kZ);  
  TGraph *grTotalEFieldPower = fp.GetTotalEFieldPeriodogram();
  grExPower->Write("grExPower");
  grEyPower->Write("grEyPower");
  grEzPower->Write("grEzPower");
  grTotalEFieldPower->Write("grTotalEFieldPower");
  TGraph *grExPowerRet = FFTtools::makePowerSpectrumPeriodogram(grExRet);//fp.GetEFieldPeriodogram(FieldPoint::Coord_t::kX, true);
  TGraph *grEyPowerRet = FFTtools::makePowerSpectrumPeriodogram(grEyRet);//fp.GetEFieldPeriodogram(FieldPoint::Coord_t::kY, true);
  TGraph *grEzPowerRet = fp.GetEFieldPeriodogram(FieldPoint::Coord_t::kZ, true);  
  TGraph *grTotalEFieldPowerRet = fp.GetTotalEFieldPeriodogram(true);
  grExPowerRet->Write("grExPowerRet");
  grEyPowerRet->Write("grEyPowerRet");
  grEzPowerRet->Write("grEzPowerRet");
  grTotalEFieldPowerRet->Write("grTotalEFieldPowerRet");

  std::cout<<std::setprecision(10);
  double sumPower1 = FFTtools::sumVoltageSquared(grEy, -1, -1)/grEy->GetN();
  double sumPower2 = FFTtools::sumVoltageSquared(grEyRet, -1, -1)/grEyRet->GetN();
  std::cout<<"Sum vsquared, sum vsquared (ret) = "<<sumPower1<<", "<<sumPower2<<std::endl;
  std::cout<<"Sum FFT power, sum FFT power (ret) = "<<FFTtools::sumPower(grEyPower)<<", "<<FFTtools::sumPower(grEyPowerRet)<<std::endl;
  
  TGraph *grEMag = fp.GetEFieldMagTimeDomain();
  TGraph *grBMag = fp.GetBFieldMagTimeDomain();
  TGraph *grSMag = fp.GetPoyntingMagTimeDomain();
  grEMag->Write("grEMag");
  grBMag->Write("grBMag");
  grSMag->Write("grSMag");
  TGraph *grEMagRet = fp.GetEFieldMagTimeDomain(true);
  TGraph *grBMagRet = fp.GetBFieldMagTimeDomain(true);
  TGraph *grSMagRet = fp.GetPoyntingMagTimeDomain(true);
  grEMagRet->Write("grEMagRet");
  grBMagRet->Write("grBMagRet");
  grSMagRet->Write("grSMagRet");

  // Transform E field magnitude
  TGraph *grEMagPower = FFTtools::makePowerSpectrumPeriodogram(grEMag);
  setGraphAttr(grEMagPower);
  grEMagPower->GetYaxis()->SetTitle("|E|^{2}");
  grEMagPower->GetXaxis()->SetTitle("Frequency [Hz]");
  grEMagPower->Write("grEMagPower");
  TGraph *grEMagPowerRet = FFTtools::makePowerSpectrumPeriodogram(grEMagRet);
  setGraphAttr(grEMagPowerRet);
  grEMagPowerRet->GetYaxis()->SetTitle("|E|^{2}");
  grEMagPowerRet->GetXaxis()->SetTitle("Frequency [Hz]");
  grEMagPowerRet->Write("grEMagPowerRet");

  std::cout<<"Power magnitude, power sum = "<<FFTtools::sumPower(grEMagPower, -1, -1)<<", "<<FFTtools::sumPower(grTotalEFieldPower, -1, -1)<<std::endl;
  
  TGraph *grDipolePower = fp.GetDipolePowerTimeDomain(true);
  grDipolePower->Write("grDipolePower");
  int retLength = grDipolePower->GetN()/2 + 1;
  double *dipolePowerX = grDipolePower->GetX();
  double *dipolePowerY = grDipolePower->GetY();
  double deltaT = dipolePowerX[1] - dipolePowerX[0];
  double deltaF = 1/(deltaT*grDipolePower->GetN());
  FFTWComplex *dipolePowerFFT = FFTtools::doFFT(grDipolePower->GetN(), dipolePowerY);  
  double *newX = new double[retLength];
  double *newY = new double[retLength];
  double tempF = 0;
  for (int i = 0; i < retLength; i++) {
    float fft = FFTtools::getAbs(dipolePowerFFT[i]);
    if (i > 0 && i < retLength - 1) fft *= 2;
    fft /= grDipolePower->GetN();
    newX[i] = tempF;
    newY[i] = fft;
    tempF += deltaF;
  }
  TGraph* grDipolePowerFFT = new TGraph(retLength, newX, newY);
  delete [] dipolePowerFFT;
  delete [] newX;
  delete [] newY;
  setGraphAttr(grDipolePowerFFT);
  grDipolePowerFFT->GetXaxis()->SetTitle("Frequency [Hz]");
  grDipolePowerFFT->GetYaxis()->SetTitle("Power [W]");
  grDipolePowerFFT->Write("grDipolePowerFFT");
  
  fout->Close();
  delete fout;
  
  return 0;
}
