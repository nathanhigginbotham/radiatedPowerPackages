// dipoleRadiationPattern.cxx
// Basic script for generating components of the electric field
 
// STL
#include <iostream>
#include <cmath>

// ROOT includes
#include "TFile.h"
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
  ROOT::Math::XYZPoint antennaPoint(0.02, 0.0, 0.0);
  ROOT::Math::XYZVector dipoleDir(0.0, 1.0, 0.0);
  // Declare the field point
  FieldPoint fp(antennaPoint, dipoleDir, "/home/sjones/work/qtnm/trajectories/90DegOnAxis.root");
  fp.GenerateFields(4.5e-7);  
  std::cout<<"Generated fields"<<std::endl;
  
  TFile *fout = new TFile("dipoleOutputFile.root", "RECREATE");
  fout->cd();
  // TGraph *grEx = fp.GetEFieldTimeDomain(FieldPoint::Coord_t::kX);
  // TGraph *grEy = fp.GetEFieldTimeDomain(FieldPoint::Coord_t::kY);
  // TGraph *grEz = fp.GetEFieldTimeDomain(FieldPoint::Coord_t::kZ);
  // TGraph *grBx = fp.GetBFieldTimeDomain(FieldPoint::Coord_t::kX);
  // TGraph *grBy = fp.GetBFieldTimeDomain(FieldPoint::Coord_t::kY);
  // TGraph *grBz = fp.GetBFieldTimeDomain(FieldPoint::Coord_t::kZ);
  // grEx->Write("grEx");
  // grEy->Write("grEy");
  // grEz->Write("grEz");
  // grBx->Write("grBx");
  // grBy->Write("grBy");
  // grBz->Write("grBz");

  // TGraph *grExPower = fp.GetEFieldPowerSpectrumNorm(FieldPoint::Coord_t::kX, false);
  // TGraph *grEyPower = fp.GetEFieldPowerSpectrumNorm(FieldPoint::Coord_t::kY, false);
  // TGraph *grEzPower = fp.GetEFieldPowerSpectrumNorm(FieldPoint::Coord_t::kY, false);
  // TGraph *grTotalEFieldPower = fp.GetTotalEFieldPowerSpectrumNorm(false);
  // grExPower->Write("grExPower");
  // grEyPower->Write("grEyPower");
  // grEzPower->Write("grEzPower");
  // grTotalEFieldPower->Write("grTotalEFieldPower");

  // TGraph *grSx = fp.GetPoyntingVecTimeDomain(FieldPoint::Coord_t::kX, false);
  // TGraph *grSy = fp.GetPoyntingVecTimeDomain(FieldPoint::Coord_t::kY, false);
  // TGraph *grSz = fp.GetPoyntingVecTimeDomain(FieldPoint::Coord_t::kZ, false);
  // grSx->Write("grSx");
  // grSy->Write("grSy");
  // grSz->Write("grSz");
  // TGraph* grSxFFT = MakePowerSpectrumNorm(grSx);
  // TGraph* grSyFFT = MakePowerSpectrumNorm(grSy);
  // TGraph* grSzFFT = MakePowerSpectrumNorm(grSz);
  // grSxFFT->Write("grSxFFT");
  // grSyFFT->Write("grSyFFT");
  // grSzFFT->Write("grSzFFT");
  
  // double totalPowerIntegral    = IntegratePowerNorm(grTotalEFieldPower);
  // std::cout<<"Total power integral, total power integral ret = "<<totalPowerIntegral<<" V^2 m^-2"<<std::endl;

  // Now get the component voltages
  
  TGraph* grDipoleVoltagePower    = fp.GetDipoleTotalVoltagePowerSpectrumNorm(false);
  TGraph* grDipoleVoltagePowerRet = fp.GetDipoleTotalVoltagePowerSpectrumNorm(true);
  grDipoleVoltagePower->Write("grDipoleVoltagePower");
  grDipoleVoltagePowerRet->Write("grDipoleVoltagePowerRet");
  
  TGraph *grSMag = fp.GetPoyntingMagTimeDomain();
  grSMag->Write("grSMag");
  TGraph *grSMagRet = fp.GetPoyntingMagTimeDomain(true);
  grSMagRet->Write("grSMagRet");
  TGraph *grDipolePowerTimeDomain = fp.GetDipolePowerTimeDomain();
  grDipolePowerTimeDomain->Write("grDipolePowerTime");

  TGraph* grDipoleLoadVoltageTime = fp.GetDipoleLoadVoltageTimeDomain();
  grDipoleLoadVoltageTime->Write("grDipoleLoadVoltageTime");
  TGraph* grDipoleLoadPowerTime = fp.GetDipoleLoadPowerTimeDomain(70.0);
  grDipoleLoadPowerTime->Write("grDipoleLoadPowerTime");
  TGraph* grDipoleLoadPeriodogram = FFTtools::makePowerSpectrumPeriodogram(grDipoleLoadVoltageTime);
  for (int n = 0; n < grDipoleLoadPeriodogram->GetN(); n++) {
    grDipoleLoadPeriodogram->SetPointY(n, grDipoleLoadPeriodogram->GetPointY(n) / 70.0);
  }
  grDipoleLoadPeriodogram->Write("grDipoleLoadPeriodogram");
  std::cout<<"Mean squared amplitude = "<<FFTtools::sumPower(grDipoleLoadPeriodogram)<<std::endl;
  
  TGraph* grDipoleLoadPower = fp.GetDipoleLoadPowerSpectrumNorm(70.0);
  grDipoleLoadPower->Write("grDipoleLoadPower");
  std::cout<<"Load power integral = "<<IntegratePowerNorm(grDipoleLoadPower)*1e15<<" fW"<<std::endl;
  
  TGraph* grDipolePowerSpectrumNorm = fp.GetDipolePowerSpectrumNorm(false);
  grDipolePowerSpectrumNorm->Write("grDipolePowerSpectrumNorm");
  std::cout<<"Power integral with proper normalisation = "<<IntegratePowerNorm(grDipolePowerSpectrumNorm)*1e15<<" fW"<<std::endl;

  
  
  fout->Close();
  delete fout;
  
  return 0;
}
