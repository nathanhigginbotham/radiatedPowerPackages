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
  std::cout<<"All done"<<std::endl;

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

  TGraph *grExPower = fp.GetEFieldPeriodogram(FieldPoint::Coord_t::kX);
  TGraph *grEyPower = fp.GetEFieldPeriodogram(FieldPoint::Coord_t::kY);
  TGraph *grEzPower = fp.GetEFieldPeriodogram(FieldPoint::Coord_t::kZ);  
  TGraph *grTotalEFieldPower = fp.GetTotalEFieldPeriodogram();
  grExPower->Write("grExPower");
  grEyPower->Write("grEyPower");
  grEzPower->Write("grEzPower");
  grTotalEFieldPower->Write("grTotalEFieldPower");
  
  fout->Close();
  delete fout;
  
  return 0;
}
