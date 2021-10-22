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
  fp.GenerateFields("/home/sjones/work/qtnm/trajectories/90DegOnAxis.root", 1e-9);  
  std::cout<<"All done"<<std::endl;

  TFile *fout = new TFile("outputFile.root", "RECREATE");
  fout->cd();
  TGraph *grEx = fp.GetEx();
  TGraph *grEy = fp.GetEy();
  TGraph *grEz = fp.GetEz();
  TGraph *grBx = fp.GetBx();
  TGraph *grBy = fp.GetBy();
  TGraph *grBz = fp.GetBz();
  grEx->Write("grEx");
  grEy->Write("grEy");
  grEz->Write("grEz");
  grBx->Write("grBx");
  grBy->Write("grBy");
  grBz->Write("grBz");

  fout->Close();
  delete fout;
  
  return 0;
}
