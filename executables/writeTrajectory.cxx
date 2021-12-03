/*
  writeTrajectory.cxx

  
*/

#include "ElectronDynamics/BorisSolver.h"
#include "ElectronDynamics/QTNMFields.h"
#include "BasicFunctions/Constants.h"

#include <unistd.h>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

using namespace rad;

int main(int argc, char *argv[])
{
  int opt;

  char* outputFile;
  double simTime;
  double simStepSize;
  char* inputFile;
  double pitchAngle = 90;
  bool energyLoss = false;
  bool hasInputFile = false;
  
  while((opt = getopt(argc, argv, ":o:t:s:i:p:e")) != -1) {
    switch(opt) {
    case 'o':
      outputFile = optarg;
      std::cout<<"Output file is "<<outputFile<<std::endl;
      break;
    case 't':
      simTime = atof(optarg);
      std::cout<<"Simulation time is "<<simTime<<std::endl;
      break;
    case 's':
      simStepSize = atof(optarg);
      std::cout<<"Simulation step size is "<<simStepSize<<std::endl;
      break;
    case 'i':
      inputFile = optarg;
      hasInputFile = true;
      std::cout<<"Input file is "<<inputFile<<std::endl;
      break;
    case 'p':
      pitchAngle = atof(optarg);
      std::cout<<"Pitch angle is "<<pitchAngle<<std::endl;
      break;
    case 'e':
      energyLoss = true;
      std::cout<<"Energy loss turned on"<<std::endl;
      break;
    case ':':
      std::cout<<"Option needs a value"<<std::endl;
      break;
    case '?':
      std::cout<<"Unknown option: "<<optopt<<std::endl;
      break;
    }
  }

  const double pitchAngleRad = pitchAngle * TMath::Pi() / 180;
  const double TElec = 18600; // eV
  const double gamma = TElec * TMath::Qe() / (ME * TMath::C()*TMath::C()) + 1;
  const double betaSq = 1 - 1 / pow(gamma, 2);
  const double V0 = sqrt(betaSq) * TMath::C();
  const TVector3 B0(0, 0, 1.0); // Background 1T field in Z direction
  double tau = 0.0;

  if (energyLoss) tau = 2 * R_E / (3 * TMath::C());

  TVector3 X0(0, -0.000460441, 0);
  TVector3 vInitial(V0 * sin(pitchAngleRad), 0, V0 * cos(pitchAngleRad));

  double simStartTime = 0;
  if (hasInputFile) {
    TFile* fin = new TFile(inputFile, "READ");
    TTree* intree = (TTree*)fin->Get("tree");
    double startX, startY, startZ;
    double startXVel, startYVel, startZVel;
    intree->SetBranchAddress("time", &simStartTime);
    intree->SetBranchAddress("xPos", &startX);
    intree->SetBranchAddress("yPos", &startY);
    intree->SetBranchAddress("zPos", &startZ);
    intree->SetBranchAddress("xVel", &startXVel);
    intree->SetBranchAddress("yVel", &startYVel);
    intree->SetBranchAddress("zVel", &startZVel);

    intree->GetEntry(0);
    double time0 = simStartTime;
    intree->GetEntry(1);
    double time1 = simStartTime;
    simStepSize = time1 - time0;

    intree->GetEntry(intree->GetEntries()-1);
    X0.SetX(startX);
    X0.SetY(startY);
    X0.SetZ(startZ);
    vInitial.SetX(startXVel);
    vInitial.SetY(startYVel);
    vInitial.SetZ(startZVel);
    
    fin->Close();
    delete fin;
  }

  // Set up the QTNM bathtub trap with coils at +/- 25cm
  double zc1 = -0.25;
  double zc2 = 0.25;
  double Rcoil = 0.03;
  double I = 2.0 * 0.0049 * Rcoil / MU0;
  BathtubField* bathtub = new BathtubField(Rcoil, I, zc1, zc2, B0);
  TVector3 maxVec(0, 0, zc1);
  std::cout<<"Max field perturbation = "<<(bathtub->evaluate_field_at_point(maxVec) - B0).Mag()<<std::endl;

  // Set up the Boris solver
  BorisSolver solver(bathtub, -TMath::Qe(), ME, tau);

  // Open the output ROOT file
  TFile* fout = new TFile(outputFile, "RECREATE");
  TTree* tree = new TTree("tree", "tree");

  double time;
  double xPos, yPos, zPos;
  double xVel, yVel, zVel;
  double xAcc, yAcc, zAcc;

  tree->Branch("time", &time, "time/D");
  tree->Branch("xPos", &xPos, "xPos/D");
  tree->Branch("yPos", &yPos, "yPos/D");
  tree->Branch("zPos", &zPos, "zPos/D");
  tree->Branch("xVel", &xVel, "xVel/D");
  tree->Branch("yVel", &yVel, "yVel/D");
  tree->Branch("zVel", &zVel, "zVel/D");
  tree->Branch("xAcc", &xAcc, "xAcc/D");
  tree->Branch("yAcc", &yAcc, "yAcc/D");
  tree->Branch("zAcc", &zAcc, "zAcc/D");

  int nTimeSteps = simTime / simStepSize;
  
  fout->Close();
  delete fout;
  
  return 0;
}