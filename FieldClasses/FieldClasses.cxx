// FieldClasses.cxx

#include <cassert>

#include "FieldClasses.h"

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TVector3.h"

rad::FieldPoint::FieldPoint() {
  Ex = new TGraph();
  Ey = new TGraph();
  Ez = new TGraph();
  Bx = new TGraph();
  By = new TGraph();
  Bz = new TGraph();
  antennaPoint = TVector3(0.0, 0.0, 0.0);
}

rad::FieldPoint::~FieldPoint() {
  delete Ex;
  delete Ey;
  delete Ez;
  delete Bx;
  delete By;
  delete Bz;
}

// Parametrised constructor
rad::FieldPoint::FieldPoint(TVector3 inputAntenna) {
  antennaPoint = inputAntenna;
}

void rad::FieldPoint::ResetFields() {
  Ex->Clear();
  Ey->Clear();
  Ez->Clear();
  Bx->Clear();
  By->Clear();
  Bz->Clear();
}

// From an input TFile generate the E and B fields for a given time
// inputFile is the input file which should have the relevant branches
// maxTime is the final time in seconds (if less than the time in the file)
void rad::FieldPoint::GenerateFields(const char* inputFile, const double maxTime) {
  ResetFields();
  
  TFile *fin = new TFile(inputFile, "READ");
  assert(fin);
  TTree *tree = (TTree*)fin->Get("tree");

  // Set variables
  double time;
  double xPos, yPos, zPos;
  double xVel, yVel, zVel;
  double xAcc, yAcc, zAcc;
  tree->SetBranchAddress("time", &time);
  tree->SetBranchAddress("xPos", &xPos);
  tree->SetBranchAddress("yPos", &yPos);
  tree->SetBranchAddress("zPos", &zPos);
  tree->SetBranchAddress("xVel", &xVel);
  tree->SetBranchAddress("yVel", &yVel);
  tree->SetBranchAddress("zVel", &zVel);
  tree->SetBranchAddress("xAcc", &xAcc);
  tree->SetBranchAddress("yAcc", &yAcc);
  tree->SetBranchAddress("zAcc", &zAcc);
  
  // Loop through the entries and get the fields at each point
  for (int e = 0; e < tree->GetEntries(); e++) {
    tree->GetEntry(e);
    if (time > maxTime) break;
    TVector3 ePos(xPos, yPos, zPos);
    TVector3 eVel(xVel, yVel, zVel);
    TVector3 eAcc(xAcc, yAcc, zAcc);
    TVector3 EField = CalcEField(antennaPoint, ePos, eVel, eAcc);
    TVector3 BField = CalcEField(antennaPoint, ePos, eVel, eAcc);

    Ex->SetPoint(Ex->GetN(), time, EField.X());
    Ey->SetPoint(Ey->GetN(), time, EField.Y());
    Ez->SetPoint(Ez->GetN(), time, EField.Z());
    Bx->SetPoint(Bx->GetN(), time, BField.X());
    By->SetPoint(By->GetN(), time, BField.Y());
    Bz->SetPoint(Bz->GetN(), time, BField.Z());
  }

  fin->Close();
  delete fin;
}
