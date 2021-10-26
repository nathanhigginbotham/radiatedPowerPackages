// FieldClasses.cxx

#include <cassert>
#include <cmath>

#include "BasicFunctions/Constants.h"
#include "FieldClasses.h"

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TVector3.h"

#include "FFTtools.h"

rad::FieldPoint::FieldPoint() {
  EField[0] = new TGraph();
  EField[1] = new TGraph();
  EField[2] = new TGraph();
  BField[0] = new TGraph();
  BField[1] = new TGraph();
  BField[2] = new TGraph();
  antennaPoint = TVector3(0.0, 0.0, 0.0);
  pos[0] = new TGraph();
  pos[1] = new TGraph();
  pos[2] = new TGraph();
  vel[0] = new TGraph();
  vel[1] = new TGraph();
  vel[2] = new TGraph();
  acc[0] = new TGraph();
  acc[1] = new TGraph();
  acc[2] = new TGraph();
  tPrime = new TGraph();
}

rad::FieldPoint::~FieldPoint() {
  delete EField[0];
  delete EField[1];
  delete EField[2];
  delete BField[0];
  delete BField[1];
  delete BField[2];

  delete pos[0];
  delete pos[1];
  delete pos[2];
  delete vel[0];
  delete vel[1];
  delete vel[2];
  delete acc[0];
  delete acc[1];
  delete acc[2];

  delete tPrime;
}

// Parametrised constructor
rad::FieldPoint::FieldPoint(TVector3 inputAntenna) {
  EField[0] = new TGraph();
  EField[1] = new TGraph();
  EField[2] = new TGraph();
  BField[0] = new TGraph();
  BField[1] = new TGraph();
  BField[2] = new TGraph();
  pos[0] = new TGraph();
  pos[1] = new TGraph();
  pos[2] = new TGraph();
  vel[0] = new TGraph();
  vel[1] = new TGraph();
  vel[2] = new TGraph();
  acc[0] = new TGraph();
  acc[1] = new TGraph();
  acc[2] = new TGraph();
  tPrime = new TGraph();
  antennaPoint = inputAntenna;
}

void rad::FieldPoint::ResetFields() {
  EField[0]->Clear();
  EField[1]->Clear();
  EField[2]->Clear();
  BField[0]->Clear();
  BField[1]->Clear();
  BField[2]->Clear();
  pos[0]->Clear();
  pos[1]->Clear();
  pos[2]->Clear();
  vel[0]->Clear();
  vel[1]->Clear();
  vel[2]->Clear();
  acc[0]->Clear();
  acc[1]->Clear();
  acc[2]->Clear();
  tPrime->Clear();
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
    TVector3 EFieldCalc = CalcEField(antennaPoint, ePos, eVel, eAcc);
    TVector3 BFieldCalc = CalcBField(antennaPoint, ePos, eVel, eAcc);

    EField[0]->SetPoint(EField[0]->GetN(), time, EFieldCalc.X());
    EField[1]->SetPoint(EField[1]->GetN(), time, EFieldCalc.Y());
    EField[2]->SetPoint(EField[2]->GetN(), time, EFieldCalc.Z());
    BField[0]->SetPoint(BField[0]->GetN(), time, BFieldCalc.X());
    BField[1]->SetPoint(BField[1]->GetN(), time, BFieldCalc.Y());
    BField[2]->SetPoint(BField[2]->GetN(), time, BFieldCalc.Z());

    pos[0]->SetPoint(pos[0]->GetN(), time, xPos);
    pos[1]->SetPoint(pos[1]->GetN(), time, yPos);
    pos[2]->SetPoint(pos[2]->GetN(), time, zPos);
    vel[0]->SetPoint(vel[0]->GetN(), time, xVel);
    vel[1]->SetPoint(vel[1]->GetN(), time, yVel);
    vel[2]->SetPoint(vel[2]->GetN(), time, zVel);
    acc[0]->SetPoint(acc[0]->GetN(), time, xAcc);
    acc[1]->SetPoint(acc[1]->GetN(), time, yAcc);
    acc[2]->SetPoint(acc[2]->GetN(), time, zAcc);

    tPrime->SetPoint(tPrime->GetN(), time, CalcRetardedTime(antennaPoint, ePos, time));
  }

  fin->Close();
  delete fin;
}

TGraph* rad::FieldPoint::GetEFieldTimeDomain(Coord_t coord) {
  TGraph* gr = 0;
  if (coord == kX) {
    gr = (TGraph*)EField[0]->Clone("grEx");
    gr->GetYaxis()->SetTitle("E_{x} [V m^{-1}]");
  }
  else if (coord == kY) {
    gr = (TGraph*)EField[1]->Clone("grEy");
    gr->GetYaxis()->SetTitle("E_{y} [V m^{-1}]");
  }
  else if (coord = kZ) {
    gr = (TGraph*)EField[2]->Clone("grEz");
    gr->GetYaxis()->SetTitle("E_{z} [V m^{-1}]");
  }
  setGraphAttr(gr);
  gr->GetXaxis()->SetTitle("Time [s]");
  
  return gr;
}

TGraph* rad::FieldPoint::GetEFieldMagTimeDomain() {
  TGraph* grMag = new TGraph();
  assert((EField[0]->GetN() == EField[1]->GetN()) && (EField[1]->GetN() == EField[2]->GetN()));
  
  for (int i = 0; i < EField[0]->GetN(); i++) {
    double mag = sqrt( pow(EField[0]->GetPointY(i), 2) + pow(EField[1]->GetPointY(i), 2) + pow(EField[2]->GetPointY(i), 2) );
    grMag->SetPoint(grMag->GetN(), EField[0]->GetPointX(i), mag);
  }
  setGraphAttr(grMag);
  grMag->GetXaxis()->SetTitle("Time [s]");
  grMag->GetYaxis()->SetTitle("|E| [V m^{-1}]");
  
  return grMag;
}

TGraph* rad::FieldPoint::GetBFieldTimeDomain(Coord_t coord) {
  TGraph* gr = 0;
  if (coord == kX) {
    gr = (TGraph*)BField[0]->Clone("grBx");
    gr->GetYaxis()->SetTitle("B_{x} [T]");
  }
  else if (coord == kY) {
    gr = (TGraph*)BField[1]->Clone("grBy");
    gr->GetYaxis()->SetTitle("B_{y} [T]");
  }
  else if (coord = kZ) {
    gr = (TGraph*)BField[2]->Clone("grBz");
    gr->GetYaxis()->SetTitle("B_{z} [T]");
  }
  setGraphAttr(gr);
  gr->GetXaxis()->SetTitle("Time [s]");
  
  return gr;
}

TGraph* rad::FieldPoint::GetBFieldMagTimeDomain() {
  TGraph* grMag = new TGraph();
  assert((BField[0]->GetN() == BField[1]->GetN()) && (BField[1]->GetN() == BField[2]->GetN()));
  
  for (int i = 0; i < BField[0]->GetN(); i++) {
    double mag = sqrt( pow(BField[0]->GetPointY(i), 2) + pow(BField[1]->GetPointY(i), 2) + pow(BField[2]->GetPointY(i), 2) );
    grMag->SetPoint(grMag->GetN(), BField[0]->GetPointX(i), mag);
  }
  setGraphAttr(grMag);
  grMag->GetXaxis()->SetTitle("Time [s]");
  grMag->GetYaxis()->SetTitle("|B| [T]");
  
  return grMag;
}

TGraph* rad::FieldPoint::GetPoyntingMagTimeDomain() {
  TGraph* grSMag = new TGraph();
  TGraph* grEMag = GetEFieldMagTimeDomain();
  TGraph* grBMag = GetBFieldMagTimeDomain();
  assert(grEMag->GetN() == grBMag->GetN());

  for (int i = 0; i < grEMag->GetN(); i++) {
    double smag = grEMag->GetPointY(i) * grBMag->GetPointY(i) / MU0;
    grSMag->SetPoint(grSMag->GetN(), grEMag->GetPointX(i), smag);
  }
  setGraphAttr(grSMag);
  grSMag->GetXaxis()->SetTitle("Time [s]");
  grSMag->GetYaxis()->SetTitle("|S| [W m^{-2}]");
  
  return grSMag;
}

/*
  Frequency domain functions
*/

TGraph* rad::FieldPoint::GetEFieldPeriodogram(Coord_t coord) {
  TGraph* grFFT = 0;
  if (coord == kX) {
    grFFT = FFTtools::makePowerSpectrumPeriodogram(EField[0]);
    grFFT->GetYaxis()->SetTitle("E_{x}^{2} [V^{2} m^{-2}]");
  }
  else if (coord == kY) {
    grFFT = FFTtools::makePowerSpectrumPeriodogram(EField[1]);
    grFFT->GetYaxis()->SetTitle("E_{y}^{2} [V^{2} m^{-2}]");
  }
  else if (coord == kZ) {
    grFFT = FFTtools::makePowerSpectrumPeriodogram(EField[2]);
    grFFT->GetYaxis()->SetTitle("E_{z}^{2} [V^{2} m^{-2}]");
  }
  setGraphAttr(grFFT);
  grFFT->GetXaxis()->SetTitle("Frequency [Hz]");
  
  return grFFT;
}
 
TGraph* rad::FieldPoint::GetTotalEFieldPeriodogram() {
  TGraph *grTotal = new TGraph();
  TGraph *grX = GetEFieldPeriodogram(kX);
  TGraph *grY = GetEFieldPeriodogram(kY);
  TGraph *grZ = GetEFieldPeriodogram(kZ);

  for (int n = 0; n < grX->GetN(); n++) {
    double tot = grX->GetPointY(n) + grY->GetPointY(n) + grZ->GetPointY(n);
    grTotal->SetPoint(n, grX->GetPointX(n), tot);
  }
  setGraphAttr(grTotal);
  grTotal->GetXaxis()->SetTitle("Frequency [Hz]");
  
  return grTotal;
}

TGraph* rad::FieldPoint::GetDipolePowerTimeDomain() {
  TGraph *grPower = new TGraph();
  TGraph *grSMag = GetPoyntingMagTimeDomain();
  TVector3 dipoleDir(0.0, 1.0, 0.0);
  
  for (int i = 0; i < grSMag->GetN(); i++) {
    TVector3 ePos(pos[0]->GetPointY(i), pos[1]->GetPointY(i), pos[2]->GetPointY(i));
    double Ae = CalcAeHertzianDipole(0.0111, dipoleDir, ePos, antennaPoint);
    grPower->SetPoint(grPower->GetN(), grSMag->GetPointX(i), grSMag->GetPointY(i) * Ae);
  }
  setGraphAttr(grPower);
  
  delete grSMag;
  return grPower;
}
