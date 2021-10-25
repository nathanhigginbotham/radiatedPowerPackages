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
  Ex = new TGraph();
  Ey = new TGraph();
  Ez = new TGraph();
  Bx = new TGraph();
  By = new TGraph();
  Bz = new TGraph();
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
}

rad::FieldPoint::~FieldPoint() {
  delete Ex;
  delete Ey;
  delete Ez;
  delete Bx;
  delete By;
  delete Bz;

  delete pos[0];
  delete pos[1];
  delete pos[2];
  delete vel[0];
  delete vel[1];
  delete vel[2];
  delete acc[0];
  delete acc[1];
  delete acc[2];
}

// Parametrised constructor
rad::FieldPoint::FieldPoint(TVector3 inputAntenna) {
  Ex = new TGraph();
  Ey = new TGraph();
  Ez = new TGraph();
  Bx = new TGraph();
  By = new TGraph();
  Bz = new TGraph();
  pos[0] = new TGraph();
  pos[1] = new TGraph();
  pos[2] = new TGraph();
  vel[0] = new TGraph();
  vel[1] = new TGraph();
  vel[2] = new TGraph();
  acc[0] = new TGraph();
  acc[1] = new TGraph();
  acc[2] = new TGraph();
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
    TVector3 BField = CalcBField(antennaPoint, ePos, eVel, eAcc);

    Ex->SetPoint(Ex->GetN(), time, EField.X());
    Ey->SetPoint(Ey->GetN(), time, EField.Y());
    Ez->SetPoint(Ez->GetN(), time, EField.Z());
    Bx->SetPoint(Bx->GetN(), time, BField.X());
    By->SetPoint(By->GetN(), time, BField.Y());
    Bz->SetPoint(Bz->GetN(), time, BField.Z());

    pos[0]->SetPoint(pos[0]->GetN(), time, xPos);
    pos[1]->SetPoint(pos[1]->GetN(), time, yPos);
    pos[2]->SetPoint(pos[2]->GetN(), time, zPos);
    vel[0]->SetPoint(vel[0]->GetN(), time, xVel);
    vel[1]->SetPoint(vel[1]->GetN(), time, yVel);
    vel[2]->SetPoint(vel[2]->GetN(), time, zVel);
    acc[0]->SetPoint(acc[0]->GetN(), time, xAcc);
    acc[1]->SetPoint(acc[1]->GetN(), time, yAcc);
    acc[2]->SetPoint(acc[2]->GetN(), time, zAcc);
  }

  fin->Close();
  delete fin;
}

TGraph* rad::FieldPoint::GetEFieldTimeDomain(Coord_t coord) {
  TGraph* gr = 0;
  if (coord == kX) {
    gr = (TGraph*)Ex->Clone("grEx");
    gr->GetYaxis()->SetTitle("E_{x} [V m^{-1}]");
  }
  else if (coord == kY) {
    gr = (TGraph*)Ey->Clone("grEy");
    gr->GetYaxis()->SetTitle("E_{y} [V m^{-1}]");
  }
  else if (coord = kZ) {
    gr = (TGraph*)Ez->Clone("grEz");
    gr->GetYaxis()->SetTitle("E_{z} [V m^{-1}]");
  }
  setGraphAttr(gr);
  gr->GetXaxis()->SetTitle("Time [s]");
  
  return gr;
}

TGraph* rad::FieldPoint::GetEFieldMagTimeDomain() {
  TGraph* grMag = new TGraph();
  assert((Ex->GetN() == Ey->GetN()) && (Ey->GetN() == Ez->GetN()));
  
  for (int i = 0; i < Ex->GetN(); i++) {
    double mag = sqrt( pow(Ex->GetPointY(i), 2) + pow(Ey->GetPointY(i), 2) + pow(Ez->GetPointY(i), 2) );
    grMag->SetPoint(grMag->GetN(), Ex->GetPointX(i), mag);
  }
  setGraphAttr(grMag);
  grMag->GetXaxis()->SetTitle("Time [s]");
  grMag->GetYaxis()->SetTitle("|E| [V m^{-1}]");
  
  return grMag;
}

TGraph* rad::FieldPoint::GetBFieldTimeDomain(Coord_t coord) {
  TGraph* gr = 0;
  if (coord == kX) {
    gr = (TGraph*)Bx->Clone("grBx");
    gr->GetYaxis()->SetTitle("B_{x} [T]");
  }
  else if (coord == kY) {
    gr = (TGraph*)By->Clone("grBy");
    gr->GetYaxis()->SetTitle("B_{y} [T]");
  }
  else if (coord = kZ) {
    gr = (TGraph*)Bz->Clone("grBz");
    gr->GetYaxis()->SetTitle("B_{z} [T]");
  }
  setGraphAttr(gr);
  gr->GetXaxis()->SetTitle("Time [s]");
  
  return gr;
}

TGraph* rad::FieldPoint::GetBFieldMagTimeDomain() {
  TGraph* grMag = new TGraph();
  assert((Bx->GetN() == By->GetN()) && (By->GetN() == Bz->GetN()));
  
  for (int i = 0; i < Bx->GetN(); i++) {
    double mag = sqrt( pow(Bx->GetPointY(i), 2) + pow(By->GetPointY(i), 2) + pow(Bz->GetPointY(i), 2) );
    grMag->SetPoint(grMag->GetN(), Bx->GetPointX(i), mag);
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
    grFFT = FFTtools::makePowerSpectrumPeriodogram(Ex);
    grFFT->GetYaxis()->SetTitle("E_{x}^{2} [V^{2} m^{-2}]");
  }
  else if (coord == kY) {
    grFFT = FFTtools::makePowerSpectrumPeriodogram(Ey);
    grFFT->GetYaxis()->SetTitle("E_{y}^{2} [V^{2} m^{-2}]");
  }
  else if (coord == kZ) {
    grFFT = FFTtools::makePowerSpectrumPeriodogram(Ez);
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
