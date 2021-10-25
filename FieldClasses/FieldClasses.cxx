// FieldClasses.cxx

#include <cassert>

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
  Ex = new TGraph();
  Ey = new TGraph();
  Ez = new TGraph();
  Bx = new TGraph();
  By = new TGraph();
  Bz = new TGraph();
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
