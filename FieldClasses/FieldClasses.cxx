// FieldClasses.cxx

#include <cassert>
#include <cmath>
#include <vector>

#include "BasicFunctions/Constants.h"
#include "BasicFunctions/BasicFunctions.h"
#include "FieldClasses/FieldClasses.h"

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TVector3.h"
#include "TSpline.h"

#include "FFTtools.h"

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

// Default constructor
rad::FieldPoint::FieldPoint() {
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
  antennaPoint = TVector3(0, 0, 0);
  dipolePolarisation = TVector3(0, 0, 0);
  
  inputFile = "";
}

// Parametrised constructor
rad::FieldPoint::FieldPoint(const TVector3 inputAntenna, const TVector3 dipoleDir,
			    TString trajectoryFilePath) {
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
  dipolePolarisation = dipoleDir;
  
  // Now check that the input file exists
  TFile* f = new TFile(trajectoryFilePath, "read");
  assert(f);
  f->Close();
  delete f;
  inputFile = trajectoryFilePath;
}

// Copy constructor
rad::FieldPoint::FieldPoint(const FieldPoint &fp) {
  inputFile = fp.inputFile;
  antennaPoint = fp.antennaPoint;
  dipolePolarisation = fp.dipolePolarisation;
  // Clone the field graphs
  for (int coord = 0; coord < 3; coord++) {
    EField[coord] = (TGraph*)fp.EField[coord]->Clone();
    BField[coord] = (TGraph*)fp.BField[coord]->Clone();
    pos[coord] = (TGraph*)fp.pos[coord]->Clone();
    vel[coord] = (TGraph*)fp.vel[coord]->Clone();
    acc[coord] = (TGraph*)fp.acc[coord]->Clone();
  }
  tPrime = (TGraph*)fp.tPrime->Clone();
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

/*
  Takes the inputted TGraph and produces the same variable plotted using the retarded time
  grOriginal is the inputted TGraph
  returns a TGraph which takes the retarded time into account
 */
TGraph* rad::FieldPoint::MakeRetardedTimeGraph(const TGraph* grOriginal) {
  TSpline3 *sptPrime = new TSpline3("sptPrime", tPrime);
  TSpline3 *spgrOriginal = new TSpline3("spgrOriginal", grOriginal);

  TGraph *grOut = new TGraph();
  for (int i = 0; i < grOriginal->GetN(); i++) {
    if (grOriginal->GetPointX(i) < 0.2e-9) continue;
    double tRet = sptPrime->Eval(grOriginal->GetPointX(i));
    grOut->SetPoint(grOut->GetN(), grOriginal->GetPointX(i), spgrOriginal->Eval(tRet));
  }
  setGraphAttr(grOut);
  grOut->GetXaxis()->SetTitle("Time [s]");
  grOut->GetYaxis()->SetTitle(grOriginal->GetYaxis()->GetTitle());
  
  delete sptPrime;
  delete spgrOriginal;
  return grOut;
}

// From an input TFile generate the E and B fields for a given time
// maxTime is the final time in seconds (if less than the time in the file)
void rad::FieldPoint::GenerateFields(const double maxTime) {
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

TGraph* rad::FieldPoint::GetEFieldTimeDomain(Coord_t coord, const bool kUseRetardedTime,
					     int firstPoint, int lastPoint) {
  TGraph* gr = 0;
  TGraph* grOut = new TGraph();
  if (coord == kX) {
    gr = (TGraph*)EField[0]->Clone("grEx");
    grOut->GetYaxis()->SetTitle("E_{x} [V m^{-1}]");
  }
  else if (coord == kY) {
    gr = (TGraph*)EField[1]->Clone("grEy");
    grOut->GetYaxis()->SetTitle("E_{y} [V m^{-1}]");
  }
  else if (coord = kZ) {
    gr = (TGraph*)EField[2]->Clone("grEz");
    grOut->GetYaxis()->SetTitle("E_{z} [V m^{-1}]");
  }

  if (firstPoint < 0) firstPoint = 0;
  if (lastPoint < 0) lastPoint = gr->GetN() - 1;
  
  setGraphAttr(grOut);
  grOut->GetXaxis()->SetTitle("Time [s]");
  if (!kUseRetardedTime) {
    for (int i = firstPoint; i <= lastPoint; i++) {
      grOut->SetPoint(grOut->GetN(), gr->GetPointX(i), gr->GetPointY(i));
    }
  }
  else {
    TGraph* grRet = MakeRetardedTimeGraph(gr);
    for (int i = firstPoint; i <= lastPoint; i++) {
      grOut->SetPoint(grOut->GetN(), grRet->GetPointX(i), grRet->GetPointY(i));
    }
    delete grRet;
  }
  delete gr;
  return grOut;
}

TGraph* rad::FieldPoint::GetEFieldMagTimeDomain(const bool kUseRetardedTime) {
  TGraph* grMag = new TGraph();
  assert((EField[0]->GetN() == EField[1]->GetN()) && (EField[1]->GetN() == EField[2]->GetN()));
  
  for (int i = 0; i < EField[0]->GetN(); i++) {
    double mag = sqrt( pow(EField[0]->GetPointY(i), 2) + pow(EField[1]->GetPointY(i), 2) + pow(EField[2]->GetPointY(i), 2) );
    grMag->SetPoint(grMag->GetN(), EField[0]->GetPointX(i), mag);
  }
  setGraphAttr(grMag);
  grMag->GetXaxis()->SetTitle("Time [s]");
  grMag->GetYaxis()->SetTitle("|E| [V m^{-1}]");

  if (!kUseRetardedTime) {
    return grMag;
  }
  else {
    TGraph* grMagRet = MakeRetardedTimeGraph(grMag);
    delete grMag;
    return grMagRet;
  }
}

TGraph* rad::FieldPoint::GetBFieldTimeDomain(Coord_t coord, const bool kUseRetardedTime) {
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

  if (!kUseRetardedTime) {
    return gr;
  }
  else {
    TGraph* grRet = MakeRetardedTimeGraph(gr);
    delete gr;
    return grRet;
  }
}

TGraph* rad::FieldPoint::GetBFieldMagTimeDomain(const bool kUseRetardedTime) {
  TGraph* grMag = new TGraph();
  assert((BField[0]->GetN() == BField[1]->GetN()) && (BField[1]->GetN() == BField[2]->GetN()));
  
  for (int i = 0; i < BField[0]->GetN(); i++) {
    double mag = sqrt( pow(BField[0]->GetPointY(i), 2) + pow(BField[1]->GetPointY(i), 2) + pow(BField[2]->GetPointY(i), 2) );
    grMag->SetPoint(grMag->GetN(), BField[0]->GetPointX(i), mag);
  }
  setGraphAttr(grMag);
  grMag->GetXaxis()->SetTitle("Time [s]");
  grMag->GetYaxis()->SetTitle("|B| [T]");

  if (!kUseRetardedTime) {
    return grMag;
  }
  else {
    TGraph* grMagRet = MakeRetardedTimeGraph(grMag);
    delete grMag;
    return grMagRet;
  }
}

TGraph* rad::FieldPoint::GetPoyntingVecTimeDomain(Coord_t coord, const bool kUseRetardedTime) {
  TGraph* grS = new TGraph();
  for (int i = 0; i < EField[0]->GetN(); i++) {
    double comp = 0;
    if (coord == kX) {
      grS->GetYaxis()->SetTitle("S_{x} [W m^{-2}]");
      comp = EField[1]->GetPointY(i)*BField[2]->GetPointY(i) - EField[2]->GetPointY(i)*BField[1]->GetPointY(i);
    }
    else if (coord == kY) {
      comp = EField[2]->GetPointY(i)*BField[0]->GetPointY(i) - EField[0]->GetPointY(i)*BField[2]->GetPointY(i);
      grS->GetYaxis()->SetTitle("S_{y} [W m^{-2}]");
    }
    else if (coord == kZ) {
      comp = EField[0]->GetPointY(i)*BField[1]->GetPointY(i) - EField[1]->GetPointY(i)*BField[0]->GetPointY(i);
      grS->GetYaxis()->SetTitle("S_{z} [W m^{-2}]");
    }
    comp /= MU0;
    grS->SetPoint(i, EField[0]->GetPointX(i), comp);
  }
  setGraphAttr(grS);
  grS->GetXaxis()->SetTitle("Time [s]");

  if (!kUseRetardedTime) {
    return grS;
  }
  else {
    TGraph* grSRet = MakeRetardedTimeGraph(grS);
    delete grS;
    return grSRet;
  }
}

TGraph* rad::FieldPoint::GetPoyntingMagTimeDomain(const bool kUseRetardedTime) {
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

  if (!kUseRetardedTime) {
    return grSMag;
  }
  else {
    TGraph* grSMagRet = MakeRetardedTimeGraph(grSMag);
    delete grSMag;
    return grSMagRet;
  }
}

TGraph* rad::FieldPoint::GetDipolePowerTimeDomain(const bool kUseRetardedTime) {
  TGraph *grPower = new TGraph();
  TGraph *grSMag = GetPoyntingMagTimeDomain();
  TVector3 dipoleDir(0.0, 1.0, 0.0);
  
  for (int i = 0; i < grSMag->GetN(); i++) {
    TVector3 ePos(pos[0]->GetPointY(i), pos[1]->GetPointY(i), pos[2]->GetPointY(i));
    double Ae = CalcAeHertzianDipole(0.0111, dipoleDir, ePos, antennaPoint);
    grPower->SetPoint(grPower->GetN(), grSMag->GetPointX(i), grSMag->GetPointY(i) * Ae);
  }
  setGraphAttr(grPower);

  if (!kUseRetardedTime) {
    delete grSMag;
    return grPower;
  }
  else {
    TGraph* grPowerRet = MakeRetardedTimeGraph(grPower);
    delete grPower;
    delete grSMag;
    return grPowerRet;
  }
}

TGraph* rad::FieldPoint::GetDipoleComponentVoltageTimeDomain(Coord_t coord, const bool kUseRetardedTime,
							     int firstPoint, int lastPoint,
							     std::vector<GaussianNoise*> noiseTerms)
{
  TGraph* gr = new TGraph();
  TGraph* grE = GetEFieldTimeDomain(coord, false, firstPoint, lastPoint);
  gr->GetXaxis()->SetTitle("Time [s]");
  gr->GetYaxis()->SetTitle(grE->GetYaxis()->GetTitle());
  
  double fs = 1.0 / (grE->GetPointX(1) - grE->GetPointX(0));
  for (int term = 0; term < noiseTerms.size(); term++) {
    (noiseTerms.at(term))->SetSampleFreq(fs);
    (noiseTerms.at(term))->SetSigma();
  }
  
  TVector3 dipoleDir(0.0, 1.0, 0.0);
  for (int i = 0; i < grE->GetN(); i++) {
    TVector3 ePos(pos[0]->GetPointY(i), pos[1]->GetPointY(i), pos[2]->GetPointY(i));
    double Al = CalcAlHertzianDipole(0.0111, dipoleDir, ePos, antennaPoint);
    double voltage = grE->GetPointY(i) * Al;
    // Now add the noise
    for (int term = 0; term < noiseTerms.size(); term++) {
      voltage += (noiseTerms.at(term))->GetNoiseVoltage();
    }
    gr->SetPoint(gr->GetN(), grE->GetPointX(i), voltage);
  }
  setGraphAttr(gr);

  delete grE;
  
  if (!kUseRetardedTime) {
    return gr;
  }
  else {
    TGraph* grRet = MakeRetardedTimeGraph(gr);
    delete gr;
    return grRet;
  }
} 

TGraph* rad::FieldPoint::GetDipoleLoadVoltageTimeDomain(const bool kUseRetardedTime,
							int firstPoint, int lastPoint,
							std::vector<GaussianNoise*> noiseTerms) {
  TGraph* grEx = GetEFieldTimeDomain(kX, kUseRetardedTime, firstPoint, lastPoint);
  TGraph* grEy = GetEFieldTimeDomain(kY, kUseRetardedTime, firstPoint, lastPoint);
  TGraph* grEz = GetEFieldTimeDomain(kZ, kUseRetardedTime, firstPoint, lastPoint);

  TGraph* gr = new TGraph();
  gr->GetXaxis()->SetTitle("Time [s]");
  gr->GetYaxis()->SetTitle("Voltage [V]");
  setGraphAttr(gr);
  
  double fs = 1.0 / (grEx->GetPointX(1) - grEx->GetPointX(0));
  for (int term = 0; term < noiseTerms.size(); term++) {
    (noiseTerms.at(term))->SetSampleFreq(fs);
    (noiseTerms.at(term))->SetSigma();
  }

  for (int i = 0; i < grEx->GetN(); i++) {
    TVector3 EField(grEx->GetPointY(i), grEy->GetPointY(i), grEz->GetPointY(i));
    double voltage = EField.Dot(dipolePolarisation) * 0.0111 / TMath::Pi();
    voltage /= 2.0; // Account for re-radiated power
    // Now add noise
    for (int term = 0; term < noiseTerms.size(); term++) {
      voltage += (noiseTerms.at(term))->GetNoiseVoltage();
    }
    gr->SetPoint(gr->GetN(), grEx->GetPointX(i), voltage);
  }

  delete grEx;
  delete grEy;
  delete grEz;
  return gr;
}

////////////////////////////////////////////////////////////
/////////////// Frequency domain functions /////////////////
////////////////////////////////////////////////////////////

TGraph* rad::FieldPoint::GetEFieldPeriodogram(Coord_t coord, const bool kUseRetardedTime) {
  TGraph* grIn = GetEFieldTimeDomain(coord, kUseRetardedTime);
  TGraph* grFFT = FFTtools::makePowerSpectrumPeriodogram(grIn);
  
  if (coord == kX) {
    grFFT->GetYaxis()->SetTitle("E_{x}^{2} [V^{2} m^{-2}]");
  }
  else if (coord == kY) {
    grFFT->GetYaxis()->SetTitle("E_{y}^{2} [V^{2} m^{-2}]");
  }
  else if (coord == kZ) {
    grFFT->GetYaxis()->SetTitle("E_{z}^{2} [V^{2} m^{-2}]");
  }
  setGraphAttr(grFFT);
  grFFT->GetXaxis()->SetTitle("Frequency [Hz]");

  delete grIn;
  return grFFT;
}

TGraph* rad::FieldPoint::GetEFieldPowerSpectrumNorm(Coord_t coord, const bool kUseRetardedTime) {
  TGraph* grIn = GetEFieldTimeDomain(coord, kUseRetardedTime);
  TGraph* grFFT = MakePowerSpectrumNorm(grIn);
  
  if (coord == kX) {
    grFFT->GetYaxis()->SetTitle("E_{x}^{2} (#Deltat)^{2} [V^{2} m^{-2} s^{2}]");
  }
  else if (coord == kY) {
    grFFT->GetYaxis()->SetTitle("E_{y}^{2} (#Deltat)^{2} [V^{2} m^{-2} s^{2}]");
  }
  else if (coord == kZ) {
    grFFT->GetYaxis()->SetTitle("E_{z}^{2} (#Deltat)^{2} [V^{2} m^{-2} s^{2}]");
  }
  setGraphAttr(grFFT);
  grFFT->GetXaxis()->SetTitle("Frequency [Hz]");

  delete grIn;
  return grFFT;
}

// Total electric field in the frequency domain
TGraph* rad::FieldPoint::GetTotalEFieldPeriodogram(const bool kUseRetardedTime) {
  TGraph *grTotal = new TGraph();
  TGraph *grX = GetEFieldPeriodogram(kX, kUseRetardedTime);
  TGraph *grY = GetEFieldPeriodogram(kY, kUseRetardedTime);
  TGraph *grZ = GetEFieldPeriodogram(kZ, kUseRetardedTime);

  for (int n = 0; n < grX->GetN(); n++) {
    double tot = grX->GetPointY(n) + grY->GetPointY(n) + grZ->GetPointY(n);
    grTotal->SetPoint(n, grX->GetPointX(n), tot);
  }
  setGraphAttr(grTotal);
  grTotal->GetXaxis()->SetTitle("Frequency [Hz]");
  
  return grTotal;
}

TGraph* rad::FieldPoint::GetTotalEFieldPowerSpectrumNorm(const bool kUseRetardedTime) {
  TGraph *grTotal = new TGraph();
  TGraph *grX = GetEFieldPowerSpectrumNorm(kX, kUseRetardedTime);
  TGraph *grY = GetEFieldPowerSpectrumNorm(kY, kUseRetardedTime);
  TGraph *grZ = GetEFieldPowerSpectrumNorm(kZ, kUseRetardedTime);

  for (int n = 0; n < grX->GetN(); n++) {
    double tot = grX->GetPointY(n) + grY->GetPointY(n) + grZ->GetPointY(n);
    grTotal->SetPoint(n, grX->GetPointX(n), tot);
  }
  setGraphAttr(grTotal);
  grTotal->GetYaxis()->SetTitle("E^{2} (#Deltat)^{2} [V^{2} m^{-2} s^{2}]");
  grTotal->GetXaxis()->SetTitle("Frequency [Hz]");

  delete grX;
  delete grY;
  delete grZ;
  return grTotal;
}

TGraph* rad::FieldPoint::GetDipoleComponentVoltagePowerSpectrumNorm(Coord_t coord, const bool kUseRetardedTime, int firstPoint, int lastPoint, std::vector<GaussianNoise*> noiseTerms) {
  TGraph* grVTime = GetDipoleComponentVoltageTimeDomain(coord, kUseRetardedTime, firstPoint, lastPoint, noiseTerms);
  TGraph* grPower = MakePowerSpectrumNorm(grVTime);
  
  if (coord == kX) {
    grPower->GetYaxis()->SetTitle("V_{x}^{2} (#Deltat)^{2} [V^{2} s^{2}]");
  }
  else if (coord == kY) {
    grPower->GetYaxis()->SetTitle("V_{y}^{2} (#Deltat)^{2} [V^{2} s^{2}]");
  }
  else if (coord == kZ) {
    grPower->GetYaxis()->SetTitle("V_{z}^{2} (#Deltat)^{2} [V^{2} s^{2}]");
  }
  setGraphAttr(grPower);
  grPower->GetXaxis()->SetTitle("Frequency [Hz]");
  
  delete grVTime;
  return grPower;
}

TGraph* rad::FieldPoint::GetDipoleTotalVoltagePowerSpectrumNorm(const bool kUseRetardedTime,
								int firstPoint, int lastPoint,
								std::vector<GaussianNoise*> noiseTerms) {
  TGraph* grTotal = new TGraph();
  TGraph *grX = GetDipoleComponentVoltagePowerSpectrumNorm(kX, kUseRetardedTime, firstPoint, lastPoint, noiseTerms);
  TGraph *grY = GetDipoleComponentVoltagePowerSpectrumNorm(kY, kUseRetardedTime, firstPoint, lastPoint, noiseTerms);
  TGraph *grZ = GetDipoleComponentVoltagePowerSpectrumNorm(kZ, kUseRetardedTime, firstPoint, lastPoint, noiseTerms);
  for (int n = 0; n < grX->GetN(); n++) {
    double tot = grX->GetPointY(n) + grY->GetPointY(n) + grZ->GetPointY(n);
    grTotal->SetPoint(n, grX->GetPointX(n), tot);
  }
  setGraphAttr(grTotal);
  grTotal->GetYaxis()->SetTitle("V^{2} (#Deltat)^{2} [V^{2} s^{2}]");
  grTotal->GetXaxis()->SetTitle("Frequency [Hz]");

  delete grX;
  delete grY;
  delete grZ;
  return grTotal;
}

TGraph* rad::FieldPoint::GetDipolePowerSpectrumNorm(const bool kUseRetardedTime,
						    int firstPoint, int lastPoint,
						    std::vector<GaussianNoise*> noiseTerms) {
  TGraph* grVoltagePower = GetDipoleTotalVoltagePowerSpectrumNorm(kUseRetardedTime, firstPoint, lastPoint, noiseTerms);
  TGraph* grDipolePower = new TGraph();

  for (int i = 0; i < grVoltagePower->GetN(); i++) {
    double powerWatts = grVoltagePower->GetPointY(i) * TMath::C() * EPSILON0;
    grDipolePower->SetPoint(i, grVoltagePower->GetPointX(i), powerWatts);
  }
  setGraphAttr(grDipolePower);
  grDipolePower->GetXaxis()->SetTitle("Frequency [Hz]");
  grDipolePower->GetYaxis()->SetTitle("Power [W]");
  
  delete grVoltagePower;
  return grDipolePower;
}

TGraph* rad::FieldPoint::GetDipoleLoadPowerSpectrumNorm(const double resistance,
							const bool kUseRetardedTime,
							int firstPoint, int lastPoint,
							std::vector<GaussianNoise*> noiseTerms) {
  TGraph* grVoltage = GetDipoleLoadVoltageTimeDomain(kUseRetardedTime, firstPoint, lastPoint, noiseTerms);
  TGraph* grPower = MakePowerSpectrumNorm(grVoltage);

  for (int i = 0; i < grPower->GetN(); i++) {
    double powerWatts = grPower->GetPointY(i) / resistance;
    grPower->SetPoint(i, grPower->GetPointX(i), powerWatts);
  }
  setGraphAttr(grPower);
  grPower->GetXaxis()->SetTitle("Frequency [Hz]");
  grPower->GetYaxis()->SetTitle("Power #times (#Delta t)^{2} [W s^{2}]");
  
  delete grVoltage;
  return grPower;
}

// Assorted useful functions
double rad::FieldPoint::GetFinalTime() {
  TFile *fin = new TFile(inputFile, "READ");
  assert(fin);
  TTree* tree = (TTree*)fin->Get("tree");
  double lastTime;
  tree->SetBranchAddress("time", &lastTime);
  tree->GetEntry(tree->GetEntries()-1);
  fin->Close();
  delete tree;
  delete fin;
  return lastTime;
}

double rad::FieldPoint::GetSampleRate() {
  TFile *fin = new TFile(inputFile, "READ");
  assert(fin);
  TTree* tree = (TTree*)fin->Get("tree");
  double time;
  tree->SetBranchAddress("time", &time);
  double time0, time1;
  tree->GetEntry(0);
  time0 = time;
  tree->GetEntry(1);
  time1 = time;
  double fs = 1.0 / (time1 - time0);
  delete tree;
  fin->Close();
  delete fin;
  return fs;
}
