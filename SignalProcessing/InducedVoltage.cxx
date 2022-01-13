// InducedVoltage.cxx

#include "SignalProcessing/InducedVoltage.h"
#include "Antennas/IAntenna.h"
#include "FieldClasses/FieldClasses.h"
#include "BasicFunctions/BasicFunctions.h"

#include "TString.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>

rad::InducedVoltage::~InducedVoltage() {
  delete grVoltage;
}

rad::InducedVoltage::InducedVoltage(TString trajectoryFilePath, IAntenna* myAntenna,
				    const bool kUseRetardedTime) {
  theFile = trajectoryFilePath;
  theAntenna = myAntenna;
  UseRetardedTime = kUseRetardedTime;
  grVoltage = new TGraph();
}

void rad::InducedVoltage::GenerateVoltage(double minTime, double maxTime) {
  FieldPoint fp(theFile, theAntenna);
  if (minTime == -1) minTime = 0.0;
  if (maxTime == -1) maxTime = fp.GetFinalTime();

  // To avoid running out of memory, generate the fields in more manageable chunks
  // Avoids having massive versions of unnecessary graphs
  const double chunkSize = 25e-6;
  double thisChunk = minTime + chunkSize;
  if (thisChunk > maxTime) thisChunk = maxTime;
  double lastChunk = minTime;
  std::cout<<"Generating voltages"<<std::endl;
  while (thisChunk <= maxTime && thisChunk != lastChunk) {
    fp.GenerateFields(lastChunk, thisChunk);
    TGraph* voltageTemp = fp.GetAntennaLoadVoltageTimeDomain(UseRetardedTime);
    
    // Account for antenna bandwidth
    // if (theAntenna->GetBandwidthUpperLimit() != DBL_MAX ||
    // 	theAntenna->GetBandwidthLowerLimit() != -DBL_MAX) {
    //   std::cout<<"Accounting for antenna bandwidth..."<<std::endl;
    //   voltageTemp = BandPassFilter(voltageTemp, theAntenna->GetBandwidthLowerLimit(), theAntenna->GetBandwidthUpperLimit());
    // }
    
    // Now write this to the main voltage graph
    std::cout<<"Writing to main voltage graph"<<std::endl;
    for (int i = 0; i < voltageTemp->GetN(); i++) {
      grVoltage->SetPoint(grVoltage->GetN(), voltageTemp->GetPointX(i), voltageTemp->GetPointY(i));
    }
    delete voltageTemp;
    lastChunk = thisChunk;
    thisChunk += chunkSize;
    if (thisChunk > maxTime) thisChunk = maxTime;
  }

}

rad::InducedVoltage::InducedVoltage(const InducedVoltage &iv) {
  grVoltage = (TGraph*)iv.grVoltage->Clone();
  theAntenna = iv.theAntenna;
  theFile = iv.theFile;
  UseRetardedTime = iv.UseRetardedTime;
}

void rad::InducedVoltage::ResetVoltage() {
  delete grVoltage;
  grVoltage = new TGraph();
}

double rad::InducedVoltage::GetFinalTime() {
  TFile *fin = new TFile(theFile, "READ");
  assert(fin);
  TTree* tree = (TTree*)fin->Get("tree");
  double lastTime;
  tree->SetBranchAddress("time", &lastTime);
  tree->GetEntry(tree->GetEntries()-1);
  delete tree;
  fin->Close();                                                                 
  delete fin;
  return lastTime;
}

double rad::InducedVoltage::GetUpperAntennaBandwidth() {
  return theAntenna->GetBandwidthUpperLimit();
}

double rad::InducedVoltage::GetLowerAntennaBandwidth() {
  return theAntenna->GetBandwidthLowerLimit();
}

void rad::InducedVoltage::ApplyAntennaBandwidth() {
  grVoltage = BandPassFilter(grVoltage, theAntenna->GetBandwidthLowerLimit(), theAntenna->GetBandwidthUpperLimit());
}

TGraph* rad::InducedVoltage::GetPowerPeriodogram(const double loadResistance) {
  TGraph* pgram = MakePowerSpectrumPeriodogram(grVoltage);
  ScaleGraph(grVoltage, 1.0 / loadResistance);
  return pgram;
}
