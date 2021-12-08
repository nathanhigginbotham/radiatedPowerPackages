// InducedVoltage.cxx

#include "SignalProcessing/InducedVoltage.h"
#include "Antennas/IAntenna.h"
#include "FieldClasses/FieldClasses.h"

#include "TString.h"
#include "TGraph.h"

#include <iostream>

rad::InducedVoltage::~InducedVoltage() {
  delete grVoltage;
}

rad::InducedVoltage::InducedVoltage(TString trajectoryFilePath, IAntenna* myAntenna,
				    double minTime, double maxTime, const bool kUseRetardedTime) {
  FieldPoint fp(trajectoryFilePath, myAntenna);
  if (minTime == -1) minTime = 0.0;
  if (maxTime == -1) maxTime = fp.GetFinalTime();

  // To avoid running out of memory, generate the fields in more manageable chunks
  // Avoids having massive versions of unnecessary graphs
  const double chunkSize = 25e-6;
  double thisChunk = minTime + chunkSize;
  double lastChunk = minTime;
  grVoltage = new TGraph();
  std::cout<<"Generating voltages"<<std::endl;
  while (thisChunk <= maxTime && thisChunk != lastChunk) {
    fp.GenerateFields(lastChunk, thisChunk);
    TGraph* voltageTemp = fp.GetAntennaLoadVoltageTimeDomain(kUseRetardedTime);
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
}
