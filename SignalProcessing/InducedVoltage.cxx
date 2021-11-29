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
				    const double maxTime, const bool kUseRetardedTime) {
  FieldPoint fp(trajectoryFilePath, myAntenna);
  fp.GenerateFields(maxTime);
  std::cout<<"Generating voltages"<<std::endl;
  grVoltage = fp.GetAntennaLoadVoltageTimeDomain(kUseRetardedTime);
}

rad::InducedVoltage::InducedVoltage(const InducedVoltage &iv) {
  grVoltage = (TGraph*)iv.grVoltage->Clone();
}
