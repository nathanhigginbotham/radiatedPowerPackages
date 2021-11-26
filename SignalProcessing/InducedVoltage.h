/*
  InducedVoltage.h
  Lighter-weight representation of the voltage induced on an antenna compared to FieldPoint
  Just contains the actual output voltage across the terminals of the chosen antenna
*/

#ifndef INDUCED_VOLTAGE_H
#define INDUCED_VOLTAGE_H

#include "Antennas/IAntenna.h"

#include "TGraph.h"

namespace rad
{
  class InducedVoltage
  {
  private:
    TGraph* grVoltage;
    
  public:
    /// Constructor for a voltage
    /// \param trajectoryFilePath Path to the file containing the electron trajectory
    /// \param myAntenna Pointer to the chosen antenna
    InducedVoltage(TString trajectoryFilePath, IAntenna* myAntenna);
    
    ~InducedVoltage();
    
  };
}


#endif
