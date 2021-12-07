/*
  InducedVoltage.h
  Lighter-weight representation of the voltage induced on an antenna compared to FieldPoint
  Just contains the actual output voltage across the terminals of the chosen antenna
*/

#ifndef INDUCED_VOLTAGE_H
#define INDUCED_VOLTAGE_H

#include "Antennas/IAntenna.h"

#include "TGraph.h"
#include "TString.h"

namespace rad
{
  class InducedVoltage
  {
  private:
    // Open circuit voltage with no signal processing performed
    TGraph* grVoltage;
    
  public:
    /// Constructor for a voltage
    /// \param trajectoryFilePath Path to the file containing the electron trajectory
    /// \param myAntenna Pointer to the chosen antenna
    /// \param minTime First time point to generate voltage from
    /// \param maxTime Maximum time period to generate voltage for
    /// \param kUseRetardedTime Boolean to select the use of the retarded time
    InducedVoltage(TString trajectoryFilePath, IAntenna* myAntenna,
		   double minTime=-1, double maxTime=-1,
		   const bool kUseRetardedTime=false);

    /// Destructor
    ~InducedVoltage();

    /// Copy constructor
    InducedVoltage(const InducedVoltage &iv);

    /// Returns the voltage graph
    TGraph* GetVoltageGraph() { return grVoltage; }
  };
}


#endif
