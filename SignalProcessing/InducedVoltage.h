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
    TString theFile;
    TGraph* grVoltage;
    IAntenna* theAntenna;
    bool UseRetardedTime;
    
  public:
    /// Constructor for a voltage
    /// \param trajectoryFilePath Path to the file containing the electron trajectory
    /// \param myAntenna Pointer to the chosen antenna
    /// \param kUseRetardedTime Boolean to select the use of the retarded time
    InducedVoltage(TString trajectoryFilePath, IAntenna* myAntenna,
		   const bool kUseRetardedTime=false);

    /// Generates the voltage between the specified times
    /// \param minTime First time point to generate voltage from
    /// \param maxTime Maximum time period to generate voltage for
    void GenerateVoltage(double minTime=-1, double maxTime=-1);
    
    /// Destructor
    ~InducedVoltage();

    /// Copy constructor
    InducedVoltage(const InducedVoltage &iv);

    /// Returns the voltage graph
    TGraph* GetVoltageGraph() { return grVoltage; }

    /// Deletes the TGraph pointer to free up memory
    void ResetVoltage();

    // Returns the input file name
    TString GetFilename() { return theFile; }

    // Returns the final time in the file
    double GetFinalTime();
  };
}


#endif
