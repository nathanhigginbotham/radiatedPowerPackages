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

#include <vector>

namespace rad
{
  class InducedVoltage
  {
  private:
    // Open circuit voltage with no signal processing performed
    TString theFile;
    TGraph* grVoltage;
    std::vector<IAntenna*> theAntennas;
    bool UseRetardedTime;
    double chunkSize;

    /// Produce a shifted version of a voltage graph according to the delay of the antenna
    TGraph* DelayVoltage(TGraph* grIn, IAntenna* ant);
    
  public:
    /// Constructor for a voltage
    /// \param trajectoryFilePath Path to the file containing the electron trajectory
    /// \param myAntenna Pointer to the chosen antenna
    /// \param kUseRetardedTime Boolean to select the use of the retarded time
    InducedVoltage(TString trajectoryFilePath, IAntenna* myAntenna,
		   const bool kUseRetardedTime=false);

    /// Constructor taking a vector of antennas as an argument
    /// \param trajectoryFilePath Path to the file containing the electron trajectory
    /// \param antennaVec Vector of antenna points to generate the voltages at
    /// \param kUseRetardedTime Boolean to select the use of the retarded time
    InducedVoltage(TString trajectoryFilePath, std::vector<IAntenna*> antennaVec,
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
    TGraph* GetVoltageGraph();

    /// Deletes the TGraph pointer to free up memory
    void ResetVoltage();

    // Returns the input file name
    TString GetFilename() { return theFile; }

    // Returns the final time in the file
    double GetFinalTime();

    // Returns the upper antenna bandwidth
    double GetUpperAntennaBandwidth();

    // Returns the lower antenna bandwidth limit
    double GetLowerAntennaBandwidth();

    // Applies the antenna bandwidth to the signal
    void ApplyAntennaBandwidth();

    /// Returns periodogram
    /// \param loadResistance The resistance of the load circuit
    TGraph* GetPowerPeriodogram(const double loadResistance);

    /// Returns the chunk size used in generating the voltage
    double GetChunkSize() { return chunkSize; }
  };
}


#endif
