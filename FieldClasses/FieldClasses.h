// FieldClasses.h
#ifndef FIELD_CLASSES_H
#define FIELD_CLASSES_H

#include "TFile.h"
#include "TGraph.h"
#include "TVector3.h"

#include "BasicFunctions/BasicFunctions.h"

namespace rad
{
  class FieldPoint
  {
  private:
    TVector3 antennaPoint;
    // Time series for field components
    TGraph* EField[3];
    TGraph* BField[3];

    TGraph* pos[3];
    TGraph* vel[3];
    TGraph* acc[3];

    TGraph* tPrime; // Relationship between time and retarded time
    
    void ResetFields();    
    
  public:
    enum Coord_t{
      kX, kY, kZ
    };
    FieldPoint();
    FieldPoint(const TVector3 inputAntenna);
    ~FieldPoint();
    
    void GenerateFields(const char* inputFile, const double maxTime);

    TGraph* GetEFieldTimeDomain(Coord_t coord);
    TGraph* GetBFieldTimeDomain(Coord_t coord);
   
    TGraph* GetEFieldMagTimeDomain();
    TGraph* GetBFieldMagTimeDomain();
    TGraph* GetPoyntingMagTimeDomain();
    
    TGraph* GetEFieldPeriodogram(Coord_t coord);
    TGraph* GetTotalEFieldPeriodogram();

    TGraph* GetDipolePowerTimeDomain();
  };

}

#endif
