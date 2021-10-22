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
    TGraph* Ex;
    TGraph* Ey;
    TGraph* Ez;
    TGraph* Bx;
    TGraph* By;
    TGraph* Bz;

    void ResetFields();    
    
  public:
    FieldPoint();
    FieldPoint(const TVector3 inputAntenna);
    ~FieldPoint();
    
    void GenerateFields(const char* inputFile, const double maxTime);

    TGraph* GetEx();
    TGraph* GetEy();
    TGraph* GetEz();
    TGraph* GetBx();
    TGraph* GetBy();
    TGraph* GetBz();
  };

}

#endif
