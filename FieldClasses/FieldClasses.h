// FieldClasses.h
#ifndef FIELD_CLASSES_H
#define FIELD_CLASSES_H

#include <vector>

#include "TFile.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TString.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"

#include "BasicFunctions/BasicFunctions.h"
#include "SignalProcessing/NoiseFunc.h"
#include "Antennas/IAntenna.h"

namespace rad
{
  class FieldPoint
  {
  private:
    TString inputFile;
    
    // Time series for field components
    TGraph* EField[3];
    TGraph* BField[3];

    // Time series of electron position
    TGraph* pos[3];

    TGraph* tPrime; // Relationship between time and retarded time
    
    void ResetFields();

    TGraph* MakeRetardedTimeGraph(const TGraph* grOriginal);

    IAntenna* myAntenna; // The chosen antenna, contains position and direction info
    
  public:
    enum Coord_t{
      kX, kY, kZ
    };

    FieldPoint(TString trajectoryFilePath, IAntenna* myAntenna);
    
    ~FieldPoint();
    
    FieldPoint(const FieldPoint &fp);
    
    void GenerateFields(const double maxTime);

    TGraph* GetPositionTimeDomain(Coord_t coord, const bool kUseRetardedTime=false,
				  int firstPoint=-1, int lastPoint=-1);
    
    TGraph* GetEFieldTimeDomain(Coord_t coord, const bool kUseRetardedTime=false,
				int firstPoint=-1, int lastPoint=-1);
    TGraph* GetBFieldTimeDomain(Coord_t coord, const bool kUseRetardedTime=false);
    TGraph* GetPoyntingVecTimeDomain(Coord_t coord, const bool kUseRetardedTime=false);
    
    TGraph* GetEFieldMagTimeDomain(const bool kUseRetardedTime=false);
    TGraph* GetBFieldMagTimeDomain(const bool kUseRetardedTime=false);
    TGraph* GetPoyntingMagTimeDomain(const bool kUseRetardedTime=false);
    
    TGraph* GetAntennaLoadVoltageTimeDomain(const bool kUseRetardedTime=false,
					    int firstPoint=-1, int lastPoint=-1);
    TGraph* GetAntennaLoadPowerTimeDomain(const double loadResistance,
					  const bool kUseRetardedTime=false,
					  int firstPoint=-1, int lastPoint=-1);

    // Frequency domain functions
    TGraph* GetEFieldPeriodogram(Coord_t coord, const bool kUseRetardedTime=false);
    TGraph* GetTotalEFieldPeriodogram(const bool kUseRetardedTime=false);

    // These have a normalisation such that they shouldn't be altered by varying
    // things such as zero padding and the sample rate
    TGraph* GetEFieldPowerSpectrumNorm(Coord_t coord, const bool kUseRetardedTime=false);
    TGraph* GetTotalEFieldPowerSpectrumNorm(const bool kUseRetardedTime=false);
    
    TGraph* GetAntennaLoadPowerSpectrumNorm(const double resistance,
					    const bool kUseRetardedTime=false,
					    int firstPoint=-1, int lastPoint=-1,
					    std::vector<GaussianNoise*> noiseTerms={});
    
    // Functions for various useful things such as the final time in a file
    double GetFinalTime();
    double GetSampleRate();
  };

}

#endif
