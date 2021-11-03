// FieldClasses.h
#ifndef FIELD_CLASSES_H
#define FIELD_CLASSES_H

#include <vector>

#include "TFile.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TString.h"

#include "BasicFunctions/BasicFunctions.h"
#include "SignalProcessing/NoiseFunc.h"

namespace rad
{
  class FieldPoint
  {
  private:
    TString inputFile;
    
    TVector3 antennaPoint;
    TVector3 dipolePolarisation;
    // Time series for field components
    TGraph* EField[3];
    TGraph* BField[3];

    TGraph* pos[3];
    TGraph* vel[3];
    TGraph* acc[3];

    TGraph* tPrime; // Relationship between time and retarded time
    
    void ResetFields();

    TGraph* MakeRetardedTimeGraph(const TGraph* grOriginal);
    
  public:
    enum Coord_t{
      kX, kY, kZ
    };
    FieldPoint();
    FieldPoint(const TVector3 inputAntenna, const TVector3 dipoleDir, TString trajectoryFilePath);
    ~FieldPoint();
    
    void GenerateFields(const double maxTime);

    TGraph* GetEFieldTimeDomain(Coord_t coord, const bool kUseRetardedTime=false,
				int firstPoint=-1, int lastPoint=-1);
    TGraph* GetBFieldTimeDomain(Coord_t coord, const bool kUseRetardedTime=false);
    TGraph* GetPoyntingVecTimeDomain(Coord_t coord, const bool kUseRetardedTime=false);
    
    TGraph* GetEFieldMagTimeDomain(const bool kUseRetardedTime=false);
    TGraph* GetBFieldMagTimeDomain(const bool kUseRetardedTime=false);
    TGraph* GetPoyntingMagTimeDomain(const bool kUseRetardedTime=false);

    TGraph* GetDipoleComponentVoltageTimeDomain(Coord_t coord, const bool kUseRetardedTime=false,
						int firstPoint=-1, int lastPoint=-1,
						std::vector<GaussianNoise*> noiseTerms={});
    TGraph* GetDipoleLoadVoltageTimeDomain(const bool kUseRetardedTime=false,
					   int firstPoint=-1, int lastPoint=-1,
					   std::vector<GaussianNoise*> noiseTerms={});
    TGraph* GetDipolePowerTimeDomain(const bool kUseRetardedTime=false);

    // Frequency domain functions
    TGraph* GetEFieldPeriodogram(Coord_t coord, const bool kUseRetardedTime=false);
    TGraph* GetTotalEFieldPeriodogram(const bool kUseRetardedTime=false);

    // These have a normalisation such that they shouldn't be altered by varying
    // things such as zero padding and the sample rate
    TGraph* GetEFieldPowerSpectrumNorm(Coord_t coord, const bool kUseRetardedTime=false);
    TGraph* GetTotalEFieldPowerSpectrumNorm(const bool kUseRetardedTime=false);

    TGraph* GetDipoleComponentVoltagePowerSpectrumNorm(Coord_t coord, const bool kUseRetardedTime=false,
						       int firstPoint=-1, int lastPoint=-1,
						       std::vector<GaussianNoise*> noiseTerms={});
    TGraph* GetDipoleTotalVoltagePowerSpectrumNorm(const bool kUseRetardedTime=false,
						   int firstPoint=-1, int lastPoint=-1,
						   std::vector<GaussianNoise*> noiseTerms={});

    // Calculate the power collected by a Hertzian dipole in the frequency domain
    TGraph* GetDipolePowerSpectrumNorm(const bool kUseRetardedTime=false,
				       int firstPoint=-1, int lastPoint=-1,
				       std::vector<GaussianNoise*> noiseTerms={});

    // Functions for various useful things such as the final time in a file
    double GetFinalTime();
    double GetSampleRate();
  };

}

#endif
