// BasicFunctions.h
#ifndef BASIC_FUNCTIONS_H
#define BASIC_FUNCTIONS_H

#include "TVector3.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "TGraph.h"
#include "FFTtools.h"
#include "TH1.h"
#include "TH2.h"

namespace rad
{

  ROOT::Math::XYZVector CalcEField(const ROOT::Math::XYZPoint fieldPoint, const ROOT::Math::XYZPoint ePosition, const ROOT::Math::XYZVector eVelocity, const ROOT::Math::XYZVector eAcceleration);

  ROOT::Math::XYZVector CalcBField(const ROOT::Math::XYZPoint fieldPoint, const ROOT::Math::XYZPoint ePosition, const ROOT::Math::XYZVector eVelocity, const ROOT::Math::XYZVector eAcceleration);
 
  ROOT::Math::XYZVector CalcPoyntingVec(const ROOT::Math::XYZPoint fieldPoint, const ROOT::Math::XYZPoint ePosition, const ROOT::Math::XYZVector eVelocity, const ROOT::Math::XYZVector eAcceleration);

  ROOT::Math::XYZVector CalcPoyntingVec(const ROOT::Math::XYZVector EField, const ROOT::Math::XYZVector BField);

  void setGraphAttr(TGraph *gr);
  
  void SetHistAttr(TH1 *h);
  
  void SetHistAttr(TH2 *h);

  double CalcAeHertzianDipole(const double wavelength, const ROOT::Math::XYZVector dipoleDir,
			      const ROOT::Math::XYZPoint ePosition,
			      const ROOT::Math::XYZPoint antennaPoint);
  
  double CalcAlHertzianDipole(const double wavelength, const ROOT::Math::XYZVector dipoleDir,
			      const ROOT::Math::XYZPoint ePosition,
			      const ROOT::Math::XYZPoint antennaPoint);
  
  double CalcRetardedTime(const ROOT::Math::XYZPoint fieldPoint, const ROOT::Math::XYZPoint ePosition, const double labTime);

  // Produces power spectrum with the desired normalisation
  TGraph* MakePowerSpectrumNorm(const TGraph* grWave);

  // Integrate the power spectrum
  double IntegratePowerNorm(const TGraph* grFFT, Int_t firstBin=-1, Int_t lastBin=-1);

  // Implements a simple band pass filter
  TGraph* BandPassFilter(const TGraph* grWave, const double minFreq, const double maxFreq);

  // PDF of the Rayleigh distribution
  double RayleighPDF(const double x, const double sigma);

  // CDF of the Rayleigh distribution
  double RayleighCDF(const double x, const double sigma);

  double RayleighPDFFunc(double *x, double *par);
  
  void AddWhiteNoiseFrequencyDomainPowerNorm(TGraph* grIn, const double Teff);
}
 
#endif
